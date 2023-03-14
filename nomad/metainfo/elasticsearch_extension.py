#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

'''
This elasticsearch extension for the Metainfo allows to define how quantities are
added to Elasticsearch indices.

This extension supports two search indices: ``entry_index`` and ``material_index``.
There are three different types of "searchable documents": ``entry_type``, ``material_type``,
``material_entry_type``. Entry documents are indexed in the entry index; material documents
in the material index. The material entry documents are nested documents in material documents.

The document types are subsets of the metainfo schema; documents have the exact same
structure as archives, but with only some of the quantities. Which quantities are in these
documents can be defined in the metainfo by using the :class:`Elasticsearch` annotation on
quantity definitions.

Entry and material entry documents start with the metainfo entry root section.
Material documents start with the ``results.material`` sub-section. Nested material entry
documents are placed under the ``entries`` key within a material document. This is the only
exception, where the material document structure deviates from the metainfo/archive structure.

A quantity can appear in multiple document types. All indexed quantities
appear by default in entry documents. If specified quantities are also put in either
the material document or a nested material entry document within the material document.
The material quantities describe the material itself
(e.g. formula, elements, system type, symmetries). These quantities are always in all
entries of the same material. Material entry quantities describe individual results
and metadata that are contributed by the entries of this material (e.g. published, embargo,
band gap, available properties). The values contributed by different entries of the same
material may vary.

Here is a small metainfo example:

.. code-block:: python

    class Entry(MSection):

        entry_id = Quantity(
            type=str,
            a_elasticsearch=Elasticsearch(material_entry_type))

        upload_create_time = Quantity(
            type=Datetime,
            a_elasticsearch=Elasticsearch())

        results = SubSection(sub_section=Results.m_def, a_elasticsearch=Elasticsearch())


    class Results(MSection):

        material = SubSection(sub_section=Material.m_def, a_elasticsearch=Elasticsearch())
        properties = SubSection(sub_section=Properties.m_def, a_elasticsearch=Elasticsearch())


    class Material(MSection):

        material_id = Quantity(
            type=str,
            a_elasticsearch=Elasticsearch(material_type))

        formula = Quantity(
            type=str,
            a_elasticsearch=[
                Elasticsearch(material_type),
                Elasticsearch(material_type, field='text', mapping='text')])


    class Properties(MSection):

        available_properties = Quantity(
            type=str, shape=['*'],
            a_elasticsearch=Elasticsearch(material_entry_type))

        band_gap = Quantity(
            type=float, unit='J',
            a_elasticsearch=Elasticsearch(material_entry_type))


The resulting indices with a single entry in them would look like this. Entry index:

.. code-block:: json

    [
        {
            "entry_id": "de54f1",
            "upload_create_time": "2021-02-01 01:23:12",
            "results": {
                "material": {
                    "material_id": "23a8bf",
                    "formula": "H2O"
                },
                "properties": {
                    "available_properties": ["dos", "bs", "band_gap", "energy_total_0"],
                    "band_gap": 0.283e-12
                }
            }
        }
    ]


And material index:

.. code-block:: json

    [
        {
            "material_id": "23a8bf",
            "formula": "H2O"
            "entries": [
                {
                    "entry_id": "de54f1",
                    "results": {
                        "properties": {
                            "available_properties": ["dos", "bs", "band_gap", "energy_total_0"],
                            "band_gap": 0.283e-12
                        }
                    }
                }
            ]
        }
    ]

You can freely define sub-sections and quantities. The only fixed structures that are
required from the metainfo are:
- the root section has an ``entry_id``
- materials are placed in ``results.material``
- the ``results.material`` sub-section has a ``material_id``
- the ``results.material`` sub-section has no property called ``entries``

This extension resolves references during indexing and basically treats referenced
sub-sections as if they were direct sub-sections.

.. autofunction:: index_entry
.. autofunction:: index_entries
.. autofunction:: create_indices


.. autoclass:: Elasticsearch
.. autoclass:: DocumentType
.. autoclass:: Index
'''


from typing import Union, Any, Dict, cast, Set, List, Callable, Tuple, DefaultDict
from collections import defaultdict
import numpy as np
import re

from nomad import config, utils

from .metainfo import (
    MSectionBound, Section, Quantity, MSection, MEnum, Datetime, Reference, DefinitionAnnotation,
    Definition, QuantityReference, Unit)


class DocumentType():
    '''
    DocumentType allows to create Elasticsearch index mappings and documents based on
    Metainfo definitions and instances. Genrally this class should not be used outside
    the elasticsearch_extension module.

    Attributes:
        root_section_def: The section definition that serves as the root for all documents.
            mapping: The elasticsearch mapping definition.
        indexed_properties: All definitions (quantities and sub sections) that are covered
            by documents of this type.
        quantities: All elasticsearch quantities that in documents of this type. A dictionary
            with full qualified name as key and :class:`Elasticsearch` annotations as
            values.
        metrics: All metrics in this document type. A dictionary with metric names as
            keys and tuples of elasticsearch metric aggregation and respective
            :class:`Elasticsearch` metainfo annotation as values.
        id_field: The quantity (and elasticsearch field) name that is used as unique
            identifier for this type of documents.
    '''
    def __init__(self, name: str, id_field: str):
        self.name = name
        self.id_field = id_field
        self.root_section_def = None
        self.mapping: Dict[str, Any] = None
        self.indexed_properties: Set[Definition] = set()
        self.nested_object_keys: List[str] = []
        self.nested_sections: List[SearchQuantity] = []
        self.quantities: Dict[str, SearchQuantity] = {}
        self.suggestions: Dict[str, Elasticsearch] = {}
        self.metrics: Dict[str, Tuple[str, SearchQuantity]] = {}

    def _reset(self):
        self.indexed_properties.clear()
        self.nested_object_keys.clear()
        self.nested_sections.clear()
        self.quantities.clear()
        self.metrics.clear()

    def create_index_doc(self, root: MSection):
        '''
        Creates an indexable document from the given archive.
        '''
        suggestions: DefaultDict = defaultdict(list)

        def transform(quantity, section, value, path):
            '''
            Custom transform function that possibly transforms the indexed
            values and also gathers the suggestion values for later storage.
            '''
            elasticsearch_annotations = quantity.m_get_annotations(Elasticsearch, as_list=True)
            for elasticsearch_annotation in elasticsearch_annotations:
                if elasticsearch_annotation.field is None:
                    if elasticsearch_annotation.suggestion:
                        # The suggestions may have a different doc_type: we
                        # don't serialize them if the doc types don't match.
                        if self != entry_type and elasticsearch_annotation.doc_type != self:
                            continue

                        # The suggestion values are saved into a temporary
                        # dictionary. The actual path of the data in the
                        # metainfo is used as a key. The suggestions will also
                        # take into account any given variants of the input.
                        transform_function = elasticsearch_annotation.value
                        variants = elasticsearch_annotation.variants
                        if transform_function is not None:
                            if variants:
                                suggestion_value = []
                                for variant in variants(value):
                                    suggestion_value.extend(transform_function(variant))
                                suggestion_value = list(set(suggestion_value))
                            else:
                                suggestion_value = transform_function(value)
                        else:
                            suggestion_value = value
                        section_path = section.m_path()[len(root.m_path()):]
                        name = elasticsearch_annotation.property_name
                        if path:
                            suggestion_path = f"{section_path}/{path}/{name}"
                        else:
                            suggestion_path = f"{section_path}/{name}"
                        suggestions[suggestion_path].extend(suggestion_value)
                    else:
                        transform_function = elasticsearch_annotation.value
                        if transform_function is not None:
                            return transform_function(section)

            return value

        def exclude(property_, section):
            if property_ not in self.indexed_properties:
                return True

            return False

        kwargs: Dict[str, Any] = dict(
            with_meta=False,
            include_defaults=True,
            include_derived=True,
            resolve_references=True,
            exclude=exclude,
            transform=transform
        )

        result = root.m_to_dict(**kwargs)

        # Add the collected suggestion values
        for path, value in suggestions.items():
            try:
                parts = path.split("/")
                section = result
                for part in parts[:-1]:
                    if part == "":
                        continue
                    try:
                        part = int(part)
                    except ValueError:
                        pass
                    section = section[part]
                section[parts[-1]] = value
            except KeyError:
                # TODO This typically happens when a suggestion is stored in a
                # referenced object that is in the metadata/results sections,
                # e.g. referenced authored that are stored in EELS DB
                # measurements
                logger = utils.get_logger(__name__, doc_type=self.name)
                logger.warn("could not add suggestion to es", path=path)

        # TODO deal with metadata
        metadata = result.get('metadata')
        if metadata is not None:
            del(result['metadata'])
            result.update(**metadata)

        return result

    def create_mapping(
            self, section_def: Section, prefix: str = None,
            auto_include_subsections: bool = False):
        '''
        Creates an Elasticsearch mapping for the given root section. It traverses all
        sub-sections to create the mapping. It will not create the mapping for nested
        documents. These have to be created manually (e.g. by :func:`create_indices`).
        Will override the existing mapping.

        Arguments:
            section_def: The section definition to create a mapping for.
            prefix: The qualified name of the section within the search index. This
                is used to create the qualified names of quantities and sub-sections.
            auto_include_subsections: Considers all sub and sub sub sections regardless
                of any annotation in the sub section definitions. By default only
                sub sections with elasticsearch annotation are traversed.
        '''
        mappings: Dict[str, Any] = {}

        if self == material_type and prefix is None:
            mappings['n_entries'] = {'type': 'integer'}

        for quantity_def in section_def.all_quantities.values():
            elasticsearch_annotations = quantity_def.m_get_annotations(Elasticsearch, as_list=True)
            for elasticsearch_annotation in elasticsearch_annotations:
                if self != entry_type and elasticsearch_annotation.doc_type != self:
                    continue

                if prefix is None:
                    qualified_name = quantity_def.name
                else:
                    qualified_name = f'{prefix}.{quantity_def.name}'
                if elasticsearch_annotation.suggestion:
                    self.suggestions[qualified_name] = elasticsearch_annotation
                is_section_reference = isinstance(quantity_def.type, Reference)
                is_section_reference &= not isinstance(quantity_def.type, QuantityReference)
                if is_section_reference:
                    # Treat referenced sections as sub-sections
                    assert quantity_def.type.target_section_def is not None
                    # TODO e.g. viewers, entry_coauthors, etc. ... should be treated as multiple inner docs
                    # assert quantity_def.is_scalar

                    reference_mapping = self.create_mapping(
                        cast(Section, quantity_def.type.target_section_def),
                        prefix=qualified_name)
                    if len(reference_mapping['properties']) > 0:
                        mappings[quantity_def.name] = reference_mapping
                else:
                    mapping = mappings.setdefault(elasticsearch_annotation.property_name, {})
                    fields = elasticsearch_annotation.fields
                    if len(fields) > 0:
                        mapping.setdefault('fields', {}).update(**fields)

                    else:
                        mapping.update(**elasticsearch_annotation.mapping)

                self.indexed_properties.add(quantity_def)
                self._register(elasticsearch_annotation, prefix)

        for sub_section_def in section_def.all_sub_sections.values():
            annotation = sub_section_def.m_get_annotations(Elasticsearch)
            if annotation is None and not auto_include_subsections:
                continue

            assert not isinstance(annotation, list), \
                'sub sections can only have one elasticsearch annotation'
            continue_with_auto_include_subsections = auto_include_subsections or (
                False if annotation is None else annotation.auto_include_subsections)

            if prefix is None:
                qualified_name = sub_section_def.name
            else:
                qualified_name = f'{prefix}.{sub_section_def.name}'

            # TODO deal with metadata
            qualified_name = re.sub(r'\.?metadata', '', qualified_name)
            qualified_name = None if qualified_name == '' else qualified_name

            sub_section_mapping = self.create_mapping(
                sub_section_def.sub_section, prefix=qualified_name,
                auto_include_subsections=continue_with_auto_include_subsections)

            nested = annotation is not None and annotation.nested
            if nested:
                sub_section_mapping['type'] = 'nested'

            if len(sub_section_mapping['properties']) > 0:
                if sub_section_def.name == 'metadata':
                    mappings.update(**sub_section_mapping['properties'])
                else:
                    mappings[sub_section_def.name] = sub_section_mapping
                self.indexed_properties.add(sub_section_def)
                if nested and qualified_name not in self.nested_object_keys:
                    self.nested_object_keys.append(qualified_name)

                    search_quantity = SearchQuantity(annotation=annotation, doc_type=self, prefix=prefix)
                    self.nested_sections.append(search_quantity)
                    self.nested_object_keys.sort(key=lambda item: len(item))

        self.mapping = dict(properties=mappings)
        return self.mapping

    def _register(self, annotation, prefix):
        search_quantity = SearchQuantity(annotation=annotation, doc_type=self, prefix=prefix)
        name = search_quantity.qualified_name

        assert name not in self.quantities or self.quantities[name] == search_quantity, \
            'Search quantity names must be unique: %s' % name

        self.quantities[name] = search_quantity

        if annotation.metrics is not None:
            for name, metric in annotation.metrics.items():
                assert name not in self.metrics, 'Metric names must be unique: %s' % name
                self.metrics[name] = (metric, search_quantity)

        if self == entry_type:
            annotation.search_quantity = search_quantity

    def __repr__(self):
        return self.name


class Index():
    '''
    Allows to access an Elasticsearch index. It forwards method calls to Python's
    Elasticsearch package for Elasticsearch document APIs like search, index, get, mget,
    bulk, etc. It adds the necessary doc_type and index parameters for you.

    Arguments:
        doc_type: The :class:`DocumentType` instance that describes the document type
            of this index.
        index_config_key: The ``nomad.config.elastic`` config key that holds the name
            for this index.

    Attributes:
        elastic_client: The used Elasticsearch Python package client.
        index_name: The name of the index in Elasticsearch.
    '''
    def __init__(self, doc_type: DocumentType, index_config_key: str):
        self.doc_type = doc_type
        self.index_config_key = index_config_key

    def __elasticsearch_operation(self, name: str, *args, **kwargs):
        if 'index' not in kwargs:
            kwargs['index'] = self.index_name

        results = getattr(self.elastic_client, name)(*args, **kwargs)
        return results

    @property
    def index_name(self):
        return getattr(config.elastic, self.index_config_key)

    @property
    def elastic_client(self):
        from nomad.infrastructure import elastic_client
        return elastic_client

    def __getattr__(self, name):
        if name not in ['get', 'index', 'mget', 'bulk', 'search']:
            return super().__getattribute__(name)

        def wrapper(*args, **kwargs):
            return self.__elasticsearch_operation(name, *args, **kwargs)

        return wrapper

    def create_index(self, upsert: bool = False):
        ''' Initially creates the index with the mapping of its document type. '''
        assert self.doc_type.mapping is not None, 'The mapping has to be created first.'
        logger = utils.get_logger(__name__, index=self.index_name)
        if not self.elastic_client.indices.exists(index=self.index_name):
            # TODO the settings emulate the path_analyzer used in the v0 elasticsearch_dsl
            # based index, configured by nomad.datamodel.datamodel::path_analyzer
            self.elastic_client.indices.create(index=self.index_name, body={
                'settings': {
                    'analysis': {
                        'analyzer': {
                            'path_analyzer': {'tokenizer': 'path_tokenizer', 'type': 'custom'}
                        },
                        'tokenizer': {
                            'path_tokenizer': {'pattern': '/', 'type': 'pattern'}
                        }
                    }
                },
                'mappings': self.doc_type.mapping
            })
            logger.info('elasticsearch index created')
        elif upsert:
            self.elastic_client.indices.put_mapping(
                index=self.index_name,
                body=self.doc_type.mapping)
            logger.info('elasticsearch index updated')
        else:
            logger.info('elasticsearch index exists')

    def delete(self):
        if self.elastic_client.indices.exists(index=self.index_name):
            self.elastic_client.indices.delete(index=self.index_name)

    def refresh(self):
        self.elastic_client.indices.refresh(index=self.index_name)


# TODO type 'doc' because it's the default used by elasticsearch_dsl and the v0 entries index.
# 'entry' would be more descriptive.
entry_type = DocumentType('doc', id_field='entry_id')
material_type = DocumentType('material', id_field='material_id')
material_entry_type = DocumentType('material_entry', id_field='entry_id')

entry_index = Index(entry_type, index_config_key='entries_index')
material_index = Index(material_type, index_config_key='materials_index')


def get_tokenizer(regex):
    '''Returns a function that tokenizes a given string using the provided
    regular epression.
    '''
    def tokenizer(value):
        tokens = []
        if value:
            tokens.append(value)
            for match in re.finditer(regex, value):
                if (match):
                    token = value[match.end():]
                    if token != '':
                        # Notice how we artificially extend the token by taking the
                        # prefix and adding it at the end. This way the token
                        # remains unique so that it will be returned by
                        # ElasticSearch when "skip_duplicates" is used in the
                        # query.
                        tokens.append(f'{token} {value[:match.end()]}')
        return tokens

    return tokenizer


tokenizer_default = get_tokenizer(r'[_\s\.\/-]+')


class Elasticsearch(DefinitionAnnotation):
    '''
    A metainfo annotation for quantity definitions. This annotation can be used multiple
    times on the same quantity (e.g. to define Elasticsearch fields with differrent mapping
    types). Each annotation will create a field in the respective elasticsearch document type.

    This annotation has to be used on all sub sections that lead to quantities that should
    be included. On sub sections an inner document mapping is applied and all other
    arguments are ignored.

    Arguments:
        doc_type: An additional document type: ``material_type`` or ``material_entry_type``.
            All quantities with this annotation are automatically placed in ``entry_type``.
        mapping: The Elasticsearch mapping for the underlying elasticsearch field. The
            default depends on the quantity type. You can provide the elasticsearch type
            name, a full dictionary with additional elasticsearch mapping parameters, or
            an elasticsearch_dsl mapping object.
        field: Allows to specify sub-field name. There has to be another annotation on the
            same quantity with the default name. The custom field name is concatenated
            to the default. This will create an additional mapping for this
            quantity. In queries this can be used like an additional field, but the
            quantity is only stored once (under the quantity name) in the source document.
        value:
            A callable that is applied to the containering section to get a value for
            this quantity when saving the section in the elastic search index. By default
            this will be the serialized quantity value.
        index:
            A boolean that indicates if this quantity should be indexed or merely be
            part of the elastic document ``_source`` without being indexed for search.
        values:
            If the quantity is used in aggregations for a fixed set of values,
            use this parameter to preset these values. On aggregation, elasticsearch
            will only return values that exist in the search results. This allows to
            create 0 statistic values and return consistent set of values. If the underlying
            quantity is an Enum, the values are determined automatically.
        default_aggregation_size:
            The of values to return by default if this quantity is used in aggregation.
            If no value is given and there are not fixed value, 10 will be used.
        metrics:
            If the quantity is used as a metric for aggregating, this has to
            be used to define a valid elasticsearch metrics aggregations, e.g.
            'sum' or 'cardinality'. It is a dictionary with metric name as key,
            and elasticsearch aggregation name as values.
        many_all:
            Multiple values can be used to search based on a property. If no operator
            is given by the user, a logical or is used, i.e., all those entries that
            have at least one the values are selected (any operator). Set many_all to true
            to alter the default behavior to a logical and (all operator). If the user
            provides an operator, the provided operator is used regardless. Usually
            many_all is only sensible for properties that can have multiple value per
            entry (e.g. elements).
        auto_include_subsections:
            If true all sub and sub sub sections are considered for search even if
            there are no elasticsearch annotations in the sub section definitions.
            By default only sub sections with elasticsearch annotation are considered
            during index mapping creation.
        nested:
            If true the section is mapped to elasticsearch nested object and all queries
            become nested queries. Only applicable to sub sections.
        suggestion:
            Controls the suggestions that are built for this field. Leave
            undefined if no suggestions are required. Can be a custom callable
            that transforms a string into a list of suggestion values, or one
            of the preset strings:
            - simple: Only the value is stored as an ES field.
            - default: The value is split into tokens using whitespace, dot and
              forward slash
        variants:
            A callable that is applied to a search value to get a list of
            alternative forms of the input. Used to augment the available
            suggestions with alternative forms.
        normalizer:
            A callable that is used to transform the search input when
            targeting this field. Note that this does not affect the way the
            value is indexed.
        es_query: The Elasticsearch query type that is used when querying for the annotated
            quantity, e.g. match, term, match_phrase. Default is 'match'.

    Attributes:
        name:
            The name of the quantity (plus additional field if set).
        search_quantity: The entry type SearchQuantity associated with this annotation.
    '''
    def __init__(
            self,
            doc_type: DocumentType = entry_type,
            mapping: Union[str, Dict[str, Any]] = None,
            field: str = None,
            es_field: str = None,
            value: Callable[[MSectionBound], Any] = None,
            index: bool = True,
            values: List[str] = None,
            default_aggregation_size: int = None,
            metrics: Dict[str, str] = None,
            many_all: bool = False,
            auto_include_subsections: bool = False,
            nested: bool = False,
            suggestion: Union[str, Callable[[MSectionBound], Any]] = None,
            variants: Union[Callable[[str], List[str]]] = None,
            normalizer: Callable[[Any], Any] = None,
            es_query: str = 'match',
            _es_field: str = None):

        # TODO remove _es_field if it is not necessary anymore to enforce a specific mapping
        # for v0 compatibility
        if suggestion:
            if doc_type != entry_type:
                raise ValueError('Suggestions should only be stored in the entry index.')
            for arg in [field, mapping, es_field, _es_field]:
                if arg is not None:
                    raise ValueError(f'You cannot modify the way suggestions are mapped or named.')
            # If no tokenizer is specified, the suggestion is stored as a field
            # that holds only the original value.
            if suggestion == 'simple':
                field = 'suggestion'
            elif suggestion == 'default':
                value = tokenizer_default
            elif callable(suggestion):
                value = suggestion
            else:
                raise ValueError(
                    'Please provide the suggestion as one of the predefined '
                    'shortcuts, False or a custom callable.'
                )

        if variants and not callable(variants):
            raise ValueError('Please provide the variants as a custom callable.')
        self.variants = variants

        if normalizer and not callable(normalizer):
            raise ValueError('Please provide the normalizer as a custom callable.')
        self.normalizer = normalizer

        self._custom_mapping = mapping
        self.field = field
        self.es_query = es_query
        self._es_field = field if _es_field is None else _es_field
        self.doc_type = doc_type
        self.value = value
        self.index = index
        self._mapping: Dict[str, Any] = None

        self.default_aggregation_size = default_aggregation_size
        self.values = values
        self.metrics = metrics
        self.many_all = many_all

        self.auto_include_subsections = auto_include_subsections
        self.nested = nested
        self.suggestion = suggestion

        self.search_quantity = None

    @property
    def values(self):
        return self._values

    @values.setter
    def values(self, value):
        self._values = value
        if self.default_aggregation_size is None and self._values is not None:
            self.default_aggregation_size = len(self._values)

    @property
    def mapping(self) -> Dict[str, Any]:
        if self._mapping is not None:
            return self._mapping

        if self.suggestion:
            from elasticsearch_dsl import Completion
            # The standard analyzer will retain numbers unlike the simple
            # analyzer which is the default.
            self._mapping = Completion(analyzer='standard').to_dict()
            return self._mapping

        if self._custom_mapping is not None:
            from elasticsearch_dsl import Field

            if isinstance(self._custom_mapping, Field):
                self._mapping = self._custom_mapping.to_dict()
            elif isinstance(self._custom_mapping, str):
                self._mapping = dict(type=self._custom_mapping)
            else:
                self._mapping = self._custom_mapping

            return self._mapping

        def compute_mapping(quantity: Quantity) -> Dict[str, Any]:
            if quantity.type == str:
                return dict(type='keyword')
            elif quantity.type in [float, np.float64]:
                return dict(type='double')
            elif quantity.type == np.float32:
                return dict(type='float')
            elif quantity.type in [np.int64]:
                return dict(type='long')
            elif quantity.type in [int, np.int32]:
                return dict(type='integer')
            elif quantity.type == bool:
                return dict(type='boolean')
            elif quantity.type == Datetime:
                return dict(type='date')
            elif quantity.type == Unit:
                return dict(type='keyword')
            elif isinstance(quantity.type, QuantityReference):
                return compute_mapping(quantity.type.target_quantity_def)
            elif isinstance(quantity.type, Reference):
                raise NotImplementedError('Resolving section references is not supported.')
            elif isinstance(quantity.type, MEnum):
                return dict(type='keyword')
            else:
                raise NotImplementedError(
                    'Quantity type %s for quantity %s is not supported.' % (quantity.type, quantity))

        self._mapping = compute_mapping(cast(Quantity, self.definition))

        if not self.index:
            self._mapping['index'] = False

        return self._mapping

    @property
    def fields(self) -> Dict[str, Any]:
        if self._es_field == '' or self._es_field is None:
            return {}

        return {
            self._es_field: self.mapping
        }

    @property
    def property_name(self) -> str:
        if self.suggestion and not self.field:
            return f'{self.definition.name}__suggestion'
        return self.definition.name

    @property
    def name(self) -> str:
        if self.field is not None:
            return f'{self.property_name}.{self.field}'
        elif self.suggestion:
            return f'{self.property_name}__suggestion'
        else:
            return self.property_name

    def __repr__(self):
        if self.definition is None:
            return super().__repr__()

        return f'Elasticsearch({self.definition})'

    def m_to_dict(self):
        if self.search_quantity:
            return self.search_quantity.qualified_name
        else:
            return self.name


class SearchQuantity():
    '''
    This is used to represent search quantities. It is different from a metainfo quantity
    because the same metainfo quantity can appear multiple times at different places in
    an archive (an search index document). A search quantity is uniquely identified by
    a qualified name that pin points its place in the sub-section hierarchy.

    Arguments:
        annotation: The elasticsearch annotation that this search quantity is based on.
        doc_type: The elasticsearch document type that this search quantity appears in.
        prefix: The prefix to build the full qualified name for this search quantity.

    Attributes:
        qualified_field:
            The full qualified name of the resulting elasticsearch field in the entry
            document type. This will be the quantity name (plus additional
            field or suggestion postfix if set) with subsection names up to the
            root of the metainfo data.
        search_field:
            The full qualified name of the field in the elasticsearch index.
        qualified_name:
            Same name as qualified_field. This will be used to address the search
            property in our APIs.
        definition: The metainfo quantity definition that this search quantity is based on
        aggregatable:
            A boolean that determines, if this quantity can be used in aggregations.
    '''
    def __init__(self, annotation: Elasticsearch, doc_type: DocumentType, prefix: str):
        self.annotation = annotation
        self.doc_type = DocumentType

        qualified_field = self.annotation.definition.name

        if prefix is not None:
            qualified_field = f'{prefix}.{qualified_field}'

        if annotation.suggestion:
            qualified_field = f'{qualified_field}__suggestion'

        self.search_field = qualified_field
        if not(annotation._es_field == '' or annotation._es_field is None):
            self.search_field = f'{qualified_field}.{annotation._es_field}'

        if annotation.field is not None:
            qualified_field = f'{qualified_field}.{annotation.field}'

        self.qualified_field = qualified_field
        self.qualified_name = qualified_field

    @property
    def definition(self):
        return self.annotation.definition

    @property
    def aggregatable(self):
        if isinstance(self.definition.type, Reference):
            return False

        field_type = self.annotation.mapping['type']
        return field_type == 'keyword' or field_type == 'boolean'

    def __repr__(self):
        if self.definition is None:
            return super().__repr__()

        return f'SearchQuantity({self.qualified_field})'

    def __getattr__(self, name):
        return getattr(self.annotation, name)


def create_indices(entry_section_def: Section = None, material_section_def: Section = None):
    '''
    Creates the mapping for all document types and creates the indices in Elasticsearch.
    The indices must not exist already. Prior created mappings will be replaced.
    '''
    if entry_section_def is None:
        from nomad.datamodel import EntryArchive
        entry_section_def = EntryArchive.m_def

    if material_section_def is None:
        from nomad.datamodel.results import Material
        material_section_def = Material.m_def

    entry_type._reset()
    material_type._reset()
    material_entry_type._reset()

    entry_type.create_mapping(entry_section_def)
    material_type.create_mapping(material_section_def, auto_include_subsections=True)
    material_entry_type.create_mapping(entry_section_def, prefix='entries')

    # Here we manually add the material_entry_type mapping as a nested field
    # inside the material index. We also need to manually specify the
    # additional nested fields that come with this: the entries + all
    # nested_object_keys from material_entry_type. Notice that we need to sort
    # the list: the API expects a list sorted by name length in ascending
    # order.
    material_entry_type.mapping['type'] = 'nested'
    material_type.mapping['properties']['entries'] = material_entry_type.mapping
    material_type.nested_object_keys += ['entries'] + material_entry_type.nested_object_keys
    material_type.nested_object_keys.sort(key=lambda item: len(item))

    entry_index.create_index(upsert=True)  # TODO update the existing v0 index
    material_index.create_index()


def delete_indices():
    entry_index.delete()
    material_index.delete()


def index_entry(entry: MSection, **kwargs):
    '''
    Upserts the given entry in the entry index. Optionally updates the materials index
    as well.
    '''
    index_entries([entry], **kwargs)


def index_entries_with_materials(entries: List, refresh: bool = False):
    index_entries(entries, refresh=refresh)
    update_materials(entries, refresh=refresh)


def index_entries(entries: List, refresh: bool = False) -> Dict[str, str]:
    '''
    Upserts the given entries in the entry index. Optionally updates the materials index
    as well. Returns a dictionary of the format {entry_id: error_message} for all entries
    that failed to index.
    '''
    rv = {}
    # split into reasonably sized problems
    if len(entries) > config.elastic.bulk_size:
        for entries_part in [entries[i:i + config.elastic.bulk_size] for i in range(0, len(entries), config.elastic.bulk_size)]:
            errors = index_entries(entries_part, refresh=refresh)
            rv.update(errors)
        return rv

    if len(entries) == 0:
        return rv

    logger = utils.get_logger('nomad.search', n_entries=len(entries))

    with utils.timer(logger, 'prepare bulk index of entries actions and docs'):
        actions_and_docs = []
        for entry in entries:
            try:
                entry_index_doc = entry_type.create_index_doc(entry)

                actions_and_docs.append(dict(index=dict(_id=entry['entry_id'])))
                actions_and_docs.append(entry_index_doc)
            except Exception as e:
                logger.error('could not create entry index doc', entry_id=entry['entry_id'], exc_info=e)

        timer_kwargs: Dict[str, Any] = {}
        try:
            import json
            timer_kwargs['size'] = len(json.dumps(actions_and_docs))
            timer_kwargs['n_actions'] = len(actions_and_docs)
        except Exception:
            pass

    with utils.timer(
        logger, 'perform bulk index of entries',
        lnr_event='failed to bulk index entries',
        **timer_kwargs
    ):
        indexing_result = entry_index.bulk(
            body=actions_and_docs, refresh=refresh,
            timeout=f'{config.elastic.bulk_timeout}s',
            request_timeout=config.elastic.bulk_timeout)
        # Extract only the errors from the indexing_result
        if indexing_result['errors']:
            for item in indexing_result['items']:
                if item['index']['status'] >= 400:
                    rv[item['index']['_id']] = str(item['index']['error'])
        return rv


def update_materials(entries: List, refresh: bool = False):
    # split into reasonably sized problems
    if len(entries) > config.elastic.bulk_size:
        for entries_part in [entries[i:i + config.elastic.bulk_size] for i in range(0, len(entries), config.elastic.bulk_size)]:
            update_materials(entries_part, refresh=refresh)
        return

    if len(entries) == 0:
        return

    logger = utils.get_logger('nomad.search', n_entries=len(entries))

    def get_material_id(entry):
        material_id = None
        try:
            material_id = entry.results.material.material_id
        except AttributeError:
            pass
        return material_id

    # Get all entry and material ids.
    entry_ids, material_ids = set(), set()
    entries_dict = {}
    for entry in entries:
        entries_dict[entry.entry_id] = entry
        entry_ids.add(entry.entry_id)
        material_id = get_material_id(entry)
        if material_id is not None:
            material_ids.add(material_id)

    logger = logger.bind(n_materials=len(material_ids))

    # Get existing materials for entries' material ids (i.e. the entry needs to be added
    # or updated).
    with utils.timer(logger, 'get existing materials', lnr_event='failed to get existing materials'):
        if material_ids:
            elasticsearch_results = material_index.mget(
                body={
                    'docs': [dict(_id=material_id) for material_id in material_ids]
                },
                request_timeout=config.elastic.bulk_timeout)
            existing_material_docs = [
                doc['_source'] for doc in elasticsearch_results['docs'] if '_source' in doc]
        else:
            existing_material_docs = []

    # Get old materials that still have one of the entries, but the material id has changed
    # (i.e. the materials where entries need to be removed due entries having different
    # materials now).
    with utils.timer(logger, 'get old materials', lnr_event='failed to get old materials'):
        elasticsearch_results = material_index.search(body={
            'size': len(entry_ids),
            'query': {
                'bool': {
                    'must': {
                        'nested': {
                            'path': 'entries',
                            'query': {
                                'terms': {
                                    'entries.entry_id': list(entry_ids)
                                }
                            }
                        }
                    },
                    'must_not': {
                        'terms': {
                            'material_id': list(material_ids)
                        }
                    }
                }
            }
        })
        old_material_docs = [hit['_source'] for hit in elasticsearch_results['hits']['hits']]

    # Compare and create the appropriate materials index actions
    # First, we go through the existing materials. The following cases need to be covered:
    # - an entry needs to be updated within its existing material (standard case)
    # - an entry needs to be added to an existing material (new entry case)
    # - there is an entry with no existing material (new material case)
    # - there is an entry that moves from one existing material to another (super rare
    #   case where an entry's material id changed within the set of other entries' material ids)
    # This n + m complexity with n=number of materials and m=number of entries

    # We create lists of bulk operations. Each list only contains enough materials to
    # have the ammount of entries in all these materials roughly match the desired bulk size.
    # Using materials as a measure might not be good enough, if a single material has
    # lots of nested entries.
    _actions_and_docs_bulks: List[List[Any]] = []
    _n_entries_in_bulk = [0]

    def add_action_or_doc(action_or_doc):
        if len(_actions_and_docs_bulks) == 0 or _n_entries_in_bulk[0] > config.elastic.bulk_size:
            _n_entries_in_bulk[0] = 0
            _actions_and_docs_bulks.append([])
        _actions_and_docs_bulks[-1].append(action_or_doc)
        if 'entries' in action_or_doc:
            _n_entries_in_bulk[0] = _n_entries_in_bulk[0] + len(action_or_doc['entries'])

    material_docs = []
    material_docs_dict = {}
    remaining_entry_ids = set(entry_ids)
    for material_doc in existing_material_docs:
        material_id = material_doc['material_id']
        material_docs_dict[material_id] = material_doc
        material_entries = material_doc['entries']
        material_entries_to_remove = []
        for index, material_entry in enumerate(material_entries):
            entry_id = material_entry['entry_id']
            entry = entries_dict.get(entry_id)
            if entry is None:
                # The entry was not changed.
                continue
            else:
                # Update the material, there might be slight changes even if it is made
                # from entry properties that are "material defining", e.g. changed external
                # material quantities like new AFLOW prototypes
                try:
                    material_doc.update(**material_type.create_index_doc(entry.results.material))
                except Exception as e:
                    logger.error('could not create material index doc', exc_info=e)

            new_material_id = get_material_id(entry)
            if new_material_id != material_id:
                # Remove the entry, it moved to another material. But the material cannot
                # run empty, because another entry had this material id.
                material_entries_to_remove.append(index)
            else:
                # Update the entry.
                try:
                    material_entries[index] = material_entry_type.create_index_doc(entry)
                except Exception as e:
                    logger.error('could not create material index doc', exc_info=e)
                remaining_entry_ids.remove(entry_id)
        for index in reversed(material_entries_to_remove):
            del(material_entries[index])

        add_action_or_doc(dict(index=dict(_id=material_id)))
        add_action_or_doc(material_doc)
        material_docs.append(material_doc)

    for entry_id in remaining_entry_ids:
        entry = entries_dict.get(entry_id)
        material_id = get_material_id(entry)
        if material_id is not None:
            material_doc = material_docs_dict.get(material_id)
            if material_doc is None:
                # The material does not yet exist. Create it.
                try:
                    material_doc = material_type.create_index_doc(entry.results.material)
                except Exception as e:
                    logger.error('could not create material index doc', exc_info=e)
                material_docs_dict[material_id] = material_doc
                add_action_or_doc(dict(create=dict(_id=material_id)))
                add_action_or_doc(material_doc)
                material_docs.append(material_doc)
            # The material does exist (now), but the entry is new.
            try:
                material_doc.setdefault('entries', []).append(material_entry_type.create_index_doc(entry))
            except Exception as e:
                logger.error('could not create material entry index doc', exc_info=e)

    # Second, we go through the old materials. The following cases need to be covered:
    # - the old materials are empty (standard case)
    # - an entry needs to be removed but the material still has entries (new material id case 1)
    # - an entry needs to be removed and the material is now "empty" (new material id case 2)
    for material_doc in old_material_docs:
        material_id = material_doc['material_id']
        material_entries = material_doc['entries']
        material_entries_to_remove = []
        for index, material_entry in enumerate(material_entries):
            entry_id = material_entry['entry_id']
            if entry_id in entry_ids:
                # The entry does not belong to this material anymore and needs to be removed.
                material_entries_to_remove.append(index)
        for index in reversed(material_entries_to_remove):
            del(material_entries[index])
        if len(material_entries) == 0:
            # The material is empty now and needs to be removed.
            add_action_or_doc(dict(delete=dict(_id=material_id)))
        else:
            # The material needs to be updated
            add_action_or_doc(dict(index=dict(_id=material_id)))
            add_action_or_doc(material_doc)
            material_docs.append(material_doc)

    # Third, we potentially cap the number of entries in a material. We ensure that only
    # a certain amounts of entries are stored with all metadata. The rest will only
    # have their entry id.
    all_n_entries_capped = 0
    all_n_entries = 0
    for material_doc in material_docs:
        material_entries = material_doc.get('entries', [])
        material_doc['n_entries'] = len(material_entries)
        if len(material_entries) > config.elastic.entries_per_material_cap:
            material_doc['entries'] = material_entries[0:config.elastic.entries_per_material_cap]

        all_n_entries_capped += len(material_entries)
        all_n_entries += material_doc['n_entries']

    # Execute the created actions in bulk.
    timer_kwargs: Dict[str, Any] = {}
    try:
        import json
        timer_kwargs['size'] = len(json.dumps(_actions_and_docs_bulks))
        timer_kwargs['n_actions'] = sum([len(bulk) for bulk in _actions_and_docs_bulks])
        timer_kwargs['n_entries'] = all_n_entries
        timer_kwargs['n_entries_capped'] = all_n_entries_capped
    except Exception:
        pass

    with utils.timer(
        logger, 'perform bulk index of materials',
        lnr_event='failed to bulk index materials',
        **timer_kwargs
    ):
        for bulk in _actions_and_docs_bulks:
            material_index.bulk(
                body=bulk, refresh=False,
                timeout=f'{config.elastic.bulk_timeout}s',
                request_timeout=config.elastic.bulk_timeout)

    if refresh:
        entry_index.refresh()
        material_index.refresh()
