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

        upload_time = Quantity(
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
            "upload_time": "2021-02-01 01:23:12",
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


from typing import Union, Any, Dict, cast, Set, List, Callable, Tuple
import numpy as np
import re

from nomad import config, utils

from .metainfo import (
    MSectionBound, Section, Quantity, MSection, MEnum, Datetime, Reference, DefinitionAnnotation,
    Definition, QuantityReference)


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
        self.nested_object_keys: List[str] = list()
        self.quantities: Dict[str, SearchQuantity] = {}
        self.metrics: Dict[str, Tuple[str, SearchQuantity]] = {}

    def _reset(self):
        self.indexed_properties.clear()
        self.nested_object_keys.clear()
        self.quantities.clear()
        self.metrics.clear()

    def create_index_doc(self, root: MSection):
        '''
        Creates an indexable document from the given archive.
        '''
        def transform(quantity, section, value):
            '''
            Custom transform function for m_to_dict that will resolve references
            on a per-quantity basis based on the annotation setup.
            '''
            elasticsearch_annotations = quantity.m_get_annotations(Elasticsearch, as_list=True)
            for elasticsearch_annotation in elasticsearch_annotations:
                if elasticsearch_annotation.field is None:
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

        for quantity_def in section_def.all_quantities.values():
            elasticsearch_annotations = quantity_def.m_get_annotations(Elasticsearch, as_list=True)
            for elasticsearch_annotation in elasticsearch_annotations:
                if self != entry_type and elasticsearch_annotation.doc_type != self:
                    continue

                is_section_reference = isinstance(quantity_def.type, Reference)
                is_section_reference &= not isinstance(quantity_def.type, QuantityReference)
                if is_section_reference:
                    # Treat referenced sections as sub-sections
                    assert quantity_def.type.target_section_def is not None
                    # TODO e.g. owners, coauthors, etc. ... should be treated as multiple inner docs
                    # assert quantity_def.is_scalar

                    if prefix is None:
                        reference_prefix = quantity_def.name
                    else:
                        reference_prefix = f'{prefix}.{quantity_def.name}'
                    reference_mapping = self.create_mapping(
                        cast(Section, quantity_def.type.target_section_def),
                        prefix=reference_prefix)
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
                'sub sections can onyl have one elasticsearch annotation'
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
        if 'doc_type' not in kwargs:
            kwargs['doc_type'] = self.doc_type.name
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
                'mappings': {
                    self.doc_type.name: self.doc_type.mapping
                }
            })
            logger.info('elasticsearch index created')
        elif upsert:
            self.elastic_client.indices.put_mapping(
                index=self.index_name,
                doc_type=self.doc_type.name,
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
            If true, a sub-field called 'suggestion' is automatically created for
            this metainfo. This stores autocompletion suggestions for this
            value.

    Attributes:
        name:
            The name of the quantity (plus additional field if set).
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
            suggestion: bool = False,
            _es_field: str = None):

        # TODO remove _es_field if it is not necessary anymore to enforce a specific mapping
        # for v0 compatibility

        self._custom_mapping = mapping
        if suggestion:
            if field is None:
                field = "suggestion"
            else:
                raise ValueError("Cannot override suggestion field name.")
        self.field = field
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
            self._mapping = Completion().to_dict()
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
            elif quantity.type in [int, np.int32]:
                return dict(type='integer')
            elif quantity.type == np.int64:
                return dict(type='long')
            elif quantity.type == bool:
                return dict(type='boolean')
            elif quantity.type == Datetime:
                return dict(type='date')
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
        return self.definition.name

    @property
    def name(self) -> str:
        if self.field is not None:
            return f'{self.property_name}.{self.field}'
        else:
            return self.property_name

    def __repr__(self):
        if self.definition is None:
            return super().__repr__()

        return f'Elasticsearch({self.definition})'


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
            document type. This will be the quantity name (plus additional field if set)
            with subsection names up to the root of the metainfo data.
        search_field:
            The full qualified name of the field in the elasticsearch index.
        qualified_name:
            Same name as qualified_field. This will be used to address the search
            property in our APIs.
        definition: The metainfo quantity definition that this search quantity is based on
        aggregateable:
            A boolean that determines, if this quantity can be used in aggregations.
    '''
    def __init__(self, annotation: Elasticsearch, doc_type: DocumentType, prefix: str):
        self.annotation = annotation
        self.doc_type = DocumentType

        qualified_field = self.annotation.definition.name

        if prefix is not None:
            qualified_field = f'{prefix}.{qualified_field}'

        if annotation.field is not None:
            qualified_field = f'{qualified_field}.{annotation.field}'

        self.qualified_field = qualified_field
        self.qualified_name = qualified_field

        self.search_field = qualified_field
        if not(annotation._es_field == '' or annotation._es_field is None):
            self.search_field = f'{qualified_field}.{annotation._es_field}'

    @property
    def definition(self):
        return self.annotation.definition

    @property
    def aggregateable(self):
        if isinstance(self.definition.type, Reference):
            return False

        return self.annotation.mapping['type'] == 'keyword'

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


def index_entry(entry: MSection, update_material: bool = False):
    '''
    Upserts the given entry in the entry index. Optionally updates the materials index
    as well.
    '''
    index_entries([entry], update_materials=update_material)


_max_entries_index_size = 1000


def index_entries(entries: List, update_materials: bool = True, refresh: bool = False):
    '''
    Upserts the given entries in the entry index. Optionally updates the materials index
    as well.
    '''
    # split into reasonably sized problems
    if len(entries) > _max_entries_index_size:
        for entries_part in [entries[i:i + _max_entries_index_size] for i in range(0, len(entries), _max_entries_index_size)]:
            index_entries(entries_part, update_materials=update_materials)
        return

    if len(entries) == 0:
        return

    # Index the entries themselves.
    actions_and_docs = []
    for entry in entries:
        actions_and_docs.append(dict(index=dict(_id=entry['entry_id'])))
        entry_index_doc = entry_type.create_index_doc(entry)
        actions_and_docs.append(entry_index_doc)

    elasticsearch_results = entry_index.bulk(body=actions_and_docs, refresh=True)

    if not update_materials:
        return

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

    # Get existing materials for entries' material ids (i.e. the entry needs to be added
    # or updated).
    if material_ids:
        elasticsearch_results = material_index.mget(body={
            'docs': [dict(_id=material_id) for material_id in material_ids]
        })
        existing_material_docs = [
            doc['_source'] for doc in elasticsearch_results['docs'] if '_source' in doc]
    else:
        existing_material_docs = []

    # Get old materials that still have one of the entries, but the material id has changed
    # (i.e. the materials where entries need to be removed due entries having different
    # materials now).
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
    actions_and_docs = []
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
                material_doc.update(**material_type.create_index_doc(entry.results.material))

            new_material_id = get_material_id(entry)
            if new_material_id != material_id:
                # Remove the entry, it moved to another material. But the material cannot
                # run empty, because another entry had this material id.
                material_entries_to_remove.append(index)
            else:
                # Update the entry.
                material_entries[index] = material_entry_type.create_index_doc(entry)
                remaining_entry_ids.remove(entry_id)
        for index in reversed(material_entries_to_remove):
            del(material_entries[index])

        actions_and_docs.append(dict(index=dict(_id=material_id)))
        actions_and_docs.append(material_doc)

    for entry_id in remaining_entry_ids:
        entry = entries_dict.get(entry_id)
        material_id = get_material_id(entry)
        if material_id is not None:
            material_doc = material_docs_dict.get(material_id)
            if material_doc is None:
                # The material does not yet exist. Create it.
                material_doc = material_type.create_index_doc(entry.results.material)
                material_docs_dict[material_id] = material_doc
                actions_and_docs.append(dict(create=dict(_id=material_id)))
                actions_and_docs.append(material_doc)
            # The material does exist (now), but the entry is new.
            material_doc.setdefault('entries', []).append(material_entry_type.create_index_doc(entry))

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
            actions_and_docs.append(dict(delete=dict(_id=material_id)))
        else:
            # The material needs to be updated
            actions_and_docs.append(dict(index=dict(_id=material_id)))
            actions_and_docs.append(material_doc)

    # Execute the created actions in bulk.
    if len(actions_and_docs) > 0:
        material_index.bulk(body=actions_and_docs, refresh=True)

    if refresh:
        entry_index.refresh()
        if update_materials:
            material_index.refresh()
