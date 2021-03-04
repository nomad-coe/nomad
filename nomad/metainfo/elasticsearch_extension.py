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


from typing import Union, Any, Dict, cast, Set, List
import numpy as np

from nomad import config

from .metainfo import (
    Section, Quantity, MSection, MEnum, Datetime, Reference, DefinitionAnnotation,
    Definition, MetainfoError, QuantityReference)


class DocumentType():
    '''
    DocumentType allows to create Elasticsearch index mappings and documents based on
    Metainfo definitions and instances. Genrally this class should not be used outside
    the elasticsearch_extension module.
    '''
    def __init__(self, name: str):
        self.name = name
        self.root_section_def = None
        self.mapping: Dict[str, Any] = None
        self.indexed_properties: Set[Definition] = set()

    def create_index_doc(self, root: MSection):
        '''
        Creates an indexable document from the given archive.
        '''
        return root.m_to_dict(
            with_meta=False, include_defaults=True, include_derived=True,
            resolve_references=True,
            exclude=lambda property_, section: property_ not in self.indexed_properties)

    def create_mapping(self, section_def: Section):
        '''
        Creates an Elasticsearch mapping for the given root section. It traverses all
        sub-sections to create the mapping. It will not create the mapping for nested
        documents. These have to be created manually (e.g. by :func:`create_indices`).
        '''
        mappings: Dict[str, Any] = {}

        for quantity_def in section_def.all_quantities.values():
            elasticsearch_annotations = quantity_def.m_get_annotations(Elasticsearch, as_list=True)
            for elasticsearch_annotation in elasticsearch_annotations:
                is_section_reference = isinstance(quantity_def.type, Reference) and not isinstance(quantity_def.type, QuantityReference)
                if not is_section_reference and self != entry_type and elasticsearch_annotation.doc_type != self:
                    continue

                if is_section_reference:
                    # Treat referenced sections as sub-sections
                    assert quantity_def.type.target_section_def is not None
                    assert quantity_def.is_scalar

                    reference_mapping = self.create_mapping(cast(Section, quantity_def.type.target_section_def))
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

        for sub_section_def in section_def.all_sub_sections.values():
            if sub_section_def.m_get_annotations(Elasticsearch) is None:
                continue

            assert not sub_section_def.repeats, 'elasticsearch fields in repeating sub sections are not supported'
            sub_section_mapping = self.create_mapping(sub_section_def.sub_section)
            if len(sub_section_mapping['properties']) > 0:
                mappings[sub_section_def.name] = sub_section_mapping
                self.indexed_properties.add(sub_section_def)

        self.mapping = dict(properties=mappings)
        return self.mapping


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

    def create_index(self):
        ''' Initially creates the index with the mapping of its document type. '''
        assert self.doc_type.mapping is not None, 'The mapping has to be created first.'
        self.elastic_client.indices.create(index=self.index_name, body={
            'mappings': {
                self.doc_type.name: self.doc_type.mapping
            }
        })


entry_type = DocumentType('entry')
material_type = DocumentType('material')
material_entry_type = DocumentType('material_entry')

entry_index = Index(entry_type, index_config_key='index_name')
material_index = Index(material_type, index_config_key='materials_index_name')


class Elasticsearch(DefinitionAnnotation):
    '''
    A metainfo annotation for quantity definitions. This annotation can be used multiple
    times on the same quantity (e.g. to define Elasticsearch fields with differrent mapping
    types).

    Arguments:
        doc_type: An additional document type: ``material_type`` or ``material_entry_type``.
            All quantities with this annotation are automatically placed in ``entry_type``.
        mapping: The Elasticsearch mapping type or full dictionary for this mapping. The
            default depends on the quantity type.
        field: Allows to specify a field name. There has to be another annotation on the
            same quantity without a field.
    '''
    def __init__(
            self,
            doc_type: DocumentType = entry_type,
            mapping: Union[str, Dict[str, Any]] = None,
            field: str = None):

        if isinstance(mapping, str):
            self._mapping = dict(type=mapping)
        else:
            self._mapping = mapping

        self._field = field
        self.doc_type = doc_type

    def _compute_mapping(self, quantity: Quantity):
        if quantity.type == str:
            return dict(type='keyword')
        elif quantity.type in [float, np.float64] and quantity.is_scalar:
            return dict(type='double')
        elif quantity.type == np.float32 and quantity.is_scalar:
            return dict(type='float')
        elif quantity.type in [int, np.int32] and quantity.is_scalar:
            return dict(type='integer')
        elif quantity.type == np.int64 and quantity.is_scalar:
            return dict(type='long')
        elif quantity.type == bool:
            return dict(type='boolean')
        elif quantity.type == Datetime:
            return dict(type='date')
        elif isinstance(quantity.type, QuantityReference):
            return self._compute_mapping(quantity.type.target_quantity_def)
        elif isinstance(quantity.type, Reference):
            raise MetainfoError('References cannot be indexed.')
        elif isinstance(quantity.type, MEnum):
            return dict(type='keyword')
        else:
            raise NotImplementedError(
                'Quantity type %s for quantity %s is not supported.' % (quantity.type, quantity))

    @property
    def mapping(self):
        if self._mapping is not None:
            return self._mapping

        return self._compute_mapping(cast(Quantity, self.definition))

    @property
    def fields(self):
        if self._field is None:
            return {}

        return {
            self._field: self.mapping
        }

    @property
    def property_name(self):
        return self.definition.name


def create_indices(entry_section_def: Section = None, material_section_def: Section = None):
    '''
    Initially creates the mapping for all document types and creates the indices in Elasticsearch.
    The indices must not exist already.
    '''
    if entry_section_def is None:
        from nomad.datamodel import EntryArchive
        entry_section_def = EntryArchive.m_def

    if material_section_def is None:
        from nomad.datamodel.encyclopedia import Material
        material_section_def = Material.m_def

    entry_type.create_mapping(entry_section_def)
    material_type.create_mapping(material_section_def)
    material_entry_type.create_mapping(entry_section_def)
    material_entry_type.mapping['type'] = 'nested'
    material_type.mapping['properties']['entries'] = material_entry_type.mapping

    entry_index.create_index()
    material_index.create_index()


def index_entry(entry: MSection, update_material: bool = False):
    '''
    Upserts the given entry in the entry index. Optionally updates the materials index
    as well.
    '''
    index_entries([entry], update_materials=update_material)


_max_entries_index_size = 1000


def index_entries(entries: List, update_materials: bool = True):
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
        actions_and_docs.append(entry_type.create_index_doc(entry))
    elasticsearch_results = entry_index.bulk(body=actions_and_docs, refresh=True)

    if not update_materials:
        return

    # Get all entry and material ids.
    entry_ids, material_ids = set(), set()
    entries_dict = {}
    for entry in entries:
        entries_dict[entry.entry_id] = entry
        entry_ids.add(entry.entry_id)
        if entry.results.material is not None:
            material_ids.add(entry.results.material.material_id)

    # Get existing materials for entries' material ids (i.e. the entry needs to be added
    # or updated).
    elasticsearch_results = material_index.mget(body={
        'docs': [dict(_id=material_id) for material_id in material_ids]
    })
    existing_material_docs = [
        doc['_source'] for doc in elasticsearch_results['docs'] if '_source' in doc]

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

            if entry.results.material.material_id != material_id:
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
        material_id = entry.results.material.material_id
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
    # - an entry needs to be removed but the material still as entries (new material id case 1)
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
