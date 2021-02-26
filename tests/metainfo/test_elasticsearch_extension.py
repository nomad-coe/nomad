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

from typing import List
import pytest

from nomad import config
from nomad.metainfo import MSection, Quantity, SubSection, Datetime
from nomad.metainfo.elasticsearch_extension import (
    Elasticsearch, create_indices, index_entry, index_entries,
    entry_type, material_type, material_entry_type, entry_index, material_index)

from ..conftest import clear_elastic_infra


class Material(MSection):

    material_id = Quantity(
        type=str,
        a_elasticsearch=Elasticsearch(material_type))

    formula = Quantity(
        type=str,
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(material_type, field='text', mapping='text')])

    springer_labels = Quantity(
        type=str, shape=['*'],
        a_elasticsearch=(Elasticsearch(material_type)))


class Properties(MSection):

    available_properties = Quantity(
        type=str, shape=['*'],
        a_elasticsearch=Elasticsearch(material_entry_type))

    band_gap = Quantity(
        type=float, unit='J',
        a_elasticsearch=Elasticsearch(material_entry_type))


class Results(MSection):

    material = SubSection(sub_section=Material.m_def)
    properties = SubSection(sub_section=Properties.m_def)


class Entry(MSection):

    entry_id = Quantity(
        type=str,
        a_elasticsearch=Elasticsearch(material_entry_type))

    upload_time = Quantity(
        type=Datetime,
        a_elasticsearch=Elasticsearch())

    results = SubSection(sub_section=Results.m_def)


def assert_mapping(mapping: dict, path: str, es_type: str, field: str = None):
    for segment in path.split('.'):
        assert 'properties' in mapping

        mapping = mapping['properties'].get(segment)

    if es_type is None:
        assert mapping is None
    else:
        assert mapping is not None
        if field is None:
            assert mapping.get('type') == es_type
        else:
            assert mapping['fields'][field].get('type') == es_type


def assert_entry_indexed(entry: Entry):
    entry_doc = entry_index.get(id=entry.entry_id)['_source']
    assert entry_doc['entry_id'] == entry.entry_id
    material_doc = material_index.get(id=entry.results.material.material_id)['_source']
    assert material_doc['material_id'] == entry.results.material.material_id
    assert any(entry_doc['entry_id'] == entry.entry_id for entry_doc in material_doc['entries'])


def assert_entries_indexed(entries: List[Entry]):
    '''
    Assert that the given entries and only the given entries and their materials are
    indexed.
    '''
    entry_docs = [
        hit['_source'] for hit in
        entry_index.search(body=dict(query=dict(match_all={})))['hits']['hits']]

    entry_ids = sorted([entry.entry_id for entry in entries])
    entry_doc_ids = sorted([entry_doc['entry_id'] for entry_doc in entry_docs])
    assert entry_doc_ids == entry_ids

    material_docs = [
        hit['_source'] for hit in
        material_index.search(body=dict(query=dict(match_all={})))['hits']['hits']]

    material_docs_based_entry_specs = sorted([
        entry['entry_id'] + '-' + material['material_id']
        for material in material_docs
        for entry in material['entries']])

    entry_specs = sorted([
        entry.entry_id + '-' + entry.results.material.material_id for entry in entries])

    assert material_docs_based_entry_specs == entry_specs

    for material_doc in material_docs:
        material = next(
            entry.results.material
            for entry in entries
            if entry.results.material.material_id == material_doc['material_id'])

        for quantity in Material.m_def.quantities:
            if material.m_is_set(quantity):
                assert material_doc[quantity.name] == getattr(material, quantity.name)
            else:
                quantity.name not in material_doc


@pytest.fixture
def example_entry():
    entry = Entry()
    entry.entry_id = 'test_entry_id'
    entry.m_create(Results).m_create(Material, material_id='test_material_id')

    return entry


@pytest.fixture
def elastic_client(elastic_function):
    yield elastic_function


@pytest.fixture(scope='module', autouse=True)
def local_elastic_infra(elastic_infra):
    # We apply a new different mapping to the test indices, we have to clear this
    # once the tests of this module have completed.
    yield
    clear_elastic_infra()


@pytest.fixture
def indices(elastic_client):
    # remove whatever the "old" infrastructure created
    elastic_client.indices.delete(index=config.elastic.index_name)
    for index in elastic_client.indices.get_alias(config.elastic.materials_index_name):
        elastic_client.indices.delete(index=index)

    return create_indices(Entry.m_def, Material.m_def)


def test_mappings(indices):
    entry_mapping, material_mapping = entry_type.mapping, material_type.mapping

    assert_mapping(entry_mapping, 'entry_id', 'keyword')
    assert_mapping(entry_mapping, 'upload_time', 'date')
    assert_mapping(entry_mapping, 'results.material.material_id', 'keyword')
    assert_mapping(entry_mapping, 'results.material.formula', 'keyword')
    assert_mapping(entry_mapping, 'results.material.formula', 'text', 'text')
    assert_mapping(entry_mapping, 'results.properties.available_properties', 'keyword')
    assert_mapping(entry_mapping, 'results.properties.band_gap', 'double')

    assert_mapping(material_mapping, 'material_id', 'keyword')
    assert_mapping(material_mapping, 'formula', 'keyword')
    assert_mapping(material_mapping, 'formula', 'text', 'text')
    assert_mapping(material_mapping, 'entries', 'nested')
    assert_mapping(material_mapping, 'entries.entry_id', 'keyword')
    assert_mapping(material_mapping, 'entries.upload_time', None)
    assert_mapping(material_mapping, 'entries.results.properties.available_properties', 'keyword')


def test_index_entry(elastic_client, indices, example_entry):
    index_entry(example_entry, update_material=True)
    assert_entry_indexed(example_entry)


# The parameters are before, to_index, after. Before and after describes what is in the index
# before and after the test. To_index describes what is to be indexed in the test.
# Both are strings with csv, where each value is a dash separated pair of numbers. The
# first number describes an entry id, the second the material id of that entry.
@pytest.mark.parametrize('before,to_index,after', [
    pytest.param('', '1-1', '1-1', id='add-first-material'),
    pytest.param('1-1', '1-1', '1-1', id='updated-material'),
    pytest.param('1-1', '2-2', '1-1, 2-2', id='added-new-material'),
    pytest.param('1-1', '2-1', '1-1, 2-1', id='added-entry-to-material'),
    pytest.param('1-1', '1-2', '1-2', id='moved-entry-between-materials-empty'),
    pytest.param('1-1, 2-1', '1-2', '2-1, 1-2', id='moved-entry-between-materials-remaining'),
    pytest.param('1-1', '1-1*', '1-1*', id='update-material-property')
])
def test_index_entries(elastic_client, indices, before, to_index, after):
    def create_entry(spec: str, material_kwargs: dict = None):
        entry_id, material_id = spec.split('-')
        changed = material_id.endswith('*')
        if changed:
            material_id = material_id[:-1]
        entry = Entry(entry_id=entry_id)
        entry.m_create(Results).m_create(
            Material, material_id=material_id,
            springer_labels=['A', 'B'] if changed else ['A'])
        return entry

    def create_entries(spec: str):
        return [
            create_entry(entry_spec.strip())
            for entry_spec in spec.split(',')
            if entry_spec.strip() != '']

    index_entries(create_entries(before), update_materials=True)
    index_entries(create_entries(to_index), update_materials=True)

    assert_entries_indexed(create_entries(after))
