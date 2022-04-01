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
import numpy as np
from elasticsearch_dsl import Keyword

from nomad import config
from nomad.utils.exampledata import ExampleData
from nomad.metainfo import MSection, Quantity, SubSection, Datetime, Unit, MEnum
from nomad.metainfo.elasticsearch_extension import (
    Elasticsearch, create_indices, index_entries_with_materials,
    entry_type, material_type, material_entry_type, entry_index, material_index)

from tests.conftest import clear_elastic_infra
from tests.app.v1.routers.common import perform_quantity_search_test


@pytest.fixture(scope="module")
def example_data_normalizers(elastic_module, raw_files_module, mongo_module, test_user, other_test_user, normalized):
    data = ExampleData(main_author=test_user)
    upload_id = "normalizer_upload"

    data.create_upload(
        upload_id=upload_id,
        published=True
    )
    data.create_entry(
        upload_id=upload_id,
        entry_id="normalizer_entry",
        material_id="normalizer_material",
        mainfile="test_content/test_entry/mainfile.json",
        results={"material": {"chemical_formula_hill": "H2O"}}
    )

    data.save()

    yield

    data.delete()
    from nomad.search import search
    assert search(query=dict(upload_id=upload_id)).pagination.total == 0


@pytest.mark.parametrize('quantity, search_str, response_str', [
    pytest.param("results.material.chemical_formula_hill", "OH2", "H2O", id='formula-reorder')
])
def test_normalizer(quantity, search_str, response_str, api_v1, example_data_normalizers):
    """Test that the normalizer specified for different annotations works
    properly in the API queries.
    """
    for resource in ["entries", "materials"]:
        perform_quantity_search_test(quantity, resource, search_str, response_str, api_v1)


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
        a_elasticsearch=(Elasticsearch(material_type, mapping=Keyword())))


class Data(MSection):
    n_points = Quantity(
        type=int,
        derived=lambda data: len(data.points[0]) if data.points is not None else 0,
        a_elasticseach=Elasticsearch(material_entry_type))

    points = Quantity(type=np.dtype(np.float64), shape=['*', '*'])

    n_series = Quantity(
        type=int,
        derived=lambda data: len(data.series) if data.series is not None else 0)

    series = Quantity(type=np.dtype(np.float64), shape=['*'])


class Dos(MSection):
    channel = Quantity(type=int, a_elasticsearch=Elasticsearch(material_entry_type))


class Properties(MSection):

    available_properties = Quantity(
        type=str, shape=['*'],
        a_elasticsearch=Elasticsearch(material_entry_type))

    band_gap = Quantity(
        type=float, unit='J',
        a_elasticsearch=Elasticsearch(material_entry_type))

    data = Quantity(type=Data, a_elasticsearch=Elasticsearch(material_entry_type))

    n_series = Quantity(type=Data.n_series, a_elasticsearch=Elasticsearch())

    dos = SubSection(sub_section=Dos.m_def, repeats=True, a_elasticsearch=Elasticsearch(nested=True))


class Results(MSection):

    material = SubSection(sub_section=Material.m_def)
    properties = SubSection(sub_section=Properties.m_def)


class User(MSection):
    user_id = Quantity(type=str, a_elasticsearch=Elasticsearch())
    name = Quantity(type=str, a_elasticsearch=Elasticsearch())


class Entry(MSection):

    entry_id = Quantity(
        type=str,
        a_elasticsearch=Elasticsearch(material_entry_type))

    upload_id = Quantity(
        type=str,
        a_elasticsearch=Elasticsearch(material_entry_type, metrics=dict(uploads='cardinality')))

    mainfile = Quantity(
        type=str,
        a_elasticsearch=Elasticsearch(index=False, value=lambda _: 'other_mainfile'))  # type: ignore

    files = Quantity(
        type=str, shape=['*'],
        a_elasticsearch=[
            Elasticsearch(_es_field='keyword'),
            Elasticsearch(mapping='text', field='path', _es_field='')])

    upload_create_time = Quantity(
        type=Datetime,
        a_elasticsearch=Elasticsearch())

    entry_create_time = Quantity(
        type=Datetime,
        a_elasticsearch=Elasticsearch())

    publish_time = Quantity(
        type=Datetime,
        a_elasticsearch=Elasticsearch())

    results = SubSection(
        sub_section=Results.m_def,
        a_eleasticsearch=Elasticsearch(auto_include_subsections=True))
    data = SubSection(sub_section=Data.m_def)

    viewers = Quantity(
        type=User, shape=['*'], a_elasticsearch=Elasticsearch())

    not_indexed = Quantity(type=str)


def assert_mapping(mapping: dict, path: str, es_type: str, field: str = None, **kwargs):
    for segment in path.split('.'):
        assert 'properties' in mapping
        mapping = mapping['properties'].get(segment)

    if es_type is None:
        assert mapping is None
    else:
        assert mapping is not None
        if field is not None:
            mapping = mapping['fields'][field]

        assert mapping.get('type') == es_type
        for key, value in kwargs.items():
            assert mapping.get(key) == value


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


@pytest.fixture(scope='module')
def indices(elastic_infra):
    # remove whatever the infrastructure created by default
    from nomad.infrastructure import elastic_client
    try:
        elastic_client.indices.delete(index=config.elastic.entries_index)
        elastic_client.indices.delete(index=config.elastic.materials_index)
    except Exception:
        pass

    create_indices(Entry.m_def, Material.m_def)
    yield
    # re-establish the default elasticsearch setup.
    clear_elastic_infra()


def test_mappings(indices):
    entry_mapping, material_mapping = entry_type.mapping, material_type.mapping

    assert_mapping(entry_mapping, 'entry_id', 'keyword')
    assert_mapping(entry_mapping, 'mainfile', 'keyword', index=False)
    assert_mapping(entry_mapping, 'upload_create_time', 'date')
    assert_mapping(entry_mapping, 'results.material.material_id', 'keyword')
    assert_mapping(entry_mapping, 'results.material.formula', 'keyword')
    assert_mapping(entry_mapping, 'results.material.formula', 'text', 'text')
    assert_mapping(entry_mapping, 'results.properties.available_properties', 'keyword')
    assert_mapping(entry_mapping, 'results.properties.band_gap', 'double')
    assert_mapping(entry_mapping, 'results.properties.data.n_points', 'integer')
    assert_mapping(entry_mapping, 'results.properties.n_series', 'integer')
    assert_mapping(entry_mapping, 'viewers.user_id', 'keyword')
    assert_mapping(entry_mapping, 'viewers.name', 'keyword')
    assert_mapping(entry_mapping, 'files', 'text')
    assert_mapping(entry_mapping, 'files', 'keyword', 'keyword')
    assert_mapping(entry_mapping, 'results.properties.dos', 'nested')
    assert_mapping(entry_mapping, 'results.properties.dos.channel', 'integer')

    assert_mapping(material_mapping, 'material_id', 'keyword')
    assert_mapping(material_mapping, 'formula', 'keyword')
    assert_mapping(material_mapping, 'formula', 'text', 'text')
    assert_mapping(material_mapping, 'entries', 'nested')
    assert_mapping(material_mapping, 'entries.entry_id', 'keyword')
    assert_mapping(material_mapping, 'entries.upload_create_time', None)
    assert_mapping(material_mapping, 'entries.results.properties.available_properties', 'keyword')
    assert_mapping(material_mapping, 'entries.results.properties.data.n_points', 'integer')
    assert_mapping(material_mapping, 'entries.results.properties.dos', 'nested')
    assert_mapping(material_mapping, 'entries.results.properties.dos.channel', 'integer')

    formula_annotations = Material.formula.m_get_annotations(Elasticsearch)
    assert entry_type.quantities.get('results.material.formula').annotation == formula_annotations[0]
    assert entry_type.quantities.get('results.material.formula.text').annotation == formula_annotations[1]

    files_annotations = Entry.files.m_get_annotations(Elasticsearch)
    assert files_annotations[0].name == 'files'
    assert files_annotations[1].name == 'files.path'
    assert entry_type.quantities.get('files').annotation == files_annotations[0]
    assert entry_type.quantities.get('files.path').annotation == files_annotations[1]

    assert entry_type.metrics['uploads'] == ('cardinality', entry_type.quantities['upload_id'])

    assert 'viewers.user_id' in entry_type.quantities
    assert Entry.viewers in entry_type.indexed_properties
    assert User.user_id in entry_type.indexed_properties
    assert Entry.viewers not in material_entry_type.indexed_properties
    assert User.user_id not in material_entry_type.indexed_properties


def test_index_docs(indices):
    user = User(user_id='test_user_id', name='Test User')
    entry = Entry(
        entry_id='test_entry_id', mainfile='test_mainfile', viewers=[user, user],
        not_indexed='value')
    data = entry.m_create(Data, points=[[0.1, 0.2], [1.1, 1.2]])
    results = entry.m_create(Results)
    results.m_create(
        Material,
        material_id='test_material_id',
        formula='H20', springer_labels=['water'])
    results.m_create(
        Properties,
        data=data, n_series=data, band_gap=1e-12, available_properties=['data', 'band_gap'])

    entry_doc = entry_type.create_index_doc(entry)
    material_entry_doc = material_entry_type.create_index_doc(entry)

    assert entry_doc == {
        'entry_id': 'test_entry_id',
        'mainfile': 'other_mainfile',
        'viewers': [
            {
                'user_id': 'test_user_id',
                'name': 'Test User'
            },
            {
                'user_id': 'test_user_id',
                'name': 'Test User'
            }
        ],
        'results': {
            'material': {
                'material_id': 'test_material_id',
                'formula': 'H20',
                'springer_labels': ['water']
            },
            'properties': {
                'available_properties': ['data', 'band_gap'],
                'band_gap': 1e-12,
                'data': {
                    'n_points': 2
                },
                'n_series': 0
            }
        }
    }

    assert material_entry_doc == {
        'entry_id': 'test_entry_id',
        'results': {
            'properties': {
                'available_properties': ['data', 'band_gap'],
                'band_gap': 1e-12,
                'data': {
                    'n_points': 2
                },
            }
        }
    }


def test_index_entry(elastic, indices, example_entry):
    index_entries_with_materials([example_entry], refresh=True)
    assert_entry_indexed(example_entry)


@pytest.mark.parametrize('metainfo_type,es_type', [
    [str, 'keyword'],
    [np.dtype(np.float32), 'float'],
    [float, 'double'],
    [np.dtype(np.float64), 'double'],
    [int, 'integer'],
    [np.dtype(np.int32), 'integer'],
    [np.dtype(np.int64), 'integer'],
    [bool, 'boolean'],
    [Datetime, 'date'],
    [Unit, 'keyword'],
    [MEnum('test'), 'keyword'],
])
def test_mapping_detection(metainfo_type, es_type):
    '''Tests that the mappings are correctly determined from the quantity type.
    '''
    a = Quantity(type=metainfo_type, a_elasticsearch=Elasticsearch())
    a.__init_metainfo__()
    assert a.m_annotations['elasticsearch'].mapping['type'] == es_type


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
def test_index_entries(elastic, indices, before, to_index, after):
    index_entries_with_materials(create_entries(before), refresh=True)
    index_entries_with_materials(create_entries(to_index), refresh=True)

    assert_entries_indexed(create_entries(after))


@pytest.mark.parametrize('cap, entries', [
    pytest.param(2, 1, id='below-cap'),
    pytest.param(2, 3, id='above-cap')
])
def test_index_materials_capped(elastic, indices, monkeypatch, cap, entries):
    monkeypatch.setattr('nomad.config.elastic.entries_per_material_cap', cap)
    index_entries_with_materials(create_entries(','.join([f'{i}-1' for i in range(1, entries + 1)])), refresh=True)

    material_docs = [
        hit['_source'] for hit in
        material_index.search(body=dict(query=dict(match_all={})))['hits']['hits']]
    for material_doc in material_docs:
        assert len(material_doc['entries']) <= cap
        assert material_doc['n_entries'] == entries
