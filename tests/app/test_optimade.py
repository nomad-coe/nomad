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

import json
import pytest

from nomad.processing import Upload
from nomad import search

from nomad.app.optimade import parse_filter, url

from tests.app.test_app import BlueprintClient
from tests.conftest import clear_elastic, clear_raw_files


@pytest.fixture(scope='session')
def api(session_client):
    return BlueprintClient(session_client, '/optimade/v1')


@pytest.fixture(scope='session')
def index_api(session_client):
    return BlueprintClient(session_client, '/optimade/index/v1')


def test_index(index_api):
    rv = index_api.get('/info')
    assert rv.status_code == 200

    rv = index_api.get('/links')
    assert rv.status_code == 200
    data = json.loads(rv.data)
    assert data['data'][0]['attributes']['base_url']['href'].endswith('optimade')


def test_get_entry(published: Upload):
    calc_id = list(published.calcs)[0].calc_id
    with published.upload_files.read_archive(calc_id) as archive:
        data = archive[calc_id]
        assert data['section_metadata']['dft']['optimade'] is not None

    search_result = search.SearchRequest().search_parameter('calc_id', calc_id).execute_paginated()['results'][0]
    assert 'dft.optimade.chemical_formula_hill' in search.flat(search_result)


def test_no_optimade(mongo, elastic, raw_files, api):
    from tests.app.utils import Upload
    upload = Upload()
    upload.create_test_structure(1, 2, 1, [], 0)
    upload.create_test_structure(2, 2, 1, [], 0, optimade=False)
    upload.create_upload_files()
    search.refresh()

    rv = api.get('/calculations')
    assert rv.status_code == 200
    data = json.loads(rv.data)

    assert data['meta']['data_returned'] == 1


@pytest.fixture(scope='module')
def example_structures(elastic_infra, mongo_infra, raw_files_infra):
    clear_elastic(elastic_infra)
    mongo_infra.drop_database('test_db')

    from tests.app.utils import Upload
    upload = Upload()
    upload.create_test_structure(1, 2, 1, [], 0)
    upload.create_test_structure(2, 2, 1, ['C'], 0)
    upload.create_test_structure(3, 2, 1, [], 1)
    upload.create_test_structure(4, 1, 1, [], 0)
    upload.create_upload_files()
    search.refresh()

    yield
    clear_elastic(elastic_infra)
    clear_raw_files()


@pytest.mark.parametrize('query, results', [
    ('nelements > 1', 4),
    ('nelements >= 2', 4),
    ('nelements > 2', 1),
    ('nelements < 4', 4),
    ('nelements < 3', 3),
    ('nelements <= 3', 4),
    ('nelements != 2', 1),
    ('1 < nelements', 4),
    ('elements HAS "H"', 4),
    ('elements HAS ALL "H", "O"', 4),
    ('elements HAS ALL "H", "C"', 1),
    ('elements HAS ANY "H", "C"', 4),
    ('elements HAS ANY "C"', 1),
    ('elements HAS ONLY "C"', 0),
    ('elements HAS ONLY "H", "O"', 3),
    ('nelements >= 2 AND elements HAS ONLY "H", "O"', 3),
    ('nelements >= 2 AND elements HAS ALL "H", "O", "C"', 1),
    ('nelements >= 2 AND NOT elements HAS ALL "H", "O", "C"', 3),
    ('NOT nelements = 2 AND elements HAS ANY "H", "O", "C"', 1),
    ('NOT nelements = 3 AND NOT elements HAS ONLY "H", "O"', 0),
    ('elements:elements_ratios HAS "H":>0.66', 2),
    ('elements:elements_ratios HAS ALL "O":>0.33', 3),
    ('elements:elements_ratios HAS ALL "O":>0.33,"O":<0.34', 2),
    ('elements IS KNOWN', 4),
    ('elements IS UNKNOWN', 0),
    ('chemical_formula_reduced = "H2O"', 2),
    ('chemical_formula_reduced CONTAINS "H2"', 3),
    ('chemical_formula_reduced CONTAINS "H"', 4),
    ('chemical_formula_reduced CONTAINS "C"', 1),
    ('chemical_formula_reduced STARTS "H2"', 3),
    ('chemical_formula_reduced STARTS WITH "H2"', 3),
    ('chemical_formula_reduced ENDS WITH "C"', 1),
    ('chemical_formula_reduced ENDS "C"', 1),
    ('chemical_formula_hill CONTAINS "1"', 0),
    ('chemical_formula_hill STARTS WITH "H" AND chemical_formula_hill ENDS WITH "O"', 3),
    ('NOT chemical_formula_descriptive ENDS WITH "1"', 4),
    ('chemical_formula_descriptive CONTAINS "C" AND NOT chemical_formula_descriptive STARTS WITH "O"', 1),
    ('NOT chemical_formula_anonymous STARTS WITH "A"', 0),
    ('chemical_formula_anonymous CONTAINS "AB2" AND chemical_formula_anonymous ENDS WITH "C"', 1),
    ('nsites >=3 AND elements LENGTH = 2', 2),
    ('elements LENGTH = 2', 3),
    ('elements LENGTH 2', 3),
    ('elements LENGTH = 3', 1),
    ('dimension_types LENGTH = 0', 3),
    ('dimension_types LENGTH = 1', 1),
    ('nelements = 2 AND dimension_types LENGTH = 1', 1),
    ('nelements = 3 AND dimension_types LENGTH = 1', 0),
    ('nelements = 3 OR dimension_types LENGTH = 1', 2),
    ('nelements > 1 OR dimension_types LENGTH = 1 AND nelements = 2', 4),
    ('(nelements > 1 OR dimension_types LENGTH = 1) AND nelements = 2', 3),
    ('NOT dimension_types LENGTH 1', 3),
    ('nelements LENGTH = 1', -1),
    ('LENGTH nelements = 1', -1),
    ('chemical_formula_anonymous starts with "A"', -1),
    ('elements HAS ONY "H", "O"', -1)
])
def test_optimade_parser(example_structures, query, results):
    if results >= 0:
        query = parse_filter(query)
        result = search.SearchRequest(query=query).execute_paginated()
        assert result['pagination']['total'] == results
    else:
        with pytest.raises(Exception):
            query = parse_filter(query)


def test_url():
    assert url('endpoint', param='value').endswith('/optimade/v1/endpoint?param=value')


def test_list_endpoint(api, example_structures):
    rv = api.get('/structures')
    assert rv.status_code == 200
    data = json.loads(rv.data)
    for entry in ['data', 'links', 'meta']:
        assert entry in data
    assert len(data['data']) == 4


def assert_eq_attrib(data, key, ref, item=None):
    if item is None:
        assert data['data']['attributes'][key] == ref
    else:
        assert data['data'][item]['attributes'][key] == ref


@pytest.mark.parametrize('limit, number, results', [
    (1, 1, 1), (1, 5, 0), (5, 1, 4)
])
def test_list_endpoint_pagination(api, example_structures, limit, number, results):
    rv = api.get('/structures?page_limit=%d&page_number=%d' % (limit, number))
    assert rv.status_code == 200
    data = json.loads(rv.data)
    assert len(data['data']) == results


@pytest.mark.parametrize('sort, order', [
    ('nelements', 1), ('-nelements', -1)
])
def test_list_endpoint_sort(api, example_structures, sort, order):
    rv = api.get('/structures?sort=%s' % sort)
    assert rv.status_code == 200
    data = json.loads(rv.data)['data']

    assert len(data) > 0
    for i, item in enumerate(data):
        if i > 0:
            if order == 1:
                assert item['attributes']['nelements'] >= data[i - 1]['attributes']['nelements']
            else:
                assert item['attributes']['nelements'] <= data[i - 1]['attributes']['nelements']


def test_list_endpoint_response_fields(api, example_structures):
    rv = api.get('/structures?response_fields=nelements,elements')
    assert rv.status_code == 200
    data = json.loads(rv.data)
    ref_elements = [['H', 'O'], ['C', 'H', 'O'], ['H', 'O'], ['H', 'O']]
    data['data'] = sorted(data['data'], key=lambda x: x['id'])
    for i in range(len(data['data'])):
        rf = sorted(list(data['data'][i]['attributes'].keys()))
        assert rf == ['elements', 'nelements']
        assert_eq_attrib(data, 'elements', ref_elements[i], i)
        assert_eq_attrib(data, 'nelements', len(ref_elements[i]), i)


def test_single_endpoint_response_fields(api, example_structures):
    rv = api.get('/structures/%s?response_fields=nelements,elements' % 'test_calc_id_1')
    assert rv.status_code == 200
    data = json.loads(rv.data)
    ref_elements = ['H', 'O']
    rf = sorted(list(data['data']['attributes'].keys()))
    assert rf == ['elements', 'nelements']
    assert_eq_attrib(data, 'elements', ref_elements)
    assert_eq_attrib(data, 'nelements', len(ref_elements))


def test_single_endpoint(api, example_structures):
    rv = api.get('/structures/%s' % 'test_calc_id_1')
    assert rv.status_code == 200
    data = json.loads(rv.data)
    for key in ['type', 'id', 'attributes']:
        assert key in data['data']
    fields = ['elements', 'nelements', 'elements_ratios',
              'chemical_formula_descriptive', 'chemical_formula_reduced',
              'chemical_formula_hill', 'chemical_formula_anonymous',
              'dimension_types', 'lattice_vectors', 'cartesian_site_positions',
              'nsites', 'species_at_sites', 'species']
    for field in fields:
        assert field in data['data']['attributes']


def test_base_info_endpoint(api):
    rv = api.get('/info')
    assert rv.status_code == 200
    data = json.loads(rv.data)
    for key in ['type', 'id', 'attributes']:
        assert key in data['data']
    assert data['data']['type'] == 'info'
    assert data['data']['id'] == '/'


@pytest.mark.parametrize('entry_type', ['calculations', 'structures'])
def test_entry_info_endpoint(api, entry_type):
    rv = api.get('/info/%s' % entry_type)
    assert rv.status_code == 200
    data = json.loads(rv.data)
    for key in ['description', 'properties', 'formats', 'output_fields_by_format']:
        assert key in data['data']

    assert '_nmd_atoms' in data['data']['properties']
    assert '_nmd_dft_system' in data['data']['properties']


def test_links_endpoint(api, example_structures):
    rv = api.get('/links')
    assert rv.status_code == 200
    data = json.loads(rv.data)
    assert data['data'][0]['attributes']['base_url']['href'].endswith('optimade/index')


def test_structures_endpoint(api, example_structures):
    rv = api.get('/structures')
    assert rv.status_code == 200
    data = json.loads(rv.data)
    assert len(data['data']) == 4
    for d in data['data']:
        for key in ['id', 'attributes']:
            assert d.get(key) is not None
        required_keys = [
            'last_modified', 'elements', 'nelements', 'elements_ratios', 'chemical_formula_descriptive',
            'chemical_formula_reduced', 'chemical_formula_anonymous', 'dimension_types',
            'cartesian_site_positions', 'nsites', 'species_at_sites', 'species', 'structure_features']
        for key in required_keys:
            assert key in d['attributes']


def test_structure_endpoint(api, example_structures):
    rv = api.get('/structures/%s' % 'test_calc_id_1')
    assert rv.status_code == 200
    data = json.loads(rv.data)
    assert data.get('data') is not None
    attr = data['data'].get('attributes')
    assert attr is not None
    assert attr.get('elements') == ['H', 'O']
    assert len(attr.get('dimension_types')) == 3


def test_calculations_endpoint(api, example_structures):
    rv = api.get('/calculations')
    assert rv.status_code == 200
    data = json.loads(rv.data)
    assert len(data['data']) == 4
    for d in data['data']:
        for key in ['id', 'attributes']:
            assert d.get(key) is not None
        required_keys = ['last_modified']
        for key in required_keys:
            assert key in d['attributes']


def test_calculation_endpoint(api, example_structures):
    rv = api.get('/calculations/%s' % 'test_calc_id_1')
    assert rv.status_code == 200
    data = json.loads(rv.data)
    assert data.get('data') is not None
    attr = data['data'].get('attributes')
    assert attr is not None
    assert len(attr) == 2


def test_nmd_properties(api, example_structures):
    rv = api.get('/structures/%s' % 'test_calc_id_1?response_fields=_nmd_atoms,_nmd_dft_system,_nmd_doesnotexist')
    assert rv.status_code == 200
    data = json.loads(rv.data)
    assert data.get('data') is not None
    attr = data['data'].get('attributes')
    assert attr is not None
    assert attr.get('_nmd_atoms') == ['H', 'O']
    assert '_nmd_dft_system' in attr
    assert '_nmd_doesnotexist' not in attr
