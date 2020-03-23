# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import json
import pytest

from nomad.processing import Upload
from nomad import search

from nomad.app.optimade import parse_filter, url

from tests.app.test_app import BlueprintClient
from tests.conftest import clear_elastic, create_test_structure


@pytest.fixture(scope='session')
def api(session_client):
    return BlueprintClient(session_client, '/optimade')


def test_get_entry(published: Upload):
    calc_id = list(published.calcs)[0].calc_id
    with published.upload_files.read_archive(calc_id) as archive:
        data = archive[calc_id]
    assert 'OptimadeEntry' in data
    search_result = search.SearchRequest().search_parameter('calc_id', calc_id).execute_paginated()['results'][0]
    assert 'dft.optimade.chemical_formula_hill' in search.flat(search_result)


def test_no_optimade(mongo, elastic, api):
    create_test_structure(1, 2, 1, [], 0)
    create_test_structure(2, 2, 1, [], 0, optimade=False)
    search.refresh()

    rv = api.get('/calculations')
    assert rv.status_code == 200
    data = json.loads(rv.data)

    assert data['meta']['data_returned'] == 1


@pytest.fixture(scope='module')
def example_structures(elastic_infra, mongo_infra):
    clear_elastic(elastic_infra)
    mongo_infra.drop_database('test_db')

    create_test_structure(1, 2, 1, [], 0)
    create_test_structure(2, 2, 1, ['C'], 0)
    create_test_structure(3, 2, 1, [], 1)
    create_test_structure(4, 1, 1, [], 0)
    search.refresh()

    yield
    clear_elastic(elastic_infra)


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
    ('nsites >=3 AND LENGTH elements = 2', 2),
    ('LENGTH elements = 2', 3),
    ('LENGTH elements = 3', 1),
    ('LENGTH dimension_types = 0', 3),
    ('LENGTH dimension_types = 1', 1),
    ('nelements = 2 AND LENGTH dimension_types = 1', 1),
    ('nelements = 3 AND LENGTH dimension_types = 1', 0),
    ('nelements = 3 OR LENGTH dimension_types = 1', 2),
    ('nelements > 1 OR LENGTH dimension_types = 1 AND nelements = 2', 4),
    ('(nelements > 1 OR LENGTH dimension_types = 1) AND nelements = 2', 3),
    ('NOT LENGTH dimension_types = 1', 3),
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
    assert url('endpoint', param='value').endswith('/optimade/endpoint?param=value')


def test_list_endpoint(api, example_structures):
    rv = api.get('/calculations')
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


def test_list_endpoint_request_fields(api, example_structures):
    rv = api.get('/calculations?request_fields=nelements,elements')
    assert rv.status_code == 200
    data = json.loads(rv.data)
    ref_elements = [['H', 'O'], ['C', 'H', 'O'], ['H', 'O'], ['H', 'O']]
    for i in range(len(data['data'])):
        rf = list(data['data'][i]['attributes'].keys())
        rf.sort()
        assert rf == ['elements', 'nelements']
        assert_eq_attrib(data, 'elements', ref_elements[i], i)
        assert_eq_attrib(data, 'nelements', len(ref_elements[i]), i)


def test_single_endpoint_request_fields(api, example_structures):
    rv = api.get('/calculations/%s?request_fields=nelements,elements' % 'test_calc_id_1')
    assert rv.status_code == 200
    data = json.loads(rv.data)
    ref_elements = ['H', 'O']
    rf = list(data['data']['attributes'].keys())
    assert rf == ['elements', 'nelements']
    assert_eq_attrib(data, 'elements', ref_elements)
    assert_eq_attrib(data, 'nelements', len(ref_elements))


def test_single_endpoint(api, example_structures):
    rv = api.get('/calculations/%s' % 'test_calc_id_1')
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


def test_calculation_info_endpoint(api):
    rv = api.get('/info/calculations')
    assert rv.status_code == 200
    data = json.loads(rv.data)
    for key in ['description', 'properties', 'formats', 'output_fields_by_format']:
        assert key in data['data']
