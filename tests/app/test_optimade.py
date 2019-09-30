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

from typing import List
import json
import numpy as np
import pytest

from nomad.processing import Upload
from nomad import search
from nomad.parsing import LocalBackend
from nomad.datamodel import CalcWithMetadata

from nomad.app.optimade import parse_filter, url

from tests.app.test_app import BlueprintClient
from tests.test_normalizing import run_normalize
from tests.conftest import clear_elastic


@pytest.fixture(scope='function')
def api(client):
    return BlueprintClient(client, '/optimade')


def test_get_entry(published: Upload):
    calc_id = list(published.calcs)[0].calc_id

    with published.upload_files.archive_file(calc_id) as f:
        data = json.load(f)

    assert 'OptimadeStructureEntry' in data
    search_result = search.SearchRequest().search_parameter('calc_id', calc_id).execute_paginated()['results'][0]
    assert 'optimade' in search_result


def create_test_structure(meta_info, id: int, h: int, o: int, extra: List[str], periodicity: int):
    atom_labels = ['H' for i in range(0, h)] + ['O' for i in range(0, o)] + extra
    test_vector = np.array([0, 0, 0])

    backend = LocalBackend(meta_info, False, True)  # type: ignore
    backend.openSection('section_run')
    backend.addValue('program_name', 'test_code')
    backend.openSection('section_system')

    backend.addArrayValues('atom_labels', np.array(atom_labels))
    backend.addArrayValues(
        'atom_positions', np.array([test_vector for i in range(0, len(atom_labels))]))
    backend.addArrayValues(
        'lattice_vectors', np.array([test_vector, test_vector, test_vector]))
    backend.addArrayValues(
        'configuration_periodic_dimensions',
        np.array([True for _ in range(0, periodicity)] + [False for _ in range(periodicity, 3)]))

    backend.closeSection('section_system', 0)
    backend.closeSection('section_run', 0)

    backend = run_normalize(backend)
    calc = CalcWithMetadata(
        upload_id='test_uload_id', calc_id='test_calc_id_%d' % id, mainfile='test_mainfile',
        published=True, with_embargo=False)
    calc.apply_domain_metadata(backend)
    search.Entry.from_calc_with_metadata(calc).save()


@pytest.fixture(scope='module')
def example_structures(meta_info, elastic_infra):
    clear_elastic(elastic_infra)
    create_test_structure(meta_info, 1, 2, 1, [], 0)
    create_test_structure(meta_info, 2, 2, 1, ['C'], 0)
    create_test_structure(meta_info, 3, 2, 1, [], 1)
    create_test_structure(meta_info, 4, 1, 1, [], 0)
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
    ('LENGTH elements = 2', 3),
    ('LENGTH elements = 3', 1),
    ('LENGTH dimension_types = 0', 3),
    ('LENGTH dimension_types = 1', 1),
    ('nelements = 2 AND LENGTH dimension_types = 1', 1),
    ('nelements = 3 AND LENGTH dimension_types = 1', 0),
    ('nelements = 3 OR LENGTH dimension_types = 1', 2),
    ('nelements > 1 OR LENGTH dimension_types = 1 AND nelements = 2', 4),
    ('(nelements > 1 OR LENGTH dimension_types = 1) AND nelements = 2', 3),
    ('NOT LENGTH dimension_types = 1', 3)
])
def test_optimade_parser(example_structures, query, results):
    query = parse_filter(query)
    result = search.SearchRequest(query=query).execute_paginated()
    assert result['pagination']['total'] == results


def test_url():
    assert url('endpoint', param='value').endswith('/optimade/endpoint?param=value')


def test_list_endpoint(api, example_structures):
    rv = api.get('/calculations')
    assert rv.status_code == 200
    data = json.loads(rv.data)
    # TODO replace with real assertions
    # print(json.dumps(data, indent=2))


def test_list_endpoint_request_fields(api, example_structures):
    rv = api.get('/calculations?request_fields=nelements,elements')
    assert rv.status_code == 200
    data = json.loads(rv.data)
    # TODO replace with real assertions
    # print(json.dumps(data, indent=2))


def test_single_endpoint(api, example_structures):
    rv = api.get('/calculations/%s' % 'test_calc_id_1')
    assert rv.status_code == 200
    data = json.loads(rv.data)
    # TODO replace with real assertions
    # print(json.dumps(data, indent=2))


def test_base_info_endpoint(api):
    rv = api.get('/info')
    assert rv.status_code == 200
    data = json.loads(rv.data)
    # TODO replace with real assertions
    # print(json.dumps(data, indent=2))


def test_calculation_info_endpoint(api):
    rv = api.get('/info/calculation')
    assert rv.status_code == 200
    data = json.loads(rv.data)
    # TODO replace with real assertions
    # print(json.dumps(data, indent=2))


# TODO test single with request_fields
# TODO test errors
# TODO test response format (deny everything but json)
