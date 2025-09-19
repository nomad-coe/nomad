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

import pytest

from nomad.app.v1.routers.apps import parse_jmespath
from nomad.config import config
from nomad.config.models.plugins import AppEntryPoint
from nomad.config.models.ui import App, Column
from nomad.metainfo.elasticsearch_extension import entry_type

BASE = config.services.api_base_path.rstrip('/') + '/apps'


def test_get_invalid_app_returns_404(client, no_warn):
    resp = client.get(f'{BASE}/entry-points/not-an-app')
    assert resp.status_code == 404
    body = resp.json()
    assert 'Could not find an app with the path' in body.get('detail', '')


def test_list_apps_returns_data_list(client, no_warn):
    from nomad.config.models.plugins import AppEntryPoint

    resp = client.get(f'{BASE}/entry-points')
    assert resp.status_code == 200
    data = resp.json().get('data', [])
    assert isinstance(data, list)

    apps_eps = [
        ep
        for ep in config.plugins.entry_points.filtered_values()
        if isinstance(ep, AppEntryPoint)
    ]
    expected_count = len(apps_eps)

    assert len(data) == expected_count

    returned_ids = {item['id'] for item in data}
    expected_ids = {ep.app.path for ep in apps_eps}
    assert returned_ids == expected_ids


def test_get_specific_app_returns_app_and_search_quantities(client, no_warn):
    resp = client.get(f'{BASE}/entry-points')
    apps = resp.json()['data']
    assert apps, 'Expected at least one app'
    app_id = apps[0]['id']

    resp2 = client.get(f'{BASE}/entry-points/{app_id}')
    assert resp2.status_code == 200
    body = resp2.json()
    assert 'app' in body and 'search_quantities' in body
    assert body['app']['path'] == app_id
    assert isinstance(body['search_quantities'], dict)


@pytest.mark.parametrize(
    'payload',
    [
        {},
        {'query': {}},
    ],
)
def test_search_quantities_missing_fields(client, no_warn, payload):
    resp = client.post(f'{BASE}/search-quantities', json=payload)
    assert resp.status_code == 422


def test_search_quantities_invalid_app_path(client, no_warn):
    payload = {
        'app_path': 'does-not-exist',
        'query': {'input': 'energy'},
        'pagination': {'page': 1, 'page_size': 5},
    }
    resp = client.post(f'{BASE}/search-quantities', json=payload)
    assert resp.status_code == 404
    body = resp.json()
    assert 'Could not find an app with the path' in body.get('detail', '')


def test_search_quantities_suggestions_and_exact_match(client, no_warn):
    all_q = list(entry_type.quantities.keys())
    assert all_q, 'No search quantities registered!'
    example = all_q[0]

    partial = example.split('.')[-1]
    resp = client.post(
        f'{BASE}/search-quantities',
        json={'query': {'input': partial}, 'pagination': {'page': 1, 'page_size': 10}},
    )
    assert resp.status_code == 200
    suggestions = resp.json()
    assert isinstance(suggestions, list) and suggestions

    resp2 = client.post(
        f'{BASE}/search-quantities',
        json={'query': {'input': example}, 'pagination': {'page': 1, 'page_size': 10}},
    )
    assert resp2.status_code == 200
    quantities = [item['quantity'] for item in resp2.json()]
    assert example in quantities


def test_search_quantities_filter_by_aggregatable(client, no_warn):
    resp = client.post(
        f'{BASE}/search-quantities',
        json={
            'query': {'aggregatable': True},
            'pagination': {'page': 1, 'page_size': 50},
        },
    )
    assert resp.status_code == 200
    out = resp.json()
    assert out, 'Expected some aggregatable quantities'
    assert all(item['aggregatable'] for item in out)


def test_search_quantities_pagination(client, no_warn):
    resp_full = client.post(
        f'{BASE}/search-quantities',
        json={'query': {}, 'pagination': {'page': 1, 'page_size': 1000}},
    )
    full = resp_full.json()
    total = len(full)
    if total < 2:
        pytest.skip('Not enough quantities for pagination')

    resp_p1 = client.post(
        f'{BASE}/search-quantities',
        json={'query': {}, 'pagination': {'page': 1, 'page_size': 1}},
    )
    resp_p2 = client.post(
        f'{BASE}/search-quantities',
        json={'query': {}, 'pagination': {'page': 2, 'page_size': 1}},
    )
    assert resp_p1.status_code == 200
    assert resp_p2.status_code == 200

    q1 = resp_p1.json()[0]['quantity']
    q2 = resp_p2.json()[0]['quantity']
    assert q1 != q2


@pytest.mark.parametrize(
    'search_quantity, error, status_code',
    [
        pytest.param(
            'results.material.chemical_formula_hill',
            None,
            200,
            id='plain',
        ),
        pytest.param(
            'missing',
            'Could not load the search quantity "missing" used in the results table column.',
            422,
            id='missing',
        ),
        pytest.param(
            'min_by(results.properties.electronic.band_gap[*], &value).type',
            None,
            200,
            id='jmespath',
        ),
        pytest.param(
            'min_by(results.properties.electronic.band_gap[*], &missing).type',
            'Could not load the search quantity "results.properties.electronic.band_gap.missing" used in the results table column.',
            422,
            id='jmespath-missing-extra',
        ),
        pytest.param(
            'data.datetime#pynxtools.nomad.schema.Root#datetime',
            None,
            200,
            id='explicit-dtype',
        ),
    ],
)
def test_app_search_quantity_validation(
    client, no_warn, monkeypatch, search_quantity, error, status_code
):
    """
    Test that search quantities used in apps are validated correctly.
    """
    # Patch the endpoint to use the given app
    monkeypatch.setattr('nomad.app.v1.routers.apps.app_cache', {})
    monkeypatch.setattr(
        'nomad.app.v1.routers.apps.app_entry_points_cache',
        {
            'test': AppEntryPoint(
                id='test',
                app=App(
                    label='test',
                    path='test',
                    category='test',
                    columns=[Column(search_quantity=search_quantity)],
                ),
            )
        },
    )

    # Try to get app, assert status code and potential error
    resp = client.get(f'{BASE}/entry-points/test')
    assert resp.status_code == status_code
    if error:
        body = resp.json()
        assert error in body.get('detail', '')


@pytest.mark.parametrize(
    'input_path, expected',
    [
        pytest.param(
            'results.material.n_elements',
            {
                'quantity': 'results.material.n_elements',
                'path': 'results.material.n_elements',
                'extras': [],
                'error': None,
                'schema': '',
            },
            id='simple-subexpression',
        ),
        pytest.param(
            'results.material.elements[0]',
            {
                'quantity': 'results.material.elements',
                'path': 'results.material.elements[0]',
                'extras': [],
                'error': None,
                'schema': '',
            },
            id='index-expression',
        ),
        pytest.param(
            'results.material.elements[0:5]',
            {
                'quantity': 'results.material.elements',
                'path': 'results.material.elements[0:5]',
                'extras': [],
                'error': None,
                'schema': '',
            },
            id='slicing',
        ),
        pytest.param(
            'results.properties.electronic.band_gap[*].value',
            {
                'quantity': 'results.properties.electronic.band_gap.value',
                'path': 'results.properties.electronic.band_gap[*].value',
                'extras': [],
                'error': None,
                'schema': '',
            },
            id='list-projection',
        ),
        pytest.param(
            'results.properties[].electronic.band_gap[].value',
            {
                'quantity': 'results.properties.electronic.band_gap.value',
                'path': 'results.properties[].electronic.band_gap[].value',
                'extras': [],
                'error': None,
                'schema': '',
            },
            id='flatten-projection',
        ),
        pytest.param(
            'map(&conversion, results.properties.catalytic.reaction.reactants[*])',
            {
                'quantity': 'results.properties.catalytic.reaction.reactants.conversion',
                'path': 'map(&conversion, results.properties.catalytic.reaction.reactants[*])',
                'extras': [],
                'error': None,
                'schema': '',
            },
            id='map',
        ),
        pytest.param(
            'min(results.properties.electronic.band_gap[*].value)',
            {
                'quantity': 'results.properties.electronic.band_gap.value',
                'path': 'min(results.properties.electronic.band_gap[*].value)',
                'extras': [],
                'error': None,
                'schema': '',
            },
            id='function-with-one-argument',
        ),
        pytest.param(
            'min_by(results.properties.electronic.band_gap[*], &value).type',
            {
                'quantity': 'results.properties.electronic.band_gap.type',
                'path': 'min_by(results.properties.electronic.band_gap[*], &value).type',
                'extras': ['results.properties.electronic.band_gap.value'],
                'error': None,
                'schema': '',
            },
            id='function-with-two-arguments',
        ),
        pytest.param(
            "results.material.topology[?label=='original'].cell.a",
            {
                'quantity': 'results.material.topology.cell.a',
                'path': "results.material.topology[?label=='original'].cell.a",
                'extras': ['results.material.topology.label'],
                'error': None,
                'schema': '',
            },
            id='filter projection',
        ),
        pytest.param(
            'results.properties.electronic.band_gap[*].value | min(@)',
            {
                'quantity': 'results.properties.electronic.band_gap.value',
                'path': 'results.properties.electronic.band_gap[*].value | min(@)',
                'extras': [],
                'error': None,
                'schema': '',
            },
            id='pipe',
        ),
        pytest.param(
            'min_by(results.properties.electronic.band_gap[*], &value).type#MySchema#int',
            {
                'quantity': 'results.properties.electronic.band_gap.type#MySchema#int',
                'path': 'min_by(results.properties.electronic.band_gap[*], &value).type',
                'extras': ['results.properties.electronic.band_gap.value#MySchema#int'],
                'error': None,
                'schema': '#MySchema#int',
            },
            id='schema-name-and-dtype-are-handled-correctly-1',
        ),
        pytest.param(
            'results.properties.electronic.band_gap[*].value#MySchema | min(@)',
            {
                'quantity': 'results.properties.electronic.band_gap.value#MySchema',
                'path': 'results.properties.electronic.band_gap[*].value | min(@)',
                'extras': [],
                'error': None,
                'schema': '#MySchema',
            },
            id='schema-name-and-dtype-are-handled-correctly-2',
        ),
        pytest.param(
            'results.material.n_elements[*',
            {
                'quantity': None,
                'path': None,
                'extras': None,
                'error': """Invalid jmespath expression: Incomplete expression:
"results.material.n_elements[*"
                              ^""",
                'schema': '',
            },
            id='syntax-error',
        ),
    ],
)
def test_parse_jmespath(input_path, expected):
    result = parse_jmespath(input_path)
    assert result == expected
