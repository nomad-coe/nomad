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
import json
from fastapi.testclient import TestClient
import httpx
from urllib.parse import urlencode
import time
import dateutil.parser

from nomad import config
from nomad.app.resources.main import app, remove_mongo
from nomad.app.resources.routers.resources import (
    aflow_prototypes_db,
    springer_materials_db,
    optimade_providers,
)


def _to_datetime(datetime_str):
    return dateutil.parser.isoparse(datetime_str).timestamp()
    # return datetime.strptime(datetime_str, '%Y-%m-%dT%H:%M:%S.%f').timestamp()


def _get_resources(api, params, is_retrieving_more, repeat=True):
    response = api.get(f'/?{urlencode(params, doseq=True)}')

    assert response.status_code == 200
    response_json = response.json()
    data = response_json['data']
    if is_retrieving_more is not None:
        assert response_json['is_retrieving_more'] == is_retrieving_more

    if not response_json['is_retrieving_more']:
        return data

    if not repeat:
        return data

    while response_json['is_retrieving_more']:
        time.sleep(1)
        response = api.get(f'/?{urlencode(params, doseq=True)}')
        assert response.status_code == 200
        response_json = response.json()
        data = response_json['data']

    return data


@pytest.fixture(scope='session')
def api():
    return TestClient(app, base_url='http://testserver/')


@pytest.fixture(scope='function')
def patched_download(monkeypatch):
    with open('tests/data/api/resources_mocked_responses.json') as f:
        responses = json.load(f)

    async def _download(session: httpx.AsyncClient, path: str) -> httpx.Response:
        async def get(path):
            response_dict = responses.get(path)
            if path is None:
                return httpx.Response(status_code=404)
            json_data = response_dict.get('json')
            response = httpx.Response(
                text=response_dict.get('text'),
                status_code=response_dict.get('status_code'),
                json=json_data if json_data else {},
            )
            return response

        response = await get(path)
        return response

    monkeypatch.setattr('nomad.app.resources.routers.resources._download', _download)


@pytest.fixture(scope='function')
def resources(mongo, monkeypatch):
    monkeypatch.setattr('nomad.config.resources.enabled', True)
    monkeypatch.setattr('nomad.config.resources.db_name', 'test_db_resources')
    remove_mongo()
    yield
    remove_mongo()


def _perform_initial_get_resources(api, params, data_length):
    data = _get_resources(api, params, is_retrieving_more=None)
    assert len(data) == data_length
    return sorted(data, key=lambda x: x['id'])


@pytest.mark.timeout(config.tests.default_timeout)
def test_initial_get_resources(api, resources, patched_download, worker):
    params = dict(
        chemical_formula_reduced='AcAg',
        wyckoff_letters=['a', 'b'],
        space_group_number=225,
        n_sites=2,
    )
    data = _perform_initial_get_resources(api, params, data_length=7)

    aflow_data = [d for d in data if d['database_name'] == aflow_prototypes_db]
    assert len(aflow_data) == 1
    assert 'Space group symbol' in aflow_data[0]['available_data']
    assert aflow_data[0]['id'] == 'AB_cF8_225_a_b'

    springer_data = [d for d in data if d['database_name'] == springer_materials_db]
    assert len(springer_data) == 0

    optimade_dbs = [provider['name'] for provider in optimade_providers.values()]
    optimade_data = [d for d in data if d['database_name'] in optimade_dbs]
    assert len(optimade_data) == 6
    assert optimade_data[0]['database_version'] == '1.0.0'
    assert optimade_data[1]['url'] == 'https://oqmd.org/materials/entry/675180'
    assert optimade_data[2]['id'] == '4815213'


@pytest.mark.timeout(config.tests.default_timeout)
def test_cached_get_resources(api, resources, patched_download, worker):
    params = dict(
        chemical_formula_reduced='Mg',
        wyckoff_letters=['a'],
        space_group_number=229,
        n_sites=1,
    )
    # do initial request
    data_initial = _perform_initial_get_resources(api, params, data_length=88)

    # repeat request
    data_repeat = _get_resources(api, params, is_retrieving_more=False)
    assert len(data_repeat) == 88
    # check if download_time is the same
    # sort data according to id
    data_repeat = sorted(data_repeat, key=lambda x: x['id'])
    for i in range(len(data_initial)):
        assert data_initial[i]['id'] == data_repeat[i]['id']
        # mongodb does not save datetime precisely
        # (https://www.mongodb.com/community/forums/t/for-date-field-dont-save-milliseconds-in-mongodb/110557)
        assert _to_datetime(data_initial[i]['download_time']) == pytest.approx(
            _to_datetime(data_repeat[i]['download_time']), 0.001
        )


@pytest.mark.timeout(config.tests.default_timeout)
def test_cache_invalidation_get_resources(
    api, resources, patched_download, worker, monkeypatch
):
    params = dict(
        chemical_formula_reduced='Mg',
        wyckoff_letters=['a'],
        space_group_number=229,
        n_sites=1,
    )
    # do initial request
    data_initial = _perform_initial_get_resources(api, params, data_length=88)

    # mimic mongo update by setting max_time_in_mongo to 0
    monkeypatch.setattr('nomad.config.resources.max_time_in_mongo', 0.0)

    # repeat request, expect that resources are downloaded again
    _get_resources(api, params, is_retrieving_more=True, repeat=False)
    monkeypatch.setattr('nomad.config.resources.max_time_in_mongo', 3600)
    data_repeat = sorted(
        _get_resources(api, params, is_retrieving_more=True, repeat=True),
        key=lambda x: x['id'],
    )
    for i in range(len(data_initial)):
        assert data_initial[i]['id'] == data_repeat[i]['id']
        assert _to_datetime(data_initial[i]['download_time']) < _to_datetime(
            data_repeat[i]['download_time']
        )
