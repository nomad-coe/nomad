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
from urllib.parse import urlencode

from nomad.datamodel import Dataset

from .test_entries import data as example_entries
from .common import assert_response

'''
These are the tests for all API operations below ``datasets``. The tests are organized
using the following type of methods: fixtures, ``perfrom_*_test``, ``assert_*``, and
``test_*``. While some ``test_*`` methods test individual API operations, some
test methods will test multiple API operations that use common aspects like
supporting queries, pagination, or the owner parameter. The test methods will use
``perform_*_test`` methods as an parameter. Similarely, the ``assert_*`` methods allow
to assert for certain aspects in the responses.
'''


@pytest.fixture(scope='function')
def data(mongo, test_user, other_test_user):
    def create_dataset(**kwargs):
        dataset = Dataset(**kwargs)
        dataset.m_get_annotations('mongo').save()

    create_dataset(
        dataset_id='dataset_1',
        user_id=test_user.user_id,
        name='test dataset 1',
        dataset_type='owned')

    create_dataset(
        dataset_id='dataset_2',
        user_id=test_user.user_id,
        name='test dataset 2',
        dataset_type='owned')

    create_dataset(
        dataset_id='dataset_listed',
        user_id=test_user.user_id,
        name='foreign test dataset',
        dataset_type='foreign')

    create_dataset(
        dataset_id='dataset_doi',
        user_id=test_user.user_id,
        name='foreign test dataset',
        dataset_type='foreign',
        doi='test_doi')


def assert_dataset(dataset, **kwargs):
    for key, value in kwargs.items():
        assert dataset[key] == value

    mongo_dataset = Dataset.m_def.a_mongo.objects(dataset_id=dataset['dataset_id']).first()
    assert mongo_dataset is not None
    for key, value in mongo_dataset.items():
        assert dataset[key] == value


@pytest.mark.parametrize('query, size, status_code', [
    pytest.param({}, 3, 200, id='empty'),
    pytest.param({'dataset_id': 'dataset_1'}, 1, 200, id='id'),
    pytest.param({'dataset_type': 'foreign'}, 1, 200, id='type'),
    pytest.param({'dataset_id': 'DOESNOTEXIST'}, 0, 200, id='id-not-exists')
])
def test_datasets(client, data, query, size, status_code):
    if len(query) == 0:
        response = client.get('datasets/')
    else:
        response = client.get('datasets/?%s' % urlencode(query, doseq=True))

    json_response = assert_response(response, status_code=status_code)

    if json_response is None:
        return

    assert len(json_response['data']) == size
    for dataset in json_response['data']:
        assert_dataset(dataset, **query)


@pytest.mark.parametrize('dataset_id, result, status_code', [
    pytest.param('dataset_1', {'dataset_id': 'dataset_1'}, 200, id='plain'),
    pytest.param('DOESNOTEXIST', None, 404, id='not-exists')
])
def test_dataset(client, data, dataset_id, result, status_code):
    response = client.get('datasets/%s' % dataset_id)

    json_response = assert_response(response, status_code=status_code)

    if json_response is None:
        return

    assert_dataset(json_response['data'], **result)


@pytest.mark.parametrize('name, dataset_type, query, entries, user, status_code', [
    pytest.param('another test dataset', 'foreign', None, None, 'test_user', 200, id='plain'),
    pytest.param('another test dataset', 'foreign', None, None, None, 401, id='no-user'),
    pytest.param('test dataset 1', 'foreign', None, None, 'test_user', 400, id='exists'),
    pytest.param('another test dataset', 'owned', None, None, 'test_user', 400, id='owned'),
    pytest.param('another test dataset', 'owned', {}, None, 'test_user', 400, id='query'),
    pytest.param('another test dataset', 'owned', None, ['id_01', 'id_02'], 'test_user', 400, id='entries')
])
def test_post_datasets(client, data, example_entries, test_user_auth, name, dataset_type, query, entries, user, status_code):
    dataset = {'name': name, 'dataset_type': dataset_type}
    if query is not None:
        dataset['query'] = query
    if entries is not None:
        dataset['entries'] = entries
    auth = None
    if user == 'test_user':
        auth = test_user_auth
    response = client.post(
        'datasets/', headers=auth, json=dataset)

    json_response = assert_response(response, status_code=status_code)
    if json_response is None:
        return

    dataset = json_response['data']
    assert_dataset(dataset, user_id=test_user.user_id, name=name, dataset_type=dataset_type)
    assert Dataset.a_mongo.objects().count() == 3

    if query is not None or entries is not None:
        assert json_response['data']['entries'] is not None
        assert len(json_response['data']['entries']) > 0


@pytest.mark.parametrize('dataset_id, user, status_code', [
    pytest.param('dataset_listed', 'test_user', 200, id='plain'),
    pytest.param('dataset_listed', None, 401, id='no-user'),
    pytest.param('dataset_listed', 'other_test_user', 401, id='wrong-user'),
    pytest.param('DOESNOTEXIST', 'test_user', 404, id='does-not-exist'),
    pytest.param('dataset_1', 'test_user', 400, id='owned'),
    pytest.param('dataset_doi', 'test_user', 400, id='with-doi'),
])
def test_delete_dataset(client, data, test_user_auth, other_test_user_auth, dataset_id, user, status_code):
    auth = None
    if user == 'test_user':
        auth = test_user_auth
    if user == 'other_test_user':
        auth = other_test_user_auth
    response = client.delete(
        'datasets/%s' % dataset_id, headers=auth)

    json_response = assert_response(response, status_code=status_code)
    if json_response is None:
        return

    assert Dataset.a_mongo.objects().count() == 2


@pytest.mark.parametrize('dataset_id, user, status_code', [
    pytest.param('dataset_1', 'test_user', 200, id='plain'),
    pytest.param('dataset_1', None, 401, id='no-user'),
    pytest.param('dataset_1', 'other_test_user', 401, id='wrong-user'),
    pytest.param('dataset_doi', 'test_user', 400, id='with-doi')
])
def test_assign_doi_dataset(client, data, test_user_auth, other_test_user_auth, dataset_id, user, status_code):
    auth = None
    if user == 'test_user':
        auth = test_user_auth
    if user == 'other_test_user':
        auth = other_test_user_auth
    response = client.post(
        'datasets/%s/doi' % dataset_id, headers=auth)

    json_response = assert_response(response, status_code=status_code)
    if json_response is None:
        return

    dataset = json_response['data']
    assert_dataset(dataset, user_id=test_user.user_id)
    assert dataset.doi is not None
