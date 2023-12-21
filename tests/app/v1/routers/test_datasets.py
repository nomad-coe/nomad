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
from urllib.parse import urlencode
from datetime import datetime

from nomad.datamodel import Dataset
from nomad import processing
from nomad.search import search
from nomad.app.v1.models import Query, Any_
from nomad.utils.exampledata import ExampleData

from tests.conftest import admin_user_id

from tests.conftest import example_data  # pylint: disable=unused-import
from .common import assert_response

"""
These are the tests for all API operations below ``datasets``. The tests are organized
using the following type of methods: fixtures, ``perfrom_*_test``, ``assert_*``, and
``test_*``. While some ``test_*`` methods test individual API operations, some
test methods will test multiple API operations that use common aspects like
supporting queries, pagination, or the owner parameter. The test methods will use
``perform_*_test`` methods as an parameter. Similarely, the ``assert_*`` methods allow
to assert for certain aspects in the responses.
"""


def create_dataset(**kwargs):
    dataset = Dataset(
        dataset_create_time=datetime.now(),
        dataset_modified_time=datetime.now(),
        **kwargs,
    )
    dataset.m_get_annotations('mongo').save()
    return dataset


@pytest.fixture(scope='function')
def data(elastic, raw_files, mongo, test_user, other_test_user):
    data = ExampleData(main_author=test_user)
    data.create_upload(upload_id='upload_1', published=True)
    data.create_entry(
        upload_id='upload_1',
        entry_id='entry_1',
        mainfile='test_content/1/mainfile.json',
        datasets=[
            create_dataset(
                dataset_id='dataset_1',
                user_id=test_user.user_id,
                dataset_name='test dataset 1',
                dataset_type='owned',
            ),
            create_dataset(
                dataset_id='dataset_2',
                user_id=test_user.user_id,
                dataset_name='test dataset 2',
                dataset_type='owned',
            ),
        ],
    )

    data.create_entry(
        upload_id='upload_1',
        entry_id='entry_2',
        mainfile='test_content/2/mainfile.json',
        datasets=[
            create_dataset(
                dataset_id='dataset_listed',
                user_id=test_user.user_id,
                dataset_name='foreign test dataset',
                dataset_type='foreign',
            ),
            create_dataset(
                dataset_id='dataset_doi',
                user_id=test_user.user_id,
                dataset_name='foreign test dataset',
                dataset_type='foreign',
                doi='test_doi',
            ),
        ],
    )
    data.create_upload(upload_id='other_data', published=True)
    for i in range(1, 4):
        data.create_entry(
            upload_id='other_data',
            entry_id='id_%02d' % i,
            mainfile='test_content/%02d/mainfile.json' % i,
        )

    data.save(with_files=False)

    return data


def assert_pagination(pagination):
    page = pagination['page']
    total = pagination['total']
    page_size = pagination['page_size']
    assert page >= 1
    assert ('page_after_value' in pagination) == (page != 1)
    assert ('next_page_after_value' in pagination) == (total > page * page_size)
    assert 'page_url' in pagination
    assert 'first_page_url' in pagination
    if 'next_page_after_value' in pagination:
        assert 'next_page_url' in pagination
    if page > 1:
        assert 'prev_page_url' in pagination


def assert_dataset(
    dataset,
    query: Query = None,
    entries: List[str] = None,
    n_entries: int = -1,
    **kwargs,
):
    for key, value in kwargs.items():
        if key == 'prefix':
            assert dataset['dataset_name'].startswith(value)
        else:
            assert dataset[key] == value

    dataset_id = dataset['dataset_id']

    mongo_dataset = Dataset.m_def.a_mongo.objects(dataset_id=dataset_id).first()
    assert mongo_dataset is not None
    for quantity in Dataset.m_def.quantities:  # pylint: disable=not-an-iterable
        if quantity in [Dataset.query, Dataset.entries]:
            if dataset.get('dataset_type') == 'owned':
                quantity.name not in dataset
            else:
                quantity.name in dataset
        elif quantity in [Dataset.pid, Dataset.doi]:
            assert quantity.name not in dataset or dataset[quantity.name] is not None
        else:
            assert quantity.name in dataset
            assert dataset[quantity.name] is not None

    if entries is not None:
        search_results = search(
            owner='user',
            query={'entry_id': Any_(any=entries)},
            user_id=dataset['user_id'],
        )
        n_entries = search_results.pagination.total
        assert n_entries <= len(entries)
    if query is not None:
        search_results = search(owner='user', query=query, user_id=dataset['user_id'])
        n_entries = search_results.pagination.total

    if n_entries == -1:
        return

    search_results = search(owner='public', query={'datasets.dataset_id': dataset_id})

    expected_n_entries = n_entries if dataset['dataset_type'] == 'owned' else 0
    assert search_results.pagination.total == expected_n_entries
    assert processing.Entry.objects(datasets=dataset_id).count() == expected_n_entries


def assert_dataset_deleted(dataset_id):
    mongo_dataset = Dataset.m_def.a_mongo.objects(dataset_id=dataset_id).first()
    assert mongo_dataset is None

    search_results = search(
        owner='admin', query={'datasets.dataset_id': dataset_id}, user_id=admin_user_id
    )
    assert search_results.pagination.total == 0
    assert processing.Entry.objects(datasets=dataset_id).count() == 0


@pytest.mark.parametrize(
    'query, size, status_code',
    [
        pytest.param({}, 4, 200, id='empty'),
        pytest.param({'dataset_id': 'dataset_1'}, 1, 200, id='id'),
        pytest.param({'dataset_name': 'test dataset 1'}, 1, 200, id='dataset_name'),
        pytest.param({'prefix': 'test dat'}, 2, 200, id='prefix'),
        pytest.param({'dataset_type': 'foreign'}, 2, 200, id='type'),
        pytest.param({'doi': 'test_doi'}, 1, 200, id='doi'),
        pytest.param({'dataset_id': 'DOESNOTEXIST'}, 0, 200, id='id-not-exists'),
    ],
)
def test_datasets(client, data, query, size, status_code):
    url = 'datasets/'
    if len(query) > 0:
        url += '?' + urlencode(query, doseq=True)
    response = client.get(url)

    assert_response(response, status_code=status_code)
    if status_code != 200:
        return

    json_response = response.json()
    assert len(json_response['data']) == size
    for dataset in json_response['data']:
        assert_dataset(dataset, **query)
    assert_pagination(json_response['pagination'])


@pytest.mark.parametrize(
    'dataset_id, result, status_code',
    [
        pytest.param('dataset_1', {'dataset_id': 'dataset_1'}, 200, id='plain'),
        pytest.param('DOESNOTEXIST', None, 404, id='not-exists'),
    ],
)
def test_dataset(client, data, dataset_id, result, status_code):
    response = client.get('datasets/%s' % dataset_id)

    assert_response(response, status_code=status_code)
    if status_code != 200:
        return

    assert_dataset(response.json()['data'], **result)


@pytest.mark.parametrize(
    'dataset_name, dataset_type, query, entries, user, status_code',
    [
        pytest.param(
            'another test dataset', 'foreign', None, None, 'test_user', 200, id='plain'
        ),
        pytest.param(
            'another test dataset', 'foreign', None, None, None, 401, id='no-user'
        ),
        pytest.param(
            'test dataset 1', 'foreign', None, None, 'test_user', 400, id='exists'
        ),
        pytest.param(
            'another test dataset', 'owned', None, None, 'test_user', 200, id='owned'
        ),
        pytest.param(
            'another test dataset',
            'owned',
            {},
            None,
            'test_user',
            200,
            id='owned-owner-query',
        ),
        pytest.param(
            'another test dataset',
            'owned',
            {},
            None,
            'other_test_user',
            200,
            id='owned-non-owner-query',
        ),
        pytest.param(
            'another test dataset',
            'owned',
            None,
            ['id_01', 'id_02'],
            'test_user',
            200,
            id='owned-owner-entries',
        ),
        pytest.param(
            'another test dataset',
            'owned',
            None,
            ['id_01', 'id_02'],
            'other_test_user',
            200,
            id='owned-non-owner-entries',
        ),
        pytest.param(
            'another test dataset',
            'foreign',
            {},
            None,
            'test_user',
            200,
            id='foreign-query',
        ),
        pytest.param(
            'another test dataset',
            'foreign',
            None,
            ['id_01', 'id_02'],
            'test_user',
            200,
            id='foreign-entries',
        ),
    ],
)
def test_post_datasets(
    client,
    data,
    example_data,
    test_user,
    test_user_auth,
    other_test_user,
    other_test_user_auth,
    dataset_name,
    dataset_type,
    query,
    entries,
    user,
    status_code,
):
    dataset = {'dataset_name': dataset_name, 'dataset_type': dataset_type}
    if query is not None:
        dataset['query'] = query
    if entries is not None:
        dataset['entries'] = entries
    auth = None
    if user == 'test_user':
        auth = test_user_auth
        user = test_user
    elif user == 'other_test_user':
        auth = other_test_user_auth
        user = other_test_user
    response = client.post('datasets/', headers=auth, json=dataset)

    assert_response(response, status_code=status_code)
    if status_code != 200:
        return

    json_response = response.json()
    dataset = json_response['data']
    assert_dataset(
        dataset,
        query=query,
        entries=entries,
        user_id=user.user_id,
        dataset_name=dataset_name,
        dataset_type=dataset_type,
    )
    assert Dataset.m_def.a_mongo.objects().count() == 5


@pytest.mark.parametrize(
    'dataset_id, user, status_code',
    [
        pytest.param('dataset_listed', 'test_user', 200, id='plain'),
        pytest.param('dataset_listed', None, 401, id='no-user'),
        pytest.param('dataset_listed', 'other_test_user', 401, id='wrong-user'),
        pytest.param('DOESNOTEXIST', 'test_user', 404, id='does-not-exist'),
        pytest.param('dataset_doi', 'test_user', 400, id='with-doi'),
    ],
)
def test_delete_dataset(
    client, data, test_user_auth, other_test_user_auth, dataset_id, user, status_code
):
    auth = None
    if user == 'test_user':
        auth = test_user_auth
    if user == 'other_test_user':
        auth = other_test_user_auth
    response = client.delete('datasets/%s' % dataset_id, headers=auth)

    assert_response(response, status_code=status_code)
    if status_code != 200:
        return

    assert Dataset.m_def.a_mongo.objects().count() == 3
    assert_dataset_deleted(dataset_id)


@pytest.mark.parametrize(
    'dataset_id, user, status_code',
    [
        pytest.param('dataset_1', 'test_user', 200, id='plain'),
        pytest.param('dataset_1', None, 401, id='no-user'),
        pytest.param('dataset_1', 'other_test_user', 401, id='wrong-user'),
        pytest.param('dataset_doi', 'test_user', 400, id='with-doi'),
        pytest.param('unpublished', 'test_user', 400, id='unpublished'),
        pytest.param('empty', 'test_user', 400, id='empty'),
    ],
)
def test_assign_doi_dataset(
    client,
    data,
    test_user,
    test_user_auth,
    other_test_user_auth,
    dataset_id,
    user,
    status_code,
):
    more_data = ExampleData(main_author=test_user)
    more_data.create_upload(upload_id='unpublished', published=False)
    more_data.create_entry(
        upload_id='unpublished',
        entry_id='unpublished',
        mainfile='test_content/1/mainfile.json',
        datasets=[
            create_dataset(
                dataset_id='unpublished',
                user_id=test_user.user_id,
                dataset_name='test unpublished entries',
                dataset_type='owned',
            )
        ],
    )

    create_dataset(
        dataset_id='empty',
        user_id=test_user.user_id,
        dataset_name='test empty dataset',
        dataset_type='owned',
    )

    more_data.save(with_files=False)

    auth = None
    if user == 'test_user':
        auth = test_user_auth
    if user == 'other_test_user':
        auth = other_test_user_auth
    response = client.post('datasets/%s/action/doi' % dataset_id, headers=auth)

    assert_response(response, status_code=status_code)
    if status_code != 200:
        return

    json_response = response.json()
    dataset = json_response['data']
    assert_dataset(dataset, user_id=test_user.user_id)
    assert dataset['doi'] is not None
