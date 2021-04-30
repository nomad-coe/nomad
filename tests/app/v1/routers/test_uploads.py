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
import os
import time
from typing import List, Dict, Any, Iterable
from tests.utils import build_url
from tests.test_files import assert_upload_files
from tests.search import assert_search_upload
from tests.app.v1.routers.common import assert_response
from nomad import config, files, infrastructure
from nomad.processing import Upload, Calc, SUCCESS
from nomad.files import UploadFiles, PublicUploadFiles
from nomad.app.v1.routers.auth import generate_upload_token
from nomad.datamodel import EntryMetadata

'''
These are the tests for all API operations below ``uploads``. The tests are organized
using the following type of methods: fixtures, ``perfrom_*_test``, ``assert_*``, and
``test_*``. While some ``test_*`` methods test individual API operations, some
test methods will test multiple API operations that use common aspects like
supporting queries, pagination, or the owner parameter. The test methods will use
``perform_*_test`` methods as an parameter. Similarely, the ``assert_*`` methods allow
to assert for certain aspects in the responses.
'''


def perform_get(client, base_url, user_auth=None, accept='application/json', **query_args):
    headers = {'Accept': accept}
    if user_auth:
        headers.update(user_auth)
    response = client.get(build_url(base_url, query_args), headers=headers)
    return response


def perform_post_upload(
        client, mode, file, user_auth=None, token=None, accept='application/json', **query_args):
    ''' Posts a new upload. '''
    headers = {'Accept': accept}
    if user_auth:
        headers.update(user_auth)
    if mode == 'local_path':
        query_args.update(local_path=file)
    if token:
        query_args.update(token=token)
    url = build_url('uploads', query_args)
    if mode == 'multipart':
        with open(file, 'rb') as f:
            response = client.post(
                url, files={'file': f}, headers=headers)
    elif mode == 'stream':
        with open(file, 'rb') as f:
            response = client.post(url, data=f.read(), headers=headers)
    elif mode == 'local_path':
        response = client.post(url, headers=headers)
    else:
        assert False, 'Invalid value for mode provided'

    return response


def perform_post_upload_action(client, user_auth, upload_id, action, **query_args):
    return client.post(
        build_url(f'uploads/{upload_id}/action/{action}', query_args), headers=user_auth)


def assert_upload(response_json, **kwargs):
    data = response_json['data']
    assert 'upload_id' in response_json
    assert 'upload_id' in data
    assert 'create_time' in data

    for key, value in kwargs.items():
        assert data.get(key, None) == value
    return data


def assert_upload_does_not_exist(client, upload_id: str, user_auth):
    block_until_completed(client, upload_id, user_auth)

    response = perform_get(client, 'uploads/{upload_id}', user_auth)
    assert_response(response, 404)

    assert Upload.objects(upload_id=upload_id).first() is None
    assert Calc.objects(upload_id=upload_id).count() is 0

    mongo_db = infrastructure.mongo_client[config.mongo.db_name]
    mongo_collection = mongo_db['archive']
    assert mongo_collection.find({}).count() == 0

    upload_files = UploadFiles.get(upload_id)
    assert upload_files is None or isinstance(upload_files, PublicUploadFiles)


def assert_processing(client, upload_id, user_auth, check_search=True, check_files=True, published=False):
    response_data = block_until_completed(client, upload_id, user_auth)

    assert len(response_data['tasks']) == 4
    assert response_data['tasks_status'] == SUCCESS
    assert response_data['current_task'] == 'cleanup'
    assert not response_data['process_running']

    response = perform_get(client, f'uploads/{upload_id}/entries', user_auth)
    assert_response(response, 200)
    response_json = response.json()
    response_data = response_json['data']
    for entry in response_json['data']:
        assert entry['tasks_status'] == SUCCESS
        assert entry['current_task'] == 'archiving'
        assert len(entry['tasks']) == 3
        assert response_json['pagination']['total'] < response_json['pagination']['page_size']

    entries = get_upload_entries_metadata(response_data)
    if check_files:
        expected_file_class = files.PublicUploadFiles if published else files.StagingUploadFiles
        assert_upload_files(upload_id, entries, expected_file_class)
    if check_search:
        assert_search_upload(entries, additional_keys=['atoms', 'dft.system'])


def assert_gets_published(client, upload_id, user_auth, from_oasis=False, **query_args):
    with_embargo = query_args.get('with_embargo', True)
    embargo_length = query_args.get('embargo_length', 36)

    block_until_completed(client, upload_id, user_auth)

    upload_proc = Upload.objects(upload_id=upload_id).first()
    assert upload_proc is not None
    assert upload_proc.published is True
    assert upload_proc.from_oasis == from_oasis
    if with_embargo:
        assert upload_proc.embargo_length == embargo_length

    with upload_proc.entries_metadata() as entries:
        for entry in entries:
            assert entry.with_embargo == with_embargo

    assert_upload_files(upload_id, entries, files.PublicUploadFiles, published=True)


def assert_entry(entry, **kwargs):
    ''' Checks the content of a returned entry dictionary. '''
    assert 'upload_id' in entry
    assert 'entry_id' in entry
    assert 'calc_id' not in entry
    assert 'create_time' in entry
    assert not entry['process_running'] and not entry['tasks_running']
    for key, value in kwargs.items():
        assert entry.get(key, None) == value


def assert_pagination(pagination, expected_pagination):
    ''' Checks that the contents of `paginaion` matches what is expected. '''
    for key, value in expected_pagination.items():
        if value is None:
            assert key not in pagination, f'No value expected for {key}, got {pagination[key]}'
        elif value is Any:
            assert pagination.get(key) is not None, f'Value expected for {key}, got None'
        else:
            assert pagination.get(key) == value, f'For {key} we expecte {value}, but got {pagination.get(key)}'


def block_until_completed(client, upload_id: str, user_auth):
    ''' Blocks until the processing of the given upload is finished. '''
    start_time = time.time()
    while time.time() - start_time < config.tests.default_timeout:
        time.sleep(0.1)
        response = client.get('uploads/%s' % upload_id, headers=user_auth)
        if response.status_code == 200:
            response_json = response.json()
            assert_upload(response_json)
            response_data = response_json['data']
            if not response_data['process_running'] and not response_data['tasks_running']:
                return response_data
        elif response.status_code == 404:
            return None
        else:
            raise Exception(
                'unexpected status code while blocking for upload processing: %s' %
                str(response.status_code))
    raise Exception('Timed out while waiting for upload processing to finish')


def get_upload_entries_metadata(entries: List[Dict[str, Any]]) -> Iterable[EntryMetadata]:
    ''' Create a iterable of :class:`EntryMetadata` from a API upload json record. '''
    return [
        EntryMetadata(domain='dft', calc_id=entry['entry_id'], mainfile=entry['mainfile'])
        for entry in entries]


@pytest.mark.parametrize('kwargs', [
    pytest.param(
        dict(
            expected_upload_ids=['id_embargo', 'id_unpublished', 'id_published', 'id_processing', 'id_empty'],
            expected_pagination={
                'total': 5, 'page': 1, 'page_after_value': None, 'next_page_after_value': None,
                'page_url': Any, 'next_page_url': None, 'prev_page_url': None, 'first_page_url': Any}
        ), id='no-args'),
    pytest.param(
        dict(
            user='other_test_user',
            expected_upload_ids=[],
        ), id='other_test_user'),
    pytest.param(
        dict(
            query_params={'is_processing': True},
            expected_upload_ids=['id_processing'],
        ), id='filter-is_processing-True'),
    pytest.param(
        dict(
            query_params={'is_processing': False},
            expected_upload_ids=['id_embargo', 'id_unpublished', 'id_published', 'id_empty'],
        ), id='filter-is_processing-False'),
    pytest.param(
        dict(
            query_params={'is_published': True},
            expected_upload_ids=['id_embargo', 'id_published'],
        ), id='filter-is_published-True'),
    pytest.param(
        dict(
            query_params={'is_published': False},
            expected_upload_ids=['id_unpublished', 'id_processing', 'id_empty'],
        ), id='filter-is_published-False'),
    pytest.param(
        dict(
            query_params={'upload_id': 'id_published'},
            expected_upload_ids=['id_published'],
        ), id='filter-upload_id-single'),
    pytest.param(
        dict(
            query_params={'upload_id': ['id_published', 'id_embargo']},
            expected_upload_ids=['id_embargo', 'id_published'],
        ), id='filter-upload_id-multiple'),
    pytest.param(
        dict(
            query_params={'upload_name': 'name_published'},
            expected_upload_ids=['id_published'],
        ), id='filter-upload_name-single'),
    pytest.param(
        dict(
            query_params={'upload_name': ['name_published', 'name_embargo']},
            expected_upload_ids=['id_embargo', 'id_published'],
        ), id='filter-upload_name-multiple'),
    pytest.param(
        dict(
            query_params={'page_size': 2},
            expected_upload_ids=['id_embargo', 'id_unpublished'],
            expected_pagination={
                'total': 5, 'page': 1, 'page_after_value': None, 'next_page_after_value': '1',
                'page_url': Any, 'next_page_url': Any, 'prev_page_url': None, 'first_page_url': Any}
        ), id='pag-page-1'),
    pytest.param(
        dict(
            query_params={'page_size': 2, 'page': 2},
            expected_upload_ids=['id_published', 'id_processing'],
            expected_pagination={
                'total': 5, 'page': 2, 'page_after_value': '1', 'next_page_after_value': '3',
                'page_url': Any, 'next_page_url': Any, 'prev_page_url': Any, 'first_page_url': Any}
        ), id='pag-page-2'),
    pytest.param(
        dict(
            query_params={'page_size': 2, 'page': 3},
            expected_upload_ids=['id_empty'],
            expected_pagination={
                'total': 5, 'page': 3, 'page_after_value': '3', 'next_page_after_value': None,
                'page_url': Any, 'next_page_url': None, 'prev_page_url': Any, 'first_page_url': Any}
        ), id='pag-page-3'),
    pytest.param(
        dict(
            query_params={'page_size': 2, 'page': 4},
            expected_status_code=400
        ), id='pag-page-out-of-range'),
    pytest.param(
        dict(
            query_params={'page_size': 2, 'order': 'desc'},
            expected_upload_ids=['id_empty', 'id_processing'],
            expected_pagination={
                'total': 5, 'page': 1, 'page_after_value': None, 'next_page_after_value': '1',
                'page_url': Any, 'next_page_url': Any, 'prev_page_url': None, 'first_page_url': Any}
        ), id='pag-page-order-desc'),
    pytest.param(
        dict(
            query_params={'order_by': 'upload_id'},
            expected_status_code=422
        ), id='pag-invalid-order_by')])
def test_get_uploads(
        client, mongo_module, test_user_auth, other_test_user_auth, admin_user_auth, example_data, kwargs):
    ''' Makes a get request to uploads in various different ways. '''
    # Extract kwargs
    user = kwargs.get('user', 'test_user')
    query_params = kwargs.get('query_params', {})
    expected_status_code = kwargs.get('expected_status_code', 200)
    expected_upload_ids = kwargs.get('expected_upload_ids', None)
    expected_pagination = kwargs.get('expected_pagination', {})
    user_auth = {
        'test_user': test_user_auth,
        'other_test_user': other_test_user_auth,
        'admin_user': admin_user_auth}[user]
    # Api call
    response = perform_get(client, 'uploads', user_auth=user_auth, **query_params)
    # Verify result
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        response_json = response.json()
        response_data = response_json['data']

        if expected_upload_ids is not None:
            assert len(response_data) == len(expected_upload_ids), (
                f'Wrong number of records returned, expected {len(expected_upload_ids)}, got {len(response_data)}')
            found_upload_ids = [upload['upload_id'] for upload in response_data]
            assert expected_upload_ids == found_upload_ids, (
                f'Wrong upload is list returned. Expected {repr(expected_upload_ids)}, got {repr(found_upload_ids)}.')

        assert_pagination(response_json['pagination'], expected_pagination)


@pytest.mark.parametrize('user, upload_id, expected_status_code', [
    pytest.param('test_user', 'id_unpublished', 200, id='valid-upload_id'),
    pytest.param('test_user', 'silly_value', 404, id='invalid-upload_id'),
    pytest.param('other_test_user', 'id_unpublished', 401, id='no-access'),
    pytest.param('admin_user', 'id_unpublished', 200, id='admin-access')])
def test_get_upload(
        client, mongo_module, test_user_auth, other_test_user_auth, admin_user_auth, example_data,
        user, upload_id, expected_status_code):
    ''' Tests the endpoint for getting an upload by upload_id. '''
    user_auth = {
        'test_user': test_user_auth,
        'other_test_user': other_test_user_auth,
        'admin_user': admin_user_auth}[user]
    response = perform_get(client, f'uploads/{upload_id}', user_auth)
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        assert_upload(response.json())


@pytest.mark.parametrize('kwargs', [
    pytest.param(
        dict(
            expected_data_len=2,
            expected_response={'processing_successful': 2, 'processing_failed': 0},
            expected_pagination={
                'total': 2, 'page': 1, 'page_after_value': None, 'next_page_after_value': None,
                'page_url': Any, 'next_page_url': None, 'prev_page_url': None, 'first_page_url': Any}),
        id='no-args'),
    pytest.param(
        dict(
            user='other_test_user',
            expected_status_code=401),
        id='no-access'),
    pytest.param(
        dict(
            user='admin_user',
            expected_data_len=2),
        id='admin-access'),
    pytest.param(
        dict(
            upload_id='silly_value',
            expected_status_code=404),
        id='invalid-upload_id'),
    pytest.param(
        dict(
            query_args={'page_size': 1},
            expected_data_len=1,
            expected_response={'processing_successful': 2, 'processing_failed': 0},
            expected_pagination={
                'total': 2, 'page': 1, 'page_after_value': None, 'next_page_after_value': '0', 'order_by': 'mainfile',
                'page_url': Any, 'next_page_url': Any, 'prev_page_url': None, 'first_page_url': Any}),
        id='pag-page-1'),
    pytest.param(
        dict(
            query_args={'page_size': 1, 'page': 1},
            expected_data_len=1,
            expected_response={'processing_successful': 2, 'processing_failed': 0},
            expected_pagination={
                'total': 2, 'page': 1, 'page_after_value': None, 'next_page_after_value': '0', 'order_by': 'mainfile',
                'page_url': Any, 'next_page_url': Any, 'prev_page_url': None, 'first_page_url': Any}),
        id='pag-page-1-by-page'),
    pytest.param(
        dict(
            query_args={'page_size': 1, 'page': 2},
            expected_data_len=1,
            expected_response={'processing_successful': 2, 'processing_failed': 0},
            expected_pagination={
                'total': 2, 'page': 2, 'page_after_value': '0', 'next_page_after_value': None, 'order_by': 'mainfile',
                'page_url': Any, 'next_page_url': None, 'prev_page_url': Any, 'first_page_url': Any}),
        id='pag-page-2-by-page'),
    pytest.param(
        dict(
            query_args={'page_size': 1, 'page_after_value': '0'},
            expected_data_len=1,
            expected_response={'processing_successful': 2, 'processing_failed': 0},
            expected_pagination={
                'total': 2, 'page': 2, 'page_after_value': '0', 'next_page_after_value': None, 'order_by': 'mainfile',
                'page_url': Any, 'next_page_url': None, 'prev_page_url': Any, 'first_page_url': Any}),
        id='pag-page-2-by-page_after_value'),
    pytest.param(
        dict(
            query_args={'page_size': 0},
            expected_data_len=0,
            expected_response={'processing_successful': 2, 'processing_failed': 0},
            expected_pagination={
                'total': 2, 'page': 1, 'page_after_value': None, 'next_page_after_value': None, 'order_by': 'mainfile',
                'page_url': Any, 'next_page_url': None, 'prev_page_url': None, 'first_page_url': None}),
        id='pag-page_size-zero'),
    pytest.param(
        dict(
            query_args={'page_size': 1, 'page': 3},
            expected_status_code=400),
        id='pag-out-of-rage-page'),
    pytest.param(
        dict(
            query_args={'page_size': 1, 'page_after_value': '1'},
            expected_status_code=400),
        id='pag-out-of-rage-page_after_value'),
    pytest.param(
        dict(
            query_args={'page_size': 1, 'order_by': 'parser'},
            expected_data_len=1,
            expected_response={'processing_successful': 2, 'processing_failed': 0},
            expected_pagination={
                'total': 2, 'page': 1, 'page_after_value': None, 'next_page_after_value': '0', 'order_by': 'parser',
                'page_url': Any, 'next_page_url': Any, 'prev_page_url': None, 'first_page_url': Any}),
        id='pag-order_by-parser'),
    pytest.param(
        dict(
            query_args={'page_size': 1, 'order_by': 'calc_id'},
            expected_status_code=422),
        id='pag-order_by-illegal'),
    pytest.param(
        dict(
            query_args={'page_size': 1, 'page': 2, 'page_after_value': '0'},
            expected_status_code=422),
        id='pag-overspecified')])
def test_get_upload_entries(
        client, mongo_module, test_user_auth, other_test_user_auth, admin_user_auth, example_data,
        kwargs):
    '''
    Fetches the entries for a specific upload, by calling uploads/{upload_id}/entries,
    with the provided query paramters, and checks the result.
    '''
    upload_id = kwargs.get('upload_id', 'id_embargo')
    user = kwargs.get('user', 'test_user')
    query_args = kwargs.get('query_args', {})
    expected_status_code = kwargs.get('expected_status_code', 200)
    expected_data_len = kwargs.get('expected_data_len', 2)
    expected_response = kwargs.get('expected_response', {})
    expected_pagination = kwargs.get('expected_pagination', {})
    user_auth = {
        'test_user': test_user_auth,
        'other_test_user': other_test_user_auth,
        'admin_user': admin_user_auth}[user]

    response = perform_get(client, f'uploads/{upload_id}/entries', user_auth, **query_args)
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        response_json = response.json()
        response_data = response_json['data']

        if expected_data_len is not None:
            assert len(response_data) == expected_data_len

        for entry in response_data:
            assert_entry(entry)

        for key, value in expected_response.items():
            assert response_json.get(key, None) == value

        pagination = response_json['pagination']
        assert_pagination(pagination, expected_pagination)


@pytest.mark.parametrize('upload_id, entry_id, user, expected_status_code', [
    pytest.param('id_embargo', 'id_embargo', 'test_user', 200, id='ok'),
    pytest.param('id_embargo', 'id_embargo', 'other_test_user', 401, id='no-access'),
    pytest.param('id_embargo', 'id_embargo', 'admin_user', 200, id='admin-access'),
    pytest.param('silly_value', 'id_embargo', 'test_user', 404, id='invalid-upload_id'),
    pytest.param('id_embargo', 'silly_value', 'test_user', 404, id='invalid-entry_id')])
def test_get_upload_entry(
        client, mongo_module, test_user_auth, other_test_user_auth, admin_user_auth, example_data,
        upload_id, entry_id, user, expected_status_code):
    '''
    Fetches an entry via a call to uploads/{upload_id}/entries/{entry_id} and checks it.
    '''
    user_auth = {
        'test_user': test_user_auth,
        'other_test_user': other_test_user_auth,
        'admin_user': admin_user_auth}[user]
    response = perform_get(client, f'uploads/{upload_id}/entries/{entry_id}', user_auth)
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        response_json = response.json()
        assert response_json['entry_id'] == entry_id
        response_data = response_json['data']
        assert_entry(response_data)


@pytest.mark.parametrize('user, upload_id, path, accept, expected_status_code, expected_mime_type, expected_content', [
    pytest.param(
        'test_user', 'id_unpublished', 'test_content/test_entry/1.aux', '*',
        200, 'application/octet-stream', 'content', id='unpublished-file'),
    pytest.param(
        'test_user', 'id_unpublished', 'test_content/test_entry/', 'application/json',
        200, 'application/json', ['1.aux', '2.aux', '3.aux', '4.aux', 'mainfile.json'], id='unpublished-dir-json'),
    pytest.param(
        'test_user', 'id_unpublished', 'test_content/test_entry/', '*',
        200, 'text/html; charset=utf-8', ['1.aux', '2.aux', '3.aux', '4.aux', 'mainfile.json'], id='unpublished-dir-html'),
    pytest.param(
        'test_user', 'id_unpublished', '', 'application/json',
        200, 'application/json', ['test_content/'], id='unpublished-dir-json-root'),
    pytest.param(
        'other_test_user', 'id_unpublished', 'test_content/test_entry/1.aux', '*',
        401, None, None, id='unpublished-file-unauthorized'),
    pytest.param(
        'admin_user', 'id_unpublished', 'test_content/test_entry/1.aux', '*',
        200, 'application/octet-stream', 'content', id='unpublished-file-admin-auth'),
    pytest.param(
        'test_user', 'id_published', 'test_content/subdir/test_entry_01/1.aux', '*',
        200, 'application/octet-stream', 'content', id='published-file'),
    pytest.param(
        'test_user', 'id_published', 'test_content/subdir/test_entry_01', 'application/json',
        200, 'application/json', ['1.aux', '2.aux', '3.aux', '4.aux', 'mainfile.json'], id='published-dir-json'),
    pytest.param(
        'test_user', 'id_published', 'test_content/subdir/test_entry_01', '*',
        200, 'text/html; charset=utf-8', ['1.aux', '2.aux', '3.aux', '4.aux', 'mainfile.json'], id='published-dir-html'),
    pytest.param(
        'test_user', 'id_published', '', 'application/json',
        200, 'application/json', ['test_content/'], id='published-dir-json-root'),
    pytest.param(
        'other_test_user', 'id_published', 'test_content/subdir/test_entry_01/1.aux', '*',
        401, None, None, id='published-file-unauthorized'),
    pytest.param(
        'admin_user', 'id_published', 'test_content/subdir/test_entry_01/1.aux', '*',
        200, 'application/octet-stream', 'content', id='published-file-admin-auth'),
])
def test_get_upload_raw_path(
        client, example_data, test_user_auth, other_test_user_auth, admin_user_auth,
        user, upload_id, path, accept, expected_status_code, expected_mime_type, expected_content):
    user_auth = {
        'test_user': test_user_auth,
        'other_test_user': other_test_user_auth,
        'admin_user': admin_user_auth}[user]
    response = perform_get(client, f'uploads/{upload_id}/raw/{path}', user_auth=user_auth, accept=accept)
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        mime_type = response.headers.get('content-type')
        assert mime_type == expected_mime_type
        if mime_type == 'application/octet-stream':
            if expected_content:
                assert expected_content in response.text, 'Expected content not found'
        elif mime_type == 'application/json':
            data = response.json()
            assert data['path'] == (path.rstrip('/') or '.')
            if expected_content:
                assert data['content'] == expected_content, 'Incorrect list of files returned'
        elif mime_type == 'text/html':
            assert response.text, 'No response text'
            if expected_content:
                for name in expected_content:
                    assert name in response.text


@pytest.mark.parametrize('mode, name, user, use_upload_token, empty, publish_directly, test_limit, accept_json, expected_status_code', [
    pytest.param('multipart', 'test_name', 'test_user', False, False, None, False, True, 200, id='multipart'),
    pytest.param('multipart', None, 'test_user', False, False, None, False, True, 200, id='multipart-no-name'),
    pytest.param('stream', 'test_name', 'test_user', False, False, None, False, True, 200, id='stream'),
    pytest.param('stream', None, 'test_user', False, False, None, False, True, 200, id='stream-no-name'),
    pytest.param('stream', None, 'test_user', False, False, None, False, False, 200, id='stream-no-accept-json'),
    pytest.param('multipart', None, 'invalid', False, False, None, False, True, 401, id='multipart-no-name-invalid-cred'),
    pytest.param('stream', None, 'invalid', False, False, None, False, True, 401, id='stream-no-name-invalid-cred'),
    pytest.param('multipart', 'test_name', 'test_user', True, False, None, False, True, 200, id='multipart-token'),
    pytest.param('stream', 'test_name', 'test_user', True, False, None, False, True, 200, id='stream-token'),
    pytest.param('multipart', 'test_name', 'invalid', True, False, None, False, True, 401, id='multipart-token-invalid-cred'),
    pytest.param('stream', 'test_name', 'invalid', True, False, None, False, True, 401, id='stream-token-invalid-cred'),
    pytest.param('local_path', None, 'admin_user', False, False, None, False, True, 200, id='local_path'),
    pytest.param('multipart', None, None, False, False, None, False, True, 401, id='not-logged-in-multipart'),
    pytest.param('stream', None, None, False, False, None, False, True, 401, id='not-logged-in-stream'),
    pytest.param('local_path', None, None, False, False, None, False, True, 401, id='not-logged-in-local_path'),
    pytest.param('local_path', None, 'test_user', False, False, None, False, True, 401, id='not-admin-local_path'),
    pytest.param('stream', 'test_name', 'test_user', False, False, True, False, True, 200, id='publish_directly'),
    pytest.param('stream', 'test_name', 'test_user', False, True, True, False, True, 200, id='publish_directly-empty'),
    pytest.param('stream', 'test_name', 'test_user', False, False, None, True, True, 400, id='upload-limit-exceeded')])
def test_post_upload(
        client, mongo, proc_infra, monkeypatch, test_user, admin_user, test_user_auth, admin_user_auth,
        empty_upload, non_empty_example_upload,
        mode, name, user, use_upload_token, empty, publish_directly, test_limit, accept_json, expected_status_code):
    '''
    Posts an upload, with different arguments.
    '''
    if user == 'test_user':
        user_auth = test_user_auth
        token = generate_upload_token(test_user)
    elif user == 'admin_user':
        user_auth = admin_user_auth
        token = generate_upload_token(admin_user)
    elif user == 'invalid':
        user_auth = {'Authorization': 'Bearer JUST-MADE-IT-UP'}
        token = 'invalid.token'
    else:
        user_auth = None
        token = None
    # Use either token or bearer token for the post operation
    user_auth_post = user_auth
    if use_upload_token:
        user_auth_post = None
    else:
        token = None

    if empty:
        upload_file = empty_upload
    else:
        upload_file = non_empty_example_upload

    if test_limit:
        monkeypatch.setattr('nomad.config.services.upload_limit', 0)

    accept = 'application/json' if accept_json else '*'

    response = perform_post_upload(
        client, mode, upload_file, user_auth_post, token, accept=accept,
        name=name, publish_directly=publish_directly)

    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        if accept_json:
            response_json = response.json()
            upload_id = response_json['upload_id']
            expected_name = name
            if not expected_name and mode in ('multipart', 'local_path'):
                expected_name = os.path.basename(non_empty_example_upload)
            assert_upload(response_json, name=expected_name)
            if mode == 'local_path':
                assert response_json['data']['upload_path'] == non_empty_example_upload

            assert_processing(client, upload_id, user_auth, published=(publish_directly and not empty))

            if publish_directly:
                upload_proc = Upload.objects(upload_id=upload_id).first()
                if empty:
                    assert not upload_proc.published
                else:
                    assert_gets_published(client, upload_id, test_user_auth, with_embargo=False)
        else:
            assert 'Thanks for uploading' in response.text


@pytest.mark.parametrize('user, oasis_uploader, oasis_upload_id, oasis_deployment_id, expected_status_code', [
    pytest.param('test_user', 'test_user', 'oasis_upload_id', 'an_id', 200, id='ok'),
    pytest.param('test_user', 'test_user', 'id_unpublished_w', 'an_id', 400, id='dulicate'),
    pytest.param('test_user', None, 'oasis_upload_id', 'an_id', 400, id='missing-oasis_uploader_id'),
    pytest.param('test_user', 'test_user', None, 'an_id', 400, id='missing-oasis_upload_id'),
    pytest.param('test_user', 'test_user', 'oasis_upload_id', None, 400, id='missing-oasis_deployment_id'),
    pytest.param('other_test_user', 'test_user', 'oasis_upload_id', 'an_id', 401, id='not-oasis-admin')])
def test_post_upload_oasis(
        client, mongo, proc_infra, oasis_example_upload, example_data_writeable,
        test_user, other_test_user, test_user_auth, other_test_user_auth,
        user, oasis_uploader, oasis_upload_id, oasis_deployment_id, expected_status_code):

    auth_dict = {
        'test_user': (test_user.user_id, test_user_auth),
        'other_test_user': (other_test_user.user_id, other_test_user_auth),
        None: (None, None)}
    __, user_auth = auth_dict[user]
    oasis_uploader_id, __ = auth_dict[oasis_uploader]

    response = perform_post_upload(
        client, 'stream', oasis_example_upload, user_auth,
        oasis_upload_id=oasis_upload_id,
        oasis_uploader_id=oasis_uploader_id,
        oasis_deployment_id=oasis_deployment_id)

    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        response_json = response.json()
        upload_id = response_json['upload_id']
        assert upload_id == oasis_upload_id
        assert_upload(response_json)
        assert_processing(client, upload_id, user_auth, published=True, check_search=False)
        assert_gets_published(client, upload_id, user_auth, from_oasis=True, with_embargo=False)


@pytest.mark.parametrize('kwargs', [
    pytest.param(
        dict(
            expected_status_code=200),
        id='no-args'),
    pytest.param(
        dict(
            query_args={'with_embargo': True, 'embargo_length': 12},
            expected_status_code=200),
        id='non-standard-embargo'),
    pytest.param(
        dict(
            query_args={'embargo_length': 24},
            expected_status_code=200),
        id='non-standard-embargo-length-only'),
    pytest.param(
        dict(
            query_args={'embargo_length': 100},
            expected_status_code=400),
        id='illegal-embargo-length'),
    pytest.param(
        dict(
            query_args={'with_embargo': False},
            expected_status_code=200),
        id='no-embargo'),
    pytest.param(
        dict(
            upload_id='id_empty_w',
            expected_status_code=400),
        id='empty'),
    pytest.param(
        dict(
            upload_id='id_processing_w',
            expected_status_code=400),
        id='processing'),
    pytest.param(
        dict(
            upload_id='id_published_w',
            expected_status_code=401),
        id='already-published'),
    pytest.param(
        dict(
            user='other_test_user',
            expected_status_code=401),
        id='not-my-upload')])
def test_post_upload_action_publish(
        client, proc_infra, example_data_writeable, test_user_auth, other_test_user_auth, admin_user_auth,
        kwargs):
    ''' Tests the publish action with various arguments. '''
    upload_id = kwargs.get('upload_id', 'id_unpublished_w')
    query_args = kwargs.get('query_args', {})
    expected_status_code = kwargs.get('expected_status_code', 200)
    user = kwargs.get('user', 'test_user')
    user_auth = {
        'test_user': test_user_auth,
        'other_test_user': other_test_user_auth,
        'admin_user': admin_user_auth}[user]

    response = perform_post_upload_action(client, user_auth, upload_id, 'publish', **query_args)

    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        upload = assert_upload(response.json())
        assert upload['current_process'] == 'publish_upload'
        assert upload['process_running']

        assert_gets_published(client, upload_id, user_auth, **query_args)


@pytest.mark.parametrize('upload_id, user, expected_status_code', [
    pytest.param('examples_template', 'test_user', 200, id='ok'),
    pytest.param('examples_template', 'other_test_user', 401, id='no-access'),
    pytest.param('id_unpublished_w', 'test_user', 400, id='not-published'),
    pytest.param('id_processing_w', 'test_user', 400, id='already-processing'),
    pytest.param('silly_value', 'test_user', 404, id='invalid-upload_id')])
def test_post_upload_action_reprocess(
        client, monkeypatch, example_data_writeable, published, test_user_auth, other_test_user_auth,
        upload_id, user, expected_status_code):
    monkeypatch.setattr('nomad.config.meta.version', 're_process_test_version')
    monkeypatch.setattr('nomad.config.meta.commit', 're_process_test_commit')
    user_auth = {
        'test_user': test_user_auth,
        'other_test_user': other_test_user_auth}[user]

    response = perform_post_upload_action(client, user_auth, upload_id, 're-process')
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        assert_processing(client, upload_id, test_user_auth, check_files=False, published=True)


@pytest.mark.parametrize('upload_id, user, expected_status_code', [
    pytest.param('id_unpublished_w', 'test_user', 200, id='delete-own'),
    pytest.param('id_unpublished_w', 'other_test_user', 401, id='delete-others-not-admin'),
    pytest.param('id_unpublished_w', 'admin_user', 200, id='delete-others-admin'),
    pytest.param('id_published_w', 'test_user', 401, id='delete-own-published'),
    pytest.param('id_published_w', 'admin_user', 200, id='delete-others-published-admin'),
    pytest.param('silly_value', 'test_user', 404, id='invalid-upload_id')])
def test_delete_upload(
        client, proc_infra, example_data_writeable,
        test_user_auth, other_test_user_auth, admin_user_auth,
        upload_id, user, expected_status_code):
    ''' Uploads a file, and then tries to delete it, with different parameters and users. '''
    user_auth = {
        'test_user': test_user_auth,
        'other_test_user': other_test_user_auth,
        'admin_user': admin_user_auth
    }[user]

    response = client.delete(f'uploads/{upload_id}', headers=user_auth)
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        assert_upload_does_not_exist(client, upload_id, test_user_auth)


@pytest.mark.parametrize('authorized, expected_status_code', [
    pytest.param(True, 200, id='ok'),
    pytest.param(False, 401, id='not-authorized')])
def test_get_command_examples(client, test_user_auth, authorized, expected_status_code):
    response = perform_get(
        client, 'uploads/command-examples', user_auth=test_user_auth if authorized else None)
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        data = response.json()
        for k in (
                'upload_url', 'upload_command', 'upload_command_with_name',
                'upload_progress_command', 'upload_command_form', 'upload_tar_command'):
            assert k in data
        assert '/api/v1/uploads' in data['upload_command']
