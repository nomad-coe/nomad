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
import io
import os
import requests
import time
import zipfile
from datetime import datetime
from typing import List, Dict, Any, Iterable
from tests.utils import build_url, set_upload_entry_metadata

from tests.test_files import (
    example_file_mainfile_different_atoms, example_file_vasp_with_binary, example_file_aux,
    example_file_unparsable, example_file_corrupt_zip, empty_file,
    assert_upload_files)
from tests.test_search import assert_search_upload
from tests.processing.test_edit_metadata import (
    assert_metadata_edited, all_coauthor_metadata, all_admin_metadata)
from tests.app.v1.routers.common import assert_response, assert_browser_download_headers
from nomad import config, files, infrastructure
from nomad.config.models import BundleImportSettings
from nomad.processing import Upload, Entry, ProcessStatus
from nomad.files import UploadFiles, StagingUploadFiles, PublicUploadFiles
from nomad.bundles import BundleExporter
from nomad.datamodel import EntryMetadata

from .test_entries import assert_archive_response

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
    headers = {}
    if accept:
        headers['Accept'] = accept
    if user_auth:
        headers.update(user_auth)

    response = client.get(build_url(base_url, query_args), headers=headers)
    return response


def perform_post_put_file(
        client, action, url, mode, file_paths, user_auth=None, token=None, accept='application/json',
        **query_args):
    ''' Posts or puts a file. '''
    if type(file_paths) == str:
        file_paths = [file_paths]
    headers = {'Accept': accept}
    if user_auth:
        headers.update(user_auth)
    if mode == 'local_path':
        assert len(file_paths) == 1
        query_args.update(local_path=file_paths[0])
    if token:
        query_args.update(token=token)
    url = build_url(url, query_args)

    if action == 'POST':
        func = client.post
    elif action == 'PUT':
        func = client.put
    else:
        assert False, f'Invalid action provided: {action}'

    if not file_paths:
        response = func(url, data='', headers=headers)
    else:
        if mode == 'multipart':
            if len(file_paths) == 1:
                with open(file_paths[0], 'rb') as f:
                    response = func(url, files={'file': f}, headers=headers)
            else:
                files_list = []
                open_files = []
                try:
                    for file_path in file_paths:
                        filename = os.path.basename(file_path)
                        f = open(file_path, 'rb')
                        open_files.append(f)
                        files_list.append(('file', (filename, f)))
                    response = func(url, files=files_list, headers=headers)
                finally:
                    for f in open_files:
                        f.close()
        elif mode == 'stream':
            assert len(file_paths) == 1
            with open(file_paths[0], 'rb') as f:
                response = func(url, data=f.read(), headers=headers)
        elif mode == 'local_path':
            response = func(url, headers=headers)
        else:
            assert False, f'Invalid value for mode provided {mode}'

    return response


def perform_post_upload_action(client, user_auth, upload_id, action, json=None, **query_args):
    return client.post(
        build_url(f'uploads/{upload_id}/action/{action}', query_args), headers=user_auth, json=json)


def assert_file_upload_and_processing(
        client, action, url, mode, user, test_auth_dict, upload_id,
        source_paths, target_path, query_args, accept_json, use_upload_token,
        expected_status_code, expected_process_status, expected_mainfiles, published,
        all_entries_should_succeed):
    '''
    Uploads a file, using the given action (POST or PUT), url, query arguments, and checks
    the results.
    '''
    source_paths = source_paths or []
    if type(source_paths) == str:
        source_paths = [source_paths]
    user_auth, token = test_auth_dict[user]
    # Use either token or bearer token for the post operation (never both)
    user_auth_action = user_auth
    if use_upload_token:
        user_auth_action = None
    else:
        token = None
    accept = 'application/json' if accept_json else '*'
    processed_response_data = None
    response = perform_post_put_file(
        client, action, url, mode, source_paths, user_auth_action, token, accept, **query_args)

    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        if accept_json:
            response_json = response.json()
            upload_id = response_json['upload_id']
            if expected_process_status:
                assert response_json['data']['process_status'] == expected_process_status
            assert_upload(response_json)
        else:
            assert 'Thanks for uploading' in response.text
            if not upload_id:
                return None, None

        if example_file_corrupt_zip in source_paths:
            processed_response_data = assert_processing_fails(client, upload_id, user_auth)
        else:
            processed_response_data = assert_processing(
                client, upload_id, user_auth, published=published,
                all_entries_should_succeed=all_entries_should_succeed)

            # Check that files got copied as expected
            for source_path in source_paths:
                upload_files = files.UploadFiles.get(upload_id)
                file_name = os.path.basename(source_path)
                if zipfile.is_zipfile(source_path):
                    with open(source_path, 'rb') as f:
                        zf = zipfile.ZipFile(f)
                        for path in zf.namelist():
                            if not path.endswith('/'):
                                target_path_full = os.path.join(target_path, path)
                                assert upload_files.raw_path_exists(target_path_full)
                                assert upload_files.raw_path_is_file(target_path_full)
                else:
                    if mode == 'stream':
                        # Must specify file_name
                        file_name = query_args['file_name']
                    target_path_full = os.path.join(target_path, file_name)
                    assert upload_files.raw_path_exists(target_path_full)
                    assert upload_files.raw_path_is_file(target_path_full)
                    assert upload_files.raw_file_size(target_path_full) == os.stat(source_path).st_size

        assert_expected_mainfiles(upload_id, expected_mainfiles)
    return response, processed_response_data


def assert_expected_mainfiles(upload_id, expected_mainfiles):
    if expected_mainfiles is not None:
        entries = [e.mainfile for e in Entry.objects(upload_id=upload_id)]
        assert set(entries) == set(expected_mainfiles), 'Wrong entries found'
        for entry in Entry.objects(upload_id=upload_id):
            if type(expected_mainfiles) != dict or expected_mainfiles[entry.mainfile]:
                assert entry.process_status == ProcessStatus.SUCCESS
            else:
                assert entry.process_status == ProcessStatus.FAILURE


def assert_upload(response_json, **kwargs):
    data = response_json['data']
    assert 'upload_id' in response_json
    assert 'upload_id' in data
    assert 'upload_create_time' in data
    assert 'main_author' in data
    assert 'coauthors' in data
    assert 'reviewers' in data
    assert 'viewers' in data
    assert 'writers' in data
    assert 'published' in data
    assert 'with_embargo' in data
    assert 'embargo_length' in data
    assert 'license' in data
    assert (data['embargo_length'] > 0) == data['with_embargo']
    if data['published']:
        assert 'publish_time' in data

    for key, value in kwargs.items():
        assert data.get(key, None) == value
    return data


def assert_upload_does_not_exist(client, upload_id: str, user_auth):
    block_until_completed(client, upload_id, user_auth)

    response = perform_get(client, 'uploads/{upload_id}', user_auth)
    assert_response(response, 404)

    assert Upload.objects(upload_id=upload_id).first() is None
    assert Entry.objects(upload_id=upload_id).count() is 0

    mongo_db = infrastructure.mongo_client[config.mongo.db_name]
    mongo_collection = mongo_db['archive']
    assert mongo_collection.find({}).count() == 0

    upload_files = UploadFiles.get(upload_id)
    assert upload_files is None or isinstance(upload_files, PublicUploadFiles)


def assert_processing(
        client, upload_id, user_auth, check_search=True, check_files=True, published=False,
        all_entries_should_succeed=True):
    response_data = block_until_completed(client, upload_id, user_auth)

    assert response_data['process_status'] in (ProcessStatus.SUCCESS, ProcessStatus.READY)
    assert not response_data['process_running']

    response_entries = perform_get(client, f'uploads/{upload_id}/entries', user_auth)
    assert_response(response_entries, 200)
    response_entries_json = response_entries.json()
    response_entries_data = response_entries_json['data']
    all_entries_succesful = True

    for entry in response_entries_data:
        entry_succeeded = entry['process_status'] == ProcessStatus.SUCCESS
        if not entry_succeeded:
            all_entries_succesful = False
            if all_entries_should_succeed:
                assert False, 'One or more entries failed to process'
        pagination = response_entries_json['pagination']
        assert pagination['total'] < pagination['page_size']

    entries = get_upload_entries_metadata(response_entries_data)
    if check_files:
        expected_file_class = files.PublicUploadFiles if published else files.StagingUploadFiles
        assert_upload_files(upload_id, entries, expected_file_class)
    if check_search and all_entries_succesful:
        assert_search_upload(entries, additional_keys=['results.material.elements', 'results.method.simulation.program_name'], upload_id=upload_id)
    return response_data


def assert_processing_fails(client, upload_id, user_auth):
    response_data = block_until_completed(client, upload_id, user_auth)

    assert response_data['process_status'] == ProcessStatus.FAILURE
    return response_data


def assert_gets_published(
        client, upload_id, user_auth, from_oasis=False, current_embargo_length=0, **query_args):
    embargo_length = query_args.get('embargo_length', current_embargo_length)

    block_until_completed(client, upload_id, user_auth)

    upload_proc = Upload.objects(upload_id=upload_id).first()
    assert upload_proc is not None
    assert upload_proc.published is True
    assert upload_proc.from_oasis == from_oasis
    assert upload_proc.embargo_length == embargo_length

    with upload_proc.entries_metadata() as entries:
        for entry in entries:
            assert entry.with_embargo == (embargo_length > 0)

    assert_upload_files(upload_id, entries, files.PublicUploadFiles, published=True)


def assert_entry(entry, **kwargs):
    ''' Checks the content of a returned entry dictionary. '''
    assert 'upload_id' in entry
    assert 'entry_id' in entry
    assert 'entry_create_time' in entry
    assert not entry['process_running']
    for key, value in kwargs.items():
        assert entry.get(key, None) == value
    assert 'entry_metadata' in entry


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
            if not response_data['process_running']:
                return response_data
        elif response.status_code == 404:
            return None
        else:
            raise Exception(
                'unexpected status code while blocking for upload processing: %s' %
                str(response.status_code))
    raise Exception('Timed out while waiting for upload processing to finish')


def get_upload_entries_metadata(entries: List[Dict[str, Any]]) -> Iterable[EntryMetadata]:
    '''
    Create a iterable of :class:`EntryMetadata` from a API upload json record, plus a
    with_embargo flag fetched from mongodb.
    '''
    return [
        EntryMetadata(
            domain='dft', entry_id=entry['entry_id'], mainfile=entry['mainfile'],
            with_embargo=Upload.get(entry['upload_id']).with_embargo)
        for entry in entries]


@pytest.mark.parametrize('kwargs', [
    pytest.param(
        dict(
            expected_upload_ids=[
                'id_embargo', 'id_embargo_w_coauthor', 'id_embargo_w_reviewer',
                'id_unpublished', 'id_unpublished_w_coauthor', 'id_unpublished_w_reviewer',
                'id_published', 'id_child_entries', 'id_processing', 'id_empty'],
            expected_pagination={
                'total': 10, 'page': 1, 'page_after_value': None, 'next_page_after_value': None,
                'page_url': Any, 'next_page_url': None, 'prev_page_url': None, 'first_page_url': Any}
        ), id='no-args'),
    pytest.param(
        dict(
            user='other_test_user',
            expected_upload_ids=[
                'id_embargo_w_coauthor', 'id_embargo_w_reviewer', 'id_unpublished_w_coauthor',
                'id_unpublished_w_reviewer']
        ), id='other_test_user'),
    pytest.param(
        dict(
            user=None,
            expected_status_code=401
        ), id='no-credentials'),
    pytest.param(
        dict(
            user='invalid',
            expected_status_code=401
        ), id='invalid-credentials'),
    pytest.param(
        dict(
            query_params={'is_processing': True},
            expected_upload_ids=['id_processing'],
        ), id='filter-is_processing-True'),
    pytest.param(
        dict(
            query_params={'is_processing': False},
            expected_upload_ids=[
                'id_embargo', 'id_embargo_w_coauthor', 'id_embargo_w_reviewer',
                'id_unpublished', 'id_unpublished_w_coauthor', 'id_unpublished_w_reviewer',
                'id_published', 'id_child_entries', 'id_empty'],
        ), id='filter-is_processing-False'),
    pytest.param(
        dict(
            query_params={'is_published': True},
            expected_upload_ids=['id_embargo', 'id_embargo_w_coauthor', 'id_embargo_w_reviewer', 'id_published'],
        ), id='filter-is_published-True'),
    pytest.param(
        dict(
            query_params={'is_published': False},
            expected_upload_ids=[
                'id_unpublished', 'id_unpublished_w_coauthor', 'id_unpublished_w_reviewer',
                'id_child_entries', 'id_processing', 'id_empty'],
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
            expected_upload_ids=['id_embargo', 'id_embargo_w_coauthor'],
            expected_pagination={
                'total': 10, 'page': 1, 'page_after_value': None, 'next_page_after_value': '1',
                'page_url': Any, 'next_page_url': Any, 'prev_page_url': None, 'first_page_url': Any}
        ), id='pag-page-1'),
    pytest.param(
        dict(
            query_params={'page_size': 2, 'page': 2},
            expected_upload_ids=['id_embargo_w_reviewer', 'id_unpublished'],
            expected_pagination={
                'total': 10, 'page': 2, 'page_after_value': '1', 'next_page_after_value': '3',
                'page_url': Any, 'next_page_url': Any, 'prev_page_url': Any, 'first_page_url': Any}
        ), id='pag-page-2'),
    pytest.param(
        dict(
            query_params={'page_size': 4, 'page': 3},
            expected_upload_ids=['id_processing', 'id_empty'],
            expected_pagination={
                'total': 10, 'page': 3, 'page_after_value': '7', 'next_page_after_value': None,
                'page_url': Any, 'next_page_url': None, 'prev_page_url': Any, 'first_page_url': Any}
        ), id='pag-page-3'),
    pytest.param(
        dict(
            query_params={'page_size': 5, 'page': 3},
            expected_status_code=400
        ), id='pag-page-out-of-range'),
    pytest.param(
        dict(
            query_params={'page_size': 2, 'order': 'desc'},
            expected_upload_ids=['id_empty', 'id_processing'],
            expected_pagination={
                'total': 10, 'page': 1, 'page_after_value': None, 'next_page_after_value': '1',
                'page_url': Any, 'next_page_url': Any, 'prev_page_url': None, 'first_page_url': Any}
        ), id='pag-page-order-desc'),
    pytest.param(
        dict(
            query_params={'order_by': 'upload_id'},
            expected_status_code=422
        ), id='pag-invalid-order_by')])
def test_get_uploads(
        client, mongo_module, test_auth_dict, example_data, kwargs):
    ''' Makes a get request to uploads in various different ways. '''
    # Extract kwargs
    user = kwargs.get('user', 'test_user')
    query_params = kwargs.get('query_params', {})
    expected_status_code = kwargs.get('expected_status_code', 200)
    expected_upload_ids = kwargs.get('expected_upload_ids', None)
    expected_pagination = kwargs.get('expected_pagination', {})
    user_auth, __token = test_auth_dict[user]
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
    pytest.param('test_user', 'id_child_entries', 200, id='valid-upload_id-w-child-entries'),
    pytest.param('test_user', 'silly_value', 404, id='invalid-upload_id'),
    pytest.param(None, 'id_unpublished', 401, id='no-credentials'),
    pytest.param('invalid', 'id_unpublished', 401, id='invalid-credentials'),
    pytest.param('other_test_user', 'id_unpublished', 401, id='no-access'),
    pytest.param('admin_user', 'id_unpublished', 200, id='admin-access')])
def test_get_upload(
        client, mongo_module, test_auth_dict, example_data,
        user, upload_id, expected_status_code):
    ''' Tests the endpoint for getting an upload by upload_id. '''
    user_auth, __token = test_auth_dict[user]
    response = perform_get(client, f'uploads/{upload_id}', user_auth)
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        assert_upload(response.json())


@pytest.mark.parametrize('kwargs', [
    pytest.param(
        dict(
            expected_data_len=1,
            expected_response={'processing_successful': 1, 'processing_failed': 0},
            expected_pagination={
                'total': 1, 'page': 1, 'page_after_value': None, 'next_page_after_value': None,
                'page_url': Any, 'next_page_url': None, 'prev_page_url': None, 'first_page_url': Any}),
        id='no-args'),
    pytest.param(
        dict(
            upload_id='id_child_entries',
            expected_data_len=3,
            expected_response={'processing_successful': 3, 'processing_failed': 0},
            expected_pagination={
                'total': 3, 'page': 1, 'page_after_value': None, 'next_page_after_value': None,
                'page_url': Any, 'next_page_url': None, 'prev_page_url': None, 'first_page_url': Any}),
        id='upload-w-child-entries'),
    pytest.param(
        dict(
            user=None,
            expected_status_code=401),
        id='no-credentials'),
    pytest.param(
        dict(
            user='invalid',
            expected_status_code=401),
        id='invalid-credentials'),
    pytest.param(
        dict(
            user='other_test_user',
            expected_status_code=401),
        id='no-access'),
    pytest.param(
        dict(
            user='admin_user',
            expected_data_len=1),
        id='admin-access'),
    pytest.param(
        dict(
            upload_id='silly_value',
            expected_status_code=404),
        id='invalid-upload_id'),
    pytest.param(
        dict(
            upload_id='id_published',
            query_args={'page_size': 5},
            expected_data_len=5,
            expected_response={'processing_successful': 23, 'processing_failed': 0},
            expected_pagination={
                'total': 23, 'page': 1, 'page_after_value': None, 'next_page_after_value': '4', 'order_by': 'mainfile',
                'page_url': Any, 'next_page_url': Any, 'prev_page_url': None, 'first_page_url': Any}),
        id='pag-page-1'),
    pytest.param(
        dict(
            upload_id='id_published',
            query_args={'page_size': 5, 'page': 1},
            expected_data_len=5,
            expected_response={'processing_successful': 23, 'processing_failed': 0},
            expected_pagination={
                'total': 23, 'page': 1, 'page_after_value': None, 'next_page_after_value': '4', 'order_by': 'mainfile',
                'page_url': Any, 'next_page_url': Any, 'prev_page_url': None, 'first_page_url': Any}),
        id='pag-page-1-by-page'),
    pytest.param(
        dict(
            upload_id='id_published',
            query_args={'page_size': 10, 'page': 3},
            expected_data_len=3,
            expected_response={'processing_successful': 23, 'processing_failed': 0},
            expected_pagination={
                'total': 23, 'page': 3, 'page_after_value': '19', 'next_page_after_value': None, 'order_by': 'mainfile',
                'page_url': Any, 'next_page_url': None, 'prev_page_url': Any, 'first_page_url': Any}),
        id='pag-page-3-by-page'),
    pytest.param(
        dict(
            upload_id='id_published',
            query_args={'page_size': 10, 'page_after_value': '19'},
            expected_data_len=3,
            expected_response={'processing_successful': 23, 'processing_failed': 0},
            expected_pagination={
                'total': 23, 'page': 3, 'page_after_value': '19', 'next_page_after_value': None, 'order_by': 'mainfile',
                'page_url': Any, 'next_page_url': None, 'prev_page_url': Any, 'first_page_url': Any}),
        id='pag-page-3-by-page_after_value'),
    pytest.param(
        dict(
            upload_id='id_published',
            query_args={'page_size': 0},
            expected_data_len=0,
            expected_response={'processing_successful': 23, 'processing_failed': 0},
            expected_pagination={
                'total': 23, 'page': 1, 'page_after_value': None, 'next_page_after_value': None, 'order_by': 'mainfile',
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
            upload_id='id_published',
            query_args={'page_size': 1, 'order_by': 'parser_name'},
            expected_data_len=1,
            expected_response={'processing_successful': 23, 'processing_failed': 0},
            expected_pagination={
                'total': 23, 'page': 1, 'page_after_value': None, 'next_page_after_value': '0', 'order_by': 'parser_name',
                'page_url': Any, 'next_page_url': Any, 'prev_page_url': None, 'first_page_url': Any}),
        id='pag-order_by-parser_name'),
    pytest.param(
        dict(
            query_args={'page_size': 1, 'order_by': 'entry_id'},
            expected_status_code=422),
        id='pag-order_by-illegal'),
    pytest.param(
        dict(
            query_args={'page_size': 1, 'page': 2, 'page_after_value': '0'},
            expected_status_code=422),
        id='pag-overspecified')])
def test_get_upload_entries(
        client, mongo_module, test_auth_dict, example_data,
        kwargs):
    '''
    Fetches the entries for a specific upload, by calling uploads/{upload_id}/entries,
    with the provided query paramters, and checks the result.
    '''
    upload_id = kwargs.get('upload_id', 'id_embargo')
    user = kwargs.get('user', 'test_user')
    query_args = kwargs.get('query_args', {})
    expected_status_code = kwargs.get('expected_status_code', 200)
    expected_data_len = kwargs.get('expected_data_len', 1)
    expected_response = kwargs.get('expected_response', {})
    expected_pagination = kwargs.get('expected_pagination', {})
    user_auth, __token = test_auth_dict[user]

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
    pytest.param('id_embargo', 'id_embargo_1', 'test_user', 200, id='ok'),
    pytest.param('id_child_entries', 'id_child_entries_child1', 'test_user', 200, id='child-entry'),
    pytest.param('id_embargo', 'id_embargo_1', None, 401, id='no-credentials'),
    pytest.param('id_embargo', 'id_embargo_1', 'invalid', 401, id='invalid-credentials'),
    pytest.param('id_embargo', 'id_embargo_1', 'other_test_user', 401, id='no-access'),
    pytest.param('id_embargo', 'id_embargo_1', 'admin_user', 200, id='admin-access'),
    pytest.param('silly_value', 'id_embargo_1', 'test_user', 404, id='invalid-upload_id'),
    pytest.param('id_embargo', 'silly_value', 'test_user', 404, id='invalid-entry_id')])
def test_get_upload_entry(
        client, mongo_module, test_auth_dict, example_data,
        upload_id, entry_id, user, expected_status_code):
    '''
    Fetches an entry via a call to uploads/{upload_id}/entries/{entry_id} and checks it.
    '''
    user_auth, __token = test_auth_dict[user]
    response = perform_get(client, f'uploads/{upload_id}/entries/{entry_id}', user_auth)
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        response_json = response.json()
        assert response_json['entry_id'] == entry_id
        response_data = response_json['data']
        assert_entry(response_data)


@pytest.mark.parametrize('args, expected_status_code, expected_mime_type, expected_content', [
    pytest.param(dict(
        user='test_user', upload_id='id_unpublished', path='test_content/id_unpublished_1/1.aux'),
        200, 'text/plain; charset=utf-8', 'content', id='unpublished-file'),
    pytest.param(dict(
        user='test_user', upload_id='id_unpublished', path='test_content/id_unpublished_1/1.aux',
        ignore_mime_type=True),
        200, 'application/octet-stream', 'content', id='unpublished-file-ignore_mime_type'),
    pytest.param(dict(
        user='other_test_user', upload_id='id_unpublished', path='test_content/id_unpublished_1/1.aux'),
        401, None, None, id='unpublished-file-unauthorized'),
    pytest.param(dict(
        user='admin_user', upload_id='id_unpublished', path='test_content/id_unpublished_1/1.aux'),
        200, 'text/plain; charset=utf-8', 'content', id='unpublished-file-admin-auth'),
    pytest.param(dict(
        user='test_user', upload_id='id_published', path='test_content/subdir/test_entry_01/mainfile.json'),
        200, 'text/plain; charset=utf-8', 'method', id='published-file'),
    pytest.param(dict(
        user='test_user', upload_id='id_published', path='test_content/subdir/test_entry_01/mainfile.json',
        ignore_mime_type=True),
        200, 'application/octet-stream', 'method', id='published-file-ignore_mime_type'),
    pytest.param(dict(
        user='admin_user', upload_id='id_published', path='test_content/subdir/test_entry_01/1.aux'),
        200, 'text/plain; charset=utf-8', 'content', id='published-file-admin-auth'),
    pytest.param(dict(
        user='test_user', upload_id='id_unpublished', path='test_content/id_unpublished_1/1.aux',
        compress=True),
        200, 'application/zip', 'content',
        id='unpublished-file-compressed'),
    pytest.param(dict(
        user='test_user', upload_id='id_unpublished', path='test_content/id_unpublished_1/',
        compress=True),
        200, 'application/zip', ['1.aux', '2.aux', '3.aux', '4.aux', 'mainfile.json'],
        id='unpublished-dir-compressed'),
    pytest.param(dict(
        user='test_user', upload_id='id_unpublished', path='', compress=True),
        200, 'application/zip',
        ['test_content', 'test_content/id_unpublished_1/1.aux', 'test_content/id_unpublished_1/mainfile.json'],
        id='unpublished-dir-compressed-root'),
    pytest.param(dict(
        user='test_user', upload_id='id_published', path='test_content/subdir/test_entry_01/1.aux',
        compress=True),
        200, 'application/zip', 'content', id='published-file-compressed'),
    pytest.param(dict(
        user='test_user', upload_id='id_published', path='test_content/subdir/test_entry_01',
        compress=True),
        200, 'application/zip', ['1.aux', '2.aux', '3.aux', '4.aux', 'mainfile.json'],
        id='published-dir-compressed'),
    pytest.param(dict(
        user='test_user', upload_id='id_published', path='', compress=True),
        200, 'application/zip', ['test_content', 'test_content/subdir/test_entry_01/1.aux'],
        id='published-dir-compressed-root'),
    pytest.param(dict(
        user='test_user', upload_id='silly_value', path='test_content/subdir/test_entry_01/1.aux',
        compress=True),
        404, None, None, id='bad-upload-id'),
    pytest.param(dict(
        user='test_user', upload_id='id_published', path='test_content/silly_name', compress=True),
        404, None, None, id='bad-path'),
    pytest.param(dict(
        user='test_user', upload_id='id_unpublished', path='test_content/id_unpublished_1/1.aux',
        offset=2),
        200, 'application/octet-stream', 'ntent\n', id='unpublished-file-offset'),
    pytest.param(dict(
        user='test_user', upload_id='id_unpublished', path='test_content/id_unpublished_1/1.aux',
        offset=2, length=4),
        200, 'application/octet-stream', 'nten', id='unpublished-file-offset-and-length'),
    pytest.param(dict(
        user='test_user', upload_id='id_published', path='test_content/subdir/test_entry_01/1.aux',
        offset=2),
        200, 'application/octet-stream', 'ntent\n', id='published-file-offset'),
    pytest.param(dict(
        user='test_user', upload_id='id_published', path='test_content/subdir/test_entry_01/1.aux',
        offset=2, length=4),
        200, 'application/octet-stream', 'nten', id='published-file-offset-and-length'),
    pytest.param(dict(
        user='test_user', upload_id='id_published', path='test_content/subdir/test_entry_01/1.aux',
        offset=-3),
        400, None, None, id='invalid-offset'),
    pytest.param(dict(
        user='test_user', upload_id='id_published', path='test_content/subdir/test_entry_01/1.aux',
        offset=3, length=-3),
        400, None, None, id='invalid-length'),
    pytest.param(dict(
        user=None, upload_id='id_unpublished', path='test_content/id_unpublished_1/1.aux'),
        401, None, None, id='no-credentials'),
    pytest.param(dict(
        user='invalid', upload_id='id_unpublished', path='test_content/id_unpublished_1/1.aux'),
        401, None, None, id='invalid-credentials'),
    pytest.param(dict(
        user='other_test_user', upload_id='id_unpublished', path='test_content/id_unpublished_1/1.aux'),
        401, None, None, id='no-access'),
    pytest.param(dict(
        user='admin_user', upload_id='id_unpublished', path='test_content/id_unpublished_1/1.aux'),
        200, 'text/plain; charset=utf-8', 'content', id='admin-access')])
def test_get_upload_raw_path(
        client, example_data, test_auth_dict,
        args, expected_status_code, expected_mime_type, expected_content):
    user = args['user']
    upload_id = args['upload_id']
    path = args['path']
    accept = args.get('accept', None)
    compress = args.get('compress', None)
    re_pattern = args.get('re_pattern', None)
    offset = args.get('offset', None)
    length = args.get('length', None)
    ignore_mime_type = args.get('ignore_mime_type', None)
    user_auth, __token = test_auth_dict[user]
    query_args = dict(
        ignore_mime_type=ignore_mime_type,
        compress=compress,
        re_pattern=re_pattern,
        offset=offset,
        length=length)

    response = perform_get(
        client, f'uploads/{upload_id}/raw/{path}', user_auth=user_auth, accept=accept,
        **query_args)

    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        mime_type = response.headers.get('Content-Type')
        if not path:
            expected_filename = upload_id + '.zip'
        else:
            expected_filename = os.path.basename(path.rstrip('/')) + ('.zip' if mime_type == 'application/zip' else '')
        assert_browser_download_headers(response, expected_mime_type, expected_filename)
        if mime_type == 'application/zip':
            if expected_content:
                with zipfile.ZipFile(io.BytesIO(response.content)) as zip_file:
                    if type(expected_content) == str:
                        # Single file - check content
                        with zip_file.open(os.path.basename(path), 'r') as f:
                            file_content = f.read()
                            assert expected_content.encode() in file_content
                    else:
                        assert type(expected_content) == list
                        # Directory - check content
                        zip_paths = zip_file.namelist()
                        # Check: only root elements specified in expected_content are allowed
                        for zip_path in zip_paths:
                            first_path_element = zip_path.split(os.path.sep)[0]
                            assert first_path_element in expected_content, f'Unexpected entry found in the zip root folder: {first_path_element}'
                        # Check: all elements specified in expected_content must exist
                        for expected_path in expected_content:
                            found = False
                            for zip_path in zip_paths:
                                if zip_path == expected_path or zip_path.startswith(expected_path + os.path.sep):
                                    found = True
                                    break
                            assert found, f'Missing expected path in zip file: {expected_path}'
        else:
            if expected_content:
                if offset is not None:
                    assert response.text == expected_content, 'Wrong content (offset and length)'
                else:
                    assert expected_content in response.text, 'Expected content not found'


@pytest.mark.parametrize('user, upload_id, path, query_args, expected_status_code, expected_content, expected_file_metadata, expected_pagination', [
    pytest.param(
        'test_user', 'id_published', 'test_content/subdir/silly_value', {},
        404, None, None, None,
        id='bad-path'),
    pytest.param(
        'test_user', 'id_published', 'test_content/subdir/test_entry_01', {},
        200, ['1.aux', '2.aux', '3.aux', '4.aux', 'mainfile.json'], None, {'total': 5},
        id='published-dir'),
    pytest.param(
        'test_user', 'id_published', 'test_content/subdir/test_entry_01', {'include_entry_info': True},
        200, ['1.aux', '2.aux', '3.aux', '4.aux', 'mainfile.json'], None, {'total': 5},
        id='published-dir-include_entry_info'),
    pytest.param(
        'test_user', 'id_published', 'test_content/subdir/test_entry_01', {'include_entry_info': True, 'page_size': 2, 'page': 3},
        200, ['mainfile.json'], None, {'total': 5},
        id='published-dir-include_entry_info-page3'),
    pytest.param(
        'test_user', 'id_published', '', {},
        200, ['test_content'], None, {'total': 1},
        id='published-dir-root'),
    pytest.param(
        'test_user', 'id_unpublished', 'test_content/id_unpublished_1/', {},
        200, ['1.aux', '2.aux', '3.aux', '4.aux', 'mainfile.json'], None, {'total': 5},
        id='unpublished-dir'),
    pytest.param(
        'test_user', 'id_unpublished', 'test_content/id_unpublished_1/', {'page_size': 3, 'page': 1},
        200, ['1.aux', '2.aux', '3.aux'], None, {'total': 5, 'next_page_after_value': '2'},
        id='unpublished-dir-page1'),
    pytest.param(
        'test_user', 'id_unpublished', 'test_content/id_unpublished_1/', {'page_size': 2, 'page': 4},
        400, None, None, None,
        id='unpublished-dir-page-out-of-range'),
    pytest.param(
        'test_user', 'id_unpublished', 'test_content/id_unpublished_1/', {'include_entry_info': True},
        200, ['1.aux', '2.aux', '3.aux', '4.aux', 'mainfile.json'], None, {'total': 5},
        id='unpublished-dir-include_entry_info'),
    pytest.param(
        'test_user', 'id_child_entries', 'test_content', {'include_entry_info': True},
        200, ['1.aux', '2.aux', '3.aux', '4.aux', 'mainfile_w_children.json'], None, {'total': 5},
        id='dir-child-entries-include_entry_info'),
    pytest.param(
        'test_user', 'id_unpublished', '', {},
        200, ['test_content'], None, {'total': 1},
        id='unpublished-dir-root'),
    pytest.param(
        'test_user', 'id_unpublished', 'test_content/id_unpublished_1/2.aux', {'include_entry_info': True},
        200, None, {'name': '2.aux', 'size': 8, 'entry_id': None, 'parser_name': None}, None,
        id='unpublished-aux-file'),
    pytest.param(
        'test_user', 'id_published', 'test_content/subdir/test_entry_01/mainfile.json', {'include_entry_info': True},
        200, None, {'name': 'mainfile.json', 'size': 3237, 'entry_id': 'id_01', 'parser_name': 'parsers/vasp'}, None,
        id='published-main-file'),
    pytest.param(
        'other_test_user', 'id_unpublished', 'test_content/id_unpublished_1', {},
        401, None, None, None,
        id='unpublished-no-access'),
    pytest.param(
        'other_test_user', 'id_embargo', 'test_content/id_embargo_1', {},
        401, None, None, None,
        id='embargoed-no-access'),
    pytest.param(
        'other_test_user', 'id_embargo_w_coauthor', 'test_content/id_embargo_w_coauthor_1', {},
        200, ['1.aux', '2.aux', '3.aux', '4.aux', 'mainfile.json'], None, {'total': 5},
        id='embargoed-coauthor-access')])
def test_get_upload_rawdir_path(
        client, example_data, test_auth_dict,
        user, upload_id, path, query_args,
        expected_status_code, expected_content, expected_file_metadata, expected_pagination):
    user_auth, __token = test_auth_dict[user]

    response = perform_get(
        client, f'uploads/{upload_id}/rawdir/{path}', user_auth=user_auth, **query_args)

    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        data = response.json()
        assert data['path'] == (path.rstrip('/') or '')
        if expected_content is not None:
            dir_content_returned = data['directory_metadata']['content']
            assert [d['name'] for d in dir_content_returned] == expected_content, 'Incorrect list of files returned'
            for d in dir_content_returned:
                if query_args.get('include_entry_info'):
                    assert (d.get('entry_id') is not None) == ('mainfile' in d['name'])
                    assert (d.get('parser_name') is not None) == ('mainfile' in d['name'])
                else:
                    assert 'entry_id' not in d and 'parser_name' not in d
        elif expected_file_metadata is not None:
            file_metadata_returned = data['file_metadata']
            for k, v in expected_file_metadata.items():
                if v is None:
                    assert k not in file_metadata_returned
                else:
                    assert file_metadata_returned.get(k) == v
        if expected_pagination is None:
            assert 'pagination' not in data
        else:
            pagination_returned = data['pagination']
            for k, v in expected_pagination.items():
                if v is None:
                    assert k not in pagination_returned
                else:
                    assert pagination_returned.get(k) == v


@pytest.mark.parametrize('upload_id, mainfile, user, status_code', [
    pytest.param('id_published', 'test_content/subdir/test_entry_01/mainfile.json', None, 200, id='published'),
    pytest.param('id_published', 'test_content/doesnotexist.json', None, 404, id='bad-mainfile'),
    pytest.param('id_doesnotexist', 'test_content/subdir/test_entry_01/mainfile.json', None, 404, id='bad-upload-id'),
    pytest.param('id_unpublished', 'test_content/id_unpublished_1/mainfile.json', None, 401, id='unpublished'),
    pytest.param('id_unpublished', 'test_content/id_unpublished_1/mainfile.json', 'test_user', 200, id='auth'),
    pytest.param('id_child_entries', 'test_content/mainfile_w_children.json', 'test_user', 200, id='entry-w-child-entries')
])
def test_get_upload_entry_archive_mainfile(
    client, example_data, test_auth_dict,
    upload_id: str, mainfile: str, user: str, status_code: int
):
    user_auth, _ = test_auth_dict[user]
    response = client.get(f'uploads/{upload_id}/archive/mainfile/{mainfile}', headers=user_auth)
    assert_response(response, status_code)
    if status_code == 200:
        assert_archive_response(response.json())


@pytest.mark.parametrize('upload_id, entry_id, user, status_code', [
    pytest.param('id_published', 'id_01', None, 200, id='published'),
    pytest.param('id_published', 'doesnotexist', None, 404, id='bad-entry-id'),
    pytest.param('id_doesnotexist', 'id_01', None, 404, id='bad-upload-id'),
    pytest.param('id_unpublished', 'id_unpublished_1', None, 401, id='unpublished'),
    pytest.param('id_unpublished', 'id_unpublished_1', 'test_user', 200, id='auth'),
    pytest.param('id_child_entries', 'id_child_entries_child1', 'test_user', 200, id='child-entry')
])
def test_get_upload_entry_archive(
    client, example_data, test_auth_dict,
    upload_id: str, entry_id: str, user: str, status_code: int
):
    user_auth, _ = test_auth_dict[user]
    response = client.get(f'uploads/{upload_id}/archive/{entry_id}', headers=user_auth)
    assert_response(response, status_code)
    if status_code == 200:
        assert_archive_response(response.json())


@pytest.mark.parametrize('mode, user, upload_id, source_paths, target_path, query_args, accept_json, use_upload_token, expected_status_code, expected_mainfiles', [
    pytest.param(
        'stream', None, 'examples_template', example_file_aux, '', {'file_name': 'blah.aux'},
        True, False, 401, None, id='no-credentials'),
    pytest.param(
        'stream', None, 'examples_template', example_file_aux, '', {'file_name': 1},
        True, False, 401, None, id='filename-not-str'),
    pytest.param(
        'stream', 'invalid', 'examples_template', example_file_aux, '', {'file_name': 'blah.aux'},
        True, False, 401, None, id='invalid-credentials'),
    pytest.param(
        'stream', 'invalid', 'examples_template', example_file_aux, '', {'file_name': 'blah.aux'},
        True, True, 401, None, id='invalid-credentials-token'),
    pytest.param(
        'multipart', 'admin_user', 'id_published_w', example_file_aux, '', {},
        True, False, 401, None, id='published'),
    pytest.param(
        'multipart', 'admin_user', 'id_processing_w', example_file_aux, '', {},
        True, False, 400, None, id='processing'),
    pytest.param(
        'multipart', 'other_test_user', 'silly_value', example_file_aux, '', {},
        True, False, 404, None, id='bad-upload_id'),
    pytest.param(
        'multipart', 'other_test_user', 'examples_template', example_file_aux, '', {},
        True, False, 401, None, id='no-access-to-upload'),
    pytest.param(
        'multipart', 'test_user', 'examples_template', None, '', {},
        True, False, 400, None, id='no-file'),
    pytest.param(
        'local_path', 'test_user', 'examples_template', example_file_aux, '', {},
        True, False, 401, None, id='local_path-not-admin'),
    pytest.param(
        'stream', 'test_user', 'examples_template', example_file_aux, '', {},
        True, False, 400, None, id='stream-no-file_name'),
    pytest.param(
        'stream', 'test_user', 'id_unpublished_w', example_file_aux, 'test_content/test_embargo_entry',
        {'file_name': 'mainfile.json', 'overwrite_if_exists': False},
        True, False, 409, None, id='cannot-overwrite-existing'),
    pytest.param(
        'stream', 'test_user', 'id_unpublished_w', None, 'test_content/test_embargo_entry',
        {
            'file_name': '2.aux',
            'copy_or_move_source_path': 'test_content/test_embargo_entry/1.aux',
            'copy_or_move': 'copy'
        },
        True, False, 409, None, id='copy-file-to-rawdir-already-exists'),
    pytest.param(
        'multipart', 'test_user', 'examples_template', example_file_aux, '', {},
        True, False, 200, ['examples_template/template.json'], id='multipart'),
    pytest.param(
        'stream', 'test_user', 'examples_template', None, '',
        {
            'file_name': 'template.json',
            'copy_or_move_source_path': 'examples_template/template.json',
            'copy_or_move': 'copy'
        },
        True, False, 200, {'template.json': True, 'examples_template/template.json': True}, id='copy-file-to-rawdir'),
    pytest.param(
        'stream', 'test_user', 'examples_template', None, '',
        {
            'file_name': 'template_2.json',
            'copy_or_move_source_path': 'examples_template/template.json',
            'copy_or_move': 'copy'
        },
        True, False, 200, {'examples_template/template.json': True}, id='copy-with-rename-file-to-rawdir'),
    pytest.param(
        'stream', 'test_user', 'examples_template', None, '',
        {
            'file_name': 'template.json',
            'copy_or_move_source_path': 'examples_template/template.json',
            'copy_or_move': 'move'
        },
        True, False, 200, {'template.json': True}, id='move-file-to-rawdir'),
    pytest.param(
        'stream', 'test_user', 'examples_template', None, '',
        {
            'file_name': 'template_2.json',
            'copy_or_move_source_path': 'examples_template/template.json',
            'copy_or_move': 'move'
        },
        True, False, 200, None, id='move-with-rename-file-to-rawdir'),
    pytest.param(
        'stream', 'test_user', 'examples_template', example_file_aux, '', {'file_name': 'blah.aux'},
        True, False, 200, ['examples_template/template.json'], id='stream'),
    pytest.param(
        'local_path', 'admin_user', 'examples_template', example_file_aux, '', {},
        True, False, 200, ['examples_template/template.json'], id='local_path'),
    pytest.param(
        'stream', 'test_user', 'examples_template', example_file_aux, '', {'file_name': 'blah.aux'},
        True, True, 200, ['examples_template/template.json'], id='token-auth'),
    pytest.param(
        'stream', 'test_user', 'examples_template', example_file_aux, 'dir1/dir2/dir3', {'file_name': 'blah.aux'},
        True, False, 200, ['examples_template/template.json'], id='file-to-subfolder'),
    pytest.param(
        'stream', 'test_user', 'examples_template', example_file_vasp_with_binary, 'dir1/dir2', {'file_name': 'tmp.zip'},
        True, False, 200, [
            'examples_template/template.json',
            'dir1/dir2/examples_vasp/xml/Si.xml',
            'dir1/dir2/examples_vasp/xml/perovskite.xml.gz'], id='zip-to-subfolder'),
    pytest.param(
        'stream', 'test_user', 'examples_template', example_file_aux, 'examples_template', {'file_name': 'template.json'},
        True, False, 200, {'examples_template/template.json': False}, id='overwrite-and-destroy-old-mainfile'),
    pytest.param(
        'stream', 'test_user', 'examples_template', example_file_vasp_with_binary, '', {'file_name': 'tmp.zip'},
        True, False, 200, [
            'examples_template/template.json',
            'examples_vasp/xml/Si.xml',
            'examples_vasp/xml/perovskite.xml.gz'], id='unzip-and-add-new-mainfiles'),
    pytest.param(
        'stream', 'test_user', 'examples_template', example_file_corrupt_zip, '', {'file_name': 'tmp.zip'},
        True, False, 400, None, id='bad-zip'),
    pytest.param(
        'multipart', 'test_user', 'examples_template', example_file_aux, 'examples_template',
        {'wait_for_processing': True},
        True, False, 200, ['examples_template/template.json'], id='wait_for_processing-auxfile-add'),
    pytest.param(
        'multipart', 'test_user', 'examples_template', example_file_mainfile_different_atoms, 'dir1/dir2',
        {'wait_for_processing': True},
        True, False, 200, [
            'examples_template/template.json',
            'dir1/dir2/template.json'], id='wait_for_processing-mainfile-add'),
    pytest.param(
        'multipart', 'test_user', 'examples_template', example_file_mainfile_different_atoms, 'dir1/dir2',
        {'wait_for_processing': True, 'include_archive': True},
        True, False, 200, [
            'examples_template/template.json',
            'dir1/dir2/template.json'], id='wait_for_processing-mainfile-add-include_archive'),
    pytest.param(
        'multipart', 'test_user', 'examples_template', example_file_mainfile_different_atoms, 'examples_template',
        {'wait_for_processing': True},
        True, False, 200, ['examples_template/template.json'], id='wait_for_processing-mainfile-overwrite'),
    pytest.param(
        'multipart', 'test_user', 'examples_template', example_file_unparsable, 'examples_template',
        {'wait_for_processing': True},
        True, False, 200, {'examples_template/template.json': False}, id='wait_for_processing-mainfile-overwrite-destroy'),
    pytest.param(
        'multipart', 'test_user', 'examples_template', example_file_vasp_with_binary, 'examples_template',
        {'wait_for_processing': True},
        True, False, 400, None, id='wait_for_processing-zipfile'),
    pytest.param(
        'multipart', 'test_user', 'examples_template', [example_file_vasp_with_binary, example_file_aux], 'dir1', {'file_name': 'tmp.zip'},
        True, False, 200, [
            'examples_template/template.json',
            'dir1/examples_vasp/xml/Si.xml',
            'dir1/examples_vasp/xml/perovskite.xml.gz'], id='upload-multiple-vasp-and-aux'),
    pytest.param(
        'multipart', 'test_user', 'examples_template', [example_file_aux, example_file_corrupt_zip], 'dir1', {'file_name': 'tmp.zip'},
        True, False, 400, ['examples_template/template.json'], id='upload-multiple-one-corrupted-zip')])
def test_put_upload_raw_path(
        client, proc_infra, non_empty_processed, example_data_writeable, test_auth_dict,
        mode, user, upload_id, source_paths, target_path, query_args, accept_json, use_upload_token,
        expected_status_code, expected_mainfiles):
    action = 'PUT'
    url = f'uploads/{upload_id}/raw/{target_path}'
    published = False
    all_entries_should_succeed = not (type(expected_mainfiles) == dict and False in expected_mainfiles.values())
    expected_process_status = ProcessStatus.SUCCESS if 'wait_for_processing' in query_args else None

    response, _ = assert_file_upload_and_processing(
        client, action, url, mode, user, test_auth_dict, upload_id,
        source_paths, target_path, query_args, accept_json, use_upload_token,
        expected_status_code, expected_process_status, expected_mainfiles, published,
        all_entries_should_succeed)

    if response.status_code == 200 and accept_json:
        response_json = response.json()
        processing = response_json['processing']
        if 'wait_for_processing' in query_args:
            assert processing
            assert processing['upload_id'] == upload_id
            assert processing['path'] == os.path.join(target_path, os.path.basename(source_paths))
            if source_paths == example_file_aux:
                # Not a mainfile
                for k in ('entry_id', 'parser_name', 'entry', 'archive'):
                    assert processing[k] is None
            else:
                # Mainfile was added
                if source_paths == example_file_unparsable:
                    expected_entry_process_status = ProcessStatus.FAILURE
                else:
                    expected_entry_process_status = ProcessStatus.SUCCESS
                assert processing['entry_id'] is not None
                assert processing['parser_name'] is not None
                assert processing['entry']['process_status'] == expected_entry_process_status
                assert (processing['archive'] is None) == (not query_args.get('include_archive'))
        else:
            assert not processing


@pytest.mark.parametrize('mode, user, expected_status_code', [
    pytest.param(
        'multipart', 'test_user', 409, id='conflict_in_concurrent_editing')])
def test_editing_raw_file(
        client, proc_infra, non_empty_processed, example_data_writeable, test_auth_dict,
        mode, user, expected_status_code):
    upload_id = 'examples_template'
    target_path = 'examples_template'
    action = 'PUT'
    url = f'uploads/{upload_id}/raw/{target_path}'
    path = 'examples_template/template.json'
    user_auth, __token = test_auth_dict[user]

    # Get an existing upload with entries
    response = perform_get(client, f'uploads/{upload_id}/archive/mainfile/{path}', user_auth=user_auth)
    assert response.status_code == 200
    response_json = response.json()
    entry_hash = response_json['data']['archive']['metadata']['entry_hash']

    # First edit
    query_args = {'file_name': 'example.json', 'wait_for_processing': True,
                  'include_archive': True, 'entry_hash': entry_hash}
    response, _ = assert_file_upload_and_processing(
        client, action, url, mode, user, test_auth_dict, upload_id,
        example_file_mainfile_different_atoms, target_path, query_args, True, False,
        200, ProcessStatus.SUCCESS, ['examples_template/template.json'], False, True)

    # Second edit on the same entry
    query_args = {'file_name': 'example.json', 'wait_for_processing': True, 'entry_hash': entry_hash}
    response, _ = assert_file_upload_and_processing(
        client, action, url, mode, user, test_auth_dict, upload_id,
        example_file_mainfile_different_atoms, target_path, query_args, True, False,
        expected_status_code, None, None, None, None)

    # Get the updated archive
    response = perform_get(client, f'uploads/{upload_id}/archive/mainfile/{path}', user_auth=user_auth)
    assert response.status_code == 200
    response_json = response.json()
    entry_hash = response_json['data']['archive']['metadata']['entry_hash']

    # Edit with the updated hash code
    query_args = {'file_name': 'example.json', 'wait_for_processing': True,
                  'include_archive': True, 'entry_hash': entry_hash}
    response, _ = assert_file_upload_and_processing(
        client, action, url, mode, user, test_auth_dict, upload_id,
        example_file_mainfile_different_atoms, target_path, query_args, True, False,
        200, ProcessStatus.SUCCESS, ['examples_template/template.json'], False, True)

    # Get the updated archive
    response = perform_get(client, f'uploads/{upload_id}/archive/mainfile/{path}', user_auth=user_auth)
    assert response.status_code == 200
    response_json = response.json()
    entry_hash = response_json['data']['archive']['metadata']['entry_hash']

    # Somebody deletes the file faster
    response = client.delete(f'uploads/{upload_id}/raw/{path}', headers=user_auth)
    assert response.status_code == 200
    time.sleep(1)

    # Editing the file which was deleted by someone else
    query_args = {'file_name': 'example.json', 'wait_for_processing': True, 'entry_hash': entry_hash}
    response, _ = assert_file_upload_and_processing(
        client, action, url, mode, user, test_auth_dict, upload_id,
        example_file_mainfile_different_atoms, target_path, query_args, True, False,
        expected_status_code, None, None, None, None)


@pytest.mark.parametrize('user, upload_id, path, expected_status_code', [
    pytest.param(
        'test_user', 'id_published_w', 'test_content/newdir', 401, id='published'),
    pytest.param(
        None, 'id_unpublished_w', 'test_content/newdir', 401, id='no-credentials'),
    pytest.param(
        'other_test_user', 'id_unpublished_w', 'test_content/newdir', 401, id='no-access'),
    pytest.param(
        'admin_user', 'id_unpublished_w', 'test_content/newdir', 200, id='admin-access'),
    pytest.param(
        'test_user', 'id_unpublished_w', 'test_content/test_embargo_entry/newdir', 200, id='ok'),
    pytest.param(
        'test_user', 'id_unpublished_w', 'test_content/chars?! "\'@#$%&\\()[]{}=+`´^~*,.;:|<>', 200, id='special-chars'),
    pytest.param(
        'test_user', 'id_unpublished_w', 'test_content/test_embargo_entry/mainfile.json/newdir', 400, id='bad-path')])
def test_post_upload_raw_create_dir_path(
        client, proc_infra, example_data_writeable, test_auth_dict,
        user, upload_id, path, expected_status_code):
    user_auth, _token = test_auth_dict[user]
    response = client.post(f'uploads/{upload_id}/raw-create-dir/{requests.utils.quote(path)}', headers=user_auth)
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        upload = Upload.get(upload_id)
        assert upload.upload_files.raw_path_exists(path)
        assert not upload.upload_files.raw_path_is_file(path)


@pytest.mark.parametrize('user, upload_id, path, use_upload_token, expected_status_code, expected_mainfiles', [
    pytest.param(
        'test_user', 'examples_template', 'examples_template/1.aux', False,
        200, ['examples_template/template.json'], id='delete-aux-file'),
    pytest.param(
        'test_user', 'examples_template', 'examples_template/template.json', False,
        200, [], id='delete-main-file'),
    pytest.param(
        'test_user', 'examples_template', '', False,
        200, [], id='delete-root'),
    pytest.param(
        'test_user', 'examples_template', 'examples_template', False,
        200, [], id='delete-subfolder'),
    pytest.param(
        'test_user', 'examples_template', 'examples_template/1.aux', True,
        200, ['examples_template/template.json'], id='delete-token-access'),
    pytest.param(
        'admin_user', 'examples_template', 'examples_template/1.aux', False,
        200, ['examples_template/template.json'], id='delete-admin-access'),
    pytest.param(
        'other_test_user', 'examples_template', 'examples_template/1.aux', False,
        401, None, id='no-access'),
    pytest.param(
        None, 'examples_template', 'examples_template/1.aux', False,
        401, None, id='no-credentials'),
    pytest.param(
        'invalid', 'examples_template', 'examples_template/1.aux', False,
        401, None, id='invalid-credentials'),
    pytest.param(
        'invalid', 'examples_template', 'examples_template/1.aux', True,
        401, None, id='invalid-credentials-token'),
    pytest.param(
        'test_user', 'id_published_w', 'examples_template/1.aux', False,
        401, None, id='published'),
    pytest.param(
        'test_user', 'id_processing_w', 'examples_template/1.aux', False,
        400, None, id='processing')])
def test_delete_upload_raw_path(
        client, proc_infra, non_empty_processed, example_data_writeable, test_auth_dict,
        user, upload_id, path, use_upload_token, expected_status_code, expected_mainfiles):
    user_auth, token = test_auth_dict[user]
    # Use either token or bearer token for the post operation (never both)
    user_auth_action = user_auth
    if use_upload_token:
        user_auth_action = None
    else:
        token = None
    if upload_id == 'id_processing_w':
        # Ensure file exists (otherwise we get 404, which is not what we want to test)
        upload_files = StagingUploadFiles(upload_id)
        upload_files.add_rawfiles('tests/data/proc/examples_template/1.aux', 'examples_template')
    query_args = dict(token=token)
    response = client.delete(build_url(f'uploads/{upload_id}/raw/{path}', query_args), headers=user_auth_action)
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        assert_processing(client, upload_id, user_auth)
        # Check that path to remove has disappeared
        upload_files = StagingUploadFiles(upload_id)
        if path == '':
            # Deleting the root folder = the folder itself should be emptied, but not deleted.
            assert not list(upload_files.raw_directory_list(''))
        else:
            # Deleting a file or folder within the raw folder - it should disappear.
            assert not upload_files.raw_path_exists(path)

        assert_expected_mainfiles(upload_id, expected_mainfiles)


@pytest.mark.parametrize('user, upload_id, kwargs', [
    pytest.param(
        'test_user', 'id_unpublished_w', dict(
            metadata=all_coauthor_metadata),
        id='edit-all'),
    pytest.param(
        'test_user', 'id_published_w', dict(
            metadata=dict(embargo_length=0)), id='lift-embargo'),
    pytest.param(
        'admin_user', 'id_published_w', dict(
            metadata=all_admin_metadata),
        id='protected-admin'),
    pytest.param(
        'test_user', 'id_unpublished_w', dict(
            metadata=dict(main_author='lhofstadter'),
            expected_error_loc=('metadata', 'main_author')),
        id='protected-not-admin'),
    pytest.param(
        'test_user', 'silly_value', dict(
            metadata=dict(upload_name='test_name'),
            expected_error_loc=('upload_id',)),
        id='bad-upload_id'),
    pytest.param(
        'admin_user', 'id_published_w', dict(
            metadata=dict(upload_name='test_name')),
        id='published-admin'),
    pytest.param(
        'test_user', 'id_published_w', dict(
            metadata=dict(upload_name='test_name')),
        id='published-not-admin'),
    pytest.param(
        None, 'id_unpublished_w', dict(
            metadata=dict(upload_name='test_name'),
            expected_status_code=401),
        id='no-credentials'),
    pytest.param(
        'invalid', 'id_unpublished_w', dict(
            metadata=dict(upload_name='test_name'),
            expected_status_code=401),
        id='invalid-credentials'),
    pytest.param(
        'other_test_user', 'id_unpublished_w', dict(
            metadata=dict(upload_name='test_name'),
            expected_error_loc=('metadata', 'upload_name')),
        id='no-access'),
    pytest.param(
        'other_test_user', 'id_unpublished_w', dict(
            metadata=dict(upload_name='test_name'),
            add_coauthor=True),
        id='coauthor-access'),
    pytest.param(
        'test_user', 'id_empty_w', dict(
            metadata=dict(upload_name='test_name')),
        id='empty-upload-ok'),
    pytest.param(
        'test_user', 'id_unpublished_w', dict(
            query={'and': [{'upload_create_time:gt': '2021-01-01'}, {'published': False}]},
            owner='user',
            metadata=dict(comment='a test comment')),
        id='query-ok'),
    pytest.param(
        'test_user', 'id_unpublished_w', dict(
            query={'and': [{'upload_create_time:gt': '2021-01-01'}, {'published': False}]},
            owner='user',
            metadata=dict(upload_name='a test name'),
            expected_error_loc=('metadata', 'upload_name')),
        id='query-cannot-edit-upload-data'),
    pytest.param(
        'test_user', 'id_unpublished_w', dict(
            query={'upload_create_time:lt': '2021-01-01'},
            owner='user',
            metadata=dict(comment='a test comment'),
            expected_error_loc=('query',)),
        id='query-no-results')])
def test_post_upload_edit(
        client, proc_infra, example_data_writeable, example_datasets, test_auth_dict, test_users_dict,
        user, upload_id, kwargs):
    '''
    Note, since the endpoint basically just forwards the request to
    `MetadataEditRequestHandler.edit_metadata`, we only do very simple verification here,
    the more extensive testnig is done in `tests.processing.test_edit_metadata`.
    '''
    user_auth, _token = test_auth_dict[user]
    user = test_users_dict.get(user)
    query = kwargs.get('query')
    owner = kwargs.get('owner')
    metadata = kwargs.get('metadata')
    entries = kwargs.get('entries')
    entries_key = kwargs.get('entries_key')
    verify_only = kwargs.get('verify_only', False)
    expected_error_loc = kwargs.get('expected_error_loc')
    expected_status_code = kwargs.get('expected_status_code')
    affected_upload_ids = kwargs.get('affected_upload_ids', [upload_id])
    expected_metadata = kwargs.get('expected_metadata', metadata)

    add_coauthor = kwargs.get('add_coauthor', False)
    if add_coauthor:
        upload = Upload.get(upload_id)
        upload.edit_upload_metadata(
            edit_request_json={'metadata': {'coauthors': user.user_id}}, user_id=upload.main_author)
        upload.block_until_complete()

    edit_request_json = dict(
        query=query, owner=owner, metadata=metadata, entries=entries, entries_key=entries_key,
        verify_only=verify_only)
    url = f'uploads/{upload_id}/edit'
    edit_start = datetime.utcnow().isoformat()[0:22]
    response = client.post(url, headers=user_auth, json=edit_request_json)
    if expected_error_loc:
        assert_response(response, 422)
        error_locs = [tuple(d['loc']) for d in response.json()['detail']]
        assert expected_error_loc in error_locs
    elif expected_status_code not in (None, 200):
        assert_response(response, expected_status_code)
    else:
        assert_response(response, 200)
        assert_metadata_edited(
            user, upload_id, query, metadata, entries, entries_key, verify_only,
            expected_metadata, affected_upload_ids, edit_start)


@pytest.mark.parametrize('mode, source_paths, query_args, user, use_upload_token, test_limit, accept_json, expected_status_code', [
    pytest.param('multipart', example_file_vasp_with_binary, dict(upload_name='test_name'), 'test_user', False, False, True, 200, id='multipart'),
    pytest.param('multipart', example_file_vasp_with_binary, dict(), 'test_user', False, False, True, 200, id='multipart-no-name'),
    pytest.param('multipart', example_file_vasp_with_binary, dict(upload_name='test_name'), 'test_user', True, False, True, 200, id='multipart-token'),
    pytest.param('stream', example_file_vasp_with_binary, dict(embargo_length=0, upload_name='test_name'), 'test_user', False, False, True, 200, id='stream-no-embargo'),
    pytest.param('stream', example_file_vasp_with_binary, dict(embargo_length=7), 'test_user', False, False, True, 200, id='stream-no-name-embargoed'),
    pytest.param('stream', example_file_vasp_with_binary, dict(embargo_length=37), 'test_user', False, False, True, 400, id='stream-invalid-embargo'),
    pytest.param('stream', example_file_vasp_with_binary, dict(upload_name='test_name'), 'test_user', True, False, True, 200, id='stream-token'),
    pytest.param('local_path', example_file_vasp_with_binary, dict(), 'admin_user', False, False, True, 200, id='local_path'),
    pytest.param('local_path', example_file_vasp_with_binary, dict(), 'test_user', False, False, True, 401, id='local_path-not-admin'),
    pytest.param('stream', example_file_vasp_with_binary, dict(), 'test_user', False, False, False, 200, id='no-accept-json'),
    pytest.param('multipart', example_file_vasp_with_binary, dict(), None, False, False, True, 401, id='no-credentials'),
    pytest.param('multipart', example_file_vasp_with_binary, dict(), 'invalid', False, False, True, 401, id='invalid-credentials'),
    pytest.param('multipart', example_file_vasp_with_binary, dict(), 'invalid', True, False, True, 401, id='invalid-credentials-token'),
    pytest.param('stream', [], dict(upload_name='test_name'), 'test_user', False, False, True, 200, id='no-file'),
    pytest.param('stream', example_file_aux, dict(file_name='1.aux'), 'test_user', False, False, True, 200, id='stream-non-zip-file'),
    pytest.param('stream', example_file_aux, dict(), 'test_user', False, False, True, 400, id='stream-non-zip-file-no-file_name'),
    pytest.param('stream', example_file_vasp_with_binary, dict(upload_name='test_name', publish_directly=True), 'test_user', False, False, True, 200, id='publish_directly'),
    pytest.param('stream', empty_file, dict(upload_name='test_name', publish_directly=True), 'test_user', False, False, True, 200, id='publish_directly-empty'),
    pytest.param('stream', example_file_vasp_with_binary, dict(upload_name='test_name'), 'test_user', False, True, True, 400, id='upload-limit-exceeded'),
    pytest.param('multipart', example_file_corrupt_zip, dict(), 'test_user', False, False, True, 200, id='bad-zip'),
    pytest.param('multipart', [example_file_aux, example_file_mainfile_different_atoms], dict(), 'test_user', False, False, True, 200, id='upload-multiple-files'),
    pytest.param('multipart', [example_file_aux, example_file_corrupt_zip], dict(), 'test_user', False, False, True, 200, id='upload-multiple-files-one-corrupt')])
def test_post_upload(
        client, mongo, proc_infra, monkeypatch, test_auth_dict,
        empty_upload, non_empty_example_upload,
        mode, source_paths, query_args, user, use_upload_token, test_limit, accept_json,
        expected_status_code):
    '''
    Posts an upload, with different arguments.
    '''
    if type(source_paths) == str:
        source_paths = [source_paths]
    if test_limit:
        monkeypatch.setattr('nomad.config.services.upload_limit', 0)

    action = 'POST'
    url = 'uploads'
    published = (query_args.get('publish_directly') and not source_paths == [empty_file])
    all_entries_should_succeed = True
    target_path = ''
    expected_mainfiles = None
    upload_id = None  # Not determined yet
    expected_process_status = None

    _, processed_response_data = assert_file_upload_and_processing(
        client, action, url, mode, user, test_auth_dict, upload_id,
        source_paths, target_path, query_args, accept_json, use_upload_token,
        expected_status_code, expected_process_status, expected_mainfiles, published,
        all_entries_should_succeed)

    if expected_status_code == 200 and processed_response_data:
        expected_upload_name = query_args.get('upload_name')
        if not expected_upload_name:
            if mode in ('multipart', 'local_path') and len(source_paths) == 1:
                expected_upload_name = os.path.basename(source_paths[0])
            elif mode == 'stream':
                expected_upload_name = query_args.get('file_name')

        assert processed_response_data.get('upload_name') == expected_upload_name

    if query_args.get('publish_directly'):
        upload_id = processed_response_data['upload_id']
        upload_proc = Upload.objects(upload_id=upload_id).first()
        if source_paths == [empty_file]:
            assert not upload_proc.published
        else:
            assert_gets_published(client, upload_id, test_auth_dict['test_user'][0], **query_args)


@pytest.mark.parametrize('kwargs', [
    pytest.param(
        dict(
            expected_status_code=200),
        id='no-args'),
    pytest.param(
        dict(
            query_args={'embargo_length': 12},
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
            query_args={'embargo_length': 0},
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
            expected_status_code=400),
        id='already-published'),
    pytest.param(
        dict(
            user=None,
            expected_status_code=401),
        id='no-credentials'),
    pytest.param(
        dict(
            user='invalid',
            expected_status_code=401),
        id='invalid-credentials'),
    pytest.param(
        dict(
            user='other_test_user',
            expected_status_code=401),
        id='no-access')])
def test_post_upload_action_publish(
        client, proc_infra, example_data_writeable,
        test_auth_dict, kwargs):
    ''' Tests the publish action with various arguments. '''
    upload_id = kwargs.get('upload_id', 'id_unpublished_w')
    query_args = kwargs.get('query_args', {})
    expected_status_code = kwargs.get('expected_status_code', 200)
    user = kwargs.get('user', 'test_user')
    user_auth, __token = test_auth_dict[user]

    response = perform_post_upload_action(client, user_auth, upload_id, 'publish', **query_args)

    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        upload = assert_upload(response.json())
        assert upload['current_process'] == 'publish_upload'
        assert upload['process_running']

        assert_gets_published(client, upload_id, user_auth, current_embargo_length=12, **query_args)


@pytest.mark.parametrize('import_settings, query_args', [
    pytest.param(
        BundleImportSettings(include_archive_files=False, trigger_processing=True),
        dict(embargo_length=0),
        id='trigger-processing'),
    pytest.param(
        BundleImportSettings(include_archive_files=True, trigger_processing=False),
        dict(embargo_length=28),
        id='no-processing')
])
def test_post_upload_action_publish_to_central_nomad(
        client, proc_infra, monkeypatch, oasis_publishable_upload,
        test_users_dict, test_auth_dict, import_settings, query_args):
    ''' Tests the publish action with to_central_nomad=True. '''
    upload_id, suffix = oasis_publishable_upload
    query_args['to_central_nomad'] = True
    embargo_length = query_args.get('embargo_length')
    expected_status_code = 200
    user = 'test_user'
    user_auth, __token = test_auth_dict[user]
    old_upload = Upload.get(upload_id)

    import_settings = config.bundle_import.default_settings.customize(import_settings)
    monkeypatch.setattr('nomad.config.bundle_import.default_settings', import_settings)
    monkeypatch.setattr('nomad.config.bundle_import.allow_bundles_from_oasis', True)

    # Finally, invoke the method to publish to central nomad
    response = perform_post_upload_action(client, user_auth, upload_id, 'publish', **query_args)

    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        upload = assert_upload(response.json())
        assert upload['current_process'] == 'publish_externally'
        assert upload['process_running']
        assert_processing(client, upload_id, user_auth, published=old_upload.published)
        assert_processing(client, upload_id + suffix, user_auth, published=old_upload.published)

        old_upload = Upload.get(upload_id)
        new_upload = Upload.get(upload_id + suffix)
        assert len(old_upload.successful_entries) == len(new_upload.successful_entries) == 1
        if embargo_length is None:
            embargo_length = old_upload.embargo_length
        old_entry = old_upload.successful_entries[0]
        new_entry = new_upload.successful_entries[0]
        old_entry_metadata_dict = old_entry.full_entry_metadata(old_upload).m_to_dict()
        new_entry_metadata_dict = new_entry.full_entry_metadata(new_upload).m_to_dict()
        for k, v in old_entry_metadata_dict.items():
            if k == 'with_embargo':
                assert new_entry_metadata_dict[k] == (embargo_length > 0)
            elif k not in (
                    'upload_id', 'entry_id', 'upload_create_time', 'entry_create_time',
                    'last_processing_time', 'publish_time', 'embargo_length',
                    'n_quantities', 'quantities'):  # TODO: n_quantities and quantities update problem?
                assert new_entry_metadata_dict[k] == v, f'Metadata not matching: {k}'
        assert new_entry.datasets == ['dataset_id']
        assert old_upload.published_to[0] == config.oasis.central_nomad_deployment_url
        assert new_upload.from_oasis and new_upload.oasis_deployment_url
        assert new_upload.embargo_length == embargo_length
        assert old_upload.upload_files.access == 'restricted' if old_upload.with_embargo else 'public'
        assert new_upload.upload_files.access == 'restricted' if new_upload.with_embargo else 'public'


@pytest.mark.parametrize('upload_id, publish, user, expected_status_code', [
    pytest.param('examples_template', True, 'admin_user', 200, id='published-admin'),
    pytest.param('examples_template', True, 'test_user', 401, id='published-not-admin'),
    pytest.param('examples_template', False, 'test_user', 200, id='not-published'),
    pytest.param('examples_template', False, None, 401, id='no-credentials'),
    pytest.param('examples_template', False, 'invalid', 401, id='invalid-credentials'),
    pytest.param('examples_template', False, 'other_test_user', 401, id='no-access'),
    pytest.param('id_processing_w', False, 'test_user', 400, id='already-processing'),
    pytest.param('silly_value', False, 'test_user', 404, id='invalid-upload_id')])
def test_post_upload_action_process(
        client, mongo, proc_infra, monkeypatch, example_data_writeable,
        non_empty_processed, internal_example_user_metadata, test_auth_dict,
        upload_id, publish, user, expected_status_code):

    if publish:
        set_upload_entry_metadata(non_empty_processed, internal_example_user_metadata)
        non_empty_processed.publish_upload()
        try:
            non_empty_processed.block_until_complete(interval=.01)
        except Exception:
            pass

    monkeypatch.setattr('nomad.config.meta.version', 're_process_test_version')
    monkeypatch.setattr('nomad.config.meta.commit', 're_process_test_commit')
    user_auth, __token = test_auth_dict[user]

    response = perform_post_upload_action(client, user_auth, upload_id, 'process')
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        assert_processing(client, upload_id, test_auth_dict['test_user'][0], check_files=False, published=True)


@pytest.mark.parametrize('upload_id, user, owner, query, include_parent_folders, expected_status_code, expect_exists, expect_not_exists', [
    pytest.param(
        'id_published_w', 'admin_user', None, {'entry_id': 'id_published_w_entry'}, False,
        401, [], [], id='published-admin'),
    pytest.param(
        'id_unpublished_w', 'other_test_user', None, {'entry_id': 'id_unpublished_w_entry'}, False,
        401, [], [], id='unpublished-no-access'),
    pytest.param(
        'id_unpublished_w', 'test_user', None, None, False,
        400, [], [], id='no-query'),
    pytest.param(
        'id_unpublished_w', 'test_user', None, {'entry_id': ['id_unpublished_w_entry', 'silly']}, False,
        200, ['test_content/test_embargo_entry/1.aux'], ['test_content/test_embargo_entry/mainfile.json'], id='ok'),
    pytest.param(
        'id_unpublished_w', 'test_user', None, {'entry_id': 'id_unpublished_w_entry'}, True,
        200, ['test_content'], ['test_content/test_embargo_entry'], id='ok-delete-folder'),
    pytest.param(
        'id_unpublished_w', 'admin_user', 'admin', {'entry_id': 'id_unpublished_w_entry'}, False,
        200, ['test_content/test_embargo_entry/1.aux'], ['test_content/test_embargo_entry/mainfile.json'], id='ok-admin-access')])
def test_post_upload_action_delete_entry_files(
        client, mongo, proc_infra, example_data_writeable, test_auth_dict,
        upload_id, user, owner, query, include_parent_folders,
        expected_status_code, expect_exists, expect_not_exists):
    user_auth, __token = test_auth_dict[user]
    json = {}
    if include_parent_folders is not None:
        json.update(include_parent_folders=include_parent_folders)
    if owner is not None:
        json.update(owner=owner)
    if query is not None:
        json.update(query=query)

    response = perform_post_upload_action(client, user_auth, upload_id, 'delete-entry-files', json=json)
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        upload = Upload.get(upload_id)
        upload.block_until_complete()
        for path in expect_exists or []:
            assert upload.upload_files.raw_path_exists(path), f'Missing expected path: {path}'
        for path in expect_not_exists or []:
            assert not upload.upload_files.raw_path_exists(path), f'Expected path not to exist: {path}'


@pytest.mark.parametrize('upload_id, user, preprocess, expected_status_code', [
    pytest.param('id_published_w', 'test_user', None, 200, id='ok'),
    pytest.param('id_published_w', 'other_test_user', None, 401, id='no-access'),
    pytest.param('id_published_w', 'other_test_user', 'make-coauthor', 200, id='ok-coauthor'),
    pytest.param('id_published_w', None, None, 401, id='no-credentials'),
    pytest.param('id_published_w', 'invalid', None, 401, id='invalid-credentials'),
    pytest.param('id_unpublished_w', 'test_user', None, 400, id='not-published'),
    pytest.param('id_published_w', 'test_user', 'lift', 400, id='already-lifted')])
def test_post_upload_action_lift_embargo(
        client, proc_infra, example_data_writeable, test_auth_dict, test_users_dict,
        upload_id, user, preprocess, expected_status_code):

    user_auth, __token = test_auth_dict[user]
    user = test_users_dict.get(user)

    if preprocess:
        if preprocess == 'lift':
            metadata = {'embargo_length': 0}
        elif preprocess == 'make-coauthor':
            metadata = {'coauthors': user.user_id}
        upload = Upload.get(upload_id)
        upload.edit_upload_metadata(dict(metadata=metadata), config.services.admin_user_id)
        upload.block_until_complete()

    response = perform_post_upload_action(client, user_auth, upload_id, 'lift-embargo')
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        assert_metadata_edited(
            user, upload_id, None, None, None, None, False,
            {'embargo_length': 0}, [upload_id], None)


@pytest.mark.parametrize('upload_id, user, expected_status_code', [
    pytest.param('id_unpublished_w', 'test_user', 200, id='delete-own'),
    pytest.param('id_unpublished_w', 'other_test_user', 401, id='delete-others-not-admin'),
    pytest.param('id_unpublished_w', 'admin_user', 200, id='delete-others-admin'),
    pytest.param('id_published_w', 'test_user', 401, id='delete-own-published'),
    pytest.param('id_published_w', 'admin_user', 200, id='delete-others-published-admin'),
    pytest.param('silly_value', 'test_user', 404, id='invalid-upload_id'),
    pytest.param('id_unpublished_w', None, 401, id='no-credentials'),
    pytest.param('id_unpublished_w', 'invalid', 401, id='invalid-credentials')])
def test_delete_upload(
        client, proc_infra, example_data_writeable, test_auth_dict,
        upload_id, user, expected_status_code):
    ''' Uploads a file, and then tries to delete it, with different parameters and users. '''
    user_auth, __token = test_auth_dict[user]

    response = client.delete(f'uploads/{upload_id}', headers=user_auth)
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        assert_upload_does_not_exist(client, upload_id, test_auth_dict['test_user'][0])


@pytest.mark.parametrize('upload_id, user, query_args, expected_status_code', [
    pytest.param(
        'id_published_w', 'test_user', dict(),
        200, id='published-owner'),
    pytest.param(
        'id_published_w', 'admin_user', dict(),
        200, id='published-admin'),
    pytest.param(
        'id_published_w', 'other_test_user', dict(),
        401, id='published-not-owner'),
    pytest.param(
        'id_published_w', 'test_user', dict(include_raw_files=False),
        200, id='published-owner-exclude-raw'),
    pytest.param(
        'id_published_w', 'test_user', dict(include_archive_files=False),
        200, id='published-owner-exclude-archive'),
    pytest.param(
        'id_unpublished_w', 'test_user', dict(),
        200, id='unpublished-owner'),
    pytest.param(
        'id_unpublished_w', 'admin_user', dict(),
        200, id='unpublished-admin'),
    pytest.param(
        'id_unpublished_w', 'other_test_user', dict(),
        401, id='unpublished-not-owner')])
def test_get_upload_bundle(
        client, proc_infra, example_data_writeable, test_auth_dict,
        upload_id, user, query_args, expected_status_code):

    include_raw_files = query_args.get('include_raw_files', True)
    include_archive_files = query_args.get('include_archive_files', True)

    url = build_url(f'uploads/{upload_id}/bundle', query_args)
    response = perform_get(client, url, user_auth=test_auth_dict[user][0])
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        with zipfile.ZipFile(io.BytesIO(response.content)) as zip_file:
            upload = Upload.get(upload_id)
            upload_files = upload.upload_files
            expected_files = set(['bundle_info.json'])
            for dirpath, __, filenames in os.walk(upload_files.os_path):
                for filename in filenames:
                    os_path = os.path.join(dirpath, filename)
                    rel_path = os.path.relpath(os_path, upload_files.os_path)
                    include = False
                    include |= rel_path.startswith('raw') and include_raw_files
                    include |= rel_path.startswith('archive') and include_archive_files
                    if include:
                        expected_files.add(rel_path)
            assert expected_files == set(zip_file.namelist())
    return


@pytest.mark.parametrize('publish, test_duplicate, user, export_args, query_args, expected_status_code', [
    pytest.param(
        True, False, 'admin_user', dict(), dict(),
        200, id='published-admin'),
    pytest.param(
        False, False, 'admin_user', dict(), dict(),
        200, id='unpublished-admin'),
    pytest.param(
        True, True, 'admin_user', dict(), dict(),
        400, id='duplicate'),
    pytest.param(
        True, False, 'other_test_user', dict(), dict(),
        401, id='not-oasis-admin'),
    pytest.param(
        True, False, None, dict(), dict(),
        401, id='no-credentials')])
def test_post_upload_bundle(
        client, proc_infra, non_empty_processed, internal_example_user_metadata, test_auth_dict,
        publish, test_duplicate, user, export_args, query_args, expected_status_code):
    # Create the bundle
    set_upload_entry_metadata(non_empty_processed, internal_example_user_metadata)
    if publish:
        non_empty_processed.publish_upload()
        non_empty_processed.block_until_complete(interval=.01)
    upload = non_empty_processed
    upload_id = upload.upload_id
    export_path = os.path.join(config.fs.tmp, 'bundle_' + upload_id)
    export_args_with_defaults = dict(
        export_as_stream=False, export_path=export_path, zipped=True, overwrite=True,
        export_settings=config.bundle_export.default_settings)
    export_args_with_defaults.update(export_args)
    BundleExporter(upload, **export_args_with_defaults).export_bundle()

    if not test_duplicate:
        # Delete the upload so we can import the bundle without id collisions
        upload.delete_upload_local()
    # Finally, import the bundle
    user_auth, __token = test_auth_dict[user]
    response = perform_post_put_file(
        client, 'POST', 'uploads/bundle', 'stream', export_path, user_auth, **query_args)
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        assert_processing(client, upload_id, user_auth, published=publish)
        upload = Upload.get(upload_id)
        assert upload.from_oasis and upload.oasis_deployment_url
    return


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
