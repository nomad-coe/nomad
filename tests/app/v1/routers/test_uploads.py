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
from typing import List, Dict, Any
from tests.utils import build_url
from tests.test_files import assert_upload_files
from tests.test_search import assert_search_upload
from nomad import config, files, infrastructure
from nomad.processing import Upload, Calc, SUCCESS
from nomad.files import UploadFiles, PublicUploadFiles
from nomad.app.v1.routers.auth import generate_upload_token

'''
These are the tests for all API operations below ``entries``. The tests are organized
using the following type of methods: fixtures, ``perfrom_*_test``, ``assert_*``, and
``test_*``. While some ``test_*`` methods test individual API operations, some
test methods will test multiple API operations that use common aspects like
supporting queries, pagination, or the owner parameter. The test methods will use
``perform_*_test`` methods as an parameter. Similarely, the ``assert_*`` methods allow
to assert for certain aspects in the responses.
'''


def perform_get(client, base_url, user_auth=None, **query_args):
    headers = user_auth
    response = client.get(build_url(base_url, query_args), headers=headers)
    return response


def perform_post_uploads(client, mode, file, user_auth=None, token=None, **query_args):
    ''' Posts a new upload. '''
    if mode == 'local_path':
        query_args.update(local_path=file)
    if token:
        query_args.update(token=token)
    url = build_url('uploads', query_args)
    if mode == 'multipart':
        with open(file, 'rb') as f:
            response = client.post(
                url, files={'file': f}, headers=user_auth)
    elif mode == 'stream':
        with open(file, 'rb') as f:
            response = client.post(url, data=f.read(), headers=user_auth)
    elif mode == 'local_path':
        response = client.post(url, headers=user_auth)
    else:
        assert False, 'Invalid value for mode provided'

    return response


def perform_post_uploads_id_action(client, user_auth, upload_id, action, **query_args):
    return client.post(
        build_url(f'uploads/{upload_id}/action/{action}', query_args), headers=user_auth)


def perform_delete_uploads(client, upload_id, user_auth=None, **query_args):
    headers = user_auth
    response = client.delete(build_url(f'uploads/{upload_id}', query_args), headers=headers)
    return response


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
    assert response.status_code == 404

    assert Upload.objects(upload_id=upload_id).first() is None
    assert Calc.objects(upload_id=upload_id).count() is 0

    mongo_db = infrastructure.mongo_client[config.mongo.db_name]
    mongo_collection = mongo_db['archive']
    assert mongo_collection.find({}).count() == 0

    upload_files = UploadFiles.get(upload_id)
    assert upload_files is None or isinstance(upload_files, PublicUploadFiles)


def assert_processing(client, upload_id, user_auth):
    response_data = block_until_completed(client, upload_id, user_auth)

    assert len(response_data['tasks']) == 4
    assert response_data['tasks_status'] == SUCCESS
    assert response_data['current_task'] == 'cleanup'
    assert not response_data['process_running']

    # TODO: Also check calcs, like the old api tests.


def assert_published(
        client, user_auth, upload_id, proc_infra, expected_status_code=200, **query_args):
    '''
    Attempts to publish the given upload and check that it is successful (unless failure
    is expected).
    '''
    response = client.get('uploads/%s' % upload_id, headers=user_auth)
    upload = assert_upload(response.json())

    # Api call to actually publish the upload
    response = perform_post_uploads_id_action(client, user_auth, upload_id, 'publish', **query_args)

    assert response.status_code == expected_status_code
    if expected_status_code != 200:
        return  # The publish request is expected to be denied - nothing more to check

    with_embargo = query_args.get('with_embargo', True)
    embargo_length = query_args.get('embargo_length', 36)

    upload = assert_upload(response.json())
    assert upload['current_process'] == 'publish_upload'
    assert upload['process_running']

    block_until_completed(client, upload_id, user_auth)

    upload_proc = Upload.objects(upload_id=upload_id).first()
    assert upload_proc is not None
    assert upload_proc.published is True
    if with_embargo:
        assert upload_proc.embargo_length == embargo_length

    with upload_proc.entries_metadata() as entries:
        for entry in entries:
            assert entry.with_embargo == with_embargo

    assert_upload_files(upload_id, entries, files.PublicUploadFiles, published=True)
    assert_search_upload(entries, additional_keys=['with_embargo'], published=True)


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


@pytest.fixture(scope='function')
def slow_processing(monkeypatch):
    ''' Slow down processing to mitigate race conditions. '''
    old_cleanup = Upload.cleanup

    def slow_cleanup(self):
        time.sleep(2)
        old_cleanup(self)

    monkeypatch.setattr('nomad.processing.data.Upload.cleanup', slow_cleanup)
    yield True
    monkeypatch.setattr('nomad.processing.data.Upload.cleanup', old_cleanup)


def test_get_uploads_empty_list(client, mongo, test_user_auth):
    ''' Gets user's uploads, without having submitted anything -> empty list. '''
    response = perform_get(client, 'uploads', test_user_auth)
    assert response.status_code == 200
    assert len(response.json()['data']) == 0


def test_get_uploads_query(client, mongo, proc_infra, slow_processing, test_user_auth, non_empty_example_upload):
    ''' Tests various ways of getting the uppload with different filtering. '''
    upload_id_to_name = {}

    def assert_get_query(test_case: str, expected_result: tuple, **query_params):
        response = perform_get(client, 'uploads', user_auth=test_user_auth, **query_params)
        assert response.status_code == 200
        expected = set(expected_result)
        data: List[Dict[str, Any]] = response.json()['data']
        for upload in data:
            upload_id = upload['upload_id']
            assert upload_id in expected, (
                f'Test case {test_case} failed - got unexpected upload: {upload_id_to_name.get(upload_id, upload_id)}')
            expected.remove(upload_id)
        assert not expected, (
            f'Test case {test_case} failed - did not find {upload_id_to_name[list(expected)[0]]}')

    # Upload #1 - published
    response = perform_post_uploads(client, 'stream', non_empty_example_upload, test_user_auth, name='name1')
    assert response.status_code == 200
    upload_id_1 = response.json()['upload_id']
    upload_id_to_name[upload_id_1] = '#1'
    assert_processing(client, upload_id_1, test_user_auth)
    assert_published(client, test_user_auth, upload_id_1, proc_infra)

    # Upload #2 - wait for processing to finish, but do not publish
    response = perform_post_uploads(client, 'stream', non_empty_example_upload, test_user_auth, name='name2')
    assert response.status_code == 200
    upload_id_2 = response.json()['upload_id']
    upload_id_to_name[upload_id_2] = '#2'
    assert_processing(client, upload_id_2, test_user_auth)

    # Upload #3 - do NOT wait for processing to finish
    response = perform_post_uploads(client, 'stream', non_empty_example_upload, test_user_auth, name='name3')
    assert response.status_code == 200
    upload_id_3 = response.json()['upload_id']
    upload_id_to_name[upload_id_3] = '#3'

    # Test filter: UploadQuery.processing
    assert_get_query('is_processing_true', [upload_id_3], is_processing=True)
    assert_get_query('is_processing_false', [upload_id_1, upload_id_2], is_processing=False)
    assert_get_query('is_processing_unset', [upload_id_1, upload_id_2, upload_id_3])
    # Test filter: published/staging
    assert_get_query('is_published_False', [upload_id_2, upload_id_3], is_published=False)
    assert_get_query('is_published_True', [upload_id_1], is_published=True)
    assert_get_query('is_publushed_unset', [upload_id_1, upload_id_2, upload_id_3])
    # Test filter: upload_id
    assert_get_query('upload_id_single', [upload_id_1], upload_id=upload_id_1)
    assert_get_query('upload_id_multiple', [upload_id_1, upload_id_3], upload_id=[upload_id_1, upload_id_3])
    # Test filter: upload_name
    assert_get_query('upload_name_single', [upload_id_1], upload_name='name1')
    assert_get_query('upload_name_multiple', [upload_id_1, upload_id_3], upload_name=['name1', 'name3'])


def test_get_uploads_id_invalid(client, mongo, test_user_auth):
    response = perform_get(client, 'uploads/1234567890', test_user_auth)
    assert response.status_code == 404


@pytest.mark.parametrize('mode, name, user, use_upload_token, expected_status_code', [
    pytest.param('multipart', 'test_name', 'test_user', False, 200, id='post-multipart'),
    pytest.param('multipart', None, 'test_user', False, 200, id='post-multipart-no-name'),
    pytest.param('stream', 'test_name', 'test_user', False, 200, id='post-stream'),
    pytest.param('stream', None, 'test_user', False, 200, id='post-stream-no-name'),
    pytest.param('multipart', None, 'invalid', False, 401, id='post-multipart-no-name-invalid-cred'),
    pytest.param('stream', None, 'invalid', False, 401, id='post-stream-no-name-invalid-cred'),
    pytest.param('multipart', 'test_name', 'test_user', True, 200, id='post-multipart-token'),
    pytest.param('stream', 'test_name', 'test_user', True, 200, id='post-stream-token'),
    pytest.param('multipart', 'test_name', 'invalid', True, 401, id='post-multipart-token-invalid-cred'),
    pytest.param('stream', 'test_name', 'invalid', True, 401, id='post-stream-token-invalid-cred'),
    pytest.param('local_path', None, 'admin_user', False, 200, id='post-local_path'),
    pytest.param('multipart', None, None, False, 401, id='post-not-logged-in-multipart'),
    pytest.param('stream', None, None, False, 401, id='post-not-logged-in-stream'),
    pytest.param('local_path', None, None, False, 401, id='post-not-logged-in-local_path'),
    pytest.param('local_path', None, 'test_user', False, 401, id='post-not-admin-local_path')])
def test_post_uploads(
        client, mongo, proc_infra, test_user, admin_user, test_user_auth, admin_user_auth, non_empty_example_upload,
        mode, name, user, use_upload_token, expected_status_code):
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

    response = perform_post_uploads(client, mode, non_empty_example_upload, user_auth_post, token, name=name)
    assert response.status_code == expected_status_code
    if expected_status_code == 200:
        response_json = response.json()
        upload_id = response_json['upload_id']
        expected_name = name
        if not expected_name and mode in ('multipart', 'local_path'):
            expected_name = os.path.basename(non_empty_example_upload)
        assert_upload(response_json, name=expected_name)
        if mode == 'local_path':
            assert response_json['data']['upload_path'] == non_empty_example_upload

        assert_processing(client, upload_id, user_auth)


@pytest.mark.parametrize('empty', [
    pytest.param(False, id='non-empty'),
    pytest.param(True, id='empty')])
def test_post_uploads_with_publish_directly(
        client, test_user_auth, empty_upload, non_empty_example_upload, proc_infra, empty):
    ''' Posts uploads with publish_directly = True. '''
    if empty:
        file = empty_upload
    else:
        file = non_empty_example_upload
    response = perform_post_uploads(client, 'stream', file, test_user_auth, publish_directly=True)
    assert response.status_code == 200
    response_json = response.json()
    upload_id = response_json['upload_id']
    assert_upload(response_json)
    assert_processing(client, upload_id, test_user_auth)
    upload_proc = Upload.objects(upload_id=upload_id).first()
    if empty:
        assert not upload_proc.published
    else:
        assert upload_proc.published

        # TODO: verify entries, as in old api tests


@pytest.mark.parametrize('query_args, expected_status_code', [
    pytest.param({}, 200, id='no-args'),
    pytest.param(dict(with_embargo=True, embargo_length=12), 200, id='non-standard-embargo'),
    pytest.param(dict(embargo_length=24), 200, id='non-standard-embargo-length-only'),
    pytest.param(dict(embargo_length=100), 400, id='illegal-embargo-length'),
    pytest.param(dict(with_embargo=False), 200, id='no-embargo')])
def test_publish(
        client, test_user_auth, non_empty_example_upload, proc_infra,
        query_args, expected_status_code):
    ''' Tests the publish action with various arguments. '''
    response = perform_post_uploads(client, 'stream', non_empty_example_upload, test_user_auth)
    assert response.status_code == 200
    response_json = response.json()
    upload_id = response_json['upload_id']
    assert_upload(response_json)
    assert_processing(client, upload_id, test_user_auth)
    assert_published(
        client, test_user_auth, upload_id, proc_infra,
        expected_status_code=expected_status_code, **query_args)


def test_publish_empty(client, test_user_auth, empty_upload, proc_infra):
    ''' Tries to publish an empty upload (without entries). Should fail. '''
    response = perform_post_uploads(client, 'stream', empty_upload, test_user_auth)
    assert response.status_code == 200
    response_json = response.json()
    upload_id = response_json['upload_id']
    assert_upload(response_json)
    assert_processing(client, upload_id, test_user_auth)
    assert_published(client, test_user_auth, upload_id, proc_infra, expected_status_code=400)


@pytest.mark.timeout(config.tests.default_timeout)
def test_upload_limit(client, mongo, test_user, test_user_auth, proc_infra, non_empty_example_upload):
    ''' Tries to violate the limit on the number of unpublished uploads. '''
    old_upload_limit = config.services.upload_limit
    config.services.upload_limit = 5
    try:
        for _ in range(0, config.services.upload_limit):
            Upload.create(user=test_user)
        response = perform_post_uploads(client, 'stream', non_empty_example_upload, test_user_auth)
        assert response.status_code == 400
        assert Upload.user_uploads(test_user).count() == config.services.upload_limit
    finally:
        config.services.upload_limit = old_upload_limit


def test_delete_id_invalid(client, mongo, test_user_auth):
    ''' Trys to delete an invalid upload_id'''
    response = perform_delete_uploads(client, upload_id='1234567890', user_auth=test_user_auth)
    assert response.status_code == 404


@pytest.mark.parametrize('publish, delete_user, expected_status_code', [
    pytest.param(False, 'test_user', 200, id='delete-own'),
    pytest.param(False, 'other_test_user', 401, id='delete-others-not-admin'),
    pytest.param(False, 'admin_user', 200, id='delete-others-admin'),
    pytest.param(True, 'test_user', 401, id='delete-own-published'),
    pytest.param(True, 'admin_user', 200, id='delete-others-published-admin')])
def test_delete(
        client, mongo, proc_infra, non_empty_example_upload,
        test_user_auth, other_test_user_auth, admin_user_auth,
        publish, delete_user, expected_status_code):
    ''' Uploads a file, and then tries to delete it, with different parameters and users. '''
    delete_auth = {
        'test_user': test_user_auth,
        'other_test_user': other_test_user_auth,
        'admin_user': admin_user_auth
    }[delete_user]

    response = perform_post_uploads(
        client, 'multipart', non_empty_example_upload, test_user_auth)
    assert response.status_code == 200
    response_json = response.json()
    upload_id = response_json['upload_id']
    assert_upload(response_json)
    assert_processing(client, upload_id, test_user_auth)
    if publish:
        assert_published(client, test_user_auth, upload_id, proc_infra)

    response = perform_delete_uploads(client, upload_id, user_auth=delete_auth)
    assert response.status_code == expected_status_code
    if expected_status_code == 200:
        assert_upload_does_not_exist(client, upload_id, test_user_auth)
