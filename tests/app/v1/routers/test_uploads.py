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
from tests.utils import build_url
from nomad import config
from nomad.processing import SUCCESS

'''
These are the tests for all API operations below ``entries``. The tests are organized
using the following type of methods: fixtures, ``perfrom_*_test``, ``assert_*``, and
``test_*``. While some ``test_*`` methods test individual API operations, some
test methods will test multiple API operations that use common aspects like
supporting queries, pagination, or the owner parameter. The test methods will use
``perform_*_test`` methods as an parameter. Similarely, the ``assert_*`` methods allow
to assert for certain aspects in the responses.
'''


def perform_uploads_post(client, mode, file, user_auth=None, **params):
    ''' Posts a new upload. '''
    if mode == 'local_path':
        params.update(local_path=file)
    url = build_url('uploads', params)
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


def perform_get(client, base_url, user_auth=None, **params):
    headers = user_auth
    params = {}
    response = client.get(build_url(base_url, params), headers=headers)
    return response


def assert_upload(response_json, **kwargs):
    data = response_json['data']
    assert 'upload_id' in response_json
    assert 'upload_id' in data
    assert 'create_time' in data

    for key, value in kwargs.items():
        assert data.get(key, None) == value


def assert_processing(client, upload_id, user_auth):
    response_data = block_until_completed(client, upload_id, user_auth)

    assert len(response_data['tasks']) == 4
    assert response_data['tasks_status'] == SUCCESS
    assert response_data['current_task'] == 'cleanup'
    assert not response_data['process_running']

    # TODO: Also check calcs, like the old api tests.


def block_until_completed(client, upload_id: str, user_auth):
    t0 = time.time()
    while time.time() - t0 < config.tests.default_timeout:
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


@pytest.mark.parametrize('mode, name, user, expected_status_code', [
    pytest.param('multipart', 'test_name', 'test_user', 200, id='post-multipart'),
    pytest.param('multipart', None, 'test_user', 200, id='post-multipart-no-name'),
    pytest.param('stream', 'test_name', 'test_user', 200, id='post-stream'),
    pytest.param('stream', None, 'test_user', 200, id='post-stream-no-name'),
    pytest.param('local_path', None, 'admin_user', 200, id='post-local_path'),
    pytest.param('multipart', None, None, 401, id='post-not-logged-in-multipart'),
    pytest.param('stream', None, None, 401, id='post-not-logged-in-stream'),
    pytest.param('local_path', None, None, 401, id='post-not-logged-in-local_path'),
    pytest.param('local_path', None, 'test_user', 401, id='post-not-admin-local_path')])
def test_uploads_post(
        client, mongo, proc_infra, test_user_auth, admin_user_auth, example_upload,
        mode, name, user, expected_status_code):

    if user == 'test_user':
        user_auth = test_user_auth
    elif user == 'admin_user':
        user_auth = admin_user_auth
    else:
        user_auth = None

    response = perform_uploads_post(client, mode, example_upload, user_auth=user_auth, name=name)
    assert response.status_code == expected_status_code
    if expected_status_code == 200:
        response_json = response.json()
        upload_id = response_json['upload_id']
        expected_name = name
        if not expected_name and mode in ('multipart', 'local_path'):
            expected_name = os.path.basename(example_upload)
        assert_upload(response_json, name=expected_name)
        if mode == 'local_path':
            assert response_json['data']['upload_path'] == example_upload

        assert_processing(client, upload_id, user_auth)
