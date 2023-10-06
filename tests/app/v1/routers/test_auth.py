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
from urllib.parse import urlencode


def perform_get_token_test(client, http_method, status_code, username, password):
    if http_method == 'post':
        response = client.post(
            'auth/token',
            data=dict(username=username, password=password))
    else:
        response = client.get('auth/token?%s' % urlencode(
            dict(username=username, password=password)))

    assert response.status_code == status_code


@pytest.mark.parametrize('http_method', ['post', 'get'])
def test_get_token(client, test_user, http_method):
    perform_get_token_test(client, http_method, 200, test_user.username, 'password')


@pytest.mark.parametrize('http_method', ['post', 'get'])
def test_get_token_bad_credentials(client, http_method):
    perform_get_token_test(client, http_method, 401, 'bad', 'credentials')


def test_get_signature_token(client, test_user_auth):
    response = client.get('auth/signature_token', headers=test_user_auth)
    assert response.status_code == 200
    assert response.json().get('signature_token') is not None


def test_get_signature_token_unauthorized(client, invalid_user_auth):
    response = client.get('auth/signature_token', headers=None)
    assert response.status_code == 401
    response = client.get('auth/signature_token', headers=invalid_user_auth)
    assert response.status_code == 401


@pytest.mark.parametrize(
    'expires_in, status_code',
    [(0, 422), (30 * 60, 200), (2 * 60 * 60, 200), (31 * 24 * 60 * 60, 422), (None, 422)])
def test_get_app_token(client, test_user_auth, expires_in, status_code):
    response = client.get(
        'auth/app_token', headers=test_user_auth, params={'expires_in': expires_in})
    assert response.status_code == status_code
    if status_code == 200:
        assert response.json().get('app_token') is not None


def test_get_app_token_unauthorized(client, invalid_user_auth):
    response = client.get('auth/app_token', headers=None,
                          params={'expires_in': 60})
    assert response.status_code == 401
    response = client.get('auth/app_token', headers=invalid_user_auth,
                          params={'expires_in': 60})
    assert response.status_code == 401
