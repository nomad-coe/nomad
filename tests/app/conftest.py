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
from fastapi.testclient import TestClient

from nomad.app.main import app
from nomad.app.v1.routers.auth import generate_simple_token, generate_upload_token
from nomad.datamodel import User


def create_auth_headers(token: str):
    return {'Authorization': f'Bearer {token}'}


@pytest.fixture(scope='session')
def user0_auth(user0: User):
    return create_auth_headers(user0.user_id)


@pytest.fixture(scope='session')
def user1_auth(user1: User):
    return create_auth_headers(user1.user_id)


@pytest.fixture(scope='session')
def user2_auth(user2: User):
    return create_auth_headers(user2.user_id)


@pytest.fixture(scope='session')
def invalid_user_auth():
    return create_auth_headers('invalid.bearer.token')


@pytest.fixture(scope='module')
def app_token_auth(user1: User):
    app_token = generate_simple_token(user1.user_id, expires_in=3600)
    return create_auth_headers(app_token)


@pytest.fixture(scope='session')
def auth_dict(users_dict, invalid_user_auth):
    """
    Returns a dictionary of the form {user_label: (auth_headers, token)}. The key 'invalid'
    contains an example of invalid credentials, and the key None contains (None, None).
    """
    auths = {
        label: (create_auth_headers(user.user_id), generate_upload_token(user))
        for label, user in users_dict.items()
    }
    auths['invalid'] = (invalid_user_auth, 'invalid.upload.token')
    auths[None] = (None, None)
    return auths


@pytest.fixture(scope='session')
def client():
    return TestClient(app, base_url='http://testserver/')
