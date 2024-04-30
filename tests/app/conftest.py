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


def create_auth_header(token: str):
    return {'Authorization': f'Bearer {token}'}


@pytest.fixture(scope='module')
def app_token_auth(user1: User):
    app_token = generate_simple_token(user1.user_id, expires_in=3600)
    return create_auth_header(app_token)


@pytest.fixture(scope='session')
def auth_headers(users_dict):
    """Return a dict: user label -> auth header.

    The key 'invalid' contains an invalid header.
    The key 'empty' contains an empty header.
    The key None contains None.
    """
    headers = {
        label: create_auth_header(user.user_id) for label, user in users_dict.items()
    }
    headers['empty'] = {}
    headers['invalid'] = create_auth_header('invalid.bearer.token')
    headers[None] = None
    return headers


@pytest.fixture(scope='session')
def upload_tokens(users_dict):
    """Return a dict: user label -> upload token.

    The key 'invalid' contains an invalid token.
    The key 'empty' contains an empty token.
    The key None contains None.
    """
    headers = {label: generate_upload_token(user) for label, user in users_dict.items()}
    headers['empty'] = {}
    headers['invalid'] = 'invalid.upload.token'
    headers[None] = None
    return headers


@pytest.fixture(scope='session')
def client():
    return TestClient(app, base_url='http://testserver/')
