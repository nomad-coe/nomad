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
from bravado.client import SwaggerClient

from nomad.app_fastapi.flask import app as flask_app

from .bravado import FlaskTestHttpClient
from ..conftest import admin_user_auth, test_user_auth  # pylint: disable=unused-import


@pytest.fixture(scope='session')
def session_client():
    flask_app.config['TESTING'] = True
    client = flask_app.test_client()

    yield client


@pytest.fixture(scope='function')
def client(mongo, session_client):
    flask_app.config['TESTING'] = True
    client = flask_app.test_client()

    yield client


@pytest.fixture(scope='function')
def bravado(client, test_user_auth):
    http_client = FlaskTestHttpClient(client, headers=test_user_auth)
    return SwaggerClient.from_url('/api/swagger.json', http_client=http_client)


@pytest.fixture(scope='function')
def admin_user_bravado_client(client, admin_user_auth, monkeypatch):
    def create_client():
        http_client = FlaskTestHttpClient(client, headers=admin_user_auth)
        return SwaggerClient.from_url('/api/swagger.json', http_client=http_client)

    monkeypatch.setattr('nomad.cli.client.create_client', create_client)


@pytest.fixture(scope='function')
def test_user_bravado_client(client, test_user_auth, monkeypatch):
    def create_client():
        http_client = FlaskTestHttpClient(client, headers=test_user_auth)
        return SwaggerClient.from_url('/api/swagger.json', http_client=http_client)

    monkeypatch.setattr('nomad.cli.client.create_client', create_client)


@pytest.fixture(scope='function')
def oasis_central_nomad_client(client, test_user_auth, monkeypatch):
    def create_client(*args, **kwargs):
        http_client = FlaskTestHttpClient(client, headers=test_user_auth)
        return SwaggerClient.from_url('/api/swagger.json', http_client=http_client)

    monkeypatch.setattr('nomad.cli.client.client._create_client', create_client)
