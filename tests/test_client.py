# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import pytest
from bravado.client import SwaggerClient
import time

from nomad.processing import SUCCESS

from tests.test_files import example_file, create_public_upload, clear_files  # noqa pylint: disable=unused-import
from tests.test_api import client as flask_client, test_user_auth  # noqa pylint: disable=unused-import
from tests.bravado_flaks import FlaskTestHttpClient


@pytest.fixture(scope='function')
def client(flask_client, repository_db, test_user_auth):
    http_client = FlaskTestHttpClient(flask_client, headers=test_user_auth)
    return SwaggerClient.from_url('/swagger.json', http_client=http_client)


def test_get_upload_command(client):
    assert client.uploads.get_upload_command().response().result.upload_command is not None


def test_upload(client, worker):
    with open(example_file, 'rb') as f:
        upload = client.uploads.upload(file=f, name='test_upload').response().result

    while upload.tasks_running:
        upload = client.uploads.get_upload(upload_id=upload.upload_id).response().result
        time.sleep(0.1)

    assert upload.tasks_status == SUCCESS


def test_get_repo_calc(client, clear_files):
    create_public_upload('test_upload', 'pp')
    repo = client.repo.get_repo_calc(upload_id='test_upload', calc_id='0').response().result
    assert repo is not None
    assert repo['calc_id'] is not None
