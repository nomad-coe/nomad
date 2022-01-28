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

from nomad.app.flask import app as flask_app


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
