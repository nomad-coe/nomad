
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


import json

from nomad import config

from tests.utils import assert_log
from tests.app import resource  # pylint: disable=unused-import


class BlueprintClient():
    def __init__(self, app_client, blueprint_url_prefix):
        self.app_client = app_client
        self.blueprint_url_prefix = blueprint_url_prefix.strip('/')

    def _delegate(self, method, path, *args, **kwargs):
        app_client_function = getattr(self.app_client, method)
        prefixed_path = '/%s/%s' % (self.blueprint_url_prefix, path.lstrip('/'))
        return app_client_function(prefixed_path, *args, **kwargs)

    def get(self, *args, **kwargs):
        return self._delegate('get', *args, **kwargs)

    def post(self, *args, **kwargs):
        return self._delegate('post', *args, **kwargs)

    def put(self, *args, **kwargs):
        return self._delegate('put', *args, **kwargs)

    def delete(self, *args, **kwargs):
        return self._delegate('delete', *args, **kwargs)


def test_alive(client):
    rv = client.get('/alive')
    assert rv.status_code == 200


def test_internal_server_error_get(client, caplog):
    rv = client.get('/api/test/ise?test_arg=value')
    assert rv.status_code == 500
    record = assert_log(caplog, 'error', 'internal server error')
    data = json.loads(record.msg)

    assert data['blueprint'] == 'api'
    assert data['endpoint'] == 'api.test_internal_server_error_resource'
    assert data['method'] == 'GET'
    assert data['args']['test_arg'] == 'value'


def test_internal_server_error_post(client, caplog):
    rv = client.post(
        '/api/test/ise',
        content_type='application/json',
        data=json.dumps(dict(test_arg='value')))
    assert rv.status_code == 500
    record = assert_log(caplog, 'error', 'internal server error')
    data = json.loads(record.msg)

    assert data['blueprint'] == 'api'
    assert data['endpoint'] == 'api.test_internal_server_error_resource'
    assert data['method'] == 'POST'
    assert data['json']['test_arg'] == 'value'


def test_docs(client):
    rv = client.get('/docs/index.html')
    rv = client.get('/docs/introduction.html')
    assert rv.status_code == 200


def test_dist(client):
    rv = client.get('/dist/nomad-%s.tar.gz' % config.version)
    assert rv.status_code == 200
