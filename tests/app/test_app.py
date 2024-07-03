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

from nomad.config import config


def test_alive(client):
    rv = client.get('/alive')
    assert rv.status_code == 200


@pytest.mark.skipif(os.getenv('RUN_DOCS_TEST') != '1', reason='Only run in build stage')
def test_docs(client):
    rv = client.get('/docs/index.html')
    assert rv.status_code == 200
    assert (
        f'max-age={config.services.html_resource_http_max_age}, must-revalidate'
        in rv.headers['Cache-Control']
    )
    assert 'Etag' in rv.headers

    rv = client.get('/docs/assets/favicon.png')
    assert rv.status_code == 200
    assert (
        f'max-age={config.services.image_resource_http_max_age}, must-revalidate'
        in rv.headers['Cache-Control']
    )
    assert 'Etag' in rv.headers

    etag = rv.headers['Etag']
    rv = client.get('/docs/assets/favicon.png', headers={'If-None-Match': etag})
    assert rv.status_code == 304
    rv = client.get('/docs/assets/favicon.png', headers={'If-None-Match': f'W/{etag}'})
    assert rv.status_code == 304


@pytest.mark.parametrize('path', ['env.js', 'artifacts.js'])
def test_gui(client, path, monkeypatch):
    monkeypatch.setattr('nomad.app.main.GuiFiles.gui_env_data', 'env.js')
    monkeypatch.setattr('nomad.app.main.GuiFiles.gui_artifacts_data', 'artifacts.js')
    monkeypatch.setattr('nomad.app.main.GuiFiles.gui_data_etag', 'etag')

    rv = client.get(f'/gui/{path}')
    assert rv.status_code == 200
    assert rv.text == path
    assert rv.headers.get('Etag') == '"etag"'

    rv = client.get(f'/gui/{path}', headers={'if-none-match': 'etag'})
    assert rv.status_code == 304

    rv = client.get(f'/gui/{path}', headers={'if-none-match': 'W/"etag"'})
    assert rv.status_code == 304

    rv = client.get(f'/gui/{path}', headers={'if-none-match': 'different-etag'})
    assert rv.status_code == 200
