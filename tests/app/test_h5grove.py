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
import h5py
from fastapi.testclient import TestClient

from nomad.app import h5grove_app
from nomad.utils.exampledata import ExampleData
from nomad.files import StagingUploadFiles
from nomad.config import config


@pytest.fixture
def h5grove_api(raw_files_function):
    h5grove_app.h5grove_router.settings.base_dir = config.fs.staging
    return TestClient(h5grove_app.app, base_url='http://testserver/')


@pytest.fixture
def upload_id():
    return 'nexus_test_upload'


@pytest.fixture
def example_data_nxs(user1, upload_id):
    data = ExampleData(main_author=user1)
    data.create_upload(upload_id)
    data.create_entry(upload_id=upload_id, entry_id='nexus_test_entry')
    data.save(with_files=False, with_mongo=True)


@pytest.mark.parametrize(
    'user, status_code, file_name',
    [
        pytest.param('invalid', 401, 'test.h5', id='invalid-credentials'),
        pytest.param(None, 401, 'test.h5', id='no-credentials'),
        pytest.param('user1', 404, 'some_other.h5', id='file-not-found'),
        pytest.param(
            'user1', 200, 'test.h5', id='file-in-published'
        ),  # TODO: Implement fetching from published upload
        pytest.param('user1', 200, 'test.h5', id='file-in-staging'),
    ],
)
def test_h5grove(
    auth_headers,
    h5grove_api,
    upload_id,
    proc_infra,
    example_data_nxs,
    user,
    status_code,
    file_name,
):
    test_file = 'test.h5'
    file_path = f'{StagingUploadFiles(upload_id=upload_id, create=True)._raw_dir}{os.sep}{test_file}'
    h5file = h5py.File(file_path, 'w')
    h5file.create_dataset('entry', data='test')
    h5file.close()
    resp = h5grove_api.get(
        f'/data/?file={file_name}&path=/entry&upload_id={upload_id}&source=raw',
        headers=auth_headers[user],
    )
    assert resp.status_code == status_code
    if status_code == 200:
        assert resp.content == b'"test"'
