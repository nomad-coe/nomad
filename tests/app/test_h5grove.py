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
from nomad import config


@pytest.fixture
def h5grove_api(raw_files):
    h5grove_app.h5grove_router.settings.base_dir = config.fs.staging
    return TestClient(h5grove_app.app, base_url='http://testserver/')


def test_h5grove(client, h5grove_api):
    file_name = "test.h5"
    file_path = f"{config.fs.staging}{os.sep}{file_name}"
    h5file = h5py.File(file_path, "w")
    h5file.create_dataset("entry", data="test")
    h5file.close()
    assert h5grove_api.get(f"/data/?file={file_name}&path=/entry").content == b'"test"'
