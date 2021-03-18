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

import time

from nomad.processing import SUCCESS
from nomad.datamodel import EntryMetadata

from tests.test_files import example_file
from tests.search.test_v0 import create_entry


def test_get_upload_command(bravado, no_warn):
    assert bravado.uploads.get_upload_command().response().result.upload_command is not None


def test_upload(bravado, proc_infra, no_warn):
    with open(example_file, 'rb') as f:
        upload = bravado.uploads.upload(file=f, name='test_upload').response().result

    while upload.tasks_running:
        upload = bravado.uploads.get_upload(upload_id=upload.upload_id).response().result
        time.sleep(0.1)

    assert upload.tasks_status == SUCCESS


def test_get_repo_calc(bravado, proc_infra, raw_files):
    create_entry(EntryMetadata(
        domain='dft', calc_id='0', upload_id='test_upload', published=True, with_embargo=False))
    repo = bravado.repo.get_repo_calc(upload_id='test_upload', calc_id='0').response().result
    assert repo is not None
    assert repo['calc_id'] is not None
