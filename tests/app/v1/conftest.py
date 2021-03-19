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

from nomad.archive import write_partial_archive_to_mongo
from nomad.app.main import app

from tests.utils import ExampleData


@pytest.fixture(scope='session')
def client():
    return TestClient(app, base_url='http://testserver/api/v1/')


@pytest.fixture(scope='module')
def example_data(elastic_module, raw_files_module, mongo_module, test_user, other_test_user):
    '''
    Provides a couple of uploads and entries including metadata, raw-data, and
    archive files.

    23 published without embargo
    1 unpublished
    1 unpublished shared
    1 published with embargo
    1 published shared with embargo

    partial archive exists only for id_01
    raw files and archive file for id_02 are missing
    id_10, id_11 reside in the same directory
    '''

    data = ExampleData(
        uploader=test_user
    )

    # one upload with two calc published with embargo, one shared
    data.create_entry(
        upload_id='id_embargo',
        calc_id='id_embargo',
        mainfile='test_content/test_embargo_entry/mainfile.json',
        shared_with=[],
        with_embargo=True)
    data.create_entry(
        upload_id='id_embargo',
        calc_id='id_embargo_shared',
        mainfile='test_content/test_embargo_entry_shared/mainfile.json',
        shared_with=[other_test_user],
        with_embargo=True)

    # one upload with two calc in staging, one shared
    data.create_entry(
        upload_id='id_unpublished',
        calc_id='id_unpublished',
        mainfile='test_content/test_entry/mainfile.json',
        with_embargo=False,
        shared_with=[],
        published=False)
    data.create_entry(
        upload_id='id_unpublished',
        calc_id='id_unpublished_shared',
        mainfile='test_content/test_entry_shared/mainfile.json',
        shared_with=[other_test_user],
        with_embargo=False,
        published=False)

    # one upload with 23 calcs published
    for i in range(1, 24):
        entry_id = 'id_%02d' % i
        mainfile = 'test_content/subdir/test_entry_%02d/mainfile.json' % i
        if i == 11:
            mainfile = 'test_content/subdir/test_entry_10/mainfile_11.json'
        data.create_entry(
            upload_id='id_published',
            calc_id=entry_id,
            mainfile=mainfile)

        if i == 1:
            archive = data.archives[entry_id]
            write_partial_archive_to_mongo(archive)

    data.save(with_files=False)
    del(data.archives['id_02'])
    data.save(with_files=True, with_es=False, with_mongo=False)

    return data
