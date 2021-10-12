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
import math

from nomad.archive import write_partial_archive_to_mongo
from nomad.datamodel import OptimadeEntry
from nomad.processing import ProcessStatus

from tests.utils import ExampleData


@pytest.fixture(scope='session')
def client(api_v1):
    return api_v1


@pytest.fixture(scope='module')
def example_data(elastic_module, raw_files_module, mongo_module, test_user, other_test_user, normalized):
    '''
    Provides a couple of uploads and entries including metadata, raw-data, and
    archive files.

    23 entries, 6 materials published without embargo
    1 entry, 1 material unpublished
    1 entry, 1 material unpublished shared
    1 entry, 1 material published with embargo
    1 entry, 1 material published shared with embargo

    partial archive exists only for id_01
    raw files and archive file for id_02 are missing
    id_10, id_11 reside in the same directory
    '''
    data = ExampleData(
        uploader=test_user)

    # one upload with two calc published with embargo, one shared
    data.create_upload(
        upload_id='id_embargo',
        upload_name='name_embargo',
        published=True,
        embargo_length=12)
    data.create_entry(
        upload_id='id_embargo',
        calc_id='id_embargo',
        material_id='id_embargo',
        mainfile='test_content/test_embargo_entry/mainfile.json')
    data.create_upload(
        upload_id='id_embargo_shared_upload',
        upload_name='name_embargo_shared',
        published=True,
        reviewers=[other_test_user.user_id],
        embargo_length=12)
    data.create_entry(
        upload_id='id_embargo_shared_upload',
        calc_id='id_embargo_shared',
        material_id='id_embargo_shared',
        mainfile='test_content/test_embargo_entry_shared/mainfile.json')

    # one upload with two calc in staging, one shared
    data.create_upload(
        upload_id='id_unpublished',
        published=False)
    data.create_entry(
        upload_id='id_unpublished',
        calc_id='id_unpublished',
        material_id='id_unpublished',
        mainfile='test_content/test_entry/mainfile.json')
    data.create_upload(
        upload_id='id_unpublished_shared_upload',
        published=False,
        reviewers=[other_test_user.user_id])
    data.create_entry(
        upload_id='id_unpublished_shared_upload',
        calc_id='id_unpublished_shared',
        material_id='id_unpublished_shared',
        mainfile='test_content/test_entry_shared/mainfile.json')

    # one upload with 23 calcs published
    data.create_upload(
        upload_id='id_published',
        upload_name='name_published',
        published=True)
    for i in range(1, 24):
        entry_id = 'id_%02d' % i
        material_id = 'id_%02d' % (int(math.floor(i / 4)) + 1)
        mainfile = 'test_content/subdir/test_entry_%02d/mainfile.json' % i
        kwargs = dict(optimade=OptimadeEntry(nelements=2, elements=['H', 'O']))
        if i == 11:
            mainfile = 'test_content/subdir/test_entry_10/mainfile_11.json'
        if i == 1:
            kwargs['pid'] = '123'
        data.create_entry(
            upload_id='id_published',
            calc_id=entry_id,
            material_id=material_id,
            mainfile=mainfile,
            **kwargs)

        if i == 1:
            archive = data.archives[entry_id]
            write_partial_archive_to_mongo(archive)

    # one upload, no calcs, still processing
    data.create_upload(
        upload_id='id_processing',
        published=False,
        process_status=ProcessStatus.RUNNING)

    # one upload, no calcs, unpublished
    data.create_upload(
        upload_id='id_empty',
        published=False)

    data.save(with_files=False)
    del(data.archives['id_02'])
    data.save(with_files=True, with_es=False, with_mongo=False)


@pytest.fixture(scope='function')
def example_data_writeable(mongo, test_user, normalized):
    data = ExampleData(uploader=test_user)

    # one upload with one entry, published
    data.create_upload(
        upload_id='id_published_w',
        published=True,
        embargo_length=12)
    data.create_entry(
        upload_id='id_published_w',
        calc_id='id_published_w_entry',
        mainfile='test_content/test_embargo_entry/mainfile.json')

    # one upload with one entry, unpublished
    data.create_upload(
        upload_id='id_unpublished_w',
        published=False,
        embargo_length=12)
    data.create_entry(
        upload_id='id_unpublished_w',
        calc_id='id_unpublished_w_entry',
        mainfile='test_content/test_embargo_entry/mainfile.json')

    # one upload, no entries, still processing
    data.create_upload(
        upload_id='id_processing_w',
        published=False,
        process_status=ProcessStatus.RUNNING)

    # one upload, no entries, unpublished
    data.create_upload(
        upload_id='id_empty_w',
        published=False)

    data.save()

    yield

    data.delete()
