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

    id_embargo:
        1 entry, 1 material, published with embargo
    id_embargo_w_coauthor:
        1 entry, 1 material, published with embargo and coauthor
    id_embargo_w_reviewer:
        1 entry, 1 material, published with embargo and reviewer
    id_unpublished:
        1 entry, 1 material, unpublished
    id_unpublished_w_coauthor:
        1 entry, 1 material, unpublished with coauthor
    id_unpublished_w_reviewer:
        1 entry, 1 material, unpublished with reviewer
    id_published:
        23 entries, 6 materials published without embargo
        partial archive exists only for id_01
        raw files and archive file for id_02 are missing
        id_10, id_11 reside in the same directory
    id_processing:
        unpublished upload without any entries, in status processing
    id_empty:
        unpublished upload without any entries
    '''
    data = ExampleData(main_author=test_user)

    # 6 uploads with different combinations of main_type and sub_type
    for main_type in ('embargo', 'unpublished'):
        for sub_type in ('', 'w_coauthor', 'w_reviewer'):
            upload_id = 'id_' + main_type + ('_' if sub_type else '') + sub_type
            if main_type == 'embargo':
                published = True
                embargo_length = 12
                upload_name = 'name_' + upload_id[3:]
            else:
                published = False
                embargo_length = 0
                upload_name = None
            calc_id = upload_id + '_1'
            coauthors = [other_test_user.user_id] if sub_type == 'w_coauthor' else None
            reviewers = [other_test_user.user_id] if sub_type == 'w_reviewer' else None
            data.create_upload(
                upload_id=upload_id,
                upload_name=upload_name,
                coauthors=coauthors,
                reviewers=reviewers,
                published=published,
                embargo_length=embargo_length)
            data.create_entry(
                upload_id=upload_id,
                calc_id=calc_id,
                material_id=upload_id,
                mainfile=f'test_content/{calc_id}/mainfile.json')

    # one upload with 23 calcs, published, no embargo
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
    data = ExampleData(main_author=test_user)

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
