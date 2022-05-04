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

from nomad import infrastructure, files
from nomad.processing import Upload
from nomad.utils.exampledata import ExampleData
from .archives.create_archives import archive_dft_bulk


def empty():
    '''
    State published upload containing one entry
    '''
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test')
    data = ExampleData(main_author=main_author)

    data.save()


def published():
    '''
    State published upload containing one entry
    '''
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test')
    coauthors = [infrastructure.user_management.get_user(username='scooper').user_id]
    reviewers = [infrastructure.user_management.get_user(username='ttester').user_id]
    data = ExampleData(main_author=main_author, coauthors=coauthors, reviewers=reviewers)

    upload_id = 'dft_upload'
    data.create_upload(upload_id=upload_id, published=True, embargo_length=0)
    entry_id = 'dft_bulk'
    data.create_entry(
        upload_id=upload_id,
        entry_id=entry_id,
        mainfile='vasp.xml',
        entry_archive=archive_dft_bulk()
    )

    data.save()


def published_with_embargo():
    '''
    State published upload but under 3 months embargo
    '''
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test')
    coauthors = [infrastructure.user_management.get_user(username='scooper').user_id]
    reviewers = [infrastructure.user_management.get_user(username='ttester').user_id]
    data = ExampleData(main_author=main_author, coauthors=coauthors, reviewers=reviewers)

    upload_id = 'dft_upload'
    data.create_upload(upload_id=upload_id, published=True, embargo_length=3)
    entry_id = 'dft_bulk'
    data.create_entry(
        upload_id=upload_id,
        entry_id=entry_id,
        mainfile='vasp.xml',
        entry_archive=archive_dft_bulk()
    )

    data.save()


def unpublished():
    '''
    State unpublished upload containing one entry
    '''
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test')
    coauthors = [infrastructure.user_management.get_user(username='scooper').user_id]
    reviewers = [infrastructure.user_management.get_user(username='ttester').user_id]
    data = ExampleData(main_author=main_author, coauthors=coauthors, reviewers=reviewers)

    upload_id = 'dft_upload'
    data.create_upload(upload_id=upload_id, published=False, embargo_length=0)
    entry_id = 'dft_bulk'
    data.create_entry(
        upload_id=upload_id,
        entry_id=entry_id,
        mainfile='vasp.xml',
        entry_archive=archive_dft_bulk()
    )

    data.save()


def multiple_entries():
    '''
    State published upload containing multiple entries
    '''
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test')
    data = ExampleData(main_author=main_author)

    upload_id = 'dft_upload_1'
    data.create_upload(upload_id=upload_id, published=True, embargo_length=0)

    for i in range(1, 7):
        entry_id = f'dft_bulk_{i}'
        data.create_entry(
            upload_id=upload_id,
            entry_id=entry_id,
            mainfile=f'vasp_{i}.xml',
            entry_archive=archive_dft_bulk()
        )

    data.save()


def multiple_uploads():
    '''
    State published upload containing multiple entries
    '''
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test')
    data = ExampleData(main_author=main_author)

    for i in range(1, 12):
        upload_id = f'dft_upload_{i}'
        data.create_upload(upload_id=upload_id, published=(i % 2 == 0), embargo_length=0)
        entry_id = f'dft_bulk_{i}'
        data.create_entry(
            upload_id=upload_id,
            entry_id=entry_id,
            mainfile=f'vasp.xml',
            entry_archive=archive_dft_bulk()
        )

    data.save()


def maximum_unpublished():
    '''
    State published upload containing multiple entries
    '''
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test')
    data = ExampleData(main_author=main_author)

    for i in range(1, 11):
        upload_id = f'dft_upload_{i}'
        data.create_upload(upload_id=upload_id, published=False, embargo_length=0)
        entry_id = f'dft_bulk_{i}'
        data.create_entry(
            upload_id=upload_id,
            entry_id=entry_id,
            mainfile=f'vasp.xml',
            entry_archive=archive_dft_bulk()
        )

    data.save()


def _browser_test(published: bool):
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test')
    coauthors = [infrastructure.user_management.get_user(username='scooper').user_id]
    reviewers = [infrastructure.user_management.get_user(username='ttester').user_id]
    data = ExampleData(main_author=main_author, coauthors=coauthors, reviewers=reviewers)

    upload_id = 'browser_test'
    data.create_upload(upload_id=upload_id, published=published, embargo_length=12 if published else 0)
    entry_id = 'dft_bulk'
    data.create_entry(
        upload_id=upload_id,
        entry_id=entry_id,
        mainfile='test_entry/vasp.xml',
        entry_archive=archive_dft_bulk()
    )
    data.save(additional_files_path='tests/data/gui/browser_test.zip')


def browser_test_published():
    _browser_test(True)


def browser_test_unpublished():
    _browser_test(False)


def archive_browser_test():
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test').user_id
    coauthors = [infrastructure.user_management.get_user(username='scooper').user_id]
    reviewers = [infrastructure.user_management.get_user(username='ttester').user_id]
    upload = Upload(
        upload_id='archive_browser_test',
        main_author=main_author,
        coauthors=coauthors,
        reviewers=reviewers)
    upload.save()
    files.StagingUploadFiles(upload_id=upload.upload_id, create=True)
    upload.staging_upload_files.add_rawfiles('tests/data/proc/examples_vasp.zip')
    upload.staging_upload_files.add_rawfiles('examples/data/custom-schema')
    upload.staging_upload_files.add_rawfiles('tests/data/datamodel/metainfo/eln/material_library')
    upload.process_upload()
