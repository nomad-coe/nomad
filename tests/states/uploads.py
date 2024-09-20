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
from .groups import init_gui_test_groups, delete_group

default_access = {'coauthors': ['scooper'], 'reviewers': ['ttester']}
twin_access = {
    'coauthors': ['scooper'],
    'reviewers': ['ttester'],
    'coauthor_groups': ['group2'],
    'reviewer_groups': ['group3'],
}


def _build_example_data(
    main_author='test',
    coauthors=None,
    reviewers=None,
    coauthor_groups=None,
    reviewer_groups=None,
):
    """
    Helper function to set access fields for example data
    """
    get_user = infrastructure.user_management.get_user
    main_author = get_user(username=main_author)
    coauthors = [get_user(username=name).user_id for name in coauthors or []]
    reviewers = [get_user(username=name).user_id for name in reviewers or []]

    groups = init_gui_test_groups()
    coauthor_groups = coauthor_groups or []
    reviewer_groups = reviewer_groups or []
    if len(coauthor_groups) + len(reviewer_groups) > 0:
        coauthor_groups = [groups[label].group_id for label in coauthor_groups]
        reviewer_groups = [groups[label].group_id for label in reviewer_groups]

    data = ExampleData(
        main_author=main_author,
        coauthors=coauthors,
        reviewers=reviewers,
        coauthor_groups=coauthor_groups,
        reviewer_groups=reviewer_groups,
    )
    return data


def _create_vasp_upload(*, access, published, embargo_length):
    infrastructure.setup()
    data = _build_example_data(**access)
    upload_id = 'dft_upload'
    data.create_upload(
        upload_id=upload_id, published=published, embargo_length=embargo_length
    )
    data.create_entry(
        upload_id=upload_id,
        entry_id='dft_bulk',
        mainfile='vasp.xml',
        entry_archive=archive_dft_bulk(),
    )
    data.save()


def empty():
    """
    State published upload containing one entry
    """
    infrastructure.setup()
    data = _build_example_data()
    data.save()


def published():
    """
    1 upload, published, 1 entry, 1 coauthor, 1 reviewer
    """
    _create_vasp_upload(access=default_access, published=True, embargo_length=0)


def published_twin_access():
    """
    1 upload, published, 1 entry,
    1 coauthor = 1 coauthor group, 1 reviewer = 1 reviewer group
    """
    _create_vasp_upload(access=twin_access, published=True, embargo_length=0)


def published_coauthor_group():
    """
    1 upload, published, 1 entry, 1 coauthor group
    """
    access = {'coauthor_groups': ['group23']}
    _create_vasp_upload(access=access, published=True, embargo_length=0)


def published_reviewer_group():
    """
    1 upload, published, 1 entry, 1 reviewer group
    """
    access = {'reviewer_groups': ['group23']}
    _create_vasp_upload(access=access, published=True, embargo_length=0)


def published_with_embargo():
    """
    1 upload, published, embargo, 1 coauthor, 1 reviewer
    """
    _create_vasp_upload(access=default_access, published=True, embargo_length=3)


def published_with_embargo_coauthor_group():
    """
    1 upload, published, 1 entry, embargo, 1 coauthor group
    """
    access = {'coauthor_groups': ['group23']}
    _create_vasp_upload(access=access, published=True, embargo_length=3)


def published_with_embargo_reviewer_group():
    """
    1 upload, published, 1 entry, embargo, 1 reviewer group
    """
    access = {'reviewer_groups': ['group23']}
    _create_vasp_upload(access=access, published=True, embargo_length=3)


def unpublished():
    """
    1 upload, unpublished, 1 entry, 1 coauthor, 1 reviewer
    """
    _create_vasp_upload(access=default_access, published=False, embargo_length=0)


def unpublished_twin_access():
    """
    1 upload, published, embargo,
    1 coauthor = 1 coauthor group, 1 reviewer = 1 reviewer group
    """
    _create_vasp_upload(access=twin_access, published=False, embargo_length=0)


def unpublished_coauthor_group():
    """
    1 upload, unpublished, 1 entry, 1 coauthor group
    """
    access = {'coauthor_groups': ['group23']}
    _create_vasp_upload(access=access, published=False, embargo_length=0)


def unpublished_reviewer_group():
    """
    1 upload, unpublished, 1 entry, 1 reviewer group
    """
    access = {'reviewer_groups': ['group23']}
    _create_vasp_upload(access=access, published=False, embargo_length=0)


def unpublished_deleted_coauthor_group():
    """
    1 upload, unpublished, 1 entry, 1 coauthor group (deleted)
    """
    access = {'coauthor_groups': ['group23']}
    _create_vasp_upload(access=access, published=False, embargo_length=0)
    delete_group('group23')


def multiple_entries():
    """
    State published upload containing multiple entries
    """
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test')
    data = ExampleData(main_author=main_author)

    upload_id = 'dft_upload_1'
    data.create_upload(upload_id=upload_id, published=False, embargo_length=0)

    for i in range(1, 7):
        data.create_entry(
            upload_id=upload_id,
            entry_id=f'dft_bulk_{i}',
            mainfile=f'vasp_{i}.xml',
            entry_archive=archive_dft_bulk(),
        )
    data.save()


def multiple_uploads():
    """
    State published upload containing multiple entries
    """
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test')
    data = ExampleData(main_author=main_author)

    for i in range(1, 12):
        upload_id = f'dft_upload_{i}'
        data.create_upload(
            upload_id=upload_id, published=(i % 2 == 0), embargo_length=0
        )
        data.create_entry(
            upload_id=upload_id,
            entry_id=f'dft_bulk_{i}',
            mainfile='vasp.xml',
            entry_archive=archive_dft_bulk(),
        )
    data.save()


def maximum_unpublished():
    """
    State published upload containing multiple entries
    """
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test')
    data = ExampleData(main_author=main_author)

    for i in range(1, 11):
        upload_id = f'dft_upload_{i}'
        data.create_upload(upload_id=upload_id, published=False, embargo_length=0)
        data.create_entry(
            upload_id=upload_id,
            entry_id=f'dft_bulk_{i}',
            mainfile='vasp.xml',
            entry_archive=archive_dft_bulk(),
        )
    data.save()


def _browser_test(published: bool):
    infrastructure.setup()
    access = {'coauthors': ['scooper'], 'reviewers': ['ttester']}
    data = _build_example_data(**access)
    upload_id = 'browser_test'
    data.create_upload(
        upload_id=upload_id, published=published, embargo_length=12 if published else 0
    )
    data.create_entry(
        upload_id=upload_id,
        entry_id='dft_bulk',
        mainfile='test_entry/vasp.xml',
        entry_archive=archive_dft_bulk(),
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
        reviewers=reviewers,
    )
    upload.save()
    files.StagingUploadFiles(upload_id=upload.upload_id, create=True)
    upload.staging_upload_files.add_rawfiles('tests/data/proc/examples_vasp.zip')
    upload.staging_upload_files.add_rawfiles('examples/data/custom-schema')
    upload.staging_upload_files.add_rawfiles(
        'tests/data/datamodel/metainfo/eln/material_library'
    )
    upload.process_upload()
    upload.block_until_complete()
