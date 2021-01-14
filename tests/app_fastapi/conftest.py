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

from typing import List
import pytest
from fastapi.testclient import TestClient
from datetime import datetime

from nomad import infrastructure, config
from nomad.archive import write_partial_archive_to_mongo
from nomad.app_fastapi.main import app
from nomad.datamodel import EntryArchive, EntryMetadata, DFTMetadata, User


def create_auth_headers(user: User):
    return {
        'Authorization': 'Bearer %s' % user.user_id
    }


@pytest.fixture(scope='module')
def test_user_auth(test_user: User):
    return create_auth_headers(test_user)


@pytest.fixture(scope='module')
def other_test_user_auth(other_test_user: User):
    return create_auth_headers(other_test_user)


@pytest.fixture(scope='module')
def admin_user_auth(admin_user: User):
    return create_auth_headers(admin_user)


@pytest.fixture(scope='session')
def client():
    return TestClient(app, base_url='http://testserver/api/v1/')


@pytest.fixture(scope='module')
def example_data(elastic_infra, raw_files_infra, mongo_infra, test_user, other_test_user, normalized):
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
    from tests.conftest import clear_raw_files, clear_elastic_infra
    from tests.test_files import create_test_upload_files

    clear_elastic_infra()
    clear_raw_files()

    archives: List[EntryArchive] = []
    archive = EntryArchive()
    entry_metadata = archive.m_create(
        EntryMetadata,
        domain='dft',
        upload_id='upload_id_1',
        upload_time=datetime.now(),
        uploader=test_user,
        published=True,
        processed=True,
        with_embargo=False,
        atoms=['H', 'O'],
        n_atoms=2,
        parser_name='parsers/vasp')
    entry_metadata.m_create(
        DFTMetadata,
        code_name='VASP',
        xc_functional='GGA',
        system='bulk')
    entry_metadata.dft.optimade = normalized.section_metadata.dft.optimade
    archive.m_update_from_dict({
        'section_run': [{}],
        'section_workflow': {}
    })

    # one upload with two calc published with embargo, one shared
    archives.clear()
    entry_metadata.m_update(
        upload_id='id_embargo',
        calc_id='id_embargo',
        mainfile='test_content/test_embargo_entry/mainfile.json',
        shared_with=[],
        with_embargo=True)
    entry_metadata.a_elastic.index()
    archives.append(archive.m_copy(deep=True))
    entry_metadata.m_update(
        calc_id='id_embargo_shared',
        mainfile='test_content/test_embargo_entry_shared/mainfile.json',
        shared_with=[other_test_user])
    entry_metadata.a_elastic.index()
    archives.append(archive.m_copy(deep=True))
    create_test_upload_files(entry_metadata.upload_id, archives)

    # one upload with two calc in staging, one shared
    archives.clear()
    entry_metadata.m_update(
        upload_id='id_unpublished',
        calc_id='id_unpublished',
        mainfile='test_content/test_entry/mainfile.json',
        with_embargo=False,
        shared_with=[],
        published=False)
    entry_metadata.a_elastic.index()
    archives.append(archive.m_copy(deep=True))
    entry_metadata.m_update(
        calc_id='id_unpublished_shared',
        mainfile='test_content/test_entry_shared/mainfile.json',
        shared_with=[other_test_user])
    entry_metadata.a_elastic.index()
    archives.append(archive.m_copy(deep=True))
    create_test_upload_files(
        entry_metadata.upload_id, archives, published=False)

    # one upload with 23 calcs published
    archives.clear()
    for i in range(1, 24):
        mainfile = 'test_content/subdir/test_entry_%02d/mainfile.json' % i
        if i == 11:
            mainfile = 'test_content/subdir/test_entry_10/mainfile_11.json'
        entry_metadata.m_update(
            upload_id='id_published',
            calc_id='id_%02d' % i,
            mainfile=mainfile,
            with_embargo=False,
            published=True,
            shared_with=[])
        entry_metadata.a_elastic.index()
        if i != 2:
            archives.append(archive.m_copy(deep=True))
        if i == 1:
            write_partial_archive_to_mongo(archive)

    infrastructure.elastic_client.indices.refresh(index=config.elastic.index_name)
    create_test_upload_files(entry_metadata.upload_id, archives)

    yield

    clear_elastic_infra()
    clear_raw_files()
