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

from typing import List, Dict, Union
import pytest
from datetime import datetime
from fastapi.testclient import TestClient

from nomad.archive import write_partial_archive_to_mongo
from nomad.datamodel import EntryArchive, EntryMetadata
from nomad.datamodel.results import Results
from nomad.metainfo.elasticsearch_extension import index_entries
from nomad.app.main import app


@pytest.fixture(scope='session')
def client():
    return TestClient(app, base_url='http://testserver/api/v1/')


class ExampleData:
    def __init__(self, dft: dict = None, **kwargs):
        self.uploads: Dict[str, List[str]] = dict()
        self.entries: Dict[str, EntryMetadata] = dict()
        self.archives: Dict[str, EntryArchive] = dict()

        self.entry_defaults = kwargs
        self.dft_defaults = dft

    def save(self, with_files: bool = True, with_mongo: bool = True, with_es: bool = True):
        from tests.test_files import create_test_upload_files
        from nomad import processing as proc

        # save to elastic and mongo
        for entry_metadata in self.entries.values():
            if with_mongo:
                mongo_entry = proc.Calc(
                    create_time=datetime.now(),
                    calc_id=entry_metadata.calc_id,
                    upload_id=entry_metadata.upload_id,
                    mainfile=entry_metadata.mainfile)
                mongo_entry.apply_entry_metadata(entry_metadata)
                mongo_entry.save()

            if with_es:
                index_entries(
                    list(self.archives.values()), update_materials=True, refresh=True)

        # create upload files
        if with_files:
            published = True
            for upload_id, entry_ids in self.uploads.items():
                archives = []
                for entry_id in entry_ids:
                    published &= self.entries[entry_id].published
                    if entry_id in self.archives:
                        archives.append(self.archives[entry_id])

                create_test_upload_files(upload_id, archives, published=published)

    def _create_entry(
            self,
            calc_id: str, upload_id: str, mainfile: str,
            results: Union[Results, dict] = None, archive: dict = None, **kwargs):

        entry_archive = EntryArchive()
        entry_metadata = entry_archive.m_create(EntryMetadata)
        entry_metadata.m_update(
            calc_id=calc_id,
            upload_id=upload_id,
            domain='dft',
            mainfile=mainfile,
            upload_time=datetime.now(),
            published=True,
            processed=True,
            with_embargo=False,
            parser_name='parsers/vasp')
        entry_metadata.m_update(**self.entry_defaults)
        entry_metadata.m_update(**kwargs)

        if results is None:
            results = {
                'material': {
                    'material_id': 'test_material_id',
                    'elements': ['H', 'O'],
                    'nelements': 2
                },
                'method': {
                    'simulation': {
                        'program_name': 'VASP',
                        'dft': {
                            'xc_functional_type': 'GGA'
                        }
                    }
                },
                'properties': {
                    'n_calculations': 1
                }
            }

        if isinstance(results, dict):
            section_results = Results.m_from_dict(results)
        else:
            section_results = results

        assert isinstance(section_results, Results)
        entry_archive.m_add_sub_section(EntryArchive.results, section_results)

        entry_archive.m_update_from_dict(dict(
            section_run=[{}],
            section_workflow={}))
        if archive is not None:
            entry_archive.m_update(**archive)

        entry_id = entry_metadata.calc_id
        self.archives[entry_id] = entry_archive
        self.entries[entry_id] = entry_metadata
        self.uploads.setdefault(entry_metadata.upload_id, []).append(entry_id)


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
    data._create_entry(
        upload_id='id_embargo',
        calc_id='id_embargo',
        mainfile='test_content/test_embargo_entry/mainfile.json',
        shared_with=[],
        with_embargo=True)
    data._create_entry(
        upload_id='id_embargo',
        calc_id='id_embargo_shared',
        mainfile='test_content/test_embargo_entry_shared/mainfile.json',
        shared_with=[other_test_user],
        with_embargo=True)

    # one upload with two calc in staging, one shared
    data._create_entry(
        upload_id='id_unpublished',
        calc_id='id_unpublished',
        mainfile='test_content/test_entry/mainfile.json',
        with_embargo=False,
        shared_with=[],
        published=False)
    data._create_entry(
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
        data._create_entry(
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
