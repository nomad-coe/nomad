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

from nomad import config
from nomad.parsing.elabftw.elabftw import ELabFTWParser
from nomad.datamodel import EntryArchive, EntryMetadata, Context
from nomad.files import StagingUploadFiles, UploadFiles
from tests.processing.test_data import run_processing


def _assert_parsed_data(
    upload_id,
    entries,
    cls,
    parser,
    assert_fnc,
    mainfile_str,
    no_archive: bool = False,
    **kwargs,
):
    upload_files = UploadFiles.get(upload_id)
    assert upload_files is not None
    assert isinstance(upload_files, cls)

    for entry in entries:
        with upload_files.raw_file(entry.mainfile, 'rb') as f:
            f.read()

        try:
            with upload_files.read_archive(entry.entry_id) as archive:
                assert entry.entry_id in archive
            if entry.entry_name and mainfile_str in entry.entry_name:
                test_archive = EntryArchive(
                    metadata=EntryMetadata(entry_id=entry.entry_id, upload_id=upload_id)
                )
                test_archive.m_context = Context()
                mainfile = '/'.join([upload_files.os_path, 'raw', entry.mainfile])
                parser.parse(
                    mainfile, test_archive, None, child_archives={'0': test_archive}
                )

                assert_fnc(test_archive)

        except KeyError:
            assert no_archive

    upload_files.close()


def _assert_elabftw(test_archive):
    assert test_archive.data is not None
    assert test_archive.data.title == 'Test'
    assert test_archive.data.id == 'ro-crate-metadata.json'
    assert len(test_archive.data.experiment_data.experiments_links) == 1
    assert len(test_archive.data.experiment_files) == 5
    assert test_archive.data.experiment_data.experiments_links[0].title == 'JSON test '
    assert test_archive.data.experiment_data.items_links[0].title == 'Untitled'
    assert test_archive.data.experiment_data.extra_fields is not None
    for item in test_archive.data.experiment_files:
        assert item.type == 'File'
        assert item.file is not None
        assert item.id is not None


@pytest.mark.timeout(config.tests.default_timeout)
def test_elabftw_parser(raw_files, proc_infra, api_v1, test_user):
    upload = run_processing(
        ('test_upload', 'tests/data/parsers/elabftw/test.eln'), test_user
    )

    assert upload.total_entries_count == 2
    assert len(upload.successful_entries) == 2

    with upload.entries_metadata() as entries:
        _assert_parsed_data(
            upload.upload_id,
            entries,
            StagingUploadFiles,
            ELabFTWParser(),
            _assert_elabftw,
            'ro-crate-metadata.json',
            published=False,
        )
