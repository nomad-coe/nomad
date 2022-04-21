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

from nomad.datamodel.context import ServerContext
from nomad.datamodel.datamodel import EntryArchive, EntryMetadata
from nomad.parsing.parser import ArchiveParser
from nomad.processing.data import Upload

from tests.normalizing.conftest import run_normalize
from tests.test_files import create_test_upload_files


def test_schema_processing(raw_files, no_warn):
    directory = 'tests/data/datamodel'
    mainfile = 'schema.archive.json'

    # create upload with example files
    upload_files = create_test_upload_files('test_upload_id', published=False, raw_files=directory)
    upload = Upload(upload_id='test_upload_id')

    # parse
    parser = ArchiveParser()
    context = ServerContext(upload=upload)
    test_archive = EntryArchive(m_context=context, metadata=EntryMetadata())
    parser.parse(
        upload_files.raw_file_object(mainfile).os_path,
        test_archive)
    run_normalize(test_archive)

    # assert archive
    assert len(test_archive.definitions.section_definitions) == 1
    assert test_archive.metadata.entry_type == 'Schema'
