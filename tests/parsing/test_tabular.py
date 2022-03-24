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
import json

from nomad.files import StagingUploadFiles
from nomad.datamodel.datamodel import EntryArchive, EntryMetadata
from nomad.datamodel.context import ServerContext
from nomad.utils import generate_entry_id, strip
from nomad.parsing.tabular import TabularDataParser
from nomad.processing import Upload


@pytest.mark.parametrize('content', [
    pytest.param(strip('''
        header_0,header_1
        0_0,0_1
        1_0,1_1
    '''), id='simple')
])
def test_tabular(raw_files, content):
    upload_files = StagingUploadFiles(upload_id='test_upload', create=True)
    upload = Upload(upload_id='test_upload')
    mainfile = 'test.csv'
    with upload_files.raw_file(mainfile, 'wt') as f:
        f.write(content)

    parser = TabularDataParser()
    keys = parser.is_mainfile(
        upload_files.raw_file_object(mainfile).os_path,
        'text/application', bytes(), '')

    assert isinstance(keys, list)
    assert len(keys) == 2

    context = ServerContext(upload=upload)
    main_archive = EntryArchive(m_context=context, metadata=EntryMetadata(
        upload_id='test_upload',
        mainfile=mainfile,
        entry_id=generate_entry_id('test_upload', mainfile)))
    child_archives = {
        key: EntryArchive(m_context=context, metadata=EntryMetadata(
            upload_id='test_upload',
            mainfile=mainfile,
            mainfile_key=key,
            entry_id=generate_entry_id('test_upload', mainfile, key)))
        for key in keys}

    parser.parse(
        upload_files.raw_file_object(mainfile).os_path,
        main_archive, None, child_archives)

    print('# main: ', json.dumps(main_archive.m_to_dict(), indent=2))
    for key in keys:
        print(f'# {key}: ', json.dumps(child_archives[key].m_to_dict(), indent=2))
