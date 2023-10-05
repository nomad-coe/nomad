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

from nomad import files
from nomad.processing import Upload, Entry, ProcessStatus
from tests.normalizing.conftest import run_processing


def _create_upload(upload_id, user_id, file_paths=None):
    if file_paths is None:
        file_paths = []
    upload = Upload(
        upload_id=upload_id,
        main_author=user_id)
    upload.save()
    files.StagingUploadFiles(upload_id=upload.upload_id, create=True)
    for file_path in file_paths:
        upload.staging_upload_files.add_rawfiles(file_path)
    upload.process_upload()
    upload.block_until_complete()
    return upload


@pytest.mark.parametrize('mainfile, assert_xpaths', [
    pytest.param('schema.archive.yaml', [], id='schema'),
    pytest.param('sample.archive.json', ['data.processes.pvd_evaporation.time'], id='sample'),
    pytest.param('PVD-P.archive.json', [], id='instrument'),
    pytest.param('Zinc_Selenide.archive.json', [], id='chemical')
])
def test_eln(mainfile, assert_xpaths, raw_files, no_warn):
    mainfile_directory = 'examples/data/eln'
    archive = run_processing(mainfile_directory, mainfile)

    for xpath in assert_xpaths:
        assert archive.m_xpath(xpath) is not None


@pytest.mark.parametrize('mainfile, assert_xpaths', [
    pytest.param('tabular-parser-col-mode.archive.yaml', ['data.My_Quantity'], id='col_mode'),
    pytest.param('tabular-parser-row-mode.archive.yaml', ['data.My_Subsection.My_Section[4].My_Quantity'],
                 id='row_mode'),
])
def test_sample_tabular(mainfile, assert_xpaths, raw_files, no_warn):
    mainfile_directory = 'examples/data/docs'
    archive = run_processing(mainfile_directory, mainfile)

    for xpath in assert_xpaths:
        assert archive.m_xpath(xpath) is not None


@pytest.mark.parametrize('test_files,number_of_entries', [
    pytest.param([
        'examples/data/docs/tabular-parser-entry-mode.archive.yaml',
        'examples/data/docs/tabular-parser-entry-mode.xlsx'
    ], 5, id='simple_entry_mode'),
    pytest.param([
        'examples/data/docs/tabular-parser-complex.archive.yaml',
        'examples/data/docs/data_file_1.csv',
        'examples/data/docs/data_file_2.csv'
    ], 6, id='complex_entry_mode')
])
def test_sample_entry_mode(mongo, test_user, raw_files, monkeypatch, proc_infra, test_files, number_of_entries):
    upload = _create_upload('test_upload_id', test_user.user_id, test_files)
    assert upload is not None
    assert upload.processed_entries_count == number_of_entries

    for entry in Entry.objects(upload_id='test_upload_id'):
        assert entry.process_status == ProcessStatus.SUCCESS


@pytest.mark.parametrize('test_files, status', [
    pytest.param([
        'examples/data/docs/tabular-parser_1_column_current-entry_to-root.archive.yaml',
        'examples/data/docs/test.xlsx'
    ], ProcessStatus.SUCCESS, id='test_1'),
    pytest.param([
        'examples/data/docs/tabular-parser_2_column_current-entry_to-path.archive.yaml',
        'examples/data/docs/test.xlsx'
    ], ProcessStatus.SUCCESS, id='test_2'),
    pytest.param([
        'examples/data/docs/tabular-parser_3_row_current-entry_to-path.archive.yaml',
        'examples/data/docs/test.xlsx'
    ], ProcessStatus.SUCCESS, id='test_3'),
    pytest.param([
        'examples/data/docs/tabular-parser_4_column_single-new-entry_to-path.archive.yaml',
        'examples/data/docs/test.xlsx'
    ], ProcessStatus.SUCCESS, id='test_4'),
    pytest.param([
        'examples/data/docs/tabular-parser_5_row_single-new-entry_to-path.archive.yaml',
        'examples/data/docs/test.xlsx'
    ], ProcessStatus.SUCCESS, id='test_5'),
    pytest.param([
        'examples/data/docs/tabular-parser_6_row_multiple-new-entries_to-root.archive.yaml',
        'examples/data/docs/test.xlsx'
    ], ProcessStatus.SUCCESS, id='test_6'),
    pytest.param([
        'examples/data/docs/tabular-parser_7_row_multiple-new-entries_to-path.archive.yaml',
        'examples/data/docs/test.xlsx'
    ], ProcessStatus.SUCCESS, id='test_7'),
    pytest.param([
        'examples/data/docs/tabular-parser_8_row_current-entry_to-path_subsubsection.archive.yaml',
        'examples/data/docs/test.xlsx'
    ], ProcessStatus.SUCCESS, id='test_8'),
    pytest.param([
        'examples/data/docs/tabular-parser_np1_row_current-entry_to-root.archive.yaml',
        'examples/data/docs/test.xlsx'
    ], ProcessStatus.SUCCESS, id='test_9'),
    pytest.param([
        'examples/data/docs/tabular-parser_np2_column_single-new-entry_to-root.archive.yaml',
        'examples/data/docs/test.xlsx'
    ], ProcessStatus.SUCCESS, id='test_10'),
    pytest.param([
        'examples/data/docs/tabular-parser_np3_row_single-new-entry_to-root.archive.yaml',
        'examples/data/docs/test.xlsx'
    ], ProcessStatus.SUCCESS, id='test_11'),
    pytest.param([
        'examples/data/docs/tabular-parser_np4_column_multiple-new-entries_to-root.archive.yaml',
        'examples/data/docs/test.xlsx'
    ], ProcessStatus.SUCCESS, id='test_12'),
    pytest.param([
        'examples/data/docs/tabular-parser_np5_column_multiple-new-entries_to-path.archive.yaml',
        'examples/data/docs/test.xlsx'
    ], ProcessStatus.SUCCESS, id='test_13'),
])
def test_tabular_doc_examples(mongo, test_user, raw_files, monkeypatch, proc_infra, test_files, status):
    upload = _create_upload('test_upload_id', test_user.user_id, test_files)
    assert upload is not None

    for entry in Entry.objects(upload_id='test_upload_id'):
        assert entry.process_status == status
