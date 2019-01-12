# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
To run this test, a celery worker must be running. The test worker provided by
the celery pytest plugin is currently not working. It results on a timeout when
reading from the redis result backend, even though all task apperently ended successfully.
"""

from typing import Generator
import pytest
from datetime import datetime
import shutil
import os.path
import json

from nomad import utils
from nomad.files import ArchiveBasedStagingUploadFiles, UploadFiles, StagingUploadFiles
from nomad.processing import Upload, Calc
from nomad.processing.base import task as task_decorator
from nomad.repo import RepoUpload

from tests.test_files import example_file, empty_file

# import fixtures
from tests.test_files import clear_files  # pylint: disable=unused-import

example_files = [empty_file, example_file]


@pytest.fixture(scope='function', autouse=True)
def mocks_forall(mocksearch, mockmongo):
    pass


@pytest.fixture(scope='function', params=example_files)
def uploaded_id(request, clear_files) -> Generator[str, None, None]:
    example_file = request.param
    example_upload_id = os.path.basename(example_file).replace('.zip', '')
    upload_files = ArchiveBasedStagingUploadFiles(example_upload_id, create=True)
    shutil.copyfile(example_file, upload_files.upload_file_os_path)

    yield example_upload_id


@pytest.fixture
def uploaded_id_with_warning(request, clear_files) -> Generator[str, None, None]:
    example_file = 'tests/data/proc/examples_with_warning_template.zip'
    example_upload_id = os.path.basename(example_file).replace('.zip', '')
    upload_files = ArchiveBasedStagingUploadFiles(example_upload_id, create=True)
    shutil.copyfile(example_file, upload_files.upload_file_os_path)

    yield example_upload_id


def run_processing(uploaded_id: str, test_user) -> Upload:
    upload = Upload.create(upload_id=uploaded_id, user=test_user)
    upload.upload_time = datetime.now()

    assert upload.status == 'RUNNING'
    assert upload.current_task == 'uploading'

    upload.process()  # pylint: disable=E1101
    upload.block_until_complete(interval=.1)

    return upload


@pytest.fixture
def processed_upload(uploaded_id, test_user, worker, no_warn) -> Upload:
    return run_processing(uploaded_id, test_user)


def assert_processing(upload: Upload, mocksearch=None):
    assert upload.completed
    assert upload.current_task == 'cleanup'
    assert upload.upload_hash is not None
    assert len(upload.errors) == 0
    assert upload.status == 'SUCCESS'

    upload_files = UploadFiles.get(upload.upload_id, is_authorized=lambda: True)
    assert isinstance(upload_files, StagingUploadFiles)

    for calc in Calc.objects(upload_id=upload.upload_id):
        assert calc.parser is not None
        assert calc.mainfile is not None
        assert calc.status == 'SUCCESS', calc.archive_id
        calc_hash = utils.archive.calc_hash(calc.archive_id)

        with upload_files.archive_file(calc_hash) as archive_json:
            archive = json.load(archive_json)
        assert 'section_run' in archive
        assert 'section_calculation_info' in archive

        with upload_files.archive_log_file(calc_hash) as f:
            assert 'a test' in f.read()
        assert len(calc.errors) == 0

        with upload_files.raw_file(calc.mainfile) as f:
            f.read()

        if mocksearch:
            repo = mocksearch[calc.archive_id]
            assert repo is not None
            assert repo.chemical_composition is not None
            assert repo.basis_set_type is not None
            assert len(repo.aux_files) == 4


# @pytest.mark.timeout(30)
def test_processing(uploaded_id, worker, mocksearch, test_user, no_warn):
    upload = run_processing(uploaded_id, test_user)
    assert_processing(upload, mocksearch)


@pytest.mark.timeout(30)
def test_processing_with_warning(uploaded_id_with_warning, worker, test_user, mocksearch):
    upload = run_processing(uploaded_id_with_warning, test_user)
    assert_processing(upload, mocksearch)


@pytest.mark.timeout(30)
def test_process_non_existing(worker, test_user, with_error):
    upload = run_processing('__does_not_exist', test_user)

    assert upload.completed
    assert upload.current_task == 'extracting'
    assert upload.status == 'FAILURE'
    assert len(upload.errors) > 0


@pytest.mark.parametrize('task', ['extracting', 'parse_all', 'cleanup', 'parsing'])
@pytest.mark.timeout(30)
def test_task_failure(monkeypatch, uploaded_id, worker, task, test_user, with_error):
    # mock the task method to through exceptions
    if hasattr(Upload, task):
        cls = Upload
    elif hasattr(Calc, task):
        cls = Calc
    else:
        assert False

    def mock(self):
        raise Exception('fail for test')
    mock.__name__ = task
    mock = task_decorator(mock)

    monkeypatch.setattr('nomad.processing.data.%s.%s' % (cls.__name__, task), mock)

    # run the test
    upload = run_processing(uploaded_id, test_user)

    assert upload.completed

    if task != 'parsing':
        assert upload.status == 'FAILURE'
        assert upload.current_task == task
        assert len(upload.errors) > 0
    else:
        # there is an empty example with no calcs, even if past parsing_all task
        utils.get_logger(__name__).error('fake')
        if upload.total_calcs > 0:  # pylint: disable=E1101
            assert upload.status == 'SUCCESS'
            assert upload.current_task == 'cleanup'
            assert len(upload.errors) == 0
            for calc in upload.all_calcs(0, 100):  # pylint: disable=E1101
                assert calc.status == 'FAILURE'
                assert calc.current_task == 'parsing'
                assert len(calc.errors) > 0
