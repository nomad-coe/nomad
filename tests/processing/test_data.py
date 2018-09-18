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

from nomad import config, files
from nomad.processing import Upload, Calc
from nomad.processing.base import task as task_decorator
from nomad.user import me
from nomad.repo import RepoCalc

from tests.test_files import example_file, empty_file

# import fixtures
from tests.test_files import clear_files  # pylint: disable=unused-import

example_files = [empty_file, example_file]


@pytest.fixture(scope='function', autouse=True)
def mocksearch_forall(mocksearch):
    pass


@pytest.fixture(scope='function', params=example_files)
def uploaded_id(request, clear_files) -> Generator[str, None, None]:
    example_file = request.param
    example_upload_id = example_file.replace('.zip', '')
    files._client.fput_object(config.files.uploads_bucket, example_upload_id, example_file)

    yield example_upload_id


def run_processing(uploaded_id: str) -> Upload:
    upload = Upload.create(upload_id=uploaded_id, user_id=me.email)
    upload.upload_time = datetime.now()

    assert upload.status == 'RUNNING'
    assert upload.current_task == 'uploading'

    upload.process()  # pylint: disable=E1101
    upload.block_until_complete(interval=.1)

    return upload


def assert_processing(upload: Upload):
    assert upload.completed
    assert upload.current_task == 'cleanup'
    assert upload.upload_hash is not None
    assert len(upload.errors) == 0
    assert upload.status == 'SUCCESS'

    for calc in Calc.objects(upload_id=upload.upload_id):
        assert calc.parser is not None
        assert calc.mainfile is not None
        assert calc.status == 'SUCCESS', calc.archive_id
        assert len(calc.errors) == 0


@pytest.mark.timeout(30)
def test_processing(uploaded_id, worker, no_warn):
    upload = run_processing(uploaded_id)
    assert_processing(upload)


@pytest.mark.parametrize('uploaded_id', [example_files[1]], indirect=True)
def test_processing_doublets(uploaded_id, worker, one_error):

    upload = run_processing(uploaded_id)
    assert upload.status == 'SUCCESS'
    assert RepoCalc.upload_exists(upload.upload_hash)  # pylint: disable=E1101

    upload = run_processing(uploaded_id)
    assert upload.status == 'FAILURE'
    assert len(upload.errors) > 0
    assert 'already' in upload.errors[0]


@pytest.mark.timeout(30)
def test_process_non_existing(worker, one_error):
    upload = run_processing('__does_not_exist')

    assert upload.completed
    assert upload.current_task == 'extracting'
    assert upload.status == 'FAILURE'
    assert len(upload.errors) > 0


@pytest.mark.parametrize('task', ['extracting', 'parse_all', 'cleanup', 'parsing'])
@pytest.mark.timeout(30)
def test_task_failure(monkeypatch, uploaded_id, worker, task, one_error):
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
    upload = run_processing(uploaded_id)

    assert upload.completed

    if task != 'parsing':
        assert upload.status == 'FAILURE'
        assert upload.current_task == task
        assert len(upload.errors) > 0
    else:
        # there is an empty example with no calcs, even if past parsing_all task
        if upload.total_calcs > 0:  # pylint: disable=E1101
            assert upload.status == 'SUCCESS'
            assert upload.current_task == 'cleanup'
            assert len(upload.errors) == 0
            for calc in upload.all_calcs(0, 100):  # pylint: disable=E1101
                assert calc.status == 'FAILURE'
                assert calc.current_task == 'parsing'
                assert len(calc.errors) > 0
