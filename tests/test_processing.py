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

import nomad.config as config
import nomad.files as files
from nomad.processing import UploadProcessing

from tests.test_files import example_file, empty_file
# import fixtures
from tests.test_search import index  # pylint: disable=unused-import
from tests.test_files import clear_files  # pylint: disable=unused-import

example_files = [empty_file, example_file]


# disable test worker for now, see docstring above
# all further celery_* fixtures become effectivly mute.
@pytest.fixture(scope='session')
def celery_session_worker():
    return None


@pytest.fixture(scope='session')
def celery_includes():
    return ['nomad.processing']


@pytest.fixture(scope='session')
def celery_config():
    return {
        'broker_url': config.celery.broker_url,
        'result_backend': config.celery.backend_url,
        'accept_content': ['json', 'pickle'],
        'task_serializer': config.celery.serializer,
        'result_serializer': config.celery.serializer
    }


@pytest.fixture(scope='function', params=example_files)
def uploaded_id(request, clear_files) -> Generator[str, None, None]:
    example_file = request.param
    example_upload_id = example_file.replace('.zip', '')
    files._client.fput_object(config.files.uploads_bucket, example_upload_id, example_file)

    yield example_upload_id


def test_processing(uploaded_id, celery_session_worker):
    run = UploadProcessing(uploaded_id)
    run.start()

    # test that the instance can be reinstantiated from a persistable representation
    run = UploadProcessing.from_result_backend(uploaded_id, run.result_tuple)

    assert run.status in ['PENDING', 'PROGRESS']

    run.get(timeout=10)

    assert run.ready()
    assert run.task_name == 'nomad.processing.close_upload'
    assert run.upload_hash is not None
    assert run.cause is None
    assert run.status == 'SUCCESS'
    for processing_result in run.processing_results:
        for stage, (status, errors) in processing_result:
            assert errors is None or len(errors) == 0
            assert status in ['ParseSuccess', 'NormalizeSuccess', 'IndexSuccess', 'PersistenceSuccess']
    for calc_proc in run.calc_processings:
        assert 'parser' in calc_proc
        assert 'mainfile' in calc_proc
        assert 'pipeline' in calc_proc
        for stage in calc_proc['pipeline']:
            assert 'stage' in stage
            assert 'status' in stage
            assert stage['status'] in ['ParseSuccess', 'NormalizeSuccess', 'IndexSuccess', 'PersistenceSuccess']
            assert 'errors' not in stage or len(stage['errors']) == 0

    run.forget()


def test_process_non_existing(celery_session_worker):
    run = UploadProcessing('__does_not_exist')
    run.start()

    run.get(timeout=10)

    assert run.ready()
    run.forget()

    assert run.task_name == 'nomad.processing.open_upload'
    assert run.status == 'SUCCESS'
    assert run.cause is not None
    assert isinstance(run.cause, KeyError)
