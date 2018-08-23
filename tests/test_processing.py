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
import time

import nomad.config as config
import nomad.files as files
from nomad.processing import start_processing, app

from tests.test_files import example_file, empty_file
# import fixtures
from tests.test_search import index  # pylint: disable=unused-import
from tests.test_files import clear_files  # pylint: disable=unused-import

example_files = [empty_file, example_file]


@pytest.fixture(scope='session')
def celery_session_worker():
    return app


@pytest.fixture(scope='session')
def celery_includes():
    return ['nomad.processing.tasks']


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


@pytest.mark.timeout(30)
def test_processing(uploaded_id, celery_session_worker):
    upload_proc = start_processing(uploaded_id)

    upload_proc.update_from_backend()

    assert upload_proc.status in ['PENDING', 'STARTED', 'PROGRESS']

    while not upload_proc.ready():
        time.sleep(1)
        upload_proc.update_from_backend()

    assert upload_proc.ready()
    assert upload_proc.current_task_name == 'cleanup'
    assert upload_proc.upload_hash is not None
    assert len(upload_proc.errors) == 0
    assert upload_proc.status == 'SUCCESS'
    for calc_proc in upload_proc.calc_procs:
        assert calc_proc.parser_name is not None
        assert calc_proc.mainfile is not None
        assert calc_proc.upload_hash is not None
        assert calc_proc.calc_hash is not None
        assert calc_proc.archive_id is not None
        assert calc_proc.status == 'SUCCESS'
        assert len(calc_proc.errors) == 0

    upload_proc.forget()


def test_process_non_existing(celery_session_worker):
    upload_proc = start_processing('__does_not_exist')

    upload_proc.get(timeout=30)

    assert upload_proc.ready()
    upload_proc.forget()

    assert upload_proc.current_task_name == 'extracting'
    assert upload_proc.status == 'FAILURE'
    assert len(upload_proc.errors) > 0
