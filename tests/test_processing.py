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

from typing import Generator
import pytest
from minio import ResponseError

import nomad.files as files
import nomad.config as config
from nomad.processing import UploadProcessing

example_files = ['data/examples_vasp.zip', 'data/empty.zip']


@pytest.fixture(scope='function', params=example_files)
def uploaded_id(request) -> Generator[str, None, None]:
    example_file = request.param
    example_upload_id = example_file.replace('.zip', '')
    files._client.fput_object(config.s3.uploads_bucket, example_upload_id, example_file)

    yield example_upload_id

    try:
        files._client.remove_object(config.s3.uploads_bucket, example_upload_id)
    except ResponseError:
        pass


def test_processing(uploaded_id):
    run = UploadProcessing(uploaded_id)
    run.start()
    run.get(timeout=30)
    assert run.ready()

    run.forget()

    assert run.task_name == 'nomad.processing.close_upload'
    assert run.cause is None
    assert run.status == 'SUCCESS'


def test_process_non_existing():
    run = UploadProcessing('__does_not_exist')
    run.start()
    run.get(timeout=30)
    assert run.ready()
    run.forget()

    assert run.task_name == 'nomad.processing.open_upload'
    assert run.status == 'SUCCESS'
    assert run.cause is not None
    assert isinstance(run.cause, KeyError)
