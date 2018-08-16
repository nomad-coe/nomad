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
from minio.error import NoSuchBucket

import nomad.config as config
import nomad.files as files
from nomad.processing import UploadProcessing

example_files = ['data/examples_vasp.zip', 'data/empty.zip']


@pytest.fixture(scope='function', params=example_files)
def uploaded_id(request) -> Generator[str, None, None]:
    example_file = request.param
    example_upload_id = example_file.replace('.zip', '')
    files._client.fput_object(config.files.uploads_bucket, example_upload_id, example_file)

    yield example_upload_id

    try:
        # remove the created uploads
        files._client.remove_object(config.files.uploads_bucket, example_upload_id)

        # remove all the created archive files
        archive_objs = files._client.list_objects(config.files.archive_bucket, recursive=True)
        errors = files._client.remove_objects(config.files.archive_bucket, [obj.object_name for obj in archive_objs])
        # you have to walk the iterator for minio to work (?!)
        for _ in errors:
            pass
    except ResponseError:
        pass
    except NoSuchBucket:
        pass


def test_processing(uploaded_id):
    run = UploadProcessing(uploaded_id)
    run.start()
    run.get(timeout=30)

    assert run.ready()
    assert run.task_name == 'nomad.processing.close_upload'
    assert run.upload_hash is not None
    assert run.cause is None
    assert run.status == 'SUCCESS'
    for status, errors, archive_id in run.parse_results:
        assert errors is None or len(errors) == 0
        assert status == 'ParseSuccess'
        assert archive_id.startswith(run.upload_hash)

    run.forget()


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
