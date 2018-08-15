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

import pytest
from minio import ResponseError
from threading import Thread
import subprocess
import shlex
import time
from typing import Generator
import json

import nomad.files as files
import nomad.config as config

example_file = './data/examples_vasp.zip'


@pytest.fixture
def uploaded_id() -> Generator[str, None, None]:
    example_upload_id = '__test_upload_id'

    files._client.fput_object(config.files.uploads_bucket, example_upload_id, example_file)
    yield example_upload_id
    try:
        files._client.remove_object(config.files.uploads_bucket, example_upload_id)
    except ResponseError:
        pass


@pytest.fixture
def upload_id() -> Generator[str, None, None]:
    example_upload_id = '__test_upload_id'
    yield example_upload_id
    try:
        files._client.remove_object(config.files.uploads_bucket, example_upload_id)
    except ResponseError:
        pass


@pytest.fixture
def archive_id() -> Generator[str, None, None]:
    example_archive_id = '__test_archive_id'
    yield example_archive_id
    try:
        files._client.remove_object(config.files.archive_bucket, example_archive_id)
    except ResponseError:
        pass


def test_presigned_url(upload_id):
    url = files.get_presigned_upload_url(upload_id)
    assert url is not None
    assert isinstance(url, str)

    upload_url = files.get_presigned_upload_url(upload_id)
    cmd = files.create_curl_upload_cmd(upload_url).replace('<ZIPFILE>', example_file)
    subprocess.call(shlex.split(cmd))

    stat = files._client.stat_object(config.files.uploads_bucket, upload_id)
    assert stat.content_type.startswith('application/octet-steam')


def test_upload(uploaded_id: str):
    with files.Upload(uploaded_id) as upload:
        assert len(upload.filelist) == 106
        # now just try to open the first file (not directory), without error
        for filename in upload.filelist:
            if filename.endswith('.xml'):
                upload.open_file(filename).close()
                break


@pytest.mark.timeout(10)
def test_upload_notification(upload_id):
    @files.upload_put_handler
    def handle_upload_put(received_upload_id: str):
        assert upload_id == received_upload_id
        raise StopIteration

    def handle_uploads():
        handle_upload_put(received_upload_id='provided by decorator')

    handle_uploads_thread = Thread(target=handle_uploads)
    handle_uploads_thread.start()

    time.sleep(1)
    test_presigned_url(upload_id)

    handle_uploads_thread.join()


def test_metadata(uploaded_id: str):
    with files.Upload(uploaded_id) as upload:
        assert upload.metadata is not None


def test_hash(uploaded_id: str):
    with files.Upload(uploaded_id) as upload:
        hash = upload.hash()
        assert hash is not None
        assert isinstance(hash, str)

        print(hash)


def test_archive(archive_id: str):
    with files.write_archive_json(archive_id) as out:
        json.dump({'test': 'value'}, out)

    result = json.load(files.open_archive_json(archive_id))

    assert 'test' in result
    assert result['test'] == 'value'
