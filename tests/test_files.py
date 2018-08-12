from typing import Generator
import pytest
import time
from minio import ResponseError
from threading import Thread
import subprocess
import shlex

import nomad.files as files
import nomad.config as config


example_file = './data/examples_vasp.zip'
example_upload_id = '__test_upload_id'


@pytest.fixture
def uploaded_id() -> str:
    files._client.fput_object(config.s3.uploads_bucket, example_upload_id, example_file)
    return example_upload_id


@pytest.fixture(autouse=True)
def tear_down():
    yield
    try:
        files._client.remove_object(config.s3.uploads_bucket, example_upload_id)
    except ResponseError:
        pass


def test_presigned_url():
    url = files.get_presigned_upload_url(example_upload_id)
    assert url is not None
    assert isinstance(url, str)

    upload_url = files.get_presigned_upload_url(example_upload_id)
    cmd = files.create_curl_upload_cmd(upload_url).replace('<ZIPFILE>', example_file)
    subprocess.call(shlex.split(cmd))

    stat = files._client.stat_object(config.s3.uploads_bucket, example_upload_id)
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
def test_upload_notification():
    @files.upload_put_handler
    def handle_upload_put(upload_id: str):
        assert upload_id == example_upload_id
        raise StopIteration

    def handle_uploads():
        handle_upload_put(upload_id='provided by decorator')

    handle_uploads_thread = Thread(target=handle_uploads)
    handle_uploads_thread.start()

    time.sleep(1)
    test_presigned_url()

    handle_uploads_thread.join()
