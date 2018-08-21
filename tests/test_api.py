import pytest
from threading import Thread
import subprocess
import shlex
import time
import json
from mongoengine import connect
from mongoengine.connection import disconnect
from minio.error import ResponseError

from nomad import config, api, files

from tests.test_processing import example_files


@pytest.fixture
def client():
    disconnect()
    connect('users_test', host=config.mongo.host, is_mock=True)

    api.app.config['TESTING'] = True
    client = api.app.test_client()

    yield client


def assert_uploads(upload_json_str, count=0, **kwargs):
    data = json.loads(upload_json_str)
    assert isinstance(data, list)
    assert len(data) == count

    if count > 0:
        assert_upload(json.dumps(data[0]), **kwargs)


def assert_upload(upload_json_str, id=None):
    data = json.loads(upload_json_str)
    assert 'id' in data
    if id is not None:
        assert id == data['id']
    assert 'create_time' in data
    assert 'presigned_url' in data

    return data


def test_no_uploads(client):
    rv = client.get('/uploads')

    assert rv.status_code == 200
    assert_uploads(rv.data, count=0)


def test_bad_upload_id(client):
    rv = client.get('/uploads/bad_id')
    assert rv.status_code == 400


def test_not_existing_upload(client):
    rv = client.get('/uploads/123456789012123456789012')
    assert rv.status_code == 404


def test_create_upload(client):
    rv = client.post('/uploads')

    assert rv.status_code == 200
    upload_id = assert_upload(rv.data)['id']

    rv = client.get('/uploads/%s' % upload_id)
    assert rv.status_code == 200
    assert_upload(rv.data, id=upload_id)

    rv = client.get('/uploads')
    assert rv.status_code == 200
    assert_uploads(rv.data, count=1, id=upload_id)


@pytest.mark.parametrize("file", example_files)
@pytest.mark.timeout(10)
def test_upload_to_upload(client, file):
    rv = client.post('/uploads')
    assert rv.status_code == 200
    upload = assert_upload(rv.data)

    @files.upload_put_handler
    def handle_upload_put(received_upload_id: str):
        assert upload['id'] == received_upload_id
        raise StopIteration

    def handle_uploads():
        handle_upload_put(received_upload_id='provided by decorator')

    handle_uploads_thread = Thread(target=handle_uploads)
    handle_uploads_thread.start()

    time.sleep(1)
    upload_url = upload['presigned_url']
    cmd = files.create_curl_upload_cmd(upload_url).replace('<ZIPFILE>', file)
    subprocess.call(shlex.split(cmd))

    handle_uploads_thread.join()

    try:
        files._client.remove_object(config.files.uploads_bucket, upload['id'])
    except ResponseError:
        assert False


@pytest.mark.parametrize("file", example_files)
@pytest.mark.timeout(10)
def test_processing(client, file):
    handle_uploads_thread = api.start_upload_handler(quit=True)

    rv = client.post('/uploads')
    assert rv.status_code == 200
    upload = assert_upload(rv.data)

    time.sleep(1)
    upload_url = upload['presigned_url']
    cmd = files.create_curl_upload_cmd(upload_url).replace('<ZIPFILE>', file)
    subprocess.call(shlex.split(cmd))

    handle_uploads_thread.join()

    while True:
        time.sleep(1)

        rv = client.get('/uploads/%s' % upload['id'])
        assert rv.status_code == 200
        upload = assert_upload(rv.data)
        assert 'upload_time' in upload
        assert 'processing' in upload

        if upload['processing']['status'] in ['SUCCESS', 'FAILURE']:
            break

    assert upload['processing']['status'] == 'SUCCESS'

    try:
        files._client.remove_object(config.files.uploads_bucket, upload['id'])
    except ResponseError:
        assert False
