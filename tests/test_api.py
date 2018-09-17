import pytest
from threading import Thread
import subprocess
import shlex
import time
import json
from mongoengine import connect
from mongoengine.connection import disconnect
from datetime import datetime, timedelta
import base64

from nomad import config
# for convinience we test the api without path prefix
services_config = config.services._asdict()
services_config.update(api_base_path='')
config.services = config.NomadServicesConfig(**services_config)

from nomad import api, files  # noqa
from nomad.processing import Upload, handle_uploads_thread  # noqa

from tests.processing.test_data import example_files  # noqa
from tests.test_files import assert_exists  # noqa

# import fixtures
from tests.test_files import clear_files, archive_id  # noqa pylint: disable=unused-import
from tests.test_normalizing import normalized_vasp_example  # noqa pylint: disable=unused-import
from tests.test_parsing import parsed_vasp_example  # noqa pylint: disable=unused-import
from tests.test_repo import example_elastic_calc  # noqa pylint: disable=unused-import


@pytest.fixture(scope='function')
def client():
    disconnect()
    connect('users_test', host=config.mongo.host, port=config.mongo.port, is_mock=True)

    api.app.config['TESTING'] = True
    client = api.app.test_client()

    yield client
    Upload._get_collection().drop()


@pytest.fixture(scope='session')
def test_user_auth():
    return {
        'Authorization': 'Basic %s' % base64.b64encode(b'me@gmail.com:nomad').decode('utf-8')
    }


@pytest.fixture(scope='session')
def test_other_user_auth():
    return {
        'Authorization': 'Basic %s' % base64.b64encode(b'other@gmail.com:nomad').decode('utf-8')
    }


def assert_uploads(upload_json_str, count=0, **kwargs):
    data = json.loads(upload_json_str)
    assert isinstance(data, list)
    assert len(data) == count

    if count > 0:
        assert_upload(json.dumps(data[0]), **kwargs)


def assert_upload(upload_json_str, id=None, **kwargs):
    data = json.loads(upload_json_str)
    assert 'upload_id' in data
    if id is not None:
        assert id == data['upload_id']
    assert 'create_time' in data
    assert 'presigned_url' in data
    assert 'upload_command' in data

    for key, value in kwargs.items():
        assert data.get(key, None) == value

    return data


def test_no_uploads(client, test_user_auth):
    rv = client.get('/uploads', headers=test_user_auth)

    assert rv.status_code == 200
    assert_uploads(rv.data, count=0)


def test_not_existing_upload(client, test_user_auth):
    rv = client.get('/uploads/123456789012123456789012', headers=test_user_auth)
    assert rv.status_code == 404


def test_stale_upload(client, test_user_auth):
    rv = client.post(
        '/uploads',
        headers=test_user_auth,
        data=json.dumps(dict(name='test_name')),
        content_type='application/json')
    assert rv.status_code == 200
    upload_id = assert_upload(rv.data)['upload_id']

    upload = Upload.get(upload_id)
    upload.create_time = datetime.now() - timedelta(days=2)
    upload.save()

    rv = client.get('/uploads/%s' % upload_id, headers=test_user_auth)
    assert rv.status_code == 200
    assert_upload(rv.data, is_stale=True)


def test_create_upload(client, test_user_auth):
    rv = client.post('/uploads', headers=test_user_auth)

    assert rv.status_code == 200
    upload_id = assert_upload(rv.data)['upload_id']

    rv = client.get('/uploads/%s' % upload_id, headers=test_user_auth)
    assert rv.status_code == 200
    assert_upload(rv.data, id=upload_id, is_stale=False)

    rv = client.get('/uploads', headers=test_user_auth)
    assert rv.status_code == 200
    assert_uploads(rv.data, count=1, id=upload_id)


def test_create_upload_with_name(client, test_user_auth):
    rv = client.post(
        '/uploads', headers=test_user_auth,
        data=json.dumps(dict(name='test_name')), content_type='application/json')
    assert rv.status_code == 200
    upload = assert_upload(rv.data)
    assert upload['name'] == 'test_name'


def test_delete_empty_upload(client, test_user_auth):
    rv = client.post('/uploads', headers=test_user_auth)

    assert rv.status_code == 200
    upload_id = assert_upload(rv.data)['upload_id']

    rv = client.delete('/uploads/%s' % upload_id, headers=test_user_auth)
    assert rv.status_code == 200

    rv = client.get('/uploads/%s' % upload_id, headers=test_user_auth)
    assert rv.status_code == 404


@pytest.mark.parametrize("file", example_files)
@pytest.mark.timeout(30)
def test_upload_to_upload(client, file, test_user_auth):
    rv = client.post('/uploads', headers=test_user_auth)
    assert rv.status_code == 200
    upload = assert_upload(rv.data)

    @files.upload_put_handler
    def handle_upload_put(received_upload_id: str):
        assert upload['upload_id'] == received_upload_id
        raise StopIteration

    def handle_uploads():
        handle_upload_put(received_upload_id='provided by decorator')

    handle_uploads_thread = Thread(target=handle_uploads)
    handle_uploads_thread.start()

    time.sleep(0.1)
    upload_url = upload['presigned_url']
    cmd = files.create_curl_upload_cmd(upload_url).replace('<ZIPFILE>', file)
    subprocess.call(shlex.split(cmd))

    handle_uploads_thread.join()

    upload_id = upload['upload_id']
    assert_exists(config.files.uploads_bucket, upload_id)


@pytest.mark.parametrize("file", example_files)
@pytest.mark.timeout(10)
def test_processing(client, file, worker, mocksearch, test_user_auth):
    handler = handle_uploads_thread(quit=True)

    rv = client.post('/uploads', headers=test_user_auth)
    assert rv.status_code == 200
    upload = assert_upload(rv.data)

    time.sleep(0.1)
    upload_url = upload['presigned_url']
    cmd = files.create_curl_upload_cmd(upload_url).replace('<ZIPFILE>', file)
    subprocess.call(shlex.split(cmd))

    handler.join()

    while True:
        time.sleep(1)

        rv = client.get('/uploads/%s' % upload['upload_id'], headers=test_user_auth)
        assert rv.status_code == 200
        upload = assert_upload(rv.data)
        assert 'upload_time' in upload
        if upload['completed']:
            break

    assert len(upload['tasks']) == 4
    assert upload['status'] == 'SUCCESS'
    assert upload['current_task'] == 'cleanup'
    calcs = upload['calcs']['results']
    for calc in calcs:
        assert calc['status'] == 'SUCCESS'
        assert calc['current_task'] == 'archiving'
        assert len(calc['tasks']) == 3
        assert_exists(config.files.uploads_bucket, upload['upload_id'])

    if upload['calcs']['pagination']['total'] > 1:
        rv = client.get('/uploads/%s?page=2&per_page=1&order_by=status' % upload['upload_id'])
        assert rv.status_code == 200
        upload = assert_upload(rv.data)
        assert len(upload['calcs']['results']) == 1

    rv = client.post(
        '/uploads/%s' % upload['upload_id'],
        headers=test_user_auth,
        data=json.dumps(dict(operation='unstage')),
        content_type='application/json')
    assert rv.status_code == 200

    rv = client.get('/uploads', headers=test_user_auth)
    assert rv.status_code == 200
    assert_uploads(rv.data, count=0)


def test_repo_calc(client, example_elastic_calc):
    rv = client.get(
        '/repo/%s/%s' % (example_elastic_calc.upload_hash, example_elastic_calc.calc_hash))
    assert rv.status_code == 200


def test_non_existing_repo_cals(client):
    rv = client.get('/repo/doesnt/exist')
    assert rv.status_code == 404


def test_repo_calcs(client, example_elastic_calc):
    rv = client.get('/repo')
    assert rv.status_code == 200
    data = json.loads(rv.data)
    results = data.get('results', None)
    assert results is not None
    assert isinstance(results, list)
    assert len(results) >= 1


def test_repo_calcs_pagination(client, example_elastic_calc):
    rv = client.get('/repo?page=1&per_page=1')
    assert rv.status_code == 200
    data = json.loads(rv.data)
    results = data.get('results', None)
    assert results is not None
    assert isinstance(results, list)
    assert len(results) == 1


def test_repo_calcs_user(client, example_elastic_calc, test_user_auth):
    rv = client.get('/repo?owner=user', headers=test_user_auth)
    assert rv.status_code == 200
    data = json.loads(rv.data)
    results = data.get('results', None)
    assert results is not None
    assert len(results) >= 1


def test_repo_calcs_user_authrequired(client, example_elastic_calc):
    rv = client.get('/repo?owner=user')
    assert rv.status_code == 401


def test_repo_calcs_user_invisible(client, example_elastic_calc, test_other_user_auth):
    rv = client.get('/repo?owner=user', headers=test_other_user_auth)
    assert rv.status_code == 200
    data = json.loads(rv.data)
    results = data.get('results', None)
    assert results is not None
    assert len(results) == 0


def test_get_archive(client, archive_id):
    rv = client.get('/archive/%s' % archive_id)
    assert rv.status_code == 302


def test_get_non_existing_archive(client):
    rv = client.get('/archive/%s' % 'doesnt/exist')
    assert rv.status_code == 404
