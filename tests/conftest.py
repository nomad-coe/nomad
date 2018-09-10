import pytest
from mongoengine import connect
from mongoengine.connection import disconnect
import time

from nomad import config


@pytest.fixture(scope='session')
def celery_includes():
    return ['nomad.processing.base']


@pytest.fixture(scope='session')
def celery_config():
    return {
        'broker_url': config.celery.broker_url
    }


@pytest.fixture(scope='function')
def worker(celery_session_worker):
    """
    Extension of the buildin celery_session_worker fixture that adds sleep to consume
    bleeding tasks.
    Processes might be completed (and therefore the test it self) before child
    processes are finished. Therefore open task request might bleed into the next test.
    """
    yield
    time.sleep(0.2)


@pytest.fixture(scope='function', autouse=True)
def mongomock(monkeypatch):
    def mock_connect(**kwargs):
        return connect('test_db', host='mongomock://localhost')

    disconnect()
    connection = mock_connect()
    monkeypatch.setattr('nomad.processing.base.mongo_connect', mock_connect)

    from nomad.user import ensure_test_users
    ensure_test_users()

    yield

    connection.drop_database('test_db')


@pytest.fixture(scope='function')
def mocksearch(monkeypatch):
    uploads = []

    def create_from_backend(_, **kwargs):
        upload_hash = kwargs.get('upload_hash', None)
        uploads.append(upload_hash)
        return {}

    def upload_exists(upload_hash):
        return upload_hash in uploads

    monkeypatch.setattr('nomad.repo.RepoCalc.create_from_backend', create_from_backend)
    monkeypatch.setattr('nomad.repo.RepoCalc.upload_exists', upload_exists)
