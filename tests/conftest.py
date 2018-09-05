import pytest
from mongoengine import connect
from mongoengine.connection import disconnect

from nomad import config


@pytest.fixture(scope='session')
def celery_includes():
    return ['nomad.processing.base']


@pytest.fixture(scope='session')
def celery_config():
    return {
        'broker_url': config.celery.broker_url,
        'result_backend': config.celery.backend_url,
        'accept_content': ['json', 'pickle'],
        'task_serializer': config.celery.serializer,
        'result_serializer': config.celery.serializer
    }


@pytest.fixture(scope='function', autouse=True)
def mongomock(monkeypatch):
    def mock_connect(**kwargs):
        return connect('test_db', host='mongomock://localhost')

    disconnect()
    connection = mock_connect()
    monkeypatch.setattr('nomad.processing.base.mongo_connect', mock_connect)
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

    monkeypatch.setattr('nomad.search.CalcElasticDocument.create_from_backend', create_from_backend)
    monkeypatch.setattr('nomad.search.CalcElasticDocument.upload_exists', upload_exists)
