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
        'broker_url': config.celery.broker_url
    }


@pytest.fixture(scope='function')
def purged_queue(celery_app):
    """
    Purges all pending tasks of the celery app before test. This is necessary to
    remove tasks from the queue that might be 'left over' from prior tests.
    """
    celery_app.control.purge()
    yield


@pytest.fixture(scope='function')
def patched_celery(monkeypatch):
    # There is a bug in celery, which prevents to use the celery_worker for multiple
    # tests: https://github.com/celery/celery/issues/4088
    # The bug has a fix from Aug 2018, but it is not yet released (TODO).
    # We monkeypatch a similar solution here.
    def add_reader(self, fds, callback, *args):
        from kombu.utils.eventio import ERR, READ, WRITE, poll

        if self.poller is None:
            self.poller = poll()

        return self.add(fds, callback, READ | ERR, args)

    monkeypatch.setattr('kombu.asynchronous.hub.Hub.add_reader', add_reader)
    yield


@pytest.fixture(scope='function')
def worker(patched_celery, purged_queue, celery_worker):
    """
    Extension of the celery_worker fixture that ensures a clean task queue before yielding.
    """
    # This wont work with the session_worker, it will already have old/unexecuted tasks
    # taken from the queue and might resubmit them. Therefore, purging the queue won't
    # help much.
    yield


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

    monkeypatch.setattr('nomad.repo.RepoCalc.create_from_backend', create_from_backend)
    monkeypatch.setattr('nomad.repo.RepoCalc.upload_exists', upload_exists)
