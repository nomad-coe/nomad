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


@pytest.fixture(scope='session')
def purged_app(celery_session_app):
    """
    Purges all pending tasks of the celery app before test. This is necessary to
    remove tasks from the queue that might be 'left over' from prior tests.
    """
    celery_session_app.control.purge()
    yield celery_session_app


@pytest.fixture()
def patched_celery(monkeypatch):
    # There is a bug in celery, which prevents to use the celery_worker for multiple
    # tests: https://github.com/celery/celery/issues/4088
    # The bug has a fix from Aug 2018, but it is not yet released (TODO).
    # We monkeypatch a similar solution here.
    def add_reader(self, fds, callback, *args):
        from kombu.utils.eventio import ERR, READ, poll

        if self.poller is None:
            self.poller = poll()

        return self.add(fds, callback, READ | ERR, args)

    monkeypatch.setattr('kombu.asynchronous.hub.Hub.add_reader', add_reader)
    yield


@pytest.fixture(scope='session')
def celery_inspect(purged_app):
    yield purged_app.control.inspect()


@pytest.fixture()
def worker(patched_celery, celery_inspect, celery_session_worker):
    """
    Extension of the celery_session_worker fixture that ensures a clean task queue.
    """
    yield

    # wait until there no more active tasks, to leave clean worker and queues for the next
    # test.
    while True:
        empty = True
        for value in celery_inspect.active().values():
            empty = empty and len(value) == 0
        if empty:
            break


@pytest.fixture(scope='function')
def mockmongo(monkeypatch):
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


@pytest.fixture(scope='function')
def no_warn(caplog):
    yield caplog
    for record in caplog.records:
        if record.levelname in ['WARNING', 'ERROR', 'CRITICAL']:
            assert False, record.msg


@pytest.fixture(scope='function')
def one_error(caplog):
    yield caplog
    count = 0
    for record in caplog.records:
        if record.levelname in ['ERROR', 'CRITICAL']:
            count += 1
            if count > 1:
                assert False, "oo many errors"
