import pytest
import logging
from sqlalchemy.orm import Session
from mongoengine import connect
from mongoengine.connection import disconnect
from contextlib import contextmanager
from collections import namedtuple
from smtpd import SMTPServer
from threading import Lock, Thread
import asyncore
import time
import pytest

from nomad import config, infrastructure


@pytest.fixture(scope="session")
def monkeysession(request):
    from _pytest.monkeypatch import MonkeyPatch
    mpatch = MonkeyPatch()
    yield mpatch
    mpatch.undo()


@pytest.fixture(scope='session', autouse=True)
def nomad_files(monkeysession):
    monkeysession.setattr('nomad.config.fs', config.FSConfig(
        tmp='.volumes/test_fs/tmp', objects='.volumes/test_fs/objects'))


@pytest.fixture(scope='session', autouse=True)
def nomad_logging():
    config.logstash = config.logstash._replace(enabled=False)
    config.console_log_level = logging.CRITICAL
    infrastructure.setup_logging()


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

    disconnect()
    connection = connect('test_db', host='mongomock://localhost')
    monkeypatch.setattr('nomad.infrastructure.setup_mongo', lambda **kwargs: None)

    yield

    connection.drop_database('test_db')


@pytest.fixture(scope='session')
def elastic():
    infrastructure.setup_elastic()
    assert infrastructure.elastic_client is not None


@contextmanager
def create_repository_db(monkeysession=None, **kwargs):
    """
    A generator that sets up and tears down a test db and monkeypatches it to the
    respective global infrastructure variables.
    """
    db_args = dict(dbname='test_nomad_fair_repo_db')
    db_args.update(**kwargs)

    old_config = config.repository_db
    new_config = config.RepositoryDBConfig(
        old_config.host,
        old_config.port,
        db_args.get('dbname'),
        old_config.user,
        old_config.password)

    if monkeysession is not None:
        monkeysession.setattr('nomad.config.repository_db', new_config)

    connection, _ = infrastructure.sqlalchemy_repository_db(**db_args)
    assert connection is not None

    # we use a transaction around the session to rollback anything that happens within
    # test execution
    trans = connection.begin()
    db = Session(bind=connection, autocommit=True)

    old_connection, old_db = None, None
    if monkeysession is not None:
        from nomad.infrastructure import repository_db_conn, repository_db
        old_connection, old_db = repository_db_conn, repository_db
        monkeysession.setattr('nomad.infrastructure.repository_db_conn', connection)
        monkeysession.setattr('nomad.infrastructure.repository_db', db)

    yield db

    if monkeysession is not None:
        monkeysession.setattr('nomad.infrastructure.repository_db_conn', old_connection)
        monkeysession.setattr('nomad.infrastructure.repository_db', old_db)
        monkeysession.setattr('nomad.config.repository_db', old_config)

    trans.rollback()
    db.expunge_all()
    db.invalidate()
    db.close_all()

    connection.close()
    connection.engine.dispose()


@pytest.fixture(scope='module')
def repository_db(monkeysession):
    with create_repository_db(monkeysession, exists=False) as db:
        yield db


@pytest.fixture(scope='function')
def expandable_repo_db(monkeysession, repository_db):
    with create_repository_db(monkeysession, dbname='test_nomad_fair_expandable_repo_db', exists=False) as db:
        yield db


@pytest.fixture(scope='function')
def clean_repository_db(repository_db):
    # do not wonder, this will not setback the id counters
    repository_db.execute('TRUNCATE uploads CASCADE;')
    yield repository_db


@pytest.fixture(scope='module')
def test_user(repository_db):
    from nomad import coe_repo
    return coe_repo.ensure_test_user(email='sheldon.cooper@nomad-fairdi.tests.de')


@pytest.fixture(scope='module')
def other_test_user(repository_db):
    from nomad import coe_repo
    return coe_repo.ensure_test_user(email='leonard.hofstadter@nomad-fairdi.tests.de')


@pytest.fixture(scope='module')
def admin_user(repository_db):
    from nomad import coe_repo
    return coe_repo.admin_user()


@pytest.fixture(scope='function')
def no_warn(caplog):
    yield caplog
    for record in caplog.get_records(when='call'):
        if record.levelname in ['WARNING', 'ERROR', 'CRITICAL']:
            assert False, record.msg


@pytest.fixture(scope='function')
def with_error(caplog):
    yield caplog
    count = 0
    for record in caplog.get_records(when='call'):
        if record.levelname in ['ERROR', 'CRITICAL']:
            count += 1

    assert count > 0


"""
Fixture for mocked SMTP server for testing.
Based on https://gist.github.com/akheron/cf3863cdc424f08929e4cb7dc365ef23.
"""

RecordedMessage = namedtuple(
    'RecordedMessage',
    'peer envelope_from envelope_recipients data',
)


class ThreadSafeList:
    def __init__(self, *args, **kwds):
        self._items = []
        self._lock = Lock()

    def clear(self):
        with self._lock:
            self._items = []

    def add(self, item):
        with self._lock:
            self._items.append(item)

    def copy(self):
        with self._lock:
            return self._items[:]


class SMTPServerThread(Thread):
    def __init__(self, messages):
        super().__init__()
        self.messages = messages
        self.host_port = None

    def run(self):
        _messages = self.messages

        class _SMTPServer(SMTPServer):
            def process_message(self, peer, mailfrom, rcpttos, data, **kwargs):
                msg = RecordedMessage(peer, mailfrom, rcpttos, data)
                _messages.add(msg)

        self.smtp = _SMTPServer(('127.0.0.1', config.mail.port), None)
        self.host_port = self.smtp.socket.getsockname()
        try:
            asyncore.loop(1)
        except Exception:
            pass

    def close(self):
        self.smtp.close()


class SMTPServerFixture:
    def __init__(self):
        self._messages = ThreadSafeList()
        self._thread = SMTPServerThread(self._messages)
        self._thread.start()

    @property
    def host_port(self):
        '''SMTP server's listening address as a (host, port) tuple'''
        while self._thread.host_port is None:
            time.sleep(0.1)
        return self._thread.host_port

    @property
    def host(self):
        return self.host_port[0]

    @property
    def port(self):
        return self.host_port[1]

    @property
    def messages(self):
        '''A list of RecordedMessage objects'''
        return self._messages.copy()

    def clear(self):
        self._messages.clear()

    def close(self):
        self._thread.close()
        self._thread.join(1)


@pytest.fixture(scope='session')
def smtpd(request):
    fixture = SMTPServerFixture()
    request.addfinalizer(fixture.close)
    return fixture


@pytest.fixture(scope='function', autouse=True)
def mails(smtpd):
    smtpd.clear()
    yield smtpd
