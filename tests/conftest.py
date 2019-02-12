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

from typing import Tuple
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
import shutil
import os.path
import datetime
import base64
from bravado.client import SwaggerClient

from nomad import config, infrastructure, files, parsing, processing, coe_repo, api

from tests import test_parsing, test_normalizing
from tests.processing import test_data as test_processing
from tests.test_files import example_file, empty_file
from tests.bravado_flask import FlaskTestHttpClient

example_files = [empty_file, example_file]


@pytest.fixture(scope="session")
def monkeysession(request):
    from _pytest.monkeypatch import MonkeyPatch
    mpatch = MonkeyPatch()
    yield mpatch
    mpatch.undo()


@pytest.fixture(scope='session', autouse=True)
def nomad_logging():
    config.logstash = config.logstash._replace(enabled=False)
    config.console_log_level = logging.CRITICAL
    infrastructure.setup_logging()


@pytest.fixture(scope='session', autouse=True)
def raw_files_infra(monkeysession):
    monkeysession.setattr('nomad.config.fs', config.FSConfig(
        tmp='.volumes/test_fs/tmp', objects='.volumes/test_fs/objects'))


@pytest.fixture(scope='function')
def raw_files(raw_files_infra):
    """ Provides cleaned out files directory structure per function. Clears files after test. """
    try:
        yield
    finally:
        try:
            shutil.rmtree(config.fs.objects)
        except FileNotFoundError:
            pass
        try:
            shutil.rmtree(config.fs.tmp)
        except FileNotFoundError:
            pass


@pytest.fixture(scope='function')
def client(monkeysession):
    api.app.config['TESTING'] = True
    client = api.app.test_client()

    yield client


@pytest.fixture(scope='session')
def celery_includes():
    return ['nomad.processing.base']


@pytest.fixture(scope='session')
def celery_config():
    return {
        'broker_url': config.celery.broker_url
    }


@pytest.fixture(scope='session')
def worker(celery_session_worker):
    """ Provides a clean worker (no old tasks) per function. Waits for all tasks to be completed. """
    pass


@pytest.fixture(scope='function')
def mongo(monkeypatch):
    """ Provides a cleaned mocked mongo per function. """

    disconnect()
    connection = connect('test_db', host='mongomock://localhost')
    monkeypatch.setattr('nomad.infrastructure.setup_mongo', lambda **kwargs: None)

    yield

    connection.drop_database('test_db')


@pytest.fixture(scope='session')
def elastic_infra(monkeysession):
    """ Provides elastic infrastructure to the session """
    monkeysession.setattr('nomad.config.elastic', config.elastic._replace(index_name='test_nomad_fairdi_calcs'))
    return infrastructure.setup_elastic()


@pytest.fixture(scope='function')
def elastic(elastic_infra):
    """ Provides a clean elastic per function. Clears elastic before test. """
    elastic_infra.delete_by_query(
        index='test_nomad_fairdi_calcs', body=dict(query=dict(match_all={})),
        wait_for_completion=True, refresh=True)
    assert infrastructure.elastic_client is not None
    return elastic_infra


@contextmanager
def create_postgres_infra(monkeysession=None, **kwargs):
    """
    A generator that sets up and tears down a test db and monkeypatches it to the
    respective global infrastructure variables.
    """
    db_args = dict(dbname='test_nomad_fairdi_repo_db')
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
def postgres_infra(monkeysession):
    """ Provides a clean coe repository db per module """
    with create_postgres_infra(monkeysession, exists=False) as db:
        yield db


@pytest.fixture(scope='function')
def proc_infra(postgres, elastic, mongo, worker):
    """ Combines all fixtures necessary for processing (postgres, elastic, worker, files, mongo) """
    return dict(
        postgres=postgres,
        elastic=elastic)


@pytest.fixture(scope='function')
def expandable_postgres(monkeysession, postgres_infra):
    """ Provides a coe repository db that can be deleted during test """
    with create_postgres_infra(monkeysession, dbname='test_nomad_fairdi_expandable_repo_db', exists=False) as db:
        yield db


@pytest.fixture(scope='function')
def postgres(postgres_infra):
    """ Provides a clean coe repository db per function. Clears db before test. """
    # do not wonder, this will not setback the id counters
    postgres_infra.execute('TRUNCATE uploads CASCADE;')
    yield postgres_infra


@pytest.fixture(scope='module')
def test_user(postgres_infra):
    from nomad import coe_repo
    return coe_repo.ensure_test_user(email='sheldon.cooper@nomad-fairdi.tests.de')


@pytest.fixture(scope='module')
def other_test_user(postgres_infra):
    from nomad import coe_repo
    return coe_repo.ensure_test_user(email='leonard.hofstadter@nomad-fairdi.tests.de')


@pytest.fixture(scope='module')
def admin_user(postgres_infra):
    from nomad import coe_repo
    return coe_repo.admin_user()


def create_auth_headers(user):
    basic_auth_str = '%s:password' % user.email
    basic_auth_bytes = basic_auth_str.encode('utf-8')
    basic_auth_base64 = base64.b64encode(basic_auth_bytes).decode('utf-8')
    return {
        'Authorization': 'Basic %s' % basic_auth_base64
    }


@pytest.fixture(scope='module')
def test_user_auth(test_user: coe_repo.User):
    return create_auth_headers(test_user)


@pytest.fixture(scope='module')
def test_other_user_auth(other_test_user: coe_repo.User):
    return create_auth_headers(other_test_user)


@pytest.fixture(scope='module')
def admin_user_auth(admin_user: coe_repo.User):
    return create_auth_headers(admin_user)


@pytest.fixture(scope='function')
def bravado(client, postgres, test_user_auth):
    http_client = FlaskTestHttpClient(client, headers=test_user_auth)
    return SwaggerClient.from_url('/swagger.json', http_client=http_client)


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


@pytest.fixture(scope='function')
def with_warn(caplog):
    yield caplog
    count = 0
    for record in caplog.get_records(when='call'):
        if record.levelname in ['WARNING']:
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
        self.smtp = None

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
        if self.smtp is not None:
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
def mails(smtpd, monkeypatch):
    smtpd.clear()

    old_config = config.mail
    new_config = config.MailConfig(
        'localhost',
        old_config.port,
        old_config.user,
        old_config.password,
        old_config.from_address)

    monkeypatch.setattr('nomad.config.mail', new_config)
    yield smtpd


@pytest.fixture(scope='session')
def example_mainfile() -> Tuple[str, str]:
    return ('parsers/template', 'tests/data/parsers/template.json')


@pytest.fixture(scope='session', params=example_files)
def example_upload(request) -> str:
    return request.param


@pytest.fixture(scope='module')
def example_user_metadata(other_test_user, test_user) -> dict:
    return {
        'comment': 'test comment',
        'with_embargo': True,
        'references': ['http://external.ref/one', 'http://external.ref/two'],
        '_uploader': other_test_user.user_id,
        'coauthors': [test_user.user_id],
        '_upload_time': datetime.datetime.now(),
        '_pid': 256
    }


@pytest.fixture(scope='function')
def parsed(example_mainfile: Tuple[str, str]) -> parsing.LocalBackend:
    """ Provides a parsed calculation in the form of a LocalBackend. """
    parser, mainfile = example_mainfile
    return test_parsing.run_parser(parser, mainfile)


@pytest.fixture(scope='function')
def normalized(parsed: parsing.LocalBackend) -> parsing.LocalBackend:
    """ Provides a normalized calculation in the form of a LocalBackend. """
    return test_normalizing.run_normalize(parsed)


@pytest.fixture(scope='function')
def uploaded(example_upload: str, raw_files) -> str:
    """
    Provides a uploaded with uploaded example file and gives the upload_id.
    Clears files after test.
    """
    example_upload_id = os.path.basename(example_upload).replace('.zip', '')
    upload_files = files.ArchiveBasedStagingUploadFiles(example_upload_id, create=True)
    shutil.copyfile(example_upload, upload_files.upload_file_os_path)

    return example_upload_id


@pytest.mark.timeout(10)
@pytest.fixture(scope='function')
def processed(uploaded: str, test_user: coe_repo.User, proc_infra) -> processing.Upload:
    """
    Provides a processed upload. Upload was uploaded with test_user.
    """
    return test_processing.run_processing(uploaded, test_user)
