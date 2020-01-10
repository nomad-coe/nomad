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

from typing import Tuple, List
import pytest
import logging
from collections import namedtuple
from smtpd import SMTPServer
from threading import Lock, Thread
import asyncore
import time
import shutil
import os.path
import datetime
from bravado.client import SwaggerClient
from flask import request, g
import elasticsearch.exceptions
from typing import List
import numpy as np

from nomadcore.local_meta_info import loadJsonFile
import nomad_meta_info

from nomad import config, infrastructure, parsing, processing, app, search
from nomad.datamodel import User, CalcWithMetadata
from nomad.parsing import LocalBackend

from tests import test_parsing
from tests.normalizing.conftest import run_normalize
from tests.processing import test_data as test_processing
from tests.test_files import example_file, empty_file
from tests.bravado_flask import FlaskTestHttpClient

test_log_level = logging.CRITICAL
example_files = [empty_file, example_file]


@pytest.fixture(scope="session")
def monkeysession(request):
    from _pytest.monkeypatch import MonkeyPatch
    mpatch = MonkeyPatch()
    yield mpatch
    mpatch.undo()


@pytest.fixture(scope='session', autouse=True)
def nomad_logging():
    config.logstash.enabled = False
    config.console_log_level = test_log_level
    infrastructure.setup_logging()


@pytest.fixture(scope='session', autouse=True)
def raw_files_infra():
    config.fs.tmp = '.volumes/test_fs/tmp'
    config.fs.staging = '.volumes/test_fs/staging'
    config.fs.public = '.volumes/test_fs/public'
    config.fs.migration_packages = '.volumes/test_fs/migration_packages'
    config.fs.prefix_size = 2


@pytest.fixture(scope='function')
def raw_files(raw_files_infra):
    """ Provides cleaned out files directory structure per function. Clears files after test. """
    directories = [config.fs.staging, config.fs.public, config.fs.migration_packages, config.fs.tmp]
    for directory in directories:
        if not os.path.exists(directory):
            os.makedirs(directory)
    try:
        yield
    finally:
        for directory in directories:
            try:
                shutil.rmtree(directory)
            except FileNotFoundError:
                pass


@pytest.fixture(scope='session')
def session_client():
    app.app.config['TESTING'] = True
    client = app.app.test_client()

    yield client


@pytest.fixture(scope='function')
def client(mongo, session_client):
    app.app.config['TESTING'] = True
    client = app.app.test_client()

    yield client


@pytest.fixture(scope='session')
def celery_includes():
    return ['nomad.processing.base']


@pytest.fixture(scope='session')
def celery_config():
    return {
        'broker_url': config.rabbitmq_url(),
        'task_queue_max_priority': 10
    }


@pytest.fixture(scope='session')
def purged_app(celery_session_app):
    """
    Purges all pending tasks of the celery app before test. This is necessary to
    remove tasks from the queue that might be 'left over' from prior tests.
    """
    celery_session_app.control.purge()
    yield celery_session_app


@pytest.fixture(scope='session')
def celery_inspect(purged_app):
    yield purged_app.control.inspect()


# It might be necessary to make this a function scoped fixture, if old tasks keep
# 'bleeding' into successive tests.
@pytest.fixture(scope='function')
def worker(mongo, celery_session_worker, celery_inspect):
    """ Provides a clean worker (no old tasks) per function. Waits for all tasks to be completed. """
    yield

    # wait until there no more active tasks, to leave clean worker and queues for the next
    # test run.
    try:
        while True:
            empty = True
            for value in celery_inspect.active().values():
                empty = empty and len(value) == 0
            if empty:
                break
    except Exception:
        pass


@pytest.fixture(scope='session')
def mongo_infra(monkeysession):
    monkeysession.setattr('nomad.config.mongo.db_name', 'test_db')
    return infrastructure.setup_mongo()


@pytest.fixture(scope='function')
def mongo(mongo_infra):
    """ Provides a cleaned mocked mongo per function. """
    mongo_infra.drop_database('test_db')
    return mongo_infra


@pytest.fixture(scope='session')
def elastic_infra(monkeysession):
    """ Provides elastic infrastructure to the session """
    monkeysession.setattr('nomad.config.elastic.index_name', 'test_nomad_fairdi_0_6')
    try:
        return infrastructure.setup_elastic()
    except Exception:
        # try to delete index, error might be caused by changed mapping
        from elasticsearch_dsl import connections
        connections.create_connection(hosts=['%s:%d' % (config.elastic.host, config.elastic.port)]) \
            .indices.delete(index='test_nomad_fairdi_0_6')
        return infrastructure.setup_elastic()


def clear_elastic(elastic):
    try:
        elastic.delete_by_query(
            index='test_nomad_fairdi_0_6', body=dict(query=dict(match_all={})),
            wait_for_completion=True, refresh=True)
    except elasticsearch.exceptions.NotFoundError:
        # it is unclear why this happens, but it happens at least once, when all tests
        # are executed
        infrastructure.setup_elastic()


@pytest.fixture(scope='function')
def elastic(elastic_infra):
    """ Provides a clean elastic per function. Clears elastic before test. """
    clear_elastic(elastic_infra)

    assert infrastructure.elastic_client is not None
    return elastic_infra


def test_user_uuid(handle):
    return '00000000-0000-0000-0000-00000000000%d' % handle


test_users = {
    test_user_uuid(0): dict(email='admin', user_id=test_user_uuid(0)),
    test_user_uuid(1): dict(email='sheldon.cooper@nomad-coe.eu', first_name='Sheldon', last_name='Cooper', user_id=test_user_uuid(1)),
    test_user_uuid(2): dict(email='leonard.hofstadter@nomad-coe.eu', first_name='Leonard', last_name='Hofstadter', user_id=test_user_uuid(2))
}


class KeycloakMock:
    def __init__(self):
        self.id_counter = 2
        self.users = dict(**test_users)

    def authorize_flask(self, *args, **kwargs):
        if 'Authorization' in request.headers and request.headers['Authorization'].startswith('Bearer '):
            user_id = request.headers['Authorization'].split(None, 1)[1].strip()
            g.oidc_access_token = user_id
            g.user = User(**self.users[user_id])

    def add_user(self, user, *args, **kwargs):
        self.id_counter += 1
        user.user_id = test_user_uuid(self.id_counter)
        self.users[user.user_id] = dict(email=user.email, first_name=user.first_name, last_name=user.last_name, user_id=user.user_id)
        return None

    def get_user(self, user_id=None, email=None):
        if user_id is not None:
            return User(**self.users[user_id])
        elif email is not None:
            for user_id, user_values in self.users.items():
                if user_values['email'] == email:
                    return User(**user_values)
            raise KeyError('Only test user emails are recognized')
        else:
            assert False, 'no token based get_user during tests'

    def search_user(self, query):
        return [
            User(**test_user) for test_user in self.users.values()
            if query in ' '.join(test_user.values())]

    @property
    def access_token(self):
        return g.oidc_access_token


_keycloak = infrastructure.keycloak


# use a session fixture in addition to the function fixture, to ensure mocked keycloak
# before other class, module, etc. scoped function are run
@pytest.fixture(scope='session', autouse=True)
def mocked_keycloak_session(monkeysession):
    monkeysession.setattr('nomad.infrastructure.keycloak', KeycloakMock())


@pytest.fixture(scope='function', autouse=True)
def mocked_keycloak(monkeypatch):
    monkeypatch.setattr('nomad.infrastructure.keycloak', KeycloakMock())


@pytest.fixture(scope='function')
def keycloak(monkeypatch):
    monkeypatch.setattr('nomad.infrastructure.keycloak', _keycloak)


@pytest.fixture(scope='function')
def proc_infra(worker, elastic, mongo, raw_files):
    """ Combines all fixtures necessary for processing (elastic, worker, files, mongo) """
    return dict(elastic=elastic)


@pytest.fixture(scope='session')
def meta_info():
    file_dir = os.path.dirname(os.path.abspath(nomad_meta_info.__file__))
    path = os.path.join(file_dir, 'all.nomadmetainfo.json')
    meta_info, _ = loadJsonFile(path)
    return meta_info


@pytest.fixture(scope='module')
def test_user():
    return User(**test_users[test_user_uuid(1)])


@pytest.fixture(scope='module')
def other_test_user():
    return User(**test_users[test_user_uuid(2)])


@pytest.fixture(scope='module')
def admin_user():
    return User(**test_users[test_user_uuid(0)])


def create_auth_headers(user: User):
    return {
        'Authorization': 'Bearer %s' % user.user_id
    }


@pytest.fixture(scope='module')
def test_user_auth(test_user: User):
    return create_auth_headers(test_user)


@pytest.fixture(scope='module')
def other_test_user_auth(other_test_user: User):
    return create_auth_headers(other_test_user)


@pytest.fixture(scope='module')
def admin_user_auth(admin_user: User):
    return create_auth_headers(admin_user)


@pytest.fixture(scope='function')
def bravado(client, test_user_auth):
    http_client = FlaskTestHttpClient(client, headers=test_user_auth)
    return SwaggerClient.from_url('/api/swagger.json', http_client=http_client)


@pytest.fixture(scope='function')
def admin_user_bravado_client(client, admin_user_auth, monkeypatch):
    def create_client():
        http_client = FlaskTestHttpClient(client, headers=admin_user_auth)
        return SwaggerClient.from_url('/api/swagger.json', http_client=http_client)

    monkeypatch.setattr('nomad.cli.client.create_client', create_client)


@pytest.fixture(scope='function')
def test_user_bravado_client(client, test_user_auth, monkeypatch):
    def create_client():
        http_client = FlaskTestHttpClient(client, headers=test_user_auth)
        return SwaggerClient.from_url('/api/swagger.json', http_client=http_client)

    monkeypatch.setattr('nomad.cli.client.create_client', create_client)


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
def smtpd(request, monkeysession):
    # on some local machines resolving the local machine takes quit a while and
    # is irrelevant for testing
    monkeysession.setattr('socket.getfqdn', lambda *args, **kwargs: 'local.server')
    fixture = SMTPServerFixture()
    request.addfinalizer(fixture.close)
    return fixture


@pytest.fixture(scope='function', autouse=True)
def mails(smtpd, monkeypatch):
    smtpd.clear()
    monkeypatch.setattr('nomad.config.mail.enabled', True)
    monkeypatch.setattr('nomad.config.mail.host', 'localhost')
    yield smtpd


@pytest.fixture(scope='session')
def example_mainfile() -> Tuple[str, str]:
    return ('parsers/template', 'tests/data/parsers/template.json')


@pytest.fixture(scope='session', params=example_files)
def example_upload(request) -> str:
    return request.param


@pytest.fixture(scope='session')
def non_empty_example_upload():
    return example_file


@pytest.fixture(scope='session')
def empty_upload():
    return empty_file


@pytest.fixture(scope='module')
def example_user_metadata(other_test_user, test_user) -> dict:
    return {
        'comment': 'test comment',
        'with_embargo': True,
        'references': ['http://external.ref/one', 'http://external.ref/two'],
        '_uploader': other_test_user.user_id,
        'coauthors': [test_user.user_id],
        '_upload_time': datetime.datetime.utcnow(),
        '_pid': 256,
        'external_id': 'external_test_id'
    }


@pytest.fixture(scope='session')
def parsed(example_mainfile: Tuple[str, str]) -> parsing.LocalBackend:
    """ Provides a parsed calculation in the form of a LocalBackend. """
    parser, mainfile = example_mainfile
    return test_parsing.run_parser(parser, mainfile)


@pytest.fixture(scope='session')
def normalized(parsed: parsing.LocalBackend) -> parsing.LocalBackend:
    """ Provides a normalized calculation in the form of a LocalBackend. """
    return run_normalize(parsed)


@pytest.fixture(scope='function')
def uploaded(example_upload: str, raw_files) -> Tuple[str, str]:
    """
    Provides a uploaded with uploaded example file and gives the upload_id.
    Clears files after test.
    """
    example_upload_id = os.path.basename(example_upload).replace('.zip', '')
    return example_upload_id, example_upload


@pytest.fixture(scope='function')
def non_empty_uploaded(non_empty_example_upload: str, raw_files) -> Tuple[str, str]:
    example_upload_id = os.path.basename(non_empty_example_upload).replace('.zip', '')
    return example_upload_id, non_empty_example_upload


@pytest.mark.timeout(config.tests.default_timeout)
@pytest.fixture(scope='function')
def processed(uploaded: Tuple[str, str], test_user: User, proc_infra) -> processing.Upload:
    """
    Provides a processed upload. Upload was uploaded with test_user.
    """
    return test_processing.run_processing(uploaded, test_user)


@pytest.mark.timeout(config.tests.default_timeout)
@pytest.fixture(scope='function')
def processeds(non_empty_example_upload: str, test_user: User, proc_infra) -> List[processing.Upload]:
    result: List[processing.Upload] = []
    for i in range(2):
        upload_id = '%s_%d' % (os.path.basename(non_empty_example_upload).replace('.zip', ''), i)
        result.append(
            test_processing.run_processing((upload_id, non_empty_example_upload), test_user))

    return result


@pytest.mark.timeout(config.tests.default_timeout)
@pytest.fixture(scope='function')
def non_empty_processed(non_empty_uploaded: Tuple[str, str], test_user: User, proc_infra) -> processing.Upload:
    """
    Provides a processed upload. Upload was uploaded with test_user.
    """
    return test_processing.run_processing(non_empty_uploaded, test_user)


@pytest.mark.timeout(config.tests.default_timeout)
@pytest.fixture(scope='function')
def published(non_empty_processed: processing.Upload, example_user_metadata) -> processing.Upload:
    """
    Provides a processed upload. Upload was uploaded with test_user.
    """
    non_empty_processed.compress_and_set_metadata(example_user_metadata)
    non_empty_processed.publish_upload()
    try:
        non_empty_processed.block_until_complete(interval=.01)
    except Exception:
        pass

    return non_empty_processed


@pytest.mark.timeout(config.tests.default_timeout)
@pytest.fixture(scope='function')
def published_wo_user_metadata(non_empty_processed: processing.Upload) -> processing.Upload:
    """
    Provides a processed upload. Upload was uploaded with test_user.
    """
    non_empty_processed.publish_upload()
    try:
        non_empty_processed.block_until_complete(interval=.01)
    except Exception:
        pass

    return non_empty_processed


@pytest.fixture
def reset_config():
    """ Fixture that resets the log-level after test. """
    service = config.service
    log_level = config.console_log_level
    yield None
    config.service = service
    config.console_log_level = log_level
    infrastructure.setup_logging()


def create_test_structure(
        meta_info, id: int, h: int, o: int, extra: List[str], periodicity: int,
        optimade: bool = True, metadata: dict = None):
    """ Creates a calculation in Elastic and Mongodb with the given properties.

    Does require initialized :func:`elastic_infra` and :func:`mongo_infra`.

    Args:
        meta_info: A legace metainfo env.
        id: A number to create ``test_calc_id_<number>`` ids.
        h: The amount of H atoms
        o: The amount of O atoms
        extra: A list of further atoms
        periodicity: The number of dimensions to repeat the structure in
        optimade: A boolean. Iff true the entry will have optimade metadata. Default is True.
        metadata: Additional (user) metadata.
    """

    atom_labels = ['H' for i in range(0, h)] + ['O' for i in range(0, o)] + extra
    test_vector = np.array([0, 0, 0])

    backend = LocalBackend(meta_info, False, True)  # type: ignore
    backend.openSection('section_run')
    backend.addValue('program_name', 'test_code')
    backend.openSection('section_system')

    backend.addArrayValues('atom_labels', np.array(atom_labels))
    backend.addArrayValues(
        'atom_positions', np.array([test_vector for i in range(0, len(atom_labels))]))
    backend.addArrayValues(
        'lattice_vectors', np.array([test_vector, test_vector, test_vector]))
    backend.addArrayValues(
        'configuration_periodic_dimensions',
        np.array([True for _ in range(0, periodicity)] + [False for _ in range(periodicity, 3)]))

    backend.closeSection('section_system', 0)
    backend.closeSection('section_run', 0)

    backend = run_normalize(backend)
    calc = CalcWithMetadata(
        upload_id='test_uload_id', calc_id='test_calc_id_%d' % id, mainfile='test_mainfile',
        published=True, with_embargo=False)
    calc.apply_domain_metadata(backend)
    if metadata is not None:
        calc.update(**metadata)

    if not optimade:
        calc.optimade = None  # type: ignore

    proc_calc = processing.Calc.from_calc_with_metadata(calc)
    proc_calc.save()
    search_entry = search.Entry.from_calc_with_metadata(calc)
    search_entry.save()

    assert processing.Calc.objects(calc_id__in=[calc.calc_id]).count() == 1
