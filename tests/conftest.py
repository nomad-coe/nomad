#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

from typing import Tuple, List
import math
import pytest
import logging
from collections import namedtuple
from smtpd import SMTPServer
from threading import Lock, Thread
import asyncore
import time
from datetime import datetime
import shutil
import os.path
import elasticsearch.exceptions
from typing import List
import json
import logging
import warnings
import os.path
from fastapi.testclient import TestClient

from nomad import config, infrastructure, processing, utils, datamodel, files
from nomad.datamodel import User, EntryArchive, OptimadeEntry
from nomad.utils import structlogging
from nomad.archive import write_archive, read_archive, write_partial_archive_to_mongo
from nomad.processing import ProcessStatus
from nomad.app.main import app
from nomad.utils.exampledata import ExampleData

from tests.parsing import test_parsing
from tests.normalizing.conftest import run_normalize
from tests.processing import test_data as test_processing
from tests.test_files import empty_file, example_file_vasp_with_binary
from tests.utils import create_template_upload_file, set_upload_entry_metadata, build_url

test_log_level = logging.CRITICAL

elastic_test_entries_index = 'nomad_entries_v1_test'
elastic_test_materials_index = 'nomad_materials_v1_test'

indices = [elastic_test_entries_index, elastic_test_materials_index]

warnings.simplefilter("ignore")


structlogging.ConsoleFormatter.short_format = True
setattr(logging, 'Formatter', structlogging.ConsoleFormatter)


@pytest.fixture(scope='function')
def tmp():
    directory = '.volumes/test_tmp'
    if os.path.exists(directory):
        shutil.rmtree(directory)
    os.mkdir(directory)
    yield directory
    shutil.rmtree(directory)


@pytest.fixture(scope="session")
def monkeysession(request):
    from _pytest.monkeypatch import MonkeyPatch
    mpatch = MonkeyPatch()
    yield mpatch
    mpatch.undo()


@pytest.fixture(scope='session', autouse=True)
def nomad_logging(monkeysession):
    utils.set_console_log_level(test_log_level)
    monkeysession.setattr('logging.Logger.setLevel', lambda *args, **kwargs: None)
    monkeysession.setattr('logging.Handler.setLevel', lambda *args, **kwargs: None)


@pytest.fixture(scope='session', autouse=True)
def raw_files_infra():
    config.fs.tmp = '.volumes/test_fs/tmp'
    config.fs.staging = '.volumes/test_fs/staging'
    config.fs.public = '.volumes/test_fs/public'
    config.fs.staging_external = os.path.abspath(config.fs.staging)
    config.fs.public_external = os.path.abspath(config.fs.public)
    config.fs.prefix_size = 2
    clear_raw_files()


@pytest.fixture(scope='module')
def raw_files_module(raw_files_infra):
    ''' Provides cleaned out files directory structure per module. Clears files before test. '''
    clear_raw_files()


@pytest.fixture(scope='function')
def raw_files_function(raw_files_infra):
    ''' Provides cleaned out files directory structure per function. Clears files before test. '''
    clear_raw_files()


@pytest.fixture(scope='function')
def raw_files(raw_files_infra):
    ''' Provides cleaned out files directory structure per function. Clears files before test. '''
    clear_raw_files()


def clear_raw_files():
    directories = [config.fs.staging, config.fs.public, config.fs.tmp]
    for directory in directories:
        try:
            shutil.rmtree(directory)
        except FileNotFoundError:
            pass

        os.makedirs(directory)


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
    '''
    Purges all pending tasks of the celery app before test. This is necessary to
    remove tasks from the queue that might be 'left over' from prior tests.
    '''
    celery_session_app.control.purge()
    yield celery_session_app
    celery_session_app.control.purge()


@pytest.fixture(scope='session')
def celery_inspect(purged_app):
    yield purged_app.control.inspect()


# It might be necessary to make this a function scoped fixture, if old tasks keep
# 'bleeding' into successive tests.
@pytest.fixture(scope='function')
def worker(mongo, celery_session_worker, celery_inspect):
    ''' Provides a clean worker (no old tasks) per function. Waits for all tasks to be completed. '''
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
        print('Exception during worker tear down.')
        import traceback
        traceback.print_exc()


@pytest.fixture(scope='session')
def mongo_infra(monkeysession):
    monkeysession.setattr('nomad.config.mongo.db_name', 'test_db')
    # disconnecting and connecting again results in an empty database with mongomock
    monkeysession.setattr('mongoengine.disconnect', lambda *args, **kwargs: None)
    return infrastructure.setup_mongo()


def clear_mongo(mongo_infra):
    # Some test cases need to reset the database connection
    infrastructure.mongo_client.drop_database('test_db')
    return infrastructure.mongo_client


@pytest.fixture(scope='module')
def mongo_module(mongo_infra):
    ''' Provides a cleaned mocked mongo per module. '''
    return clear_mongo(mongo_infra)


@pytest.fixture(scope='function')
def mongo_function(mongo_infra):
    ''' Provides a cleaned mocked mongo per function. '''
    return clear_mongo(mongo_infra)


@pytest.fixture(scope='function')
def mongo(mongo_infra):
    ''' Provides a cleaned mocked mongo per function. '''
    print('setup mongo worker')
    return clear_mongo(mongo_infra)


@pytest.fixture(scope='session')
def elastic_infra(monkeysession):
    ''' Provides elastic infrastructure to the session '''
    monkeysession.setattr('nomad.config.elastic.entries_index', elastic_test_entries_index)
    monkeysession.setattr('nomad.config.elastic.materials_index', elastic_test_materials_index)

    # attempt to remove and recreate all indices
    return clear_elastic_infra()


def clear_elastic_infra():
    '''
    Removes and re-creates all indices and mappings.
    '''
    from elasticsearch_dsl import connections
    connection = connections.create_connection(
        hosts=['%s:%d' % (config.elastic.host, config.elastic.port)])

    for index in indices:
        try:
            connection.indices.delete(index=index)
        except Exception:
            pass

    return infrastructure.setup_elastic()


def clear_elastic(elastic_infra):
    '''
    Removes all contents from the existing indices.
    '''
    try:
        for index in indices:
            elastic_infra.delete_by_query(
                index=index, body=dict(query=dict(match_all={})),
                wait_for_completion=True, refresh=True)
    except elasticsearch.exceptions.NotFoundError:
        # Happens if a test removed indices without recreating them.
        clear_elastic_infra()

    assert infrastructure.elastic_client is not None
    return elastic_infra


@pytest.fixture(scope='module')
def elastic_module(elastic_infra):
    ''' Provides a clean elastic per module. Clears elastic before test. '''
    return clear_elastic(elastic_infra)


@pytest.fixture(scope='function')
def elastic_function(elastic_infra):
    ''' Provides a clean elastic per function. Clears elastic before test. '''
    return clear_elastic(elastic_infra)


@pytest.fixture(scope='function')
def elastic(elastic_infra):
    ''' Provides a clean elastic per function. Clears elastic before test. '''
    return clear_elastic(elastic_infra)


def test_user_uuid(handle):
    return '00000000-0000-0000-0000-00000000000%d' % handle


admin_user_id = test_user_uuid(0)

test_users = {
    test_user_uuid(0): dict(username='admin', email='admin', user_id=test_user_uuid(0)),
    test_user_uuid(1): dict(username='scooper', email='sheldon.cooper@nomad-coe.eu', first_name='Sheldon', last_name='Cooper', user_id=test_user_uuid(1), is_oasis_admin=True),
    test_user_uuid(2): dict(username='lhofstadter', email='leonard.hofstadter@nomad-fairdi.tests.de', first_name='Leonard', last_name='Hofstadter', user_id=test_user_uuid(2))
}


@pytest.fixture(scope='session', autouse=True)
def configure_admin_user_id(monkeysession):
    monkeysession.setattr('nomad.config.services.admin_user_id', admin_user_id)


class KeycloakMock:
    def __init__(self):
        self.id_counter = 2
        self.users = dict(**test_users)

    def tokenauth(self, access_token: str):
        if access_token in self.users:
            return User(**self.users[access_token])
        else:
            raise infrastructure.KeycloakError('user does not exist')

    def add_user(self, user, *args, **kwargs):
        self.id_counter += 1
        user.user_id = test_user_uuid(self.id_counter)
        user.username = (user.first_name[0] + user.last_name).lower()
        self.users[user.user_id] = dict(
            email=user.email, username=user.username, first_name=user.first_name,
            last_name=user.last_name, user_id=user.user_id)
        return None

    def get_user(self, user_id=None, username=None, email=None):
        if user_id is not None:
            return User(**self.users[user_id])
        elif username is not None:
            for user_id, user_values in self.users.items():
                if user_values['username'] == username:
                    return User(**user_values)
            raise KeyError('Only test user usernames are recognized')
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
            if query in ' '.join([str(value) for value in test_user.values()])]

    def basicauth(self, username: str, password: str) -> str:
        for user in self.users.values():
            if user['username'] == username or user['email'] == username:
                return user['user_id']

        raise infrastructure.KeycloakError()


config.keycloak.realm_name = 'fairdi_nomad_test'
config.keycloak.password = 'password'

_keycloak = infrastructure.keycloak
_user_management = infrastructure.user_management


# use a session fixture in addition to the function fixture, to ensure mocked keycloak
# before other class, module, etc. scoped function are run
@pytest.fixture(scope='session', autouse=True)
def mocked_keycloak_session(monkeysession):
    monkeysession.setattr('nomad.infrastructure.keycloak', KeycloakMock())
    monkeysession.setattr('nomad.infrastructure.user_management', KeycloakMock())


@pytest.fixture(scope='function', autouse=True)
def mocked_keycloak(monkeypatch):
    monkeypatch.setattr('nomad.infrastructure.keycloak', KeycloakMock())
    monkeypatch.setattr('nomad.infrastructure.user_management', KeycloakMock())


@pytest.fixture(scope='function')
def keycloak(monkeypatch):
    monkeypatch.setattr('nomad.infrastructure.keycloak', _keycloak)
    monkeypatch.setattr('nomad.infrastructure.user_management', _user_management)


@pytest.fixture(scope='function')
def proc_infra(worker, elastic, mongo, raw_files):
    ''' Combines all fixtures necessary for processing (elastic, worker, files, mongo) '''
    return dict(elastic=elastic)


@pytest.fixture(scope='function')
def with_oasis_user_management(monkeypatch):
    from nomad.infrastructure import OasisUserManagement
    monkeypatch.setattr('nomad.infrastructure.user_management', OasisUserManagement())
    yield
    monkeypatch.setattr('nomad.infrastructure.user_management', _user_management)


@pytest.fixture(scope='module')
def test_user():
    return User(**test_users[test_user_uuid(1)])


@pytest.fixture(scope='module')
def other_test_user():
    return User(**test_users[test_user_uuid(2)])


@pytest.fixture(scope='module')
def admin_user():
    return User(**test_users[test_user_uuid(0)])


@pytest.fixture(scope='module')
def test_users_dict(test_user, other_test_user, admin_user):
    return {
        'test_user': test_user,
        'other_test_user': other_test_user,
        'admin_user': admin_user}


@pytest.fixture(scope='function')
def no_warn(caplog):
    caplog.handler.formatter = structlogging.ConsoleFormatter()
    yield caplog
    for record in caplog.get_records(when='call'):
        if record.levelname in ['WARNING', 'ERROR', 'CRITICAL']:
            try:
                msg = structlogging.ConsoleFormatter.serialize(json.loads(record.msg))
            except Exception:
                msg = record.msg
            assert False, msg


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


'''
Fixture for mocked SMTP server for testing.
Based on https://gist.github.com/akheron/cf3863cdc424f08929e4cb7dc365ef23.
'''

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
    return ('parsers/template', 'tests/data/templates/template.json')


@pytest.fixture(scope='function', params=['empty_file', 'example_file'])
def example_upload(request, tmp) -> str:
    if request.param == 'empty_file':
        return create_template_upload_file(tmp, mainfiles=[], auxfiles=0)

    return create_template_upload_file(
        tmp, mainfiles=['tests/data/proc/templates/template.json'])


@pytest.fixture(scope='function')
def non_empty_example_upload(tmp):
    return create_template_upload_file(
        tmp, mainfiles=['tests/data/proc/templates/template.json'])


@pytest.fixture(scope='session')
def non_empty_example_upload_vasp_with_binary():
    return example_file_vasp_with_binary


@pytest.fixture(scope='session')
def empty_upload():
    return empty_file


@pytest.fixture(scope='module')
def example_user_metadata(other_test_user, test_user) -> dict:
    return {
        'comment': 'test comment',
        'references': ['http://external.ref/one', 'http://external.ref/two'],
        'entry_coauthors': [other_test_user.user_id],
        '_pid': '256',
        'external_id': 'external_test_id'
    }


@pytest.fixture(scope='module')
def internal_example_user_metadata(example_user_metadata) -> dict:
    return {
        key[1:] if key[0] == '_' else key: value
        for key, value in example_user_metadata.items()}


@pytest.fixture(scope='session')
def parsed(example_mainfile: Tuple[str, str]) -> EntryArchive:
    ''' Provides a parsed entry in the form of an EntryArchive. '''
    parser, mainfile = example_mainfile
    return test_parsing.run_singular_parser(parser, mainfile)


@pytest.fixture(scope='session')
def parsed_ems() -> EntryArchive:
    ''' Provides a parsed experiment in the form of a EntryArchive. '''
    return test_parsing.run_singular_parser('parsers/eels', 'tests/data/parsers/eels.json')


@pytest.fixture(scope='session')
def normalized(parsed: EntryArchive) -> EntryArchive:
    ''' Provides a normalized entry in the form of a EntryArchive. '''
    return run_normalize(parsed)


@pytest.fixture(scope='function')
def uploaded(example_upload: str, raw_files) -> Tuple[str, str]:
    '''
    Provides a uploaded with uploaded example file and gives the upload_id.
    Clears files after test.
    '''
    example_upload_id = os.path.basename(example_upload).replace('.zip', '')
    return example_upload_id, example_upload


@pytest.fixture(scope='function')
def non_empty_uploaded(non_empty_example_upload: str, raw_files) -> Tuple[str, str]:
    example_upload_id = os.path.basename(non_empty_example_upload).replace('.zip', '')
    return example_upload_id, non_empty_example_upload


@pytest.fixture(scope='function')
def oasis_publishable_upload(
        api_v1, proc_infra, non_empty_processed: processing.Upload, internal_example_user_metadata,
        monkeypatch, test_user):
    '''
    Creates a published upload which can be used with Upload.publish_externally. Some monkeypatching
    is done which replaces IDs when importing.
    '''
    # Create a published upload
    set_upload_entry_metadata(non_empty_processed, internal_example_user_metadata)
    non_empty_processed.publish_upload()
    non_empty_processed.block_until_complete(interval=.01)

    suffix = '_2'  # Will be added to all IDs in the mirrored upload
    upload_id = non_empty_processed.upload_id

    # Do some tricks to add suffix to the ID fields
    old_bundle_init = files.UploadBundle.__init__

    def new_bundle_init(self, *args, **kwargs):
        old_bundle_init(self, *args, **kwargs)
        # Change the id's in the bundle_info dict behind the scenes when loading the bundle
        bundle_info = self.bundle_info
        bundle_info['upload_id'] += suffix
        bundle_info['upload']['_id'] += suffix
        for entry_dict in bundle_info['entries']:
            entry_dict['_id'] = utils.generate_entry_id(
                upload_id + suffix, entry_dict['mainfile'], entry_dict.get('mainfile_key'))
            entry_dict['upload_id'] += suffix

    old_bundle_import_files = files.UploadBundle.import_upload_files

    def new_bundle_import_files(self, *args, **kwargs):
        upload_files = old_bundle_import_files(self, *args, **kwargs)
        # Overwrite the archive files with files containing the updated IDs
        archive_path = upload_files.os_path
        for file_name in os.listdir(archive_path):
            if file_name.endswith('.msg'):
                full_path = os.path.join(archive_path, file_name)
                data = read_archive(full_path)
                new_data = []
                for entry_id in data.keys():
                    archive_dict = data[entry_id].to_dict()
                    section_metadata = archive_dict['metadata']
                    section_metadata['upload_id'] += suffix
                    new_entry_id = utils.generate_entry_id(
                        section_metadata['upload_id'],
                        section_metadata['mainfile'],
                        section_metadata.get('mainfile_key'))
                    section_metadata['entry_id'] = new_entry_id
                    new_data.append((new_entry_id, archive_dict))
                write_archive(full_path, len(new_data), new_data)
        return upload_files

    monkeypatch.setattr('nomad.files.UploadBundle.__init__', new_bundle_init)
    monkeypatch.setattr('nomad.files.UploadBundle.import_upload_files', new_bundle_import_files)

    # Further monkey patching
    def new_post(url, data, params={}, **kwargs):
        return api_v1.post(
            build_url(url.lstrip('/api/v1/'), params),
            data=data.read(), **kwargs)

    monkeypatch.setattr('requests.post', new_post)
    monkeypatch.setattr('nomad.config.oasis.is_oasis', True)
    monkeypatch.setattr('nomad.config.keycloak.username', test_user.username)

    monkeypatch.setattr('nomad.config.oasis.central_nomad_api_url', '/api')

    # create a dataset to also test this aspect of oasis uploads
    entry = non_empty_processed.successful_entries[0]
    datamodel.Dataset(
        dataset_id='dataset_id', dataset_name='dataset_name',
        user_id=test_user.user_id).a_mongo.save()
    entry.datasets = ['dataset_id']
    entry.save()
    return non_empty_processed.upload_id, suffix


@pytest.mark.timeout(config.tests.default_timeout)
@pytest.fixture(scope='function')
def processed(uploaded: Tuple[str, str], test_user: User, proc_infra) -> processing.Upload:
    '''
    Provides a processed upload. Upload was uploaded with test_user.
    '''
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
    '''
    Provides a processed upload. Upload was uploaded with test_user.
    '''
    return test_processing.run_processing(non_empty_uploaded, test_user)


@pytest.mark.timeout(config.tests.default_timeout)
@pytest.fixture(scope='function')
def published(non_empty_processed: processing.Upload, internal_example_user_metadata) -> processing.Upload:
    '''
    Provides a processed published upload. Upload was uploaded with test_user and is embargoed.
    '''
    set_upload_entry_metadata(non_empty_processed, internal_example_user_metadata)
    non_empty_processed.publish_upload(embargo_length=12)
    try:
        non_empty_processed.block_until_complete(interval=.01)
    except Exception:
        pass

    return non_empty_processed


@pytest.mark.timeout(config.tests.default_timeout)
@pytest.fixture(scope='function')
def published_wo_user_metadata(non_empty_processed: processing.Upload) -> processing.Upload:
    '''
    Provides a processed upload. Upload was uploaded with test_user.
    '''
    non_empty_processed.publish_upload()
    try:
        non_empty_processed.block_until_complete(interval=.01)
    except Exception:
        pass

    return non_empty_processed


@pytest.fixture(scope='module')
def example_data(elastic_module, raw_files_module, mongo_module, test_user, other_test_user, normalized):
    '''
    Provides a couple of uploads and entries including metadata, raw-data, and
    archive files.

    id_embargo:
        1 entry, 1 material, published with embargo
    id_embargo_w_coauthor:
        1 entry, 1 material, published with embargo and coauthor
    id_embargo_w_reviewer:
        1 entry, 1 material, published with embargo and reviewer
    id_unpublished:
        1 entry, 1 material, unpublished
    id_unpublished_w_coauthor:
        1 entry, 1 material, unpublished with coauthor
    id_unpublished_w_reviewer:
        1 entry, 1 material, unpublished with reviewer
    id_published:
        23 entries, 6 materials published without embargo
        partial archive exists only for id_01
        raw files and archive file for id_02 are missing
        id_10, id_11 reside in the same directory
    id_child_entries:
        1 parent entry and 2 child entries from one mainfile, 1 material, unpublished
    id_processing:
        unpublished upload without any entries, in status processing
    id_empty:
        unpublished upload without any entries
    '''
    data = ExampleData(main_author=test_user)

    # 6 uploads with different combinations of main_type and sub_type
    for main_type in ('embargo', 'unpublished'):
        for sub_type in ('', 'w_coauthor', 'w_reviewer'):
            upload_id = 'id_' + main_type + ('_' if sub_type else '') + sub_type
            if main_type == 'embargo':
                published = True
                embargo_length = 12
                upload_name = 'name_' + upload_id[3:]
            else:
                published = False
                embargo_length = 0
                upload_name = None
            entry_id = upload_id + '_1'
            coauthors = [other_test_user.user_id] if sub_type == 'w_coauthor' else None
            reviewers = [other_test_user.user_id] if sub_type == 'w_reviewer' else None
            data.create_upload(
                upload_id=upload_id,
                upload_name=upload_name,
                coauthors=coauthors,
                reviewers=reviewers,
                published=published,
                embargo_length=embargo_length)
            data.create_entry(
                upload_id=upload_id,
                entry_id=entry_id,
                material_id=upload_id,
                mainfile=f'test_content/{entry_id}/mainfile.json')

    # one upload with 23 entries, published, no embargo
    data.create_upload(
        upload_id='id_published',
        upload_name='name_published',
        published=True)
    for i in range(1, 24):
        entry_id = 'id_%02d' % i
        material_id = 'id_%02d' % (int(math.floor(i / 4)) + 1)
        mainfile = 'test_content/subdir/test_entry_%02d/mainfile.json' % i
        kwargs = dict(optimade=OptimadeEntry(nelements=2, elements=['H', 'O']))
        if i == 11:
            mainfile = 'test_content/subdir/test_entry_10/mainfile_11.json'
        if i == 1:
            kwargs['pid'] = '123'
        data.create_entry(
            upload_id='id_published',
            entry_id=entry_id,
            material_id=material_id,
            mainfile=mainfile,
            **kwargs)

        if i == 1:
            archive = data.archives[entry_id]
            write_partial_archive_to_mongo(archive)

    # 3 entries from one mainfile, 1 material, unpublished
    upload_id = 'id_child_entries'
    data.create_upload(
        upload_id=upload_id,
        upload_name='name_child_entries',
        published=False)
    for mainfile_key in (None, 'child1', 'child2'):
        data.create_entry(
            upload_id=upload_id,
            entry_id=upload_id + '_' + (mainfile_key or 'main'),
            material_id=upload_id,
            mainfile=f'test_content/mainfile_w_children.json',
            mainfile_key=mainfile_key)

    # one upload, no entries, still processing
    data.create_upload(
        upload_id='id_processing',
        published=False,
        process_status=ProcessStatus.RUNNING)

    # one upload, no entries, unpublished
    data.create_upload(
        upload_id='id_empty',
        published=False)

    data.save(with_files=False)
    del(data.archives['id_02'])
    data.save(with_files=True, with_es=False, with_mongo=False)


@pytest.fixture(scope='function')
def example_data_writeable(mongo, test_user, normalized):
    data = ExampleData(main_author=test_user)

    # one upload with one entry, published
    data.create_upload(
        upload_id='id_published_w',
        published=True,
        embargo_length=12)
    data.create_entry(
        upload_id='id_published_w',
        entry_id='id_published_w_entry',
        mainfile='test_content/test_embargo_entry/mainfile.json')

    # one upload with one entry, unpublished
    data.create_upload(
        upload_id='id_unpublished_w',
        published=False,
        embargo_length=12)
    data.create_entry(
        upload_id='id_unpublished_w',
        entry_id='id_unpublished_w_entry',
        mainfile='test_content/test_embargo_entry/mainfile.json')

    # one upload, no entries, running a blocking processing
    data.create_upload(
        upload_id='id_processing_w',
        published=False,
        process_status=ProcessStatus.RUNNING,
        current_process='publish_upload')

    # one upload, no entries, unpublished
    data.create_upload(
        upload_id='id_empty_w',
        published=False)

    data.save()

    yield

    data.delete()


@pytest.fixture(scope='function')
def example_datasets(mongo, test_user, other_test_user):
    dataset_specs = (
        ('test_dataset_1', test_user, None),
        ('test_dataset_2', test_user, 'test_doi_2'),
        ('test_dataset_3', other_test_user, None)
    )
    datasets = []
    for dataset_name, user, doi in dataset_specs:
        now = datetime.utcnow()
        dataset = datamodel.Dataset(
            dataset_id=utils.create_uuid(),
            dataset_name=dataset_name,
            doi=doi,
            user_id=user.user_id,
            dataset_create_time=now,
            dataset_modified_time=now,
            dataset_type='owned')
        dataset.a_mongo.create()
        datasets.append(dataset)

    yield datasets

    while datasets:
        datasets.pop().a_mongo.delete()


@pytest.fixture
def reset_config():
    ''' Fixture that resets configuration. '''
    service = config.meta.service
    yield None
    config.meta.service = service
    utils.set_console_log_level(test_log_level)


@pytest.fixture
def reset_infra(mongo, elastic):
    ''' Fixture that resets infrastructure after deleting db or search index. '''
    yield None


@pytest.fixture(scope='session')
def api_v1(monkeysession):
    '''
    This fixture provides an HTTP client with Python requests interface that accesses
    the fast api. The have to provide URLs that start with out leading '/' after '.../api/v1.
    This fixture also patches the actual requests. If some code is using requests to
    connect to the NOMAD v1 at ``nomad.config.client.url``, the patch will redirect to the
    fast api under test.
    '''
    test_client = TestClient(app, base_url='http://testserver/api/v1/')

    def call_test_client(method, url, *args, **kwargs):
        url = url.replace(f'{config.client.url}/v1/', '')
        return getattr(test_client, method)(url, *args, **kwargs)

    monkeysession.setattr('requests.get', lambda *args, **kwargs: call_test_client('get', *args, **kwargs))
    monkeysession.setattr('requests.put', lambda *args, **kwargs: call_test_client('put', *args, **kwargs))
    monkeysession.setattr('requests.post', lambda *args, **kwargs: call_test_client('post', *args, **kwargs))
    monkeysession.setattr('requests.delete', lambda *args, **kwargs: call_test_client('delete', *args, **kwargs))

    def __call__(self, request):
        for user in test_users.values():
            if user['username'] == self.user or user['email'] == self.user:
                request.headers['Authorization'] = f'Bearer {user["user_id"]}'
        return request

    monkeysession.setattr('nomad.client.api.Auth.__call__', __call__)

    return test_client


@pytest.fixture(scope='session')
def client_with_api_v1(api_v1, monkeysession):
    def call_requests(method, path, *args, **kwargs):
        return getattr(api_v1, method)(path, *args, **kwargs)

    monkeysession.setattr('nomad.client.api._call_requests', call_requests)
