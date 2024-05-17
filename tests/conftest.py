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
import builtins
from pathlib import Path
from io import StringIO
import pytest
import time
import os
import socket
import json
import logging
import warnings
import tempfile

from fastapi.testclient import TestClient
import socketserver

from nomad.config import config
from nomad.config.models.plugins import Schema, add_plugin, remove_plugin

# make sure to disable logstash (the logs can interfere with the testing, especially for logtransfer)
config.logstash.enabled = False  # noqa: E402  # this must be set *before* the other modules are imported

from nomad import utils
from nomad.utils import structlogging
from nomad.app.main import app

# Set up pytest to pass control to the debugger on an exception.
if os.getenv('_PYTEST_RAISE', '0') != '0':

    @pytest.hookimpl(tryfirst=True)
    def pytest_exception_interact(call):
        raise call.excinfo.value

    @pytest.hookimpl(tryfirst=True)
    def pytest_internalerror(excinfo):
        raise excinfo.value


test_log_level = logging.CRITICAL

warnings.simplefilter('ignore')

structlogging.ConsoleFormatter.short_format = True
setattr(logging, 'Formatter', structlogging.ConsoleFormatter)

pytest_plugins = (
    'celery.contrib.pytest',
    'tests.fixtures.data',
    'tests.fixtures.groups',
    'tests.fixtures.group_uploads',
    'tests.fixtures.infrastructure',
    'tests.fixtures.mails',
    'tests.fixtures.users',
)


def pytest_addoption(parser):
    help = 'Set this < 1.0 to speed up worker cleanup. May leave tasks running.'
    parser.addoption('--celery-inspect-timeout', type=float, default=1.0, help=help)
    help = (
        'Only run tests with these fixtures and exclude ones prefixed with "!".'
        'Does not consider dynamically loaded fixtures (e.g. `request.getfixturevalue`).'
    )
    parser.addoption('--fixture-filters', nargs='+', help=help)


def filter_tests_by_fixtures(items, config):
    """Filter tests by fixture names based on CLI argument `--fixture-filters`.

    Will include tests that have all the fixtures in `--fixture-filters`
    and exclude tests that have any of the fixtures prefixed with '!'.

    Does not consider dynamically loaded fixtures (e.g. `request.getfixturevalue`)."""

    fixture_filters = config.getoption('fixture_filters')
    if not fixture_filters:
        return

    must_filters = set(f for f in fixture_filters if not f.startswith('!'))
    not_filters = set(f[1:] for f in fixture_filters if f.startswith('!'))

    selected_items = []
    deselected_items = []

    for item in items:
        fixtures = getattr(item, 'fixturenames', ())
        if must_filters.issubset(fixtures) and not_filters.isdisjoint(fixtures):
            selected_items.append(item)
        else:
            deselected_items.append(item)

    config.hook.pytest_deselected(items=deselected_items)
    items[:] = selected_items


def pytest_collection_modifyitems(items, config):
    """Manipulate the list of test items (pytest hook)."""
    filter_tests_by_fixtures(items, config)


@pytest.fixture(scope='function')
def tmp():
    parent_directory = '.volumes'
    if not os.path.isdir(parent_directory):
        os.makedirs(parent_directory, exist_ok=True)
    directory = tempfile.TemporaryDirectory(dir=parent_directory, prefix='test_tmp')
    yield directory.name
    directory.cleanup()


@pytest.fixture(scope='session')
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


@pytest.fixture(scope='function')
def plugin_schema():
    """Fixture for loading a schema plugin into the config."""
    plugin_path = str(Path(__file__, '../data/schemas').resolve())
    plugin_name = 'nomadschemaexample'
    plugin = Schema(
        name=plugin_name,
        key=plugin_name,
        package_path=plugin_path,
        python_package=plugin_name,
    )
    add_plugin(plugin)

    yield

    remove_plugin(plugin)


@pytest.fixture
def reset_config():
    """Fixture that resets configuration."""
    service = config.meta.service
    yield None
    config.meta.service = service
    utils.set_console_log_level(test_log_level)


@pytest.fixture(scope='session')
def api_v1(monkeysession, user_molds):
    """
    This fixture provides an HTTP client with Python requests interface that accesses
    the fast api. The have to provide URLs that start with out leading '/' after '.../api/v1.
    This fixture also patches the actual requests. If some code is using requests to
    connect to the NOMAD v1 at ``nomad.config.client.url``, the patch will redirect to the
    fast api under test.
    """
    test_client = TestClient(app, base_url='http://testserver/api/v1/')

    def call_test_client(method, url, *args, **kwargs):
        url = url.replace(f'{config.client.url}/v1/', '')
        return getattr(test_client, method)(url, *args, **kwargs)

    monkeysession.setattr(
        'requests.get', lambda *args, **kwargs: call_test_client('get', *args, **kwargs)
    )
    monkeysession.setattr(
        'requests.put', lambda *args, **kwargs: call_test_client('put', *args, **kwargs)
    )
    monkeysession.setattr(
        'requests.post',
        lambda *args, **kwargs: call_test_client('post', *args, **kwargs),
    )
    monkeysession.setattr(
        'requests.delete',
        lambda *args, **kwargs: call_test_client('delete', *args, **kwargs),
    )

    def __call__(self, request):
        for user in user_molds.values():
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


@pytest.fixture(scope='function')
def central_logstash_mock(monkeypatch):
    monkeypatch.setattr('nomad.config.config.logstash.enabled', True)

    class TCPServerStore(socketserver.TCPServer):
        received_content = []

        def set_request_timeout(self, timeout):
            # Note acts timeout is on the LogstashCentralHandler socket
            # this seems to behave differently to self.timeout and is
            # particularly useful to interrupt blocking "handle_request()"
            self.RequestHandlerClass.timeout = timeout

    class LogstashCentralHandler(socketserver.StreamRequestHandler):
        def handle(self):
            while True:
                try:
                    line = self.rfile.readline()
                    # print(f'received {line=}')
                except socket.timeout:
                    # print(f'server timed out')
                    line = b''  # if time out, close connection

                if line == b'':
                    # print(f'received closing for LogstashCentralHandler')
                    break

                line = line.strip()
                if len(line) > 0:
                    self.server.received_content.append(line)

    host, port = config.logstash.host, int(config.logstash.tcp_port)

    # print(f"set up mock on host={host} and port={port}")
    logstash_mock_central = TCPServerStore(
        (host, port), LogstashCentralHandler, bind_and_activate=False
    )
    logstash_mock_central.allow_reuse_address = True
    logstash_mock_central.allow_reuse_port = True
    # logstash_mock_central.timeout = 0.3

    # It is in the responsibility of the test to set a time out of the test
    while True:
        try:
            logstash_mock_central.server_bind()
            logstash_mock_central.server_activate()
        except Exception:
            time.sleep(0.001)  # try again
        else:
            break

    yield logstash_mock_central

    # make sure to close the server
    logstash_mock_central.server_close()


class MockFileManager:
    def __init__(self):
        self.files = {}
        self._open = builtins.open

    def open(self, name, mode='r', buffering=-1, **options):
        name = os.path.abspath(name)
        if mode.startswith('r') and name not in self.files:
            # We have to let some files through
            return self._open(name, mode, buffering, **options)
            # This causes stracktraces not to display
            # raise IOError(2, "No such file or directory: '%s'" % name)

        if mode.startswith('w') or (mode.startswith('a') and name not in self.files):
            buf = StringIO()
            buf.close = lambda: None
            self.files[name] = buf

        buf = self.files[name]

        if mode.startswith('r'):
            buf.seek(0)
        elif mode.startswith('a'):
            buf.seek(0)

        return buf

    def write(self, name, text):
        name = os.path.abspath(name)
        buf = StringIO(text)
        buf.close = lambda: None
        self.files[name] = buf

    def read(self, name):
        name = os.path.abspath(name)
        if name not in self.files:
            raise IOError(2, "No such file or directory: '%s'" % name)

        return self.files[name].getvalue()


@pytest.fixture
def mockopen(monkeypatch):
    manager = MockFileManager()
    monkeypatch.setattr(builtins, 'open', manager.open)
    return manager
