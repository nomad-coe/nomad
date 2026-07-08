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

import time
from collections import namedtuple

import pytest
from aiosmtpd.controller import Controller

from nomad.config import config

RecordedMessage = namedtuple(
    'RecordedMessage',
    'peer envelope_from envelope_recipients data',
)


class Handler:
    def __init__(self):
        self.messages = []

    async def handle_exception(self, exc):
        return '250 Dummy'

    async def handle_DATA(self, server, session, envelope):
        peer = session.peer
        mailfrom = envelope.mail_from
        rcpttos = envelope.rcpt_tos
        data = envelope.content
        msg = RecordedMessage(peer, mailfrom, rcpttos, data)
        self.messages.append(msg)


class SMTPServer:
    def __init__(self, port=None):
        self.host_port = None
        self.smtp = None
        self.handler = None
        self.port = port

    def run(self):
        self.handler = Handler()
        # Use provided port or fall back to config
        port = self.port if self.port is not None else config.mail.port
        self.smtp = Controller(self.handler, hostname='127.0.0.1', port=port)
        self.smtp.start()
        self.host_port = self.smtp.hostname, self.smtp.port

    def close(self):
        if self.smtp is not None:
            self.smtp.stop()


class SMTPServerFixture:
    def __init__(self, port=None):
        self.server = SMTPServer(port=port)
        self.server.run()

    @property
    def host_port(self):
        """SMTP server's listening address as a (host, port) tuple"""
        while self.server.host_port is None:
            time.sleep(0.1)
        return self.server.host_port

    @property
    def host(self):
        return self.server.host_port[0]

    @property
    def port(self):
        return self.server.host_port[1]

    @property
    def messages(self):
        """A list of RecordedMessage objects"""
        return self.server.handler.messages[:]

    def clear(self):
        self.server.handler.messages = []

    def close(self):
        self.server.close()


def get_worker_port(worker_id, base_port=None):
    """
    Calculate a unique port for each xdist worker.

    Args:
        worker_id: pytest-xdist worker identifier (e.g., 'gw0', 'gw1', 'master')
        base_port: Base port number (defaults to config.mail.port)

    Returns:
        Unique port number for the worker
    """
    if base_port is None:
        base_port = config.mail.port

    if worker_id == 'master':
        return base_port

    # Extract numeric part from worker_id (e.g., 'gw0' -> 0, 'gw1' -> 1)
    worker_num = int(worker_id.replace('gw', ''))
    return base_port + worker_num + 1


@pytest.fixture(scope='session')
def smtpd(request, monkeysession, worker_id):
    """
    SMTP server fixture with xdist-aware port assignment.

    Each xdist worker gets a unique port to avoid conflicts during parallel testing.
    """
    # on some local machines resolving the local machine takes quite a while and
    # is irrelevant for testing
    monkeysession.setattr('socket.getfqdn', lambda *args, **kwargs: 'local.server')

    # Calculate unique port for this worker
    port = get_worker_port(worker_id)

    fixture = SMTPServerFixture(port=port)
    request.addfinalizer(fixture.close)
    return fixture


@pytest.fixture(scope='function')
def mails(smtpd, monkeypatch):
    smtpd.clear()
    monkeypatch.setattr('nomad.config.mail.enabled', True)
    monkeypatch.setattr('nomad.config.mail.host', 'localhost')
    monkeypatch.setattr('nomad.config.mail.port', smtpd.port)
    yield smtpd
