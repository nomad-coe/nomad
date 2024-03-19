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
    def __init__(self):
        self.host_port = None
        self.smtp = None
        self.handler = None

    def run(self):
        self.handler = Handler()
        self.smtp = Controller(
            self.handler, hostname='127.0.0.1', port=config.mail.port
        )
        self.smtp.start()
        self.host_port = self.smtp.hostname, self.smtp.port

    def close(self):
        if self.smtp is not None:
            self.smtp.stop()


class SMTPServerFixture:
    def __init__(self):
        self.server = SMTPServer()
        self.server.run()
        _ = self.host_port

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


@pytest.fixture(scope='session')
def smtpd(request, monkeysession):
    # on some local machines resolving the local machine takes quit a while and
    # is irrelevant for testing
    monkeysession.setattr('socket.getfqdn', lambda *args, **kwargs: 'local.server')
    fixture = SMTPServerFixture()
    request.addfinalizer(fixture.close)
    return fixture


@pytest.fixture(scope='function')
def mails(smtpd, monkeypatch):
    smtpd.clear()
    monkeypatch.setattr('nomad.config.mail.enabled', True)
    monkeypatch.setattr('nomad.config.mail.host', 'localhost')
    yield smtpd
