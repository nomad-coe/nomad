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

"""
.. autofunc::nomad.utils.create_uuid
.. autofunc::nomad.utils.hash
.. autofunc::nomad.utils.timer

Logging in nomad is structured. Structured logging means that log entries contain
dictionaries with quantities related to respective events. E.g. having the code,
parser, parser version, calc_id, mainfile, etc. for all events that happen during
calculation processing. This means the :func:`get_logger` and all logger functions
take keyword arguments for structured data. Otherwise :func:`get_logger` can
be used similar to the standard *logging.getLogger*.

Depending on the configuration all logs will also be send to a central logstash.

.. autofunc::nomad.utils.get_logger
.. autofunc::nomad.utils.hash
.. autofunc::nomad.utils.create_uuid
.. autofunc::nomad.utils.timer
.. autofunc::nomad.utils.lnr
"""

from typing import List
import base64
import logging
import structlog
from structlog.processors import StackInfoRenderer, format_exc_info, TimeStamper, JSONRenderer
from structlog.stdlib import LoggerFactory
import logstash
from contextlib import contextmanager
import json
import uuid
import time
import re
from werkzeug.exceptions import HTTPException
import hashlib

default_hash_len = 28
""" Length of hashes and hash-based ids (e.g. calc, upload) in nomad. """


def hash(*args, length: int = default_hash_len) -> str:
    """ Creates a websave hash of the given length based on the repr of the given arguments. """
    hash = hashlib.sha512()
    for arg in args:
        hash.update(str(arg).encode('utf-8'))

    return make_websave(hash, length=length)


def make_websave(hash, length: int = default_hash_len) -> str:
    """ Creates a websave string for a hashlib hash object. """
    if length > 0:
        return base64.b64encode(hash.digest(), altchars=b'-_')[:length].decode('utf-8')
    else:
        return base64.b64encode(hash.digest(), altchars=b'-_')[0:-2].decode('utf-8')


def sanitize_logevent(event: str) -> str:
    """
    Prepares a log event or message for analysis in elastic stack. It removes numbers,
    list, and matrices of numbers from the event string and limits its size. The
    goal is to make it easier to define aggregations over events by using event
    strings as representatives for event classes rather than event instances (with
    concrete numbers, etc).
    """
    sanitized_event = event[:120]
    sanitized_event = re.sub(r'(\d*\.\d+|\d+(\.\d*)?)', 'X', sanitized_event)
    sanitized_event = re.sub(r'((\[|\()\s*)?X\s*(,\s*X)+(\s*(\]|\)))?', 'L', sanitized_event)
    sanitized_event = re.sub(r'((\[|\()\s*)?[XL](,\s*[XL])+(\s*(\]|\)))?', 'M', sanitized_event)
    return sanitized_event


@contextmanager
def legacy_logger(logger):
    """ Context manager that makes the given logger the logger for legacy log entries. """
    LogstashHandler.legacy_logger = logger
    try:
        yield
    finally:
        LogstashHandler.legacy_logger = None


class LogstashHandler(logstash.TCPLogstashHandler):
    """
    A log handler that emits records to logstash. It also filters logs for being
    structlog entries. All other entries are diverted to a global `legacy_logger`.
    This legacy logger is supposed to be a structlog logger that turns legacy
    records into structlog entries with reasonable binds depending on the current
    execution context (e.g. parsing/normalizing, etc.). If no legacy logger is
    set, they get emitted as usual (e.g. non nomad logs, celery, dbs, etc.)
    """

    legacy_logger = None

    def filter(self, record):
        if super().filter(record):
            is_structlog = record.msg.startswith('{') and record.msg.endswith('}')
            if is_structlog:
                return True
            else:
                if LogstashHandler.legacy_logger is None:
                    return True
                else:
                    LogstashHandler.legacy_logger.log(
                        record.levelno, sanitize_logevent(record.msg), args=record.args,
                        exc_info=record.exc_info, stack_info=record.stack_info,
                        legacy_logger=record.name)

                    return False

        return False


class LogstashFormatter(logstash.formatter.LogstashFormatterBase):

    def format(self, record):
        try:
            structlog = json.loads(record.getMessage())
        except json.JSONDecodeError:
            structlog = dict(event=record.getMessage())

        # Create message dict
        message = {
            '@timestamp': self.format_timestamp(record.created),
            '@version': '1',
            'event': structlog['event'],
            'message': structlog['event'],
            'host': self.host,
            'path': record.pathname,
            'tags': self.tags,
            'type': self.message_type,

            # Extra Fields
            'level': record.levelname,
            'logger_name': record.name,
        }

        if record.name.startswith('nomad'):
            for key, value in structlog.items():
                if key in ('event', 'stack_info', 'id', 'timestamp'):
                    continue
                elif key == 'exception':
                    message['digest'] = str(value)[-256:]
                elif key in (
                        'upload_id', 'calc_id', 'mainfile',
                        'service', 'release'):
                    key = 'nomad.%s' % key
                else:
                    key = '%s.%s' % (record.name, key)

                message[key] = value
        else:
            from nomad import config
            structlog['nomad.service'] = config.service
            structlog['nomad.release'] = config.release
            message.update(structlog)

        # Add extra fields
        message.update(self.get_extra_fields(record))

        # If exception, add debug info
        if record.exc_info:
            message.update(self.get_debug_fields(record))

        return self.serialize(message)


class ConsoleFormatter(LogstashFormatter):
    @classmethod
    def serialize(cls, message_dict):
        from io import StringIO

        logger = message_dict.pop('logger_name', None)
        event = message_dict.pop('event', None)
        level = message_dict.pop('level', None)
        exception = message_dict.pop('exception', None)
        time = message_dict.pop('@timestamp', None)
        for key in ['type', 'tags', 'stack_info', 'path', 'message', 'host', '@version', 'digest']:
            message_dict.pop(key)
        keys = list(message_dict.keys())
        keys.sort()

        out = StringIO()
        out.write('%s %s %s %s' % (
            level.ljust(8), logger.ljust(20)[:20], time.ljust(19)[:19], event))
        if exception is not None:
            out.write('\n  - exception: %s' % str(exception).replace('\n', '\n    '))

        for key in keys:
            out.write('\n  - %s: %s' % (key, str(message_dict.get(key, None))))

        return out.getvalue()


def add_logstash_handler(logger):
    logstash_handler = next((
        handler for handler in logger.handlers
        if isinstance(handler, LogstashHandler)), None)

    if logstash_handler is None:
        from nomad import config
        logstash_handler = LogstashHandler(
            config.logstash.host,
            config.logstash.tcp_port, version=1)
        logstash_handler.formatter = LogstashFormatter(tags=['nomad', config.release])
        logstash_handler.setLevel(config.logstash.level)
        logger.addHandler(logstash_handler)


def configure_logging():
    from nomad import config

    # configure structlog
    log_processors = [
        StackInfoRenderer(),
        format_exc_info,
        TimeStamper(fmt="%Y-%m-%d %H:%M.%S", utc=False),
        JSONRenderer(sort_keys=True)
    ]

    default_factory = LoggerFactory()

    def logger_factory(*args):
        logger = default_factory(*args)
        logger.setLevel(logging.DEBUG)
        return logger

    structlog.configure(
        processors=log_processors,
        logger_factory=logger_factory,
        wrapper_class=structlog.stdlib.BoundLogger)

    # configure logging in general
    logging.basicConfig(level=logging.DEBUG)
    root = logging.getLogger()
    for handler in root.handlers:
        if not isinstance(handler, LogstashHandler):
            handler.setLevel(config.console_log_level)
            handler.setFormatter(ConsoleFormatter())

    # configure logstash
    if config.logstash.enabled:
        add_logstash_handler(root)

    logger = get_logger(__name__)
    logger.info('structlog configured', with_logstash=config.logstash.enabled)

    # configure log levels
    for logger in [
            'elasticsearch',
            # 'celery.app.trace', 'celery.worker.strategy',
            'urllib3.connectionpool', 'bravado', 'bravado_core', 'swagger_spec_validator']:
        logging.getLogger(logger).setLevel(logging.WARNING)


def create_uuid() -> str:
    """ Returns a web-save base64 encoded random uuid (type 4). """
    return base64.b64encode(uuid.uuid4().bytes, altchars=b'-_').decode('utf-8')[0:-2]


def get_logger(name, **kwargs):
    """
    Returns a structlog logger that is already attached with a logstash handler.
    Use additional *kwargs* to pre-bind some values to all events.
    """
    from nomad import config

    if name.startswith('nomad.'):
        name = '.'.join(name.split('.')[:2])

    logger = structlog.get_logger(name, service=config.service, release=config.release, **kwargs)
    return logger


@contextmanager
def lnr(logger, event, **kwargs):
    """
    A context manager that Logs aNd Raises all exceptions with the given logger.

    Arguments:
        logger: The logger that should be used for logging exceptions.
        event: the log message
        **kwargs: additional properties for the structured log
    """
    try:
        yield
    except HTTPException as e:
        # ignore HTTPException as they are part of the normal flask error handling
        raise e
    except Exception as e:
        logger.error(event, exc_info=e, **kwargs)
        raise e


@contextmanager
def timer(logger, event, method='info', **kwargs):
    """
    A context manager that takes execution time and produces a log entry with said time.

    Arguments:
        logger: The logger that should be used to produce the log entry.
        event: The log message/event.
        method: The log methad that should be used. Must be a valid logger method name.
            Default is 'info'.
        **kwargs: Additional logger data that is passed to the log entry.

    Returns:
        The method yields a dictionary that can be used to add further log data.
    """
    start = time.time()

    try:
        yield kwargs
    finally:
        stop = time.time()

    logger_method = getattr(logger, 'info', None)
    if logger_method is not None:
        logger_method(event, exec_time=stop - start, **kwargs)
    else:
        logger.error('Unknown logger method %s.' % method)


class archive:
    @staticmethod
    def create(upload_id: str, calc_id: str) -> str:
        return '%s/%s' % (upload_id, calc_id)

    @staticmethod
    def items(archive_id: str) -> List[str]:
        return archive_id.split('/')

    @staticmethod
    def item(archive_id: str, index: int) -> str:
        return archive.items(archive_id)[index]

    @staticmethod
    def calc_id(archive_id: str) -> str:
        return archive.item(archive_id, 1)

    @staticmethod
    def upload_id(archive_id: str) -> str:
        return archive.item(archive_id, 0)


def to_tuple(self, *args):
    return tuple(self[arg] for arg in args)


def chunks(list, n):
    """ Chunks up the given list into parts of size n. """
    for i in range(0, len(list), n):
        yield list[i:i + n]


class POPO(dict):
    """
    A dict subclass that uses attributes as key/value pairs.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def __getattr__(self, name):
        if name in self:
            return self[name]
        else:
            raise AttributeError("No such attribute: " + name)

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        if name in self:
            del self[name]
        else:
            raise AttributeError("No such attribute: " + name)


class SleepTimeBackoff:
    """
    Provides increasingly larger sleeps. Useful when
    observing long running processes with unknown runtime.
    """

    def __init__(self, start_time: float = 0.1, max_time: float = 60):
        self.current_time = start_time
        self.max_time = max_time

    def __call__(self):
        self.sleep()

    def sleep(self):
        time.sleep(self.current_time)
        self.current_time *= 2
        self.current_time = min(self.max_time, self.current_time)
