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

'''
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
'''

import logging
import structlog
from structlog.processors import StackInfoRenderer, format_exc_info, TimeStamper, JSONRenderer
from structlog.stdlib import LoggerFactory
import logstash
from contextlib import contextmanager
import json
import re

from nomad import config


def sanitize_logevent(event: str) -> str:
    '''
    Prepares a log event or message for analysis in elastic stack. It removes numbers,
    list, and matrices of numbers from the event string and limits its size. The
    goal is to make it easier to define aggregations over events by using event
    strings as representatives for event classes rather than event instances (with
    concrete numbers, etc).
    '''
    sanitized_event = event[:120]
    sanitized_event = re.sub(r'(\d*\.\d+|\d+(\.\d*)?)', 'X', sanitized_event)
    sanitized_event = re.sub(r'((\[|\()\s*)?X\s*(,\s*X)+(\s*(\]|\)))?', 'L', sanitized_event)
    sanitized_event = re.sub(r'((\[|\()\s*)?[XL](,\s*[XL])+(\s*(\]|\)))?', 'M', sanitized_event)
    return sanitized_event


@contextmanager
def legacy_logger(logger):
    ''' Context manager that makes the given logger the logger for legacy log entries. '''
    LogstashHandler.legacy_logger = logger
    try:
        yield
    finally:
        LogstashHandler.legacy_logger = None


class LogstashHandler(logstash.TCPLogstashHandler):
    '''
    A log handler that emits records to logstash. It also filters logs for being
    structlog entries. All other entries are diverted to a global `legacy_logger`.
    This legacy logger is supposed to be a structlog logger that turns legacy
    records into structlog entries with reasonable binds depending on the current
    execution context (e.g. parsing/normalizing, etc.). If no legacy logger is
    set, they get emitted as usual (e.g. non nomad logs, celery, dbs, etc.)
    '''

    legacy_logger = None

    def filter(self, record):
        if record.name == 'gunicorn.access' and 'alive' in record.args.get('r', ''):
            return False

        if super().filter(record):
            is_structlog = False
            if isinstance(record.msg, str):
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


_gunicorn_pattern_parts = [
    r'(?P<host>\S+)',  # host %h
    r'\S+',  # indent %l (unused)
    r'(?P<user>\S+)',  # user %u
    r'\[(?P<time>.+)\]',  # time %t
    r'"(?P<request>.+)"',  # request "%r"
    r'(?P<status>[0-9]+)',  # status %>s
    r'(?P<size>\S+)',  # size %b (careful, can be '-')
    r'"(?P<referer>.*)"',  # referer "%{Referer}i"
    r'"(?P<agent>.*)"',  # user agent "%{User-agent}i"
]
_gunicorn_pattern = re.compile(r'\s+'.join(_gunicorn_pattern_parts) + r'\s*\Z')


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
            'host': self.host,
            'path': record.pathname,
            'tags': self.tags,
            'type': self.message_type,

            # Extra Fields
            'level': record.levelname,
            'logger_name': record.name,

            # Nomad specific
            'nomad.service': config.meta.service,
            'nomad.release': config.meta.release
        }

        if record.name.startswith('nomad'):
            for key, value in structlog.items():
                if key in ('event', 'stack_info', 'id', 'timestamp'):
                    continue
                elif key == 'exception':
                    message['digest'] = str(value)[-256:]
                elif key in ['upload_id', 'calc_id', 'mainfile']:
                    key = 'nomad.%s' % key
                else:
                    key = '%s.%s' % (record.name, key)

                message[key] = value
        else:
            message.update(structlog)

        # Handle gunicorn access events
        if record.name == 'gunicorn.access':
            gunicorn_message = structlog['event']
            gunicorn_record = _gunicorn_pattern.match(gunicorn_message).groupdict()

            if gunicorn_record['user'] == '-':
                gunicorn_record['user'] = None

            gunicorn_record['status'] = int(gunicorn_record['status'])

            if gunicorn_record['size'] == '-':
                gunicorn_record['size'] = 0
            else:
                gunicorn_record['size'] = int(gunicorn_record['size'])

            if gunicorn_record['referer'] == '-':
                gunicorn_record['referer'] = None

            message.update({'gunicorn.%s' % key: value for key, value in gunicorn_record.items()})
            message['event'] = gunicorn_record['request']

        # Add extra fields
        message.update(self.get_extra_fields(record))

        # If exception, add debug info
        if record.exc_info:
            message.update(self.get_debug_fields(record))

        return self.serialize(message)


class ConsoleFormatter(LogstashFormatter):

    short_format = False

    @classmethod
    def serialize(cls, message_dict):
        from io import StringIO

        logger = message_dict.pop('logger_name', 'unknown logger')
        event = message_dict.pop('event', None)
        level = message_dict.pop('level', 'UNKNOWN')
        exception = message_dict.pop('exception', None)
        time = message_dict.pop('@timestamp', '1970-01-01 12:00:00')

        for key in ['type', 'tags', 'stack_info', 'path', 'message', 'host', '@version', 'digest']:
            message_dict.pop(key, None)
        keys = list(message_dict.keys())
        keys.sort()

        out = StringIO()
        out.write('%s %s %s %s' % (
            level.ljust(8), logger.ljust(20)[:20], time.ljust(19)[:19], event))
        if exception is not None:
            out.write('\n  - exception: %s' % str(exception).replace('\n', '\n    '))

        for key in keys:
            if cls.short_format and key.startswith('nomad.'):
                print_key = key[6:]
            else:
                print_key = key
            if not cls.short_format or print_key not in ['release', 'service']:
                out.write('\n  - %s: %s' % (print_key, str(message_dict.get(key, None))))
        return out.getvalue()


def add_logstash_handler(logger):
    logstash_handler = next((
        handler for handler in logger.handlers
        if isinstance(handler, LogstashHandler)), None)

    if logstash_handler is None:
        logstash_handler = LogstashHandler(
            config.logstash.host,
            config.logstash.tcp_port, version=1)
        logstash_handler.formatter = LogstashFormatter(tags=['nomad', config.meta.release])
        logstash_handler.setLevel(config.logstash.level)
        logger.addHandler(logstash_handler)


def get_logger(name, **kwargs):
    '''
    Returns a structlog logger that is already attached with a logstash handler.
    Use additional *kwargs* to pre-bind some values to all events.
    '''
    if name.startswith('nomad.'):
        name = '.'.join(name.split('.')[:2])

    logger = structlog.get_logger(name, **kwargs)
    return logger


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

    get_logger(__name__).info(
        'setup logging',
        logstash=config.logstash.enabled,
        logstash_host=config.logstash.host,
        logstash_port=config.logstash.tcp_port,
        logstash_level=config.logstash.level)

# configure log levels
for logger in [
        'elasticsearch',
        # 'celery.app.trace', 'celery.worker.strategy',
        'urllib3.connectionpool', 'bravado', 'bravado_core', 'swagger_spec_validator']:
    logging.getLogger(logger).setLevel(logging.WARNING)
