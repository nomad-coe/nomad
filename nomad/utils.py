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
parser, parser version, calc_hash, mainfile, etc. for all events that happen during
calculation processing. This means the :func:`get_logger` and all logger functions
take keyword arguments for structured data. Otherwise :func:`get_logger` can
be used similar to the standard *logging.getLogger*.

Depending on the configuration all logs will also be send to a central logstash.

.. autofunc::nomad.utils.get_logger
"""

from typing import Union, IO, cast
import hashlib
import base64
import logging
import structlog
from structlog.processors import StackInfoRenderer, format_exc_info, TimeStamper, JSONRenderer
from structlog.stdlib import LoggerFactory
import logstash
from contextlib import contextmanager
import json
import os
import sys
import uuid
import time

from nomad import config

_service = os.environ.get('NOMAD_SERVICE', 'nomad service')


class LogstashFormatterVersion1ForStructlog(logstash.formatter.LogstashFormatterBase):

    def format(self, record):
        try:
            structlog = json.loads(record.getMessage())
        except json.JSONDecodeError:
            structlog = dict(event=record.getMessage())

        # Create message dict
        message = {
            '@timestamp': self.format_timestamp(record.created),
            '@version': '1',
            'message': structlog['event'],
            'host': self.host,
            'path': record.pathname,
            'tags': self.tags,
            'type': self.message_type,

            # Extra Fields
            'level': record.levelname,
            'logger_name': record.name,
        }

        message.update(structlog)

        # Add extra fields
        message.update(self.get_extra_fields(record))

        # If exception, add debug info
        if record.exc_info:
            message.update(self.get_debug_fields(record))

        return self.serialize(message)


def add_logstash_handler(logger):
    logstash_handler = logstash.TCPLogstashHandler(
        config.logstash.host,
        config.logstash.tcp_port, version=1)
    logstash_handler.formatter = LogstashFormatterVersion1ForStructlog(tags=['nomad', _service])
    logstash_handler.setLevel(config.logstash.level)
    logger.addHandler(logstash_handler)


_logging_is_configured = False
if not _logging_is_configured:
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
    logging.basicConfig(stream=sys.stdout)
    root = logging.getLogger()
    for handler in root.handlers:
        handler.setLevel(config.console_log_level if 'pytest' not in sys.modules else logging.CRITICAL)

    # configure logstash
    if config.logstash.enabled and 'pytest' not in sys.modules:
        add_logstash_handler(root)
        root.info('Structlog configured for logstash')

    _logging_is_configured = True


def create_uuid() -> str:
    """ Returns a web-save base64 encoded random uuid (type 4). """
    return base64.b64encode(uuid.uuid4().bytes, altchars=b'-_').decode('utf-8')[0:-2]


def hash(obj: Union[IO, str]) -> str:
    """
    Returns a web-save base64 encoded 28 long hash for the given contents.
    First 28 character of an URL safe base 64 encoded sha512 digest.
    """
    hash = hashlib.sha512()
    if getattr(obj, 'read', None) is not None:
        for data in iter(lambda: cast(IO, obj).read(65536), b''):
            hash.update(data)
    elif isinstance(obj, str):
        hash.update(obj.encode('utf-8'))

    return base64.b64encode(hash.digest(), altchars=b'-_')[0:28].decode('utf-8')


def get_logger(name, **kwargs):
    """
    Returns a structlog logger that is already attached with a logstash handler.
    Use additional *kwargs* to pre-bind some values to all events.
    """
    logger = structlog.get_logger(name, service=_service, **kwargs)
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
    except Exception as e:
        logger.error(event, exc_info=e, **kwargs)
        raise e


@contextmanager
def timer(logger, event, method='info', **kwargs):
    start = time.time()

    try:
        yield
    finally:
        stop = time.time()

    logger_method = getattr(logger, 'info', None)
    if logger_method is not None:
        logger_method(event, exec_time=stop - start, **kwargs)
    else:
        logger.error('Uknown logger method %s.' % method)
