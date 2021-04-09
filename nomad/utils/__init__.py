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
.. autofunc::nomad.utils.strip
'''

from typing import List, Iterable
from collections import OrderedDict
import base64
from contextlib import contextmanager
import json
import uuid
import time
import re
import hashlib
import sys
from datetime import timedelta
import collections
import logging
import inspect

from nomad import config

default_hash_len = 28
''' Length of hashes and hash-based ids (e.g. calc, upload) in nomad. '''

try:
    from . import structlogging
    from .structlogging import legacy_logger
    from .structlogging import configure_logging

    def get_logger(name, **kwargs):
        '''
        Returns a structlog logger that is already attached with a logstash handler.
        Use additional *kwargs* to pre-bind some values to all events.
        '''
        return structlogging.get_logger(name, **kwargs)

except ImportError:
    def get_logger(name, **kwargs):
        return ClassicLogger(name, **kwargs)

    def configure_logging(console_log_level=config.console_log_level):
        import logging
        logging.basicConfig(level=console_log_level)


class ClassicLogger:
    '''
    A logger class that emulates the structlog interface, but uses the classical
    build-in Python logging.
    '''
    def __init__(self, name, **kwargs):
        self.kwargs = kwargs
        self.logger = logging.getLogger(name)

    def bind(self, **kwargs):
        all_kwargs = dict(self.kwargs)
        all_kwargs.update(**kwargs)
        return ClassicLogger(self.logger.name, **all_kwargs)

    def __log(self, method_name, event, **kwargs):
        method = getattr(self.logger, method_name)
        all_kwargs = dict(self.kwargs)
        all_kwargs.update(**kwargs)

        message = '%s (%s)' % (
            event,
            ', '.join(['%s=%s' % (str(key), str(value)) for key, value in all_kwargs.items()])
        )
        method(message)

    def __getattr__(self, key):
        return lambda *args, **kwargs: self.__log(key, *args, **kwargs)


def set_console_log_level(level):
    root = logging.getLogger()
    for handler in root.handlers:
        if isinstance(handler, (logging.StreamHandler, logging.FileHandler)):
            handler.setLevel(level)


def decode_handle_id(handle_str: str):
    result = 0
    for c in handle_str:
        ordinal = ord(c.lower())
        if 48 <= ordinal <= 57:
            number = ordinal - 48
        elif 97 <= ordinal <= 118:
            number = ordinal - 87
        else:
            raise ValueError()

        result = result * 32 + number

    return result


def hash(*args, length: int = default_hash_len) -> str:
    ''' Creates a websafe hash of the given length based on the repr of the given arguments. '''
    hash = hashlib.sha512()
    for arg in args:
        hash.update(str(arg).encode('utf-8'))

    return make_websave(hash, length=length)


def make_websave(hash, length: int = default_hash_len) -> str:
    ''' Creates a websafe string for a hashlib hash object. '''
    if length > 0:
        return base64.b64encode(hash.digest(), altchars=b'-_')[:length].decode('utf-8')
    else:
        return base64.b64encode(hash.digest(), altchars=b'-_')[0:-2].decode('utf-8')


def base64_encode(string):
    '''
    Removes any `=` used as padding from the encoded string.
    '''
    encoded = base64.urlsafe_b64encode(string).decode('utf-8')
    return encoded.rstrip("=")


def base64_decode(string):
    '''
    Adds back in the required padding before decoding.
    '''
    padding = 4 - (len(string) % 4)
    bytes = (string + ("=" * padding)).encode('utf-8')
    return base64.urlsafe_b64decode(bytes)


def create_uuid() -> str:
    ''' Returns a web-save base64 encoded random uuid (type 4). '''
    return base64.b64encode(uuid.uuid4().bytes, altchars=b'-_').decode('utf-8')[0:-2]


def adjust_uuid_size(uuid, length: int = default_hash_len):
    ''' Adds prefixing spaces to a uuid to ensure the default uuid length. '''
    uuid = uuid.rjust(length, ' ')
    assert len(uuid) == length, 'uuids must have the right fixed size'
    return uuid


@contextmanager
def lnr(logger, event, **kwargs):
    '''
    A context manager that Logs aNd Raises all exceptions with the given logger.

    Arguments:
        logger: The logger that should be used for logging exceptions.
        event: the log message
        **kwargs: additional properties for the structured log
    '''
    try:
        yield

    except Exception as e:
        # ignore HTTPException as they are part of the normal flask error handling
        if e.__class__.__name__ == 'HTTPException':
            logger.error(event, exc_info=e, **kwargs)
        raise e


@contextmanager
def timer(logger, event, method='info', **kwargs):
    '''
    A context manager that takes execution time and produces a log entry with said time.

    Arguments:
        logger: The logger that should be used to produce the log entry.
        event: The log message/event.
        method: The log method that should be used. Must be a valid logger method name.
            Default is 'info'.
        **kwargs: Additional logger data that is passed to the log entry.

    Returns:
        The method yields a dictionary that can be used to add further log data.
    '''
    start = time.time()

    try:
        yield kwargs
    finally:
        stop = time.time()

    if logger is None:
        print(event, stop - start)
        return

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
    ''' Chunks up the given list into parts of size n. '''
    for i in range(0, len(list), n):
        yield list[i:i + n]


class SleepTimeBackoff:
    '''
    Provides increasingly larger sleeps. Useful when
    observing long running processes with unknown runtime.
    '''

    def __init__(self, start_time: float = 0.1, max_time: float = 5):
        self.current_time = start_time
        self.max_time = max_time

    def __call__(self):
        self.sleep()

    def sleep(self):
        time.sleep(self.current_time)
        self.current_time *= 2
        self.current_time = min(self.max_time, self.current_time)


class ETA:
    def __init__(self, total: int, message: str, interval: int = 1000):
        self.start = time.time()
        self.total = total
        self.count = 0
        self.interval = interval
        self.interval_count = 0
        self.message = message

    def add(self, amount: int = 1):
        self.count += amount
        interval_count = int(self.count / self.interval)
        if interval_count > self.interval_count:
            self.interval_count = interval_count
            delta_t = time.time() - self.start
            eta = delta_t * (self.total - self.count) / self.count
            eta_str = str(timedelta(seconds=eta))
            sys.stdout.write('\r' + (self.message % (self.count, self.total, eta_str)))
            sys.stdout.flush()

    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args, **kwargs):
        print('')


def common_prefix(paths):
    '''
    Computes the longest common file path prefix (with respect to '/' separated segments).
    Returns empty string is ne common prefix exists.
    '''
    common_prefix = None

    for path in paths:
        if common_prefix is None:
            common_prefix = path

        index = 0
        index_last_slash = -1
        for a, b in zip(path, common_prefix):
            if a != b:
                break
            if a == '/':
                index_last_slash = index
            index += 1

        if index_last_slash == -1:
            common_prefix = ''
            break

        common_prefix = common_prefix[:index_last_slash + 1]

    if common_prefix is None:
        common_prefix = ''

    return common_prefix


class RestrictedDict(OrderedDict):
    """Dictionary-like container with predefined set of mandatory and optional
    keys and a set of forbidden values.
    """
    def __init__(self, mandatory_keys: Iterable = None, optional_keys: Iterable = None, forbidden_values: Iterable = None, lazy: bool = True):
        """
        Args:
            mandatory_keys: Keys that have to be present.
            optional_keys: Keys that are optional.
            forbidden_values: Values that are forbidden. Only supports hashable values.
            lazy: If false, the values are checked already when inserting. If
                True, the values should be manually checked by calling the
                check()-function.
        """
        super().__init__()

        if isinstance(mandatory_keys, (list, tuple, set)):
            self._mandatory_keys = set(mandatory_keys)
        elif mandatory_keys is None:
            self._mandatory_keys = set()
        else:
            raise ValueError("Please provide the mandatory_keys as a list, tuple or set.")

        if isinstance(optional_keys, (list, tuple, set)):
            self._optional_keys = set(optional_keys)
        elif optional_keys is None:
            self._optional_keys = set()
        else:
            raise ValueError("Please provide the optional_keys as a list, tuple or set.")

        if isinstance(forbidden_values, (list, tuple, set)):
            self._forbidden_values = set(forbidden_values)
        elif forbidden_values is None:
            self._forbidden_values = set()
        else:
            raise ValueError("Please provide the forbidden_values as a list or tuple of values.")

        self._lazy = lazy

    def __setitem__(self, key, value):
        if not self._lazy:

            # Check that only the defined keys are used
            if key not in self._mandatory_keys and key not in self._optional_keys:
                raise KeyError("The key '{}' is not allowed.".format(key))

            # Check that forbidden values are not used.
            try:
                match = value in self._forbidden_values
            except TypeError:
                pass  # Unhashable value will not match
            else:
                if match:
                    raise ValueError("The value '{}' is not allowed.".format(key))

        super().__setitem__(key, value)

    def check(self, recursive=False):
        # Check that only the defined keys are used
        for key in self.keys():
            if key not in self._mandatory_keys and key not in self._optional_keys:
                raise KeyError("The key '{}' is not allowed.".format(key))

        # Check that all mandatory values are all defined
        for key in self._mandatory_keys:
            if key not in self:
                raise KeyError("The mandatory key '{}' is not present.".format(key))

        # Check that forbidden values are not used.
        for key, value in self.items():
            match = False
            try:
                match = value in self._forbidden_values
            except TypeError:
                pass  # Unhashable value will not match
            else:
                if match:
                    raise ValueError("The value '{}' is not allowed but was set for key '{}'.".format(value, key))

        # Check recursively
        if recursive:
            for value in self.values():
                if isinstance(value, RestrictedDict):
                    value.check(recursive)

    def update(self, other):
        for key, value in other.items():
            self.__setitem__(key, value)

    def hash(self) -> str:
        """Creates a hash code from the contents. Ensures consistent ordering.
        """
        hash_str = json.dumps(self, sort_keys=True)

        return hash(hash_str)


def strip(docstring):
    ''' Removes any unnecessary whitespaces from a multiline doc string or description. '''
    return inspect.cleandoc(docstring)
