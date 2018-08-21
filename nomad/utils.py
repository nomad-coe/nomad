from typing import Union, IO, cast
import hashlib
import base64
import json
import logging
from contextlib import contextmanager


def hash(obj: Union[IO, str]) -> str:
    """ First 28 character of an URL safe base 64 encoded sha512 digest. """
    hash = hashlib.sha512()
    if getattr(obj, 'read', None) is not None:
        for data in iter(lambda: cast(IO, obj).read(65536), b''):
            hash.update(data)
    elif isinstance(obj, str):
        hash.update(obj.encode('utf-8'))

    return base64.b64encode(hash.digest(), altchars=b'-_')[0:28].decode('utf-8')


class DataLogger():
    def __init__(self, logger, **kwargs):
        self._logger = logger
        self.data = kwargs

    def _prepare_msg(self, base_msg):
        return '%s %s' % (base_msg, self._format_data())

    def _format_data(self, ):
        return json.dumps(self.data)

    def debug(self, msg, *args, **kwargs):
        self._logger.debug(self._prepare_msg(msg), *args, **kwargs)

    def info(self, msg, *args, **kwargs):
        self._logger.info(self._prepare_msg(msg), *args, **kwargs)

    def warn(self, msg, *args, **kwargs):
        self._logger.warn(self._prepare_msg(msg), *args, **kwargs)

    def error(self, msg, *args, **kwargs):
        self._logger.error(self._prepare_msg(msg), *args, **kwargs)

    def crit(self, msg, *args, **kwargs):
        self._logger.crit(self._prepare_msg(msg), *args, **kwargs)

    @contextmanager
    def lnr_error(self, msg, *args, **kwargs):
        """
        Will *log and raise* with an error and the given message and args/kwargs
        on all exceptions.
        """
        try:
            yield
        except Exception as e:
            self._logger.error(
                self._prepare_msg('Exception while: %s' % msg),
                exc_info=e, *args, **kwargs)
            raise e


def get_logger(name, *args, **kwargs):
    """
    Returns a :class:`DataLogger` with the data given as kwargs.
    A data logger can be used like any other logger, but will add the data to all
    log output. Allowing more structured logging.
    """
    return DataLogger(logging.getLogger(name), *args, **kwargs)
