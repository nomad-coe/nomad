from typing import Union, IO, cast
import hashlib
import base64
import logging
import structlog
from structlog.processors import StackInfoRenderer, format_exc_info, TimeStamper, JSONRenderer
from structlog.stdlib import LoggerFactory
import logstash
from contextlib import contextmanager

from nomad import config


_logging_is_configured = False
if not _logging_is_configured:
    # basic config
    logging.basicConfig(level=logging.WARNING)

    # configure logstash
    if config.logstash.enabled:
        logstash_handler = logstash.TCPLogstashHandler(
            config.logstash.host,
            config.logstash.tcp_port, version=1)
        logstash_handler.setLevel(config.logstash.level)
        logging.getLogger().addHandler(logstash_handler)

    # configure structlog
    log_processors = [
        StackInfoRenderer(),
        format_exc_info,
        TimeStamper(fmt="%Y-%m-%d %H:%M.%S", utc=False),
        JSONRenderer(sort_keys=True)
    ]
    structlog.configure(processors=log_processors, logger_factory=LoggerFactory())
    _logging_is_configured = True


def hash(obj: Union[IO, str]) -> str:
    """ First 28 character of an URL safe base 64 encoded sha512 digest. """
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
    User additional *kwargs* to pre-bind some values.
    """
    logger = structlog.get_logger(**kwargs)
    return logger


@contextmanager
def lnr(logger, event, **kwargs):
    try:
        yield
    except Exception as e:
        logger.error(event, exc_info=e, **kwargs)
        raise e


class DataObject(dict):
    """
    A simpe data class base that allows to create javascript style objects.
    Is also json serializable, if you only but json serializable contents into it.
    """
    def __getattr__(self, name):
        try:
            return self.__getitem__(name)
        except KeyError:
            raise AttributeError

    def __setattr__(self, name, val):
        return self.__setitem__(name, val)

    def __delattr__(self, name):
        assert name != 'data'
        return self.__delitem__(name)

    def update(self, dct):
        return super().update({key: value for key, value in dct.items() if value is not None})

if __name__ == '__main__':
    logger = get_logger(__name__, test='value')
    logger.info('Hi', add='cool')