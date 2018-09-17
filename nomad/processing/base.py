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

from typing import List, cast, Any
import types
from contextlib import contextmanager
import collections
import inspect
import logging
import time
import celery
from celery import Celery, Task
from celery.signals import after_setup_task_logger, after_setup_logger, worker_process_init
from mongoengine import Document, StringField, ListField, DateTimeField, IntField, \
    ReferenceField, connect, ValidationError
from mongoengine.connection import MongoEngineConnectionError
from mongoengine.base.metaclasses import TopLevelDocumentMetaclass
from pymongo import ReturnDocument
from datetime import datetime
import sys

from nomad import config, utils
import nomad.patch  # pylint: disable=unused-import


def mongo_connect():
    return connect(db=config.mongo.users_db, host=config.mongo.host, port=config.mongo.port)

if config.logstash.enabled:
    def initialize_logstash(logger=None, loglevel=logging.DEBUG, **kwargs):
        utils.add_logstash_handler(logger)
        return logger

    after_setup_task_logger.connect(initialize_logstash)
    after_setup_logger.connect(initialize_logstash)

worker_process_init.connect(lambda **kwargs: mongo_connect())

app = Celery('nomad.processing', broker=config.celery.broker_url)

# ensure elastic and mongo connections
if 'sphinx' not in sys.modules:
    connect(db=config.mongo.users_db, host=config.mongo.host, port=config.mongo.port)


PENDING = 'PENDING'
RUNNING = 'RUNNING'
FAILURE = 'FAILURE'
SUCCESS = 'SUCCESS'


class InvalidId(Exception): pass


class ProcNotRegistered(Exception): pass


class ProcMetaclass(TopLevelDocumentMetaclass):
    def __new__(cls, name, bases, attrs):
        cls = super().__new__(cls, name, bases, attrs)

        tasks = []
        setattr(cls, 'tasks', tasks)

        for name, attr in attrs.items():
            task = getattr(attr, '__task_name', None)
            if task is not None and task not in tasks:
                tasks.append(task)

        return cls


class Proc(Document, metaclass=ProcMetaclass):
    """
    Base class for objects that are involved in processing and need persistent processing
    state.

    It solves two issues. First, distributed operation (via celery) and second keeping
    state of a chain of potentially failing processing tasks. Both are controlled via
    decorators @process and @task. Subclasses should use these decorators on their
    methods. Parameters are not supported for decorated functions. Use fields on the
    document instead.

    Processing state will be persistet at appropriate
    times and must not be persistet manually. All attributes are stored to mongodb.
    The class allows to render into a JSON serializable dict via :attr:`json_dict`.

    Possible processing states are PENDING, RUNNING, FAILURE, and SUCCESS.

    Attributes:
        current_task: the currently running or last completed task
        status: the overall status of the processing
        errors: a list of errors that happened during processing. Error fail a processing
            run
        warnings: a list of warnings that happened during processing. Warnings do not
            fail a processing run
        create_time: the time of creation (not the start of processing)
        proc_time: the time that processing completed (successfully or not)
    """

    meta: Any = {
        'abstract': True,
    }

    tasks: List[str] = None
    """ the ordered list of tasks that comprise a processing run """

    current_task = StringField(default=None)
    status = StringField(default='CREATED')

    errors = ListField(StringField())
    warnings = ListField(StringField())

    create_time = DateTimeField(required=True)
    complete_time = DateTimeField()

    _async_status = StringField(default='UNCALLED')

    @property
    def completed(self) -> bool:
        """ Returns True of the process has failed or succeeded. """
        return self.status in [SUCCESS, FAILURE]

    def get_logger(self):
        return utils.get_logger(
            __name__, current_task=self.current_task, process=self.__class__.__name__,
            status=self.status)

    @classmethod
    def create(cls, **kwargs):
        """ Factory method that must be used instead of regular constructor. """
        assert cls.tasks is not None and len(cls.tasks) > 0, \
            """ the class attribute tasks must be overwritten with an acutal list """
        assert 'status' not in kwargs, \
            """ do not set the status manually, its managed """

        kwargs.setdefault('create_time', datetime.now())
        self = cls(**kwargs)
        self.status = PENDING if self.current_task is None else RUNNING
        self.save()

        return self

    @classmethod
    def get_by_id(cls, id: str, id_field: str):
        try:
            obj = cls.objects(**{id_field: id}).first()
        except ValidationError as e:
            raise InvalidId('%s is not a valid id' % id)
        except MongoEngineConnectionError as e:
            raise e

        if obj is None:
            raise KeyError('%s with id %s does not exist' % (cls.__name__, id))

        return obj

    @classmethod
    def get(cls, obj_id):
        return cls.get_by_id(str(obj_id), 'id')

    @staticmethod
    def log(logger, log_level, msg, **kwargs):
        # TODO there seems to be a bug in structlog, cannot use logger.log
        if log_level == logging.ERROR:
            logger.error(msg, **kwargs)
        elif log_level == logging.WARNING:
            logger.warning(msg, **kwargs)
        elif log_level == logging.INFO:
            logger.info(msg, **kwargs)
        elif log_level == logging.DEBUG:
            logger.debug(msg, **kwargs)
        else:
            logger.critical(msg, **kwargs)

    def fail(self, *errors, log_level=logging.ERROR, **kwargs):
        """ Allows to fail the process. Takes strings or exceptions as args. """
        assert not self.completed, 'Cannot fail a completed process.'

        failed_with_exception = False

        self.status = FAILURE

        logger = self.get_logger(**kwargs)
        for error in errors:
            if isinstance(error, Exception):
                failed_with_exception = True
                Proc.log(logger, log_level, 'task failed with exception', exc_info=error, **kwargs)

        self.errors = [str(error) for error in errors]
        self.complete_time = datetime.now()

        if not failed_with_exception:
            errors_str = "; ".join([str(error) for error in errors])
            Proc.log(logger, log_level, 'task failed', errors=errors_str, **kwargs)

        logger.debug('process failed')

        self.save()

    def warning(self, *warnings, log_level=logging.warning, **kwargs):
        """ Allows to save warnings. Takes strings or exceptions as args. """
        assert not self.completed

        logger = self.get_logger(**kwargs)

        for warning in warnings:
            warning = str(warning)
            self.warnings.append(warning)
            Proc.log(logger, log_level, 'task with warning', warning=warning)

    def _continue_with(self, task):
        tasks = self.__class__.tasks
        assert task in tasks, 'task %s must be one of the classes tasks %s' % (task, str(tasks))
        if self.current_task is None:
            assert task == tasks[0], "process has to start with first task"
        else:
            assert tasks.index(task) == tasks.index(self.current_task) + 1, \
                "tasks must be processed in the right order"

        if self.status == FAILURE:
            return False

        if self.status == PENDING:
            assert self.current_task is None
            assert task == tasks[0]
            self.status = RUNNING
            self.current_task = task
            self.get_logger().debug('started process')
        else:
            self.current_task = task
            self.get_logger().debug('successfully completed task')

        self.save()
        return True

    def _complete(self):
        if self.status != FAILURE:
            assert self.status == RUNNING, 'Can only complete a running process.'
            self.status = SUCCESS
            self.save()
            self.get_logger().debug('completed process')

    def incr_counter(self, field, value=1, other_fields=None):
        """
        Atomically increases the given field by value and return the new value.
        Optionally return also other values from the updated object to avoid
        reloads.

        Arguments:
            field: the name of the field to increment, must be a :class:`IntField`
            value: the value to increment the field by, default is 1
            other_fields: an optional list of field names that should also be returned

        Returns:
            either the value of the updated field, or a tuple with this value and a list
            of value for the given other fields
        """
        # use a primitive but atomic pymongo call
        updated_raw = self._get_collection().find_one_and_update(
            {'_id': self.id},
            {'$inc': {field: value}},
            return_document=ReturnDocument.AFTER)

        if updated_raw is None:
            raise KeyError('object does not exist, was probaly not yet written to db')

        if other_fields is None:
            return updated_raw[field]
        else:
            return updated_raw[field], [updated_raw[field] for field in other_fields]

    def block_until_complete(self, interval=0.01):
        """
        Reloads the process constrantly until it sees a completed process. Should be
        used with care as it can block indefinetly. Just intended for testing purposes.
        """
        while not self.completed:
            time.sleep(interval)
            self.reload()

    @property
    def json_dict(self) -> dict:
        """ A json serializable dictionary representation. """
        data = {
            'tasks': getattr(self.__class__, 'tasks'),
            'current_task': self.current_task,
            'status': self.status,
            'completed': self.completed,
            'errors': self.errors,
            'warnings': self.warnings,
            'create_time': self.create_time.isoformat() if self.create_time is not None else None,
            'complete_time': self.complete_time.isoformat() if self.complete_time is not None else None,
            '_async_status': self._async_status
        }
        return {key: value for key, value in data.items() if value is not None}


def task(func):
    """
    The decorator for tasks that will be wrapped in excaption handling that will fail the process.
    The task methods of a :class:`Proc` class/document comprise a sequence
    (order of methods in class namespace) of tasks. Tasks must be executed in that order.
    Completion of the last task, will put the :class:`Proc` instance into the
    SUCCESS state. Calling the first task will put it into RUNNING state. Tasks will
    only be exectued, if the process has not yet reached FAILURE state.
    """
    def wrapper(self, *args, **kwargs):
        if self.status == 'FAILURE':
            return

        self._continue_with(func.__name__)
        try:
            func(self, *args, **kwargs)
        except Exception as e:
            self.fail(e)

        if self.__class__.tasks[-1] == self.current_task and not self.completed:
            self._complete()

    setattr(wrapper, '__task_name', func.__name__)
    wrapper.__name__ = func.__name__
    return wrapper


@app.task(bind=True, ignore_results=True, max_retries=3)
def proc_task(task, cls_name, self_id, func_attr):
    """
    The celery task that is used to execute async process functions.
    It ignores results, since all results are handled via the self document.
    It retries for 3 times with a countdown of 3 on missing 'selfs', since this
    might happen in sharded, distributed mongo setups where the object might not
    have yet been propagated and therefore apear missing.
    """
    logger = utils.get_logger(__name__, cls=cls_name, id=self_id, func=func_attr)

    logger.debug('received process function call')
    all_cls = Proc.__subclasses__()
    cls = next((cls for cls in all_cls if cls.__name__ == cls_name), None)
    if cls is None:
        logger.error('document not a subcass of Proc')
        raise ProcNotRegistered('document %s not a subclass of Proc' % cls_name)

    try:
        self = cls.get(self_id)
    except KeyError as e:
        logger.warning('called object is missing')
        raise task.retry(exc=e, countdown=3)

    func = getattr(self, func_attr, None)
    if func is None:
        logger.error('called function not a function of proc class')
        self.fail('called function %s is not a function of proc class %s' % (func_attr, cls_name))
        return

    func = getattr(func, '__process_unwrapped', None)
    if func is None:
        logger.error('called function was not decorated with @process')
        self.fail('called function %s was not decorated with @process' % (func_attr, cls_name))
        return

    try:
        self._async_status = 'RECEIVED-%s' % func.__name__
        func(self)
    except Exception as e:
        self.fail(e)


def process(func):
    """
    The decorator for process functions that will be called async via celery.
    All calls to the decorated method will result in celery task requests.
    To transfer state, the instance will be saved to the database and loading on
    the celery task worker. Process methods can call other (process) functions/methods.
    """
    def wrapper(self, *args, **kwargs):
        assert len(args) == 0 and len(kwargs) == 0, 'process functions must not have arguments'
        self._async_status = 'CALLED-%s' % func.__name__
        self.save()

        self_id = self.id.__str__()
        cls_name = self.__class__.__name__

        logger = utils.get_logger(__name__, cls=cls_name, id=self_id, func=func.__name__)
        logger.debug('calling process function')
        return proc_task.s(cls_name, self_id, func.__name__).delay()

    task = getattr(func, '__task_name', None)
    if task is not None:
        setattr(wrapper, '__task_name', task)
    wrapper.__name__ = func.__name__
    setattr(wrapper, '__process_unwrapped', func)

    return wrapper
