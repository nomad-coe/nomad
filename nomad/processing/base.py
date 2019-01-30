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

from typing import List, Any
import logging
import time
from celery import Celery
from celery.signals import after_setup_task_logger, after_setup_logger, worker_process_init
from mongoengine import Document, StringField, ListField, DateTimeField, IntField, \
    ValidationError, BooleanField
from mongoengine.connection import MongoEngineConnectionError
from mongoengine.base.metaclasses import TopLevelDocumentMetaclass
from pymongo import ReturnDocument
from datetime import datetime

from nomad import config, utils, infrastructure
import nomad.patch  # pylint: disable=unused-import


if config.logstash.enabled:
    utils.configure_logging()

    def initialize_logstash(logger=None, loglevel=logging.DEBUG, **kwargs):
        utils.add_logstash_handler(logger)
        return logger

    after_setup_task_logger.connect(initialize_logstash)
    after_setup_logger.connect(initialize_logstash)


@worker_process_init.connect
def setup(**kwargs):
    infrastructure.setup()


app = Celery('nomad.processing', broker=config.celery.broker_url)
app.conf.update(worker_hijack_root_logger=False)

CREATED = 'CREATED'
PENDING = 'PENDING'
RUNNING = 'RUNNING'
FAILURE = 'FAILURE'
SUCCESS = 'SUCCESS'

PROCESS_CALLED = 'CALLED'
PROCESS_RUNNING = 'RUNNING'
PROCESS_COMPLETED = 'COMPLETED'


class InvalidId(Exception): pass


class ProcNotRegistered(Exception): pass


class ProcessAlreadyRunning(Exception): pass


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

    Possible processing states are PENDING, RUNNING, FAILURE, and SUCCESS.

    Attributes:
        current_task: the currently running or last completed task
        tasks_status: the overall status of the processing
        errors: a list of errors that happened during processing. Error fail a processing
            run
        warnings: a list of warnings that happened during processing. Warnings do not
            fail a processing run
        create_time: the time of creation (not the start of processing)
        complete_time: the time that processing completed (successfully or not)
        current_process: the currently or last run asyncronous process
        process_status: the status of the currently or last run asyncronous process
    """

    meta: Any = {
        'abstract': True,
    }

    tasks: List[str] = None
    """ the ordered list of tasks that comprise a processing run """

    current_task = StringField(default=None)
    tasks_status = StringField(default=CREATED)
    create_time = DateTimeField(required=True)
    complete_time = DateTimeField()

    errors = ListField(StringField())
    warnings = ListField(StringField())

    current_process = StringField(default=None)
    process_status = StringField(default=None)

    @property
    def tasks_running(self) -> bool:
        """ Returns True of the process has failed or succeeded. """
        return self.tasks_status not in [SUCCESS, FAILURE]

    @property
    def process_running(self) -> bool:
        """ Returns True of an asynchrounous process is currently running. """
        return self.process_status is not None and self.process_status != PROCESS_COMPLETED

    def get_logger(self):
        return utils.get_logger(
            'nomad.processing', current_task=self.current_task, proc=self.__class__.__name__,
            current_process=self.current_process, process_status=self.process_status,
            tasks_status=self.tasks_status)

    @classmethod
    def create(cls, **kwargs):
        """ Factory method that must be used instead of regular constructor. """
        assert cls.tasks is not None and len(cls.tasks) > 0, \
            """ the class attribute tasks must be overwritten with an actual list """
        assert 'tasks_status' not in kwargs, \
            """ do not set the status manually, its managed """

        kwargs.setdefault('create_time', datetime.now())
        self = cls(**kwargs)
        self.tasks_status = PENDING if self.current_task is None else RUNNING
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
        assert self.tasks_running, 'Cannot fail a completed process.'

        failed_with_exception = False

        self.tasks_status = FAILURE

        logger = self.get_logger(**kwargs)
        for error in errors:
            if isinstance(error, Exception):
                failed_with_exception = True
                Proc.log(logger, log_level, 'task failed with exception', exc_info=error)

        self.errors = [str(error) for error in errors]
        self.complete_time = datetime.now()

        if not failed_with_exception:
            errors_str = "; ".join([str(error) for error in errors])
            Proc.log(logger, log_level, 'task failed', errors=errors_str)

        logger.info('process failed')

        self.save()

    def warning(self, *warnings, log_level=logging.WARNING, **kwargs):
        """ Allows to save warnings. Takes strings or exceptions as args. """
        assert self.tasks_running

        logger = self.get_logger(**kwargs)

        for warning in warnings:
            warning = str(warning)
            self.warnings.append(warning)
            Proc.log(logger, log_level, 'task with warning', warning=warning)

    def _continue_with(self, task):
        tasks = self.__class__.tasks
        assert task in tasks, 'task %s must be one of the classes tasks %s' % (task, str(tasks))  # pylint: disable=E1135
        if self.current_task is None:
            assert task == tasks[0], "process has to start with first task"  # pylint: disable=E1136
        else:
            assert tasks.index(task) == tasks.index(self.current_task) + 1, \
                "tasks must be processed in the right order"

        if self.tasks_status == FAILURE:
            return False

        if self.tasks_status == PENDING:
            assert self.current_task is None
            assert task == tasks[0]  # pylint: disable=E1136
            self.tasks_status = RUNNING
            self.current_task = task
            self.get_logger().info('started process')
        else:
            self.current_task = task
            self.get_logger().info('successfully completed task')

        self.save()
        return True

    def _complete(self):
        if self.tasks_status != FAILURE:
            assert self.tasks_status == RUNNING, 'Can only complete a running process, process is %s' % self.tasks_status
            self.tasks_status = SUCCESS
            self.complete_time = datetime.now()
            self.save()
            self.get_logger().info('completed process')

    def block_until_complete(self, interval=0.01):
        """
        Reloads the process constantly until it sees a completed process. Should be
        used with care as it can block indefinitely. Just intended for testing purposes.
        """
        while self.tasks_running:
            time.sleep(interval)
            self.reload()


class InvalidChordUsage(Exception): pass


class Chord(Proc):
    """
    A special Proc base class that manages a chord of child processes. It saves some
    additional state to track child processes and provides methods to control that
    state.

    It uses a counter approach with atomic updates to track the number of processed
    children.

    TODO the joined attribute is not strictly necessary and only serves debugging purposes.
    Maybe it should be removed, since it also requires another save.

    TODO it is vital that sub classes and children don't miss any calls. This might
    not be practical, because in reality processes might even fail to fail.

    TODO in the current upload processing, the join functionality is not strictly necessary.
    Nothing is done after join. We only need it to report the upload completed on API
    request. We could check the join condition on each of thise API queries.

    Attributes:
        total_children (int): the number of spawed children, -1 denotes that number was not
            saved yet
        completed_children (int): the number of completed child procs
        joined (bool): true if all children are completed and the join method was already called
    """
    total_children = IntField(default=-1)
    completed_children = IntField(default=0)
    joined = BooleanField(default=False)

    meta = {
        'abstract': True
    }

    def spwaned_childred(self, total_children=1):
        """
        Subclasses must call this method after all childred have been spawned.

        Arguments:
            total_children (int): the number of spawned children
        """
        self.total_children = total_children
        self.modify(total_children=self.total_children)
        self._check_join(children=0)

    def completed_child(self):
        """ Children must call this, when they completed processing. """
        self._check_join(children=1)

    def _check_join(self, children):
        # incr the counter and get reference values atomically
        completed_children, others = self.incr_counter(
            'completed_children', children, ['total_children', 'joined'])
        total_children, joined = others

        self.get_logger().debug(
            'check for join', total_children=total_children,
            completed_children=completed_children, joined=joined)

        # check the join condition and raise errors if chord is in bad state
        if completed_children == total_children:
            if not joined:
                self.join()
                self.joined = True
                self.modify(joined=self.joined)
                self.get_logger().debug('chord is joined')
            else:
                raise InvalidChordUsage('chord cannot be joined twice.')
        elif completed_children > total_children and total_children != -1:
            raise InvalidChordUsage('chord counter is out of limits.')

    def join(self):
        """ Subclasses might overwrite to do something after all children have completed. """
        pass

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


def task(func):
    """
    The decorator for tasks that will be wrapped in exception handling that will fail the process.
    The task methods of a :class:`Proc` class/document comprise a sequence
    (order of methods in class namespace) of tasks. Tasks must be executed in that order.
    Completion of the last task, will put the :class:`Proc` instance into the
    SUCCESS state. Calling the first task will put it into RUNNING state. Tasks will
    only be executed, if the process has not yet reached FAILURE state.
    """
    def wrapper(self, *args, **kwargs):
        if self.tasks_status == FAILURE:
            return

        self._continue_with(func.__name__)
        try:
            func(self, *args, **kwargs)
        except Exception as e:
            self.fail(e)

        if self.__class__.tasks[-1] == self.current_task and self.tasks_running:
            self._complete()

    setattr(wrapper, '__task_name', func.__name__)
    wrapper.__name__ = func.__name__
    return wrapper


def all_subclasses(cls):
    """ Helper method to calculate set of all subclasses of a given class. """
    return set(cls.__subclasses__()).union(
        [s for c in cls.__subclasses__() for s in all_subclasses(c)])


all_proc_cls = {cls.__name__: cls for cls in all_subclasses(Proc)}
""" Name dictionary for all Proc classes. """


@app.task(bind=True, ignore_results=True, max_retries=3)
def proc_task(task, cls_name, self_id, func_attr):
    """
    The celery task that is used to execute async process functions.
    It ignores results, since all results are handled via the self document.
    It retries for 3 times with a countdown of 3 on missing 'selfs', since this
    might happen in sharded, distributed mongo setups where the object might not
    have yet been propagated and therefore appear missing.
    """
    logger = utils.get_logger(__name__, cls=cls_name, id=self_id, func=func_attr)

    # get the process class
    logger.debug('received process function call')
    global all_proc_cls
    cls = all_proc_cls.get(cls_name, None)
    if cls is None:
        # refind all Proc classes, since more modules might have been imported by now
        all_proc_cls = {cls.__name__: cls for cls in all_subclasses(Proc)}
        cls = all_proc_cls.get(cls_name, None)

    if cls is None:
        logger.critical('document not a subcass of Proc')
        raise ProcNotRegistered('document %s not a subclass of Proc' % cls_name)

    # get the process instance
    try:
        try:
            self = cls.get(self_id)
        except KeyError as e:
            logger.warning('called object is missing')
            raise task.retry(exc=e, countdown=3)
    except KeyError:
        logger.critical('called object is missing, retries exeeded')

    logger = self.get_logger()

    # get the process function
    func = getattr(self, func_attr, None)
    if func is None:
        logger.error('called function not a function of proc class')
        self.fail('called function %s is not a function of proc class %s' % (func_attr, cls_name))
        self.process_status = PROCESS_COMPLETED
        self.save()
        return

    # unwrap the process decorator
    func = getattr(func, '__process_unwrapped', None)
    if func is None:
        logger.error('called function was not decorated with @process')
        self.fail('called function %s was not decorated with @process' % func_attr)
        self.process_status = PROCESS_COMPLETED
        self.save()
        return

    # call the process function
    deleted = False
    try:
        self.process_status = PROCESS_RUNNING
        deleted = func(self)
    except Exception as e:
        self.fail(e)
    finally:
        if deleted is None or not deleted:
            self.process_status = PROCESS_COMPLETED
            self.save()


def process(func):
    """
    The decorator for process functions that will be called async via celery.
    All calls to the decorated method will result in celery task requests.
    To transfer state, the instance will be saved to the database and loading on
    the celery task worker. Process methods can call other (process) functions/methods on
    other :class:`Proc` instances. Each :class:`Proc` instance can only run one
    asny process at a time.
    """
    def wrapper(self, *args, **kwargs):
        assert len(args) == 0 and len(kwargs) == 0, 'process functions must not have arguments'
        if self.process_running:
            raise ProcessAlreadyRunning

        self.current_process = func.__name__
        self.process_status = PROCESS_CALLED
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
