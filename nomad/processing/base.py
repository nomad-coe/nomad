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

from typing import List, Any, Dict
import logging
import time
import os
from celery import Celery, Task
from celery.worker.request import Request
from celery.signals import after_setup_task_logger, after_setup_logger, worker_process_init, \
    celeryd_after_setup, worker_process_shutdown
from celery.utils import worker_direct
from celery.exceptions import SoftTimeLimitExceeded
from billiard.exceptions import WorkerLostError
from mongoengine import Document, StringField, ListField, DateTimeField, ValidationError
from mongoengine.connection import ConnectionFailure
from mongoengine.base.metaclasses import TopLevelDocumentMetaclass
from datetime import datetime
import functools

from nomad import config, utils, infrastructure
import nomad.patch  # pylint: disable=unused-import


if config.logstash.enabled:
    from nomad.utils import structlogging

    def initialize_logstash(logger=None, loglevel=logging.DEBUG, **kwargs):
        structlogging.add_logstash_handler(logger)
        return logger

    after_setup_task_logger.connect(initialize_logstash)
    after_setup_logger.connect(initialize_logstash)


@worker_process_init.connect
def setup(**kwargs):
    # each subprocess is supposed disconnect connect again: https://jira.mongodb.org/browse/PYTHON-2090
    try:
        from mongoengine import disconnect
        disconnect()
    except Exception:
        pass

    infrastructure.setup()
    utils.get_logger(__name__).info(
        'celery configured with acks_late=%s' % str(config.celery.acks_late))


worker_hostname = None


@celeryd_after_setup.connect
def capture_worker_name(sender, instance, **kwargs):
    global worker_hostname
    worker_hostname = sender


@worker_process_shutdown.connect
def on_worker_process_shutdown(*args, **kwargs):
    # We need to make sure not to leave open sessions: https://jira.mongodb.org/browse/PYTHON-2090
    from mongoengine.connection import disconnect
    disconnect()


app = Celery('nomad.processing', broker=config.rabbitmq_url())
app.conf.update(worker_hijack_root_logger=False)
app.conf.update(worker_max_memory_per_child=config.celery.max_memory)
if config.celery.routing == config.CELERY_WORKER_ROUTING:
    app.conf.update(worker_direct=True)

app.conf.task_queue_max_priority = 10

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


class ProcObjectDoesNotExist(Exception): pass


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
    '''
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
    '''

    meta: Any = {
        'abstract': True,
    }

    tasks: List[str] = None
    ''' the ordered list of tasks that comprise a processing run '''

    current_task = StringField(default=None)
    tasks_status = StringField(default=CREATED)
    create_time = DateTimeField(required=True)
    complete_time = DateTimeField()

    errors = ListField(StringField())
    warnings = ListField(StringField())
    last_status_message = StringField(default=None)

    current_process = StringField(default=None)
    process_status = StringField(default=None)

    worker_hostname = StringField(default=None)
    celery_task_id = StringField(default=None)

    @property
    def tasks_running(self) -> bool:
        ''' Returns True of the process has failed or succeeded. '''
        return self.tasks_status not in [SUCCESS, FAILURE]

    @property
    def process_running(self) -> bool:
        ''' Returns True of an asynchrounous process is currently running. '''
        return self.process_status is not None and self.process_status != PROCESS_COMPLETED

    @classmethod
    def process_running_mongoengine_query(cls):
        ''' Returns a mongoengine query dict (to be used in objects) to find running processes. '''
        return dict(process_status__in=[PROCESS_CALLED, PROCESS_RUNNING])

    def get_logger(self):
        return utils.get_logger(
            'nomad.processing', task=self.current_task, proc=self.__class__.__name__,
            process=self.current_process, process_status=self.process_status,
            tasks_status=self.tasks_status)

    @classmethod
    def create(cls, **kwargs):
        ''' Factory method that must be used instead of regular constructor. '''
        assert 'tasks_status' not in kwargs, \
            ''' do not set the status manually, its managed '''

        kwargs.setdefault('create_time', datetime.utcnow())
        self = cls(**kwargs)
        if len(cls.tasks) == 0:
            self.tasks_status = SUCCESS
        else:
            self.tasks_status = PENDING if self.current_task is None else RUNNING
        self.save()

        return self

    def reset(
            self, worker_hostname: str = None, force: bool = False,
            tasks_status: str = PENDING):

        ''' Resets the task chain. Assumes there no current running process. '''
        assert not self.process_running or force

        self.current_task = None
        self.process_status = None
        self.tasks_status = tasks_status
        self.errors = []
        self.warnings = []
        self.worker_hostname = worker_hostname

    @classmethod
    def reset_pymongo_update(cls, worker_hostname: str = None):
        ''' Returns a pymongo update dict part to reset calculations. '''
        return dict(
            current_task=None, process_status=None, tasks_status=PENDING, errors=[], warnings=[],
            worker_hostname=worker_hostname)

    @classmethod
    def get_by_id(cls, id: str, id_field: str):
        try:
            obj = cls.objects(**{id_field: id}).first()
        except ValidationError as e:
            raise InvalidId('%s is not a valid id' % id)
        except ConnectionFailure as e:
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
        ''' Allows to fail the process. Takes strings or exceptions as args. '''
        assert self.process_running or self.tasks_running, 'Cannot fail a completed process.'

        failed_with_exception = False

        self.tasks_status = FAILURE

        logger = self.get_logger(**kwargs)
        self.errors = []
        for error in errors:
            if isinstance(error, Exception):
                failed_with_exception = True
                self.errors.append('%s: %s' % (error.__class__.__name__, str(error)))
                Proc.log(
                    logger, log_level, 'task failed with exception',
                    exc_info=error, error=str(error))
            else:
                self.errors.append(str(error))

        self.complete_time = datetime.utcnow()

        if not failed_with_exception:
            errors_str = "; ".join([str(error) for error in errors])
            Proc.log(logger, log_level, 'task failed', errors=errors_str)

        self.on_fail()

        logger.info('process failed')
        if len(self.errors) > 0:
            self.last_satus_message = 'ERROR: %s' % self.errors[-1]

        self.save()

    def on_fail(self):
        pass

    def warning(self, *warnings, log_level=logging.WARNING, **kwargs):
        ''' Allows to save warnings. Takes strings or exceptions as args. '''
        assert self.process_running or self.tasks_running

        logger = self.get_logger(**kwargs)

        for warning in warnings:
            warning = str(warning)
            self.warnings.append(warning)
            Proc.log(logger, log_level, 'task with warning', warning=warning)

    def _continue_with(self, task):
        tasks = self.__class__.tasks
        assert task in tasks, 'task %s must be one of the classes tasks %s' % (task, str(tasks))  # pylint: disable=E1135
        if self.current_task is None:
            assert task == tasks[0], "process has to start with first task %s" % tasks[0]  # pylint: disable=E1136
        elif tasks.index(task) <= tasks.index(self.current_task):
            # task is repeated, probably the celery task of the process was reschedule
            # due to prior worker failure
            self.current_task = task
            self.get_logger().error('task is re-run')
            self.save()
            return True
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
            self.get_logger().info('task completed successfully')

        self.save()
        return True

    def _complete(self):
        if self.tasks_status != FAILURE:
            assert self.tasks_status == RUNNING, 'Can only complete a running process, process is %s' % self.tasks_status
            self.tasks_status = SUCCESS
            self.complete_time = datetime.utcnow()
            self.on_tasks_complete()
            self.save()
            self.get_logger().info('completed process')

    def on_tasks_complete(self):
        ''' Callback that is called when the list of task are completed '''
        pass

    def on_process_complete(self, process_name):
        ''' Callback that is called when the corrent process completed '''
        pass

    def block_until_complete(self, interval=0.01):
        '''
        Reloads the process constantly until it sees a completed process with finished tasks.
        Should be used with care as it can block indefinitely. Just intended for testing
        purposes.
        '''
        while self.tasks_running or self.process_running:
            time.sleep(interval)
            self.reload()

    def block_until_process_complete(self, interval=0.01):
        '''
        Reloads the process constantly until it sees a completed process. Should be
        used with care as it can block indefinitely. Just intended for testing purposes.
        '''
        while self.process_running:
            time.sleep(interval)
            self.reload()

    @classmethod
    def process_all(cls, func, query: Dict[str, Any], exclude: List[str] = []):
        '''
        Allows to run process functions for all objects on the given query. Calling
        process functions though the func:`process` wrapper might be slow, because
        it causes a save on each call. This function will use a query based update to
        do the same for all objects at once.
        '''

        running_query = dict(cls.process_running_mongoengine_query())
        running_query.update(query)
        if cls.objects(**running_query).first() is not None:
            raise ProcessAlreadyRunning('Tried to call a processing function on an already processing process.')

        cls._get_collection().update_many(query, {'$set': dict(
            current_process=func.__name__,
            process_status=PROCESS_CALLED)})

        for obj in cls.objects(**query).exclude(*exclude):
            obj._run_process(func)

    def _run_process(self, func):
        if hasattr(func, '__process_unwrapped'):
            func = getattr(func, '__process_unwrapped')

        self_id = self.id.__str__()
        cls_name = self.__class__.__name__

        queue = None
        if config.celery.routing == config.CELERY_WORKER_ROUTING and self.worker_hostname is not None:
            queue = worker_direct(self.worker_hostname).name

        priority = config.celery.priorities.get('%s.%s' % (cls_name, func.__name__), 1)

        logger = utils.get_logger(__name__, cls=cls_name, id=self_id, func=func.__name__)
        logger.debug('calling process function', queue=queue, priority=priority)

        return proc_task.apply_async(
            args=[cls_name, self_id, func.__name__],
            queue=queue, priority=priority)

    def __str__(self):
        return 'proc celery_task_id=%s worker_hostname=%s' % (self.celery_task_id, self.worker_hostname)


def task(func):
    '''
    The decorator for tasks that will be wrapped in exception handling that will fail the process.
    The task methods of a :class:`Proc` class/document comprise a sequence
    (order of methods in class namespace) of tasks. Tasks must be executed in that order.
    Completion of the last task, will put the :class:`Proc` instance into the
    SUCCESS state. Calling the first task will put it into RUNNING state. Tasks will
    only be executed, if the process has not yet reached FAILURE state.
    '''
    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        try:
            if self.tasks_status == FAILURE:
                return

            self._continue_with(func.__name__)
            try:
                func(self, *args, **kwargs)

            except Exception as e:
                self.fail(e)

            except SystemExit:
                self.fail('unexpected system exit')

            if self.__class__.tasks[-1] == self.current_task and self.tasks_running:
                self._complete()

        except Exception as e:
            # this is very critical and an indicator that the task fail error handling
            # it self failed
            self.get_logger().critical('task wrapper failed with exception', exc_info=e)

    setattr(wrapper, '__task_name', func.__name__)
    return wrapper


def all_subclasses(cls):
    ''' Helper method to calculate set of all subclasses of a given class. '''
    return set(cls.__subclasses__()).union(
        [s for c in cls.__subclasses__() for s in all_subclasses(c)])


all_proc_cls = {cls.__name__: cls for cls in all_subclasses(Proc)}
''' Name dictionary for all Proc classes. '''


class NomadCeleryRequest(Request):
    '''
    A custom celery request class that allows to catch error in the worker main
    thread, which cannot be caught on the worker threads themselves.
    '''

    def _fail(self, event, **kwargs):
        args = self._payload[0]
        # this might be run in the worker main thread, which does not have a mongo
        # connection by default
        if infrastructure.mongo_client is None:
            infrastructure.setup_mongo()

        proc = unwarp_task(self.task, *args)
        proc.fail(event, **kwargs)
        proc.process_status = PROCESS_COMPLETED
        proc.on_process_complete(None)
        proc.save()

    def on_timeout(self, soft, timeout):
        if not soft:
            self._fail('task timeout occurred', timeout=timeout)

        super().on_timeout(soft, timeout)

    def on_failure(self, exc_info, send_failed_event=True, return_ok=False):
        if isinstance(exc_info.exception, WorkerLostError):
            infrastructure.setup()
            utils.get_logger(__name__).error(
                'detected WorkerLostError', exc_info=exc_info.exception)
            self._fail(
                'task failed due to worker lost: %s' % str(exc_info.exception),
                exc_info=exc_info)

        super().on_failure(
            exc_info,
            send_failed_event=send_failed_event,
            return_ok=return_ok
        )


class NomadCeleryTask(Task):
    Request = NomadCeleryRequest


def unwarp_task(task, cls_name, self_id, *args, **kwargs):
    '''
    Retrieves the proc object that the given task is executed on from the database.
    '''
    logger = utils.get_logger(__name__, cls=cls_name, id=self_id)

    # get the process class
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
            from nomad.app import flask
            if flask.app.config['TESTING']:
                # This only happens in tests, where it is not always avoidable that
                # tasks from old test-cases bleed over.
                raise ProcObjectDoesNotExist()
            logger.warning('called object is missing, retry')
            raise task.retry(exc=e, countdown=3)
    except KeyError:
        logger.critical('called object is missing, retries exeeded', proc_id=self_id)
        raise ProcObjectDoesNotExist()

    return self


@app.task(
    bind=True, base=NomadCeleryTask, ignore_results=True, max_retries=3,
    acks_late=config.celery.acks_late, soft_time_limit=config.celery.timeout,
    time_limit=config.celery.timeout * 2)
def proc_task(task, cls_name, self_id, func_attr):
    '''
    The celery task that is used to execute async process functions.
    It ignores results, since all results are handled via the self document.
    It retries for 3 times with a countdown of 3 on missing 'selfs', since this
    might happen in sharded, distributed mongo setups where the object might not
    have yet been propagated and therefore appear missing.
    '''
    self = unwarp_task(task, cls_name, self_id)

    logger = self.get_logger()
    logger.debug('received process function call')

    self.worker_hostname = worker_hostname
    self.celery_task_id = task.request.id

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
        self.on_process_complete(None)
        self.save()
        return

    # call the process function
    deleted = False
    try:
        self.process_status = PROCESS_RUNNING
        os.chdir(config.fs.working_directory)
        with utils.timer(logger, 'process executed on worker'):
            deleted = func(self)
    except SoftTimeLimitExceeded as e:
        logger.error('exceeded the celery task soft time limit')
        self.fail(e)
    except Exception as e:
        self.fail(e)
    except SystemExit as e:
        self.fail(e)
    finally:
        if deleted is None or not deleted:
            self.process_status = PROCESS_COMPLETED
            self.on_process_complete(func.__name__)
            self.save()


def process(func):
    '''
    The decorator for process functions that will be called async via celery.
    All calls to the decorated method will result in celery task requests.
    To transfer state, the instance will be saved to the database and loading on
    the celery task worker. Process methods can call other (process) functions/methods on
    other :class:`Proc` instances. Each :class:`Proc` instance can only run one
    any process at a time.
    '''
    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        assert len(args) == 0 and len(kwargs) == 0, 'process functions must not have arguments'
        if self.process_running:
            raise ProcessAlreadyRunning('Tried to call a processing function on an already processing process.')

        self.current_process = func.__name__
        self.process_status = PROCESS_CALLED
        self.save()

        self._run_process(func)

    task = getattr(func, '__task_name', None)
    if task is not None:
        setattr(wrapper, '__task_name', task)

    setattr(wrapper, '__process_unwrapped', func)

    return wrapper
