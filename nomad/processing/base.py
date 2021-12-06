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


class ProcessStatus:
    '''
    Class holding constants related to the possible process statuses.

    Attributes:
        READY: The process is ready to start
        PENDING: The process has been called, but still waiting for a celery worker to start running.
        RUNNING: Currently running the main process function.
        WAITING_FOR_RESULT: Waiting for the result from some other process.
        SUCCESS: The last process completed successfully.
        FAILURE: The last process completed with a fatal failure.
        DELETED: Used to signal that the process results in the deletion of the object.

        STATUSES_PROCESSING: List of statuses where the process is still incomplete (no other
            process can be started).
        STATUSES_NOT_PROCESSING: The opposite of the above - statuses from which a new
            process can be started.
    '''
    READY = 'READY'
    PENDING = 'PENDING'
    RUNNING = 'RUNNING'
    WAITING_FOR_RESULT = 'WAITING_FOR_RESULT'
    SUCCESS = 'SUCCESS'
    FAILURE = 'FAILURE'
    DELETED = 'DELETED'

    STATUSES_PROCESSING = (PENDING, RUNNING, WAITING_FOR_RESULT)
    STATUSES_NOT_PROCESSING = (READY, SUCCESS, FAILURE)
    STATUSES_COMPLETED = (SUCCESS, FAILURE)
    STATUSES_VALID_IN_DB = tuple(list(STATUSES_NOT_PROCESSING) + list(STATUSES_PROCESSING))


class InvalidId(Exception): pass


class ProcNotRegistered(Exception): pass


class ProcessAlreadyRunning(Exception): pass


class ProcObjectDoesNotExist(Exception): pass


class ProcessFailure(Exception):
    '''
    Special exception class which allows the user to control how :func:`Proc.fail` should
    be called when the exception is caught from the process function.
    '''
    def __init__(self, *errors, log_level=logging.ERROR, **kwargs):
        self._errors = errors
        self._log_level = log_level
        self._kwargs = kwargs


class Proc(Document):
    '''
    Base class for objects that are subject to processing and need persistent processing
    state. The processing state is persisted in mongo db. Possible processing statuses are
    defined by :class:`ProcessStatus`.

    To initiate a process, an object subclassing Proc must first be created using
    :func:`create`. Processes are then initiated by calling a *process function* on this
    object, which is a member function marked with the decorator @process. Calling a process
    function sets the process_state to PENDING and a celery task is created, which will be
    picked up by a worker, which sets the state to RUNNING and actually executes the
    process function.

    From process_status RUNNING, the process can transition to either SUCCESS, FAILURE or
    WAITING_FOR_RESULT, or result in the deletion of the process object itself (for example
    the process for deleting an upload).

    WAITING_FOR_RESULT means the process needs to wait for the result from spawned child
    processes. When a child process finishes, it will call `try_to_join` on the parent, and
    if it is the parent's last running child process, we will switch to the parent Proc,
    put its status to RUNNING and execute the `join` method on the parent.

    To send the process to status WAITING_FOR_RESULT, the process function must return
    normally, without errors and exceptions, and return the value `ProcessStatus.WAITING_FOR_RESULT`.

    If the process deletes the object itself, the process function should instead return
    `ProcessStatus.DELETED`. If the process function returns normally, without exception
    and without a return value, the process status will be set to SUCCESS.


    Attributes:
        errors: a list of errors that happened during processing.
            NOTE: This value is managed by the framework, do not tamper with this value.
            To fail a process, an exception should be raised.
        warnings: a list of warnings that happened during processing. Warnings do not
            fail a processing run
        last_status_message: A short, human readable message from the current process, with
            information about what the current process is doing, or information about the
            completion (successful or not) of the last process, if no process is currently
            running.
        complete_time: the time that processing completed (successfully or not).
            NOTE: This value is managed by the framework, do not tamper with this value.
        current_process: the currently or last run asyncronous process
            NOTE: This value is managed by the framework, do not tamper with this value.
        process_status: one of the values defined by :class:`ProcessStatus`.
            NOTE: This value is managed by the framework, do not tamper with this value.
    '''

    meta: Any = {
        'abstract': True,
    }

    complete_time = DateTimeField()

    errors = ListField(StringField())
    warnings = ListField(StringField())
    last_status_message = StringField(default=None)

    current_process = StringField(default=None)
    process_status = StringField(default=None)

    worker_hostname = StringField(default=None)
    celery_task_id = StringField(default=None)

    @property
    def process_running(self) -> bool:
        ''' Returns True of an asynchrounous process is currently running (or waiting to run). '''
        return self.process_status in ProcessStatus.STATUSES_PROCESSING

    @classmethod
    def process_running_mongoengine_query(cls):
        ''' Returns a mongoengine query dict (to be used in objects) to find running processes. '''
        return dict(process_status__in=ProcessStatus.STATUSES_PROCESSING)

    def get_logger(self):
        return utils.get_logger(
            'nomad.processing', proc=self.__class__.__name__,
            process=self.current_process, process_status=self.process_status)

    @classmethod
    def create(cls, **kwargs):
        ''' Factory method that must be used instead of regular constructor. '''
        assert 'process_status' not in kwargs, \
            ''' do not set the status manually, its managed '''

        self = cls(**kwargs)
        self.process_status = ProcessStatus.READY
        self.save()

        return self

    def reset(
            self, worker_hostname: str = None, force: bool = False,
            process_status: str = ProcessStatus.READY):
        ''' Resets the process status. If force is not set, there must be no currently running process. '''
        assert not self.process_running or force

        self.current_process = None
        self.process_status = process_status
        self.errors = []
        self.warnings = []
        self.worker_hostname = worker_hostname

    @classmethod
    def reset_pymongo_update(cls, worker_hostname: str = None):
        ''' Returns a pymongo update dict part to reset calculations. '''
        return dict(
            current_process=None, process_status=ProcessStatus.READY,
            errors=[], warnings=[], worker_hostname=worker_hostname)

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

    def fail(self, *errors, save=True, log_level=logging.ERROR, **kwargs):
        '''
        Used to handle when a process fails. Takes strings or exceptions as args.
        The method logs the error(s), updates `self.errors`, `self.last_status_message`,
        `self.process_status`, and calls :func:`on_fail`, and if `save` == True (default)
        it also saves the object to mongodb. The positional args define the errors. An
        error should normally be an instance of Exception, if not it will be converted to
        a string.

        NOTE, processes should NOT call this method directly, or tamper with self.errors etc.
        Rather, if something goes wrong in a process, it should raise an exception!
        '''
        assert self.process_running, 'Cannot fail a completed process.'

        failed_with_exception = False

        logger = self.get_logger(**kwargs)
        self.errors = []
        for error in errors:
            if isinstance(error, Exception):
                failed_with_exception = True
                self.errors.append('%s: %s' % (error.__class__.__name__, str(error)))
                Proc.log(
                    logger, log_level, 'process failed with exception',
                    exc_info=error, error=str(error))
            else:
                self.errors.append(str(error))

        self.complete_time = datetime.utcnow()

        if not failed_with_exception:
            errors_str = "; ".join([str(error) for error in errors])
            Proc.log(logger, log_level, 'process failed', errors=errors_str)

        try:
            self.on_fail()
        except Exception as e:
            # Oh my, nothing is going our way today
            Proc.log(logger, logging.ERROR, 'on_fail failed', exc_info=e, error=str(e))

        logger.info('process failed')
        if len(self.errors) > 0:
            self.last_status_message = f'Process {self.current_process} failed: {self.errors[-1]}'

        self.process_status = ProcessStatus.FAILURE
        if save:
            self.save()

    def warning(self, *warnings, log_level=logging.WARNING, **kwargs):
        ''' Allows to save warnings. Takes strings or exceptions as args. '''
        assert self.process_running

        logger = self.get_logger(**kwargs)

        for warning in warnings:
            warning = str(warning)
            self.warnings.append(warning)
            Proc.log(logger, log_level, 'task with warning', warning=warning)

    def set_last_status_message(self, last_status_message: str):
        ''' Sets the `last_status_message` and saves. '''
        assert self.process_running
        self.last_status_message = last_status_message
        self.save()

    def on_success(self):
        ''' To be called whenever a process is about to transition to status SUCCESS. '''
        pass

    def on_fail(self):
        ''' To be called whenever a process is about to transition to status FAILURE. '''
        pass

    def on_waiting_for_result(self):
        ''' To be called whenever a process is about to transition to status WAITING_FOR_RESULT. '''
        pass

    def block_until_complete(self, interval=0.01):
        '''
        Reloads the process constantly until it sees a completed process (FAILURE or SUCCESS).
        Should be used with care as it can block indefinitely. Just intended for testing
        purposes.
        '''
        self.reload()
        while self.process_running:
            time.sleep(interval)
            self.reload()

    def block_until_complete_or_waiting_for_result(self, interval=0.01):
        '''
        Reloads the process constantly until the process is either complete or in status WAITING_FOR_RESULT.
        Should be used with care as it can block indefinitely. Just intended for testing
        purposes.
        '''
        while self.process_status in (ProcessStatus.PENDING, ProcessStatus.RUNNING):
            time.sleep(interval)
            self.reload()

    @classmethod
    def process_all(
            cls, func, query: Dict[str, Any], exclude: List[str] = [],
            process_args: List[Any] = [], process_kwargs: Dict[str, Any] = {}):
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
            process_status=ProcessStatus.PENDING)})

        for obj in cls.objects(**query).exclude(*exclude):
            obj._run_process(func, process_args, process_kwargs)

    def _run_process(self, func, process_args, process_kwargs):
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
            args=[cls_name, self_id, func.__name__, process_args, process_kwargs],
            queue=queue, priority=priority)

    def __str__(self):
        return 'proc celery_task_id=%s worker_hostname=%s' % (self.celery_task_id, self.worker_hostname)

    def parent(self) -> 'Proc':
        '''
        When running a process marked with `is_child`, and the process completes (i.e. succeeds or
        fails), this method will be invoked by the framework to determine the object's parent Proc.
        '''
        raise NotImplementedError('`parent` not implemented')

    def child_cls(self) -> 'Proc':
        '''
        When running a process which spawns child processes and transitions to WAITING_FOR_RESULT,
        this method will be invoked to determine the "child" class when it is time to try to
        join.
        '''
        raise NotImplementedError('`child_cls` not implemented')

    def _try_to_join(self) -> bool:
        '''
        Called on the parent Proc object to join (resume) the current process.
        For the join to succeed, the following must be fullfilled:
            *   This Proc must have status `WAITING_FOR_RESULT`
            *   Every instance of the child class (defined by calling :func:`child_cls`)
                linked to this Proc (if any) must be in status FAILURE or SUCCESS.
        If the join succeeds, the `process_status` will be set to `RUNNING` and True
        will be returned. Otherwise, the method just returns False. The method is written so
        that the join should succeed once and only once.
        '''
        self.reload()
        if self.process_status != ProcessStatus.WAITING_FOR_RESULT:
            self.get_logger().debug('trying to join: not waiting for result')
            return False
        children_pending_completion = self.child_cls().objects(**{
            self.id_field: self.id,
            'process_status__nin': (ProcessStatus.SUCCESS, ProcessStatus.FAILURE)}).count()
        self.get_logger().debug('trying to join', children_pending_completion=children_pending_completion)

        if not children_pending_completion:
            # We may easily get here multiple times if multiple children finish at the same time
            # To join, we need to read and update the mongo record as a single atomic operation
            old_record = self._get_collection().find_one_and_update(
                {'_id': self.id, 'process_status': ProcessStatus.WAITING_FOR_RESULT},
                {'$set': {'process_status': ProcessStatus.RUNNING}})
            if old_record and old_record['process_status'] == ProcessStatus.WAITING_FOR_RESULT:
                # We managed to update the process_status from WAITING_FOR_RESULT to RUNNING
                # I.e. we've joined!
                self.reload()
                return True
        return False

    def join(self):
        '''
        Override, if applicable, to define what to do when joined, i.e. when the process is
        resumed after all child processes are done. Should return either None (if successful)
        or `ProcessStatus.WAITING_FOR_RESULT` if it again wants to wait for child processes.
        '''
        raise NotImplementedError('`join` not implemented')


class Queueable:
    '''
    Derive from this mark a Proc class as "queuable", meaing it can have processes marked
    with `is_queueable` = True, that will be queued if the Proc is running some other process.
    '''
    pass


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
                'process failed due to worker lost: %s' % str(exc_info.exception),
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
def proc_task(task, cls_name, self_id, func_attr, process_args, process_kwargs):
    '''
    The celery task that is used to execute async process functions.
    It retries for 3 times with a countdown of 3 on missing 'proc', since this
    might happen in sharded, distributed mongo setups where the object might not
    have yet been propagated and therefore appear missing.
    '''
    proc: Proc = unwarp_task(task, cls_name, self_id)  # The current Proc object
    try_to_join = False
    deleting = False

    logger = proc.get_logger()
    logger.debug('received process function call')

    proc.worker_hostname = worker_hostname
    proc.celery_task_id = task.request.id

    # get the process function
    func = getattr(proc, func_attr, None)
    if func is None:
        logger.error('called function not a function of proc class')
        proc.fail('called function %s is not a function of proc class %s' % (func_attr, cls_name))
        return

    # unwrap the process decorator
    call_func = getattr(func, '__process_unwrapped', None)
    is_child = getattr(func, '__is_child', False)
    is_queueable = getattr(func, '__is_queueable', False)
    if call_func is None:
        logger.error('called function was not decorated with @process')
        proc.fail('called function %s was not decorated with @process' % func_attr)
        return
    if is_queueable and not isinstance(proc, Queueable):
        proc.fail('process marked as queuable but class is not a Queueable')
        return

    # call the process function
    try:
        os.chdir(config.fs.working_directory)
        with utils.timer(logger, 'process executed on worker'):
            # Actually call the process function
            proc.process_status = ProcessStatus.RUNNING
            proc.last_status_message = 'Started process: ' + func_attr
            proc.save()
            rv = call_func(proc, *process_args, **process_kwargs)
            if proc.errors:
                # Should be impossible unless the process has called self.fail or tampered
                # with self.errors somehow, which it should not do. We will treat it essentially
                # as if it had raised an exception
                proc.fail('completed with errors but no exception, should not happen', save=False)
            elif rv is None:
                # All looks good
                proc.on_success()
                proc.process_status = ProcessStatus.SUCCESS
                proc.complete_time = datetime.utcnow()
                if proc.warnings:
                    proc.last_status_message = f'Process {func_attr} completed with warnings'
                else:
                    proc.last_status_message = f'Process {func_attr} completed successfully'
                logger.info('completed process')
            elif rv == ProcessStatus.WAITING_FOR_RESULT:
                # No errors, and the process requests to wait for other processes
                proc.on_waiting_for_result()
                proc.process_status = ProcessStatus.WAITING_FOR_RESULT
                proc.save()
                # The process could be ready to resume immediately, for example if no child
                # processes were spawned or if they all finished before the parent process.
                try_to_join = True
                logger.info('waiting for results')
            elif rv == ProcessStatus.DELETED:
                # The Proc object itself to be deleted from the database
                deleting = True
            else:
                raise ValueError('Invalid return value from process function')
    except SystemExit as e:
        proc.fail(e)
        return
    except SoftTimeLimitExceeded as e:
        logger.error('exceeded the celery task soft time limit')
        proc.fail(e, save=False)
    except ProcessFailure as e:
        # Exception with details about how to call self.fail
        proc.fail(*e._errors, save=False, log_level=e._log_level, **e._kwargs)
    except Exception as e:
        proc.fail(e, save=False)

    # The proc is done running
    if is_child and proc.process_status in ProcessStatus.STATUSES_COMPLETED:
        # Save and "switch" to the parent to try to join.
        try:
            proc.save()
            proc = proc.parent()
            logger = proc.get_logger()
            try_to_join = True
        except Exception as e:
            proc.fail(e)  # "Should not happen"
            return

    while try_to_join:
        try_to_join = False
        try:
            joined = proc._try_to_join()
            if joined:
                logger.info('joined')
                rv = proc.join()
                if rv is None:
                    # Succeeded and process is done
                    proc.on_success()
                    proc.process_status = ProcessStatus.SUCCESS
                    proc.complete_time = datetime.utcnow()
                    if proc.warnings:
                        proc.last_status_message = f'Process {proc.current_process} completed with warnings'
                    else:
                        proc.last_status_message = f'Process {proc.current_process} completed successfully'
                    logger.info('completed process')
                elif rv == ProcessStatus.WAITING_FOR_RESULT:
                    # Waiting for result again
                    proc.on_waiting_for_result()
                    proc.process_status = ProcessStatus.WAITING_FOR_RESULT
                    proc.save()
                    # Need to try to join again because, as before, it could be ready to resume immediately
                    logger.info('waiting for results')
                    try_to_join = True
        except SystemExit as e:
            proc.fail(e)
            return
        except SoftTimeLimitExceeded as e:
            logger.error('exceeded the celery task soft time limit')
            proc.fail(e, save=False)
        except ProcessFailure as e:
            # Exception with details about how to call self.fail
            proc.fail(*e._errors, save=False, log_level=e._log_level, **e._kwargs)
        except Exception as e:
            proc.fail(e, save=False)

    if isinstance(proc, Queueable) and proc.process_status in ProcessStatus.STATUSES_COMPLETED:
        # Check if something has been queued up while we were working on this process.
        # Note that the new, completed status has not yet been saved to mongo
        pass  # TODO: call the next @process
    # Finally, save, unless we are deleting, or status is WAITING_FOR_RESULT (in which case
    # it should already have been saved)
    if not deleting and proc.process_status != ProcessStatus.WAITING_FOR_RESULT:
        proc.save()


def process(is_child: bool = False, is_queueable: bool = False):
    '''
    The decorator for process functions that will be called async via celery.
    All calls to the decorated method will result in celery task requests.
    To transfer state, the instance will be saved to the database and loading on
    the celery task worker. Process methods can call other (process) functions/methods on
    other :class:`Proc` instances. Each :class:`Proc` instance can only run one process
    at a time.

    Arguments:
        is_child:
            If this is a child process, which should try to join with the parent process when done.
        is_queueable:
            If invokations of this process should be queued up if a process is already running on
            the object. Default is False, which means if a process of some kind is already
            running, invoking this process will result in a :class:`ProcessAlreadyRunning`
            exception.
    '''
    def process_decorator(func):
        @functools.wraps(func)
        def wrapper(self, *args, **kwargs):
            if self.process_running:
                if is_queueable:
                    raise NotImplementedError('TODO')
                raise ProcessAlreadyRunning('Another process is already running.')

            self.current_process = func.__name__
            self.process_status = ProcessStatus.PENDING
            self.save()

            self._run_process(func, args, kwargs)

        setattr(wrapper, '__process_unwrapped', func)
        setattr(wrapper, '__is_child', is_child)
        setattr(wrapper, '__is_queueable', is_queueable)
        return wrapper
    return process_decorator
