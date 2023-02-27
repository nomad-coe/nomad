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

from typing import Any, Tuple, List, Dict, NamedTuple
import logging
import time
import os
from collections import defaultdict
from celery import Celery, Task
from celery.worker.request import Request
from celery.signals import after_setup_task_logger, after_setup_logger, worker_process_init, \
    celeryd_after_setup, worker_process_shutdown
from celery.utils import worker_direct
from celery.exceptions import SoftTimeLimitExceeded
import billiard
from billiard.exceptions import WorkerLostError
from mongoengine import Document, StringField, ListField, DateTimeField, IntField, ValidationError
from mongoengine.connection import ConnectionFailure
from datetime import datetime
import functools

from nomad import config, utils, infrastructure
from nomad.config.models import CELERY_WORKER_ROUTING
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
if config.celery.routing == CELERY_WORKER_ROUTING:
    app.conf.update(worker_direct=True)

app.conf.task_queue_max_priority = 10
app.conf.worker_redirect_stdouts = config.process.redirect_stdouts
app.conf.worker_redirect_stdouts_level = "INFO"


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


class ProcessSyncFailure(Exception): pass


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


class ProcessFlags(NamedTuple):
    ''' Flags defined by the @process and @process_local decorators. '''
    is_blocking: bool
    clear_queue_on_failure: bool
    is_child: bool
    is_local: bool


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
        queue: A list defining queued up calls, waiting to be run. Each item is a triple
            of [func_name, args, kwargs].
            NOTE: This value is managed by the framework, do not tamper with this value.
        sync_counter: An integeger, incremented every time a "sync" operation is executed,
            to ensure state consistency and atomicity. There are three types of sync operations:
            when scheduling a process, starting a process, and completing a process.
            NOTE: This value is managed by the framework, do not tamper with this value.
    '''

    id_field: str = None
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

    queue = ListField()
    sync_counter = IntField(default=0)

    @property
    def process_running(self) -> bool:
        '''
        Returns True of an asynchrounous process is currently running (or waiting to run).
        NOTE, the return value will reflect the state when the object was last updated from
        mongo, not neccesarily the current state.
        '''
        return self.process_status in ProcessStatus.STATUSES_PROCESSING

    @classmethod
    def process_running_mongoengine_query(cls):
        ''' Returns a mongoengine query dict (to be used in objects) to find running processes. '''
        return dict(process_status__in=ProcessStatus.STATUSES_PROCESSING)

    @property
    def current_process_flags(self) -> ProcessFlags:
        if not self.current_process:
            return None
        return process_flags[self.__class__.__name__][self.current_process]

    @property
    def queue_blocked(self) -> bool:
        '''
        If the queue is blocked (i.e. no new @process can be invoked on this object).
        NOTE, the return value will reflect the state when the object was last updated from
        mongo, not neccesarily the current state.
        '''
        if self.process_status in ProcessStatus.STATUSES_PROCESSING:
            # Check current process
            if self.current_process_flags.is_blocking:
                return True
            # Check queued processes
            for item in self.queue:
                func_name = item[0]
                if process_flags[self.__class__.__name__][func_name].is_blocking:
                    return True
        return False

    def get_logger(self):
        process = billiard.current_process()  # pylint: disable=no-member
        worker_id = getattr(process, '_nomad_id', None)
        if worker_id is None:
            worker_id = utils.create_uuid()
            setattr(process, '_nomad_id', worker_id)

        return utils.get_logger(
            'nomad.processing', proc=self.__class__.__name__,
            process=self.current_process, process_status=self.process_status,
            process_worker_id=worker_id)

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
            self, force: bool = False, worker_hostname: str = None,
            process_status: str = ProcessStatus.READY,
            errors: List[str] = [], clear_queue: bool = True):
        '''
        Resets the process status. If force is not set, there must be no currently running process.
        NOTE, use this with care! This should normally only be used manually, to fix processes
        that are "stuck" in status processing, for example if the worker has died etc.
        '''
        assert not self.process_running or force

        self.current_process = None
        self.process_status = process_status
        self.errors = errors
        self.warnings = []
        self.worker_hostname = worker_hostname
        if clear_queue:
            self.queue = []

    @classmethod
    def reset_pymongo_update(
            cls, worker_hostname: str = None,
            process_status=ProcessStatus.READY,
            errors: List[str] = [], clear_queue: bool = True):
        '''
        Returns a pymongo update dict part to reset a Proc.
        NOTE, use this with care! This should normally only be used manually, to fix processes
        that are "stuck" in status processing, for example if the worker has died etc.
        '''
        rv = dict(
            current_process=None, process_status=process_status,
            errors=errors, warnings=[], worker_hostname=worker_hostname)
        if clear_queue:
            rv['queue'] = []
        return rv

    @classmethod
    def get_by_id(cls, id: str, id_field: str):
        try:
            obj = cls.objects(**{id_field: id}).first()
        except ValidationError:
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

    def fail(self, *errors, log_level=logging.ERROR, complete=True, **kwargs):
        '''
        Used to handle when a process fails. Takes strings or exceptions as args.
        The method logs the error(s), updates `self.errors`, `self.last_status_message`,
        `self.process_status`, and calls :func:`on_fail`, and if `complete` == True (default)
        it also saves the object to mongodb. The positional args define the errors. An
        error should normally be an instance of Exception, if not it will be converted to
        a string.

        NOTE, processes should NOT call this method directly, or tamper with self.errors etc.
        Rather, if something goes wrong in a process, it should raise an exception!
        '''
        assert self.process_running, 'Cannot fail a completed process.'

        failed_with_exception = False

        # Log the error
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

        if not failed_with_exception:
            errors_str = "; ".join([str(error) for error in errors])
            Proc.log(logger, log_level, 'process failed', errors=errors_str)

        self.complete_time = datetime.utcnow()

        try:
            self.on_fail()
        except Exception as e:
            # Oh my, nothing is going our way today
            Proc.log(logger, logging.ERROR, 'on_fail failed', exc_info=e, error=str(e))

        logger.info('process failed')
        if len(self.errors) > 0:
            self.last_status_message = f'Process {self.current_process} failed: {self.errors[-1]}'

        self.process_status = ProcessStatus.FAILURE
        if complete:
            self._sync_complete_process(force_clear_queue_on_failure=True)

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

    def _send_to_worker(self, func_name, *args, **kwargs):
        ''' Invokes a celery task, which will prompt a worker to pick up this Proc object. '''
        self_id = self.id.__str__()
        cls_name = self.__class__.__name__

        queue = None
        if config.celery.routing == CELERY_WORKER_ROUTING and self.worker_hostname is not None:
            queue = worker_direct(self.worker_hostname).name

        priority = config.celery.priorities.get('%s.%s' % (cls_name, func_name), 1)

        logger = utils.get_logger(__name__, cls=cls_name, id=self_id, func=func_name)
        logger.debug('calling process function', queue=queue, priority=priority)

        return proc_task.apply_async(
            args=[cls_name, self_id, func_name, args, kwargs],
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
            *   No instance of the child class (defined by calling :func:`child_cls`)
                linked to this Proc (if any) must be processing.
        If the join succeeds, the `process_status` will be set to `RUNNING` and True
        will be returned. Otherwise, the method just returns False. The method is written so
        that the join should succeed once and only once.
        '''
        self.reload()
        if self.process_status != ProcessStatus.WAITING_FOR_RESULT:
            self.get_logger().debug('trying to join: not waiting for result')
            return False
        children_processing = self.child_cls().objects(**{
            self.id_field: self.id,
            'process_status__in': ProcessStatus.STATUSES_PROCESSING}).count()
        self.get_logger().debug('trying to join', children_processing=children_processing)

        if not children_processing:
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

    def _sync_schedule_process(self, func_name: str, *args, **kwargs) -> bool:
        '''
        Used to schedule a call to the @process function named `func_name` with the provided
        `args` and `kwargs`. The `args` and `kwargs` need to be json serializable. If this
        object was previously not processing, we set the state to PENDING. We return True
        if the task should be sent to a worker, False otherwise (if a process is already
        running, we should not send the task to the worker right away, instead it should
        be sent to a worker when the preceding process finishes).

        The call will fail and raise a :class:`ProcessAlreadyRunning` if a *blocking process*
        is running or has been added to the queue (blocking processes prevent further requests
        to be queued until they have been completed).

        This is the first of three *sync operations*. These should be atomic and occur in
        sequence, each call fully seeing the relevant changes of the previous operation.
        Because of propagation delays and possible race conditions, we use `sync_counter` to
        detect update collisions. Such collisions should be very unusual, but if they occur
        we reload self and try again, up to 3 times (which *should* be enough to guarantee
        success with virtually absolute certainty). Because we may need to reload, an object
        calling a sync operation should have no unsaved changes.
        '''
        try_counter = 0
        while True:
            if self.queue_blocked:
                raise ProcessAlreadyRunning('A blocking process is running or waiting to run')
            prev_process_running = self.process_running
            mongo_update = {
                '$set': {'sync_counter': self.sync_counter + 1}}
            if prev_process_running:
                # Something else is running, add to queue
                if self.queue:
                    mongo_update['$push'] = {'queue': [func_name, args, kwargs]}
                else:
                    mongo_update['$set'].update(queue=[[func_name, args, kwargs]])
            else:
                # Nothing is running
                mongo_update['$set'].update(
                    process_status=ProcessStatus.PENDING,
                    current_process=func_name,
                    last_status_message='Pending: ' + func_name)
            # Try to update self atomically. Will fail if someone else has managed to write
            # a sync op in between.
            old_record = self._get_collection().find_one_and_update(
                {
                    '$and': [
                        {'_id': self.id},
                        {
                            '$or': [
                                {'sync_counter': self.sync_counter},
                                {'sync_counter': {'$exists': False}}
                            ]
                        }
                    ]
                }, mongo_update)
            try_counter += 1
            if old_record and old_record.get('sync_counter', 0) == self.sync_counter:
                # We have successfully scheduled the process!
                self.reload()
                return not prev_process_running
            # Someone else must have written a sync op (ticked up the sync_counter) in between
            if try_counter >= 3:
                # Three failed attempts - should be virtually impossible!
                raise ProcessSyncFailure('Failed to schedule process too many times - should not happen')
            # Otherwise, sleep, reload, and try again
            time.sleep(0.1)
            self.reload()

    def _sync_start_local_process(self, func_name: str):
        '''
        Used to start a *local* process. If successful, the status transitions to RUNNING
        atomically. The call will fail and raise a :class:`ProcessAlreadyRunning` if any
        other process is currently running.

        This is one of three *sync operations*. See :func:`_sync_schedule_process`
        for more info.
        '''
        try_counter = 0
        while True:
            if self.process_running:
                raise ProcessAlreadyRunning('Another process is running or waiting to run')
            mongo_update = {
                '$set': dict(
                    sync_counter=self.sync_counter + 1,
                    process_status=ProcessStatus.RUNNING,
                    current_process=func_name,
                    last_status_message='Started: ' + func_name,
                    worker_hostname=None,
                    celery_task_id=None,
                    errors=[],
                    warnings=[])}
            # Try to update self atomically. Will fail if someone else has managed to write
            # a sync op in between.
            old_record = self._get_collection().find_one_and_update(
                {
                    '$and': [
                        {'_id': self.id},
                        {
                            '$or': [
                                {'sync_counter': self.sync_counter},
                                {'sync_counter': {'$exists': False}}
                            ]
                        }
                    ]
                }, mongo_update)
            try_counter += 1
            if old_record and old_record.get('sync_counter') == self.sync_counter:
                # We have successfully started the process!
                self.reload()
                return
            # Someone else must have written a sync op (ticked up the sync_counter) in between
            if try_counter >= 3:
                # Three failed attempts - should be virtually impossible!
                raise ProcessSyncFailure('Failed to start local process too many times - should not happen')
            # Otherwise, sleep, reload, and try again
            time.sleep(0.1)
            self.reload()

    def _sync_complete_process(self, force_clear_queue_on_failure=False) -> Tuple[str, List[Any], Dict[str, Any]]:
        '''
        Used to complete a process (when done, successful or not). Returns a triple
        containing information about the next process to run (if any), of the
        form (func_name, args, kwargs).

        There are 3 possibilities:
            1)  There is something in the queue, and the current process was successful
                -> We set status to PENDING and return the next process.
            2)  There is something in the queue, and the current process FAILED:
                -> Behaviour depends on the process decorator flag `clear_queue_on_failure`
                   and the parameter `force_clear_queue_on_failure`:
                        If either is True: We clear the queue, set status to FAILURE and return None.
                        Otherwise: We set status to PENDING and return the next process.
            3)  There is nothing in the queue:
                -> We set the status to the provided value and return None

        This is one of three *sync operations*. See :func:`_sync_schedule_process`
        for more info.
        '''
        assert self.process_status in ProcessStatus.STATUSES_COMPLETED
        # As a safety precausion, save all updates made to the object except the status
        # (We want it to have status RUNNING until the atomic read/write finishes)
        process_status = self.process_status
        self.process_status = ProcessStatus.RUNNING
        self.save()
        self.process_status = process_status
        clear_queue_on_failure = (
            force_clear_queue_on_failure or self.current_process_flags.clear_queue_on_failure)
        try_counter = 0
        while True:
            next_process = None
            mongo_update = {'$set': {'sync_counter': self.sync_counter + 1}}
            if self.queue:
                # Something in the queue
                if not clear_queue_on_failure or process_status == ProcessStatus.SUCCESS:
                    # Move on to the next process
                    next_process = self.queue[0]
                    next_func_name = next_process[0]
                    mongo_update['$pop'] = {'queue': -1}  # pops the first element
                    mongo_update['$set'].update(
                        process_status=ProcessStatus.PENDING,
                        last_status_message='Pending: ' + next_func_name,
                        current_process=next_func_name)
                else:
                    # Failed and clear_queue_on_failure is set to True - clear the queue
                    mongo_update['$set'].update(process_status=process_status, queue=[])
            else:
                mongo_update['$set'].update(process_status=process_status)
            # Try to update self atomically. Will fail if someone else has managed to write
            # a sync op in between.
            old_record = self._get_collection().find_one_and_update(
                {'_id': self.id, 'sync_counter': self.sync_counter}, mongo_update)
            try_counter += 1
            if old_record and old_record.get('sync_counter') == self.sync_counter:
                # We have successfully completed the process
                return next_process
            # Someone else must have written a sync op (ticked up the sync_counter) in between
            if try_counter >= 3:
                # Three failed attempts - should be virtually impossible!
                raise ProcessSyncFailure('Failed to complete process too many times - should not happen')
            # Make another attempt
            time.sleep(0.1)
            self.reload()


def all_subclasses(cls):
    ''' Helper method to calculate set of all subclasses of a given class. '''
    return set(cls.__subclasses__()).union(
        [s for c in cls.__subclasses__() for s in all_subclasses(c)])


all_proc_cls = {cls.__name__: cls for cls in all_subclasses(Proc)}
''' Name dictionary for all Proc classes. '''

process_flags: Dict[str, Dict[str, ProcessFlags]] = defaultdict(dict)
''' { <Proc class name>: { <process func name>: ProcessFlags } } '''


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
        logger.critical('document not a subclass of Proc')
        raise ProcNotRegistered('document %s not a subclass of Proc' % cls_name)

    # get the process instance
    try:
        try:
            self = cls.get(self_id)
        except KeyError as e:
            logger.warning('called object is missing, retry')
            raise task.retry(exc=e, countdown=3)
    except KeyError:
        logger.critical('called object is missing, retries exceeded', proc_id=self_id)
        raise ProcObjectDoesNotExist()

    return self


@app.task(
    bind=True, base=NomadCeleryTask, ignore_results=True, max_retries=3,
    acks_late=config.celery.acks_late, soft_time_limit=config.celery.timeout,
    time_limit=config.celery.timeout * 2)
def proc_task(task, cls_name, self_id, func_name, args, kwargs):
    '''
    The celery task that is used to execute async process functions.
    It retries for 3 times with a countdown of 3 in case of propagation problems, since this
    might happen in sharded, distributed mongo setups where the updates might not
    have yet propagated to everyone.
    '''
    # Obtain the Proc object. Raises exception to make celery retry if object has not propagated.
    proc: Proc = unwarp_task(task, cls_name, self_id)
    logger = proc.get_logger()
    logger.debug('Executing celery task')

    if '_meta_label' in kwargs:
        config.meta.label = kwargs['_meta_label']
        del(kwargs['_meta_label'])

    try_to_join = False
    deleting = False

    # get the process function
    func = getattr(proc, func_name, None)
    if func is None:  # "Should not happen"
        logger.error('called function not a function of proc class')
        proc.fail('called function %s is not a function of proc class %s' % (func_name, cls_name))
        return

    # unwrap the process decorator
    unwrapped_func = getattr(func, '__process_unwrapped', None)
    if unwrapped_func is None:  # "Should not happen"
        logger.error('called function was not decorated with @process')
        proc.fail('called function %s was not decorated with @process' % func_name)
        return

    # call the process function
    try:
        is_child = process_flags[cls_name][func_name].is_child
        os.chdir(config.fs.working_directory)
        with utils.timer(logger, 'process executed on worker', log_memory=True):
            # Set state to RUNNING
            proc.process_status = ProcessStatus.RUNNING
            proc.last_status_message = 'Started: ' + func_name
            proc.worker_hostname = worker_hostname
            proc.celery_task_id = task.request.id
            proc.errors = []
            proc.warnings = []
            proc.save()
            # Actually call the process function
            rv = unwrapped_func(proc, *args, **kwargs)
            if proc.errors:
                # Should be impossible unless the process has tampered with self.errors, which
                # it should not do. We will treat it essentially as if it had raised an exception
                proc.fail('completed with errors but no exception, should not happen', complete=False)
            elif rv is None:
                # All looks good
                proc.on_success()
                proc.process_status = ProcessStatus.SUCCESS
                proc.complete_time = datetime.utcnow()
                if proc.warnings:
                    proc.last_status_message = f'Process {func_name} completed with warnings'
                else:
                    proc.last_status_message = f'Process {func_name} completed successfully'
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
        proc.fail(e, complete=False)
    except ProcessFailure as e:
        # Exception with details about how to call self.fail
        proc.fail(*e._errors, log_level=e._log_level, complete=False, **e._kwargs)
    except Exception as e:
        proc.fail(e, complete=False)

    # The proc is done running
    if is_child and proc.process_status in ProcessStatus.STATUSES_COMPLETED:
        try:
            next_process = proc._sync_complete_process()
            if next_process:
                # More jobs in the queue
                func_name, args, kwargs = next_process
                proc._send_to_worker(func_name, *args, **kwargs)
                return
            # Processing finished (successful or not)
            # Switch to the parent to try to join.
            proc = proc.parent()
            logger = proc.get_logger()
            try_to_join = True
        except Exception as e:  # "Should not happen"
            proc.fail(e)
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
            proc.fail(e, complete=False)
        except ProcessFailure as e:
            # Exception with details about how to call self.fail
            proc.fail(*e._errors, complete=False, log_level=e._log_level, **e._kwargs)
        except Exception as e:
            proc.fail(e, complete=False)

    if not deleting and proc.process_status in ProcessStatus.STATUSES_COMPLETED:
        # We are about to transition from RUNNING to completed (FAILURE or SUCCESS)
        # But, if something is queued up we should actually go to PENDING instead, and
        # trigger celery again
        try:
            next_process = proc._sync_complete_process()
            if next_process:
                func_name, args, kwargs = next_process
                proc._send_to_worker(func_name, *args, **kwargs)
        except Exception as e:
            proc.fail(e)


def process(is_blocking: bool = False, clear_queue_on_failure: bool = True, is_child: bool = False):
    '''
    The decorator for process functions that will be queued up and executed async via celery.
    To transfer state, the instance will be saved to the database and loading on
    the celery task worker. Process methods can call other (process) functions/methods on
    other :class:`Proc` instances. Each :class:`Proc` instance can only run one process
    at a time. The Proc object should not have any unsaved changes when a process is invoked.

    Arguments:
        is_blocking:
            If True, this is a *blocking process*. After a blocking process has been scheduled,
            no other processes can be scheduled, until this process has finished. Attempts to
            invoke another process will in this case result in an exception.
        clear_queue_on_failure:
            If True and the process fails, we'll clear the queue, thus skipping any pending jobs.
        is_child:
            If this is a child process, which should try to join with the parent process when done
            (= when the queue is empty).
    '''
    def process_decorator(func):
        # Determine canonical class name
        cls_name, func_name = func.__qualname__.split('.')

        process_flags[cls_name][func_name] = ProcessFlags(
            is_blocking,
            clear_queue_on_failure,
            is_child,
            is_local=False)

        @functools.wraps(func)
        def wrapper(self: Proc, *args, **kwargs):
            kwargs['_meta_label'] = config.meta.label
            send_to_worker = self._sync_schedule_process(func_name, *args, **kwargs)
            if send_to_worker:
                try:
                    self._send_to_worker(func_name, *args, **kwargs)
                except Exception as e:
                    self.fail(e)
                    raise

        setattr(wrapper, '__process_unwrapped', func)
        return wrapper
    return process_decorator


def process_local(func):
    '''
    The decorator for functions that process locally. These work similarly to functions
    marked with the `@process` decorator, but they are executed directly, in the current
    thread, not via celery. Consequently, they can only be started if no other process is
    running. They are also implicitly blocking, i.e. while running, no other process
    (local or celery-based) can be started or scheduled on the same object. It can also not
    spawn child processes and wait for them using the WAITING_FOR_RESULT mechanism, or itself
    be a child process (as this means joining with a parent process when done).

    If successful, a local process can return a value to the caller (unlike celery processes).
    If unsuccessful, an exception will be raised (note that the usual process handling is
    always applied, i.e. we set self.process_status, self.errors etc. accordingly). The Proc
    object should not have any unsaved changes when a local process is invoked.
    '''
    # Determine canonical class name
    cls_name, func_name = func.__qualname__.split('.')

    process_flags[cls_name][func_name] = ProcessFlags(
        is_blocking=True,
        clear_queue_on_failure=False,  # Not relevant, since local processes are always blocking
        is_child=False,
        is_local=True)

    def wrapper(self: Proc, *args, **kwargs):
        logger = self.get_logger()
        logger.debug('Executing local process')
        self._sync_start_local_process(func_name)

        try:
            os.chdir(config.fs.working_directory)
            with utils.timer(logger, 'process executed locally', log_memory=True):
                # Actually call the process function
                rv = func(self, *args, **kwargs)
                if self.errors:
                    # Should be impossible unless the process has tampered with self.errors, which
                    # it should not do. We will treat it essentially as if it had raised an exception
                    raise RuntimeError('completed with errors but no exception, should not happen')
                # All looks good
                self.on_success()
                self.process_status = ProcessStatus.SUCCESS
                self.complete_time = datetime.utcnow()
                if self.warnings:
                    self.last_status_message = f'Process {func_name} completed with warnings'
                else:
                    self.last_status_message = f'Process {func_name} completed successfully'
                logger.info('completed process')
                return rv
        except SystemExit as e:
            self.fail(e, complete=False)
            raise
        except ProcessFailure as e:
            # Exception with details about how to call self.fail
            self.fail(*e._errors, log_level=e._log_level, complete=False, **e._kwargs)
            raise
        except Exception as e:
            self.fail(e, complete=False)
            raise
        finally:
            self._sync_complete_process()  # Queue should be empty, so nothing more to do

    return wrapper
