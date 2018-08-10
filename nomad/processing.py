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
This modules allows to (1) run a celery worker that can perform all processing
task, (2) allows to start, manage, read status of processing runs
(i.e. celery canvas, i.e. workflow of tasks), (3) contains the implementation of said
celery tasks.

We make use of celery. It is a widely popular module for running complex distributed
task workflows on a variety of task and result backends. We are using the popular
rabbitmq, redis combination. Rabbitmq allows for a very scalable distribution of tasks in
clouds and clusters, while redis provides a more centralized, reliable temporary storage
for task stati and results.

.. autoclass:: nomad.processing.ProcessRun
"""

from celery import Celery, chord, group, chain, subtask
from celery.result import result_from_tuple
from celery.signals import after_setup_task_logger, after_setup_logger
from celery.utils.log import get_task_logger
import logging
import logstash
import time
import sys
import json

import nomad.config as config
import nomad.files as files
from nomad.dependencies import parsers, parser_dict

# The legacy nomad code uses a logger called 'nomad'. We do not want that this
# logger becomes a child of this logger due to its module name starting with 'nomad.'
logger = get_task_logger(__name__.replace('nomad', 'nomad-xt'))
logger.setLevel(logging.DEBUG)

if config.logstash.enabled:
    def initialize_logstash(logger=None, loglevel=logging.INFO, **kwargs):
        handler = logstash.TCPLogstashHandler(
            config.logstash.host, config.logstash.tcp_port,
            tags=['celery'], message_type='celery', version=1)
        handler.setLevel(loglevel)
        logger.addHandler(handler)
        return logger

    after_setup_task_logger.connect(initialize_logstash)
    after_setup_logger.connect(initialize_logstash)


broker_url = 'pyamqp://%s:%s@%s//' % (
    config.celery.rabbit_user, config.celery.rabbit_password, config.celery.rabbit_host)
backend_url = 'redis://%s/0' % config.celery.redis_host
app = Celery('nomad.processing', backend=backend_url, broker=broker_url)
app.conf.update(
    accept_content=['pickle'],
    task_serializer='pickle',
    result_serializer='pickle',
)


@app.task(bind=True)
def open_upload(task, state):
    if not state.continue_with(task):
        return state

    try:
        upload = files.upload(state.upload_id)
        upload.open()
    except files.UploadError as e:
        logger.debug('Could not open upload %s: %s' % (state.upload_id, e))
        return state.fail(e)

    try:
        state.parse_specs = list()
        for filename in upload.filelist:
            for parser in parsers:
                if parser.is_mainfile(upload, filename):
                    parse_spec = (parser.name, upload.get_path(filename))
                    state.parse_specs.append(parse_spec)
    except files.UploadError as e:
        logger.warning('Could find parse specs in open upload %s: %s' % (state.upload_id, e))
        return state.fail(e)

    return state


@app.task(bind=True)
def close_upload(task, parse_results, state):
    if not state.continue_with(task):
        return state

    try:
        upload = files.upload(state.upload_id)
    except KeyError as e:
        logger.warning('No upload %s' % state.upload_id)
        return state.fail(e)

    upload.close()

    return state


@app.task(bind=True)
def distributed_parse(task, state, close_upload):
    if not state.continue_with(task):
        chord([])(close_upload.clone(args=(state,)))

    parses = group(parse.s(parse_spec) for parse_spec in state.parse_specs)
    chord(parses)(close_upload.clone(args=(state,)))


@app.task()
def parse(parse_spec):
    parser, mainfile = parse_spec

    logger.debug('Start %s for %s.' % parse_spec)
    try:
        parser_dict[parser].run(mainfile)
    except ValueError as e:
        logger.warning('%s stopped on %s/%s: %s' % (parse_spec + (e,)))
        return e
    except Exception as e:
        logger.warning('%s stopped on %s/%s: %s' % (parse_spec + (e,)), exc_info=e)
        return e

    return True  # TODO some other value?


class ProcessState():
    """
    JSON serializable state of a pending, running, or completed :class:`ProcessRun`.
    Instances are used to pass data from task to task within a process workflow.
    Instances are also used to represent state to clients via :func:`ProcessRun.status`.

    Attributes:
        upload_id: The *upload_id* of the :class:`ProcessRun`.
        parse_specs: A list of already identified parse_specs, or None.
        parse_results: A list of completed (failed or successful) parse results.
        current_task: The name of the current task of the process run.
    """

    def __init__(self, upload_id):
        self.upload_id = upload_id
        self.parse_specs = None
        self.parse_results = None

        self.status = 'PENDING'
        self.task_name = None
        self.task_id = None
        self.cause = None

    def fail(self, e):
        self.cause = e
        self.status = 'FAILURE'
        return self

    def continue_with(self, task):
        assert self.status != 'SUCCESS', 'Cannot continue on completed workflow.'

        if self.status == 'FAILURE':
            return False
        else:
            self.status = 'STARTED'
            self.task_name = task.name
            self.task_id = task.request.id
            return True

    def to_json(self):
        return json.dumps(self, indent=4)


class ProcessRun():
    """
    Represents the processing of an uploaded file. It allows to start and manage a
    processing run, retrieve status information and results.

    It is serializable (JSON, pickle). Iternaly stores
    :class:`~celery.results.AsyncResults` instance in serialized *tuple* form.
    We use the serialized form to allow serialization (i.e. storage). Keep in mind
    that the sheer `task_id` is not enough, because it does not contain
    the parent tasks. See `third comment <https://github.com/celery/celery/issues/1328>`_
    for details.

    Warning:
        You have to call :func:`forget` eventually to free all resources and the celery
        results backend.

        Anyhow, results will be deleted after 1 day, depending on `configuration
        <http://docs.celeryproject.org/en/latest/userguide/configuration.html#result-expires>`_.

    Arguments:
        upload_id: The id of the uploaded file in the object storage,
                   see also :mod:`nomad.files`.
    """
    def __init__(self, upload_id):
        self._start_state = ProcessState(upload_id)
        self.result_tuple = None

    def start(self):
        """ Initiates the processing tasks via celery canvas. """
        assert not self.is_started, 'Cannot start a started or used run.'

        finalize = close_upload.s()
        # Keep the results of the last task is the workflow.
        # The last task is started by another task, therefore it
        # is not the end of the main task chain.
        finalize_result = finalize.freeze()

        main_chain = open_upload.s(self._start_state) | distributed_parse.s(finalize)

        # start the main chain
        main_chain_result = main_chain.delay()

        # Create a singular result tree. This might not be the right way to do it.
        finalize_result.parent = main_chain_result

        # Keep the result as tuple to keep self object pickable
        self.result_tuple = finalize_result.as_tuple()

    @property
    def async_result(self):
        """ The celery async_result in its regular usable, but not serializable form. """
        return result_from_tuple(self.result_tuple)

    @property
    def is_started(self):
        """ True, if the task is started. """
        return self.result_tuple is not None

    def status(self):
        """
        Extract the current state from the various tasks involved in upload processing.

        Returns: JSON-style python object with various task information.
        """

        assert self.is_started, 'Run is not yet started.'

        async_result = self.async_result
        while async_result is not None:
            if async_result.ready():
                async_result.result.status = async_result.status
                return async_result.result
            else:
                async_result = async_result.parent
        return self._start_state

    def forget(self):
        """ Forget the results of a completed run; free all resources in the results backend. """
        assert self.ready(), 'Run is not completed.'

        async_result = self.async_result
        while async_result is not None:
            async_result.forget()
            async_result = async_result.parent

    def ready(self):
        """ Returns: True if the task has been executed. """
        assert self.is_started, 'Run is not yet started.'

        return self.async_result.ready()

    def get(self, *args, **kwargs):
        """
        Blocks until the processing has finished. Forwards args, kwargs to
        *celery.result.get* for timeouts, etc.

        Returns: The task result as :func:`status`.
        """
        assert self.is_started, 'Run is not yet started.'

        self.async_result.get(*args, **kwargs)
        return self.status()
