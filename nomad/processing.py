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


@app.task()
def open_upload(upload_id):
    try:
        upload = files.upload(upload_id)
        upload.open()
        logger.debug('Upload %s opened successfully' % upload_id)

        parse_specs = list()
        for filename in upload.filelist:
            for parser in parsers:
                if parser.is_mainfile(upload, filename):
                    parse_spec = (parser.name, upload.get_path(filename))
                    parse_specs.append(parse_spec)

        return parse_specs
    except files.UploadError as e:
        logger.debug('Could not open upload %s: %s' % (upload_id, e))
        return e
    except Exception as e:
        logger.error('Could not open upload %s: %s' % (upload_id, e), exc_info=e)
        return e


@app.task()
def close_upload(parse_results, upload_id):
    try:
        upload = files.upload(upload_id)
    except KeyError as e:
        logger.warning('No upload %s' % upload_id)
        return e

    upload.close()

    logger.debug('Upload %s closed successfully.' % upload_id)

    return parse_results


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


@app.task()
def dmap(it, task, next):
    """
    Map a given list result to clones of a task and call the next task after all
    clones have completed.

    Args:
        it: An AsyncResult that provides an iterable.
        task: A partial signature that is used to run tasks with it elements as arg.
        next: A partial task signature that is run after all task clones with an
                iterable results of the task clone results as arg.
    """
    tasks = group(task.clone([arg, ]) for arg in it)
    return chord(tasks)(subtask(next))


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
        self.upload_id = upload_id
        self.result_tuple = None

    def start(self):
        """ Initiates the processing tasks via celery canvas. """
        assert not self.is_started, 'Cannot start a started or used run.'

        finalize = close_upload.s(self.upload_id)
        # Keep the results of the last task is the workflow.
        # The last task is started by another task, therefore it
        # is not the end of the main task chain.
        finalize_result = finalize.freeze()

        main_chain = chain(open_upload.s(self.upload_id), dmap.s(parse.s(), finalize))

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

        close_result = async_result
        parses_result = close_result.parent
        open_result = parses_result.parent

        return {
            'open': open_result.state,
            'parses': parses_result.state,
            'close': close_result.state
        }

    def forget(self):
        """ Forget the results of a completed run; free all resources in the results backend. """
        assert self.ready(), 'Run is not completed.'
        self.async_result.forget()

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
