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
This modules allows to (1) run a celery worker that can perform all upload processing
tasks, (2) allows to start, manage, get status of upload processings
(i.e. celery canvas, i.e. workflow of tasks), (3) contains the implementation of said
celery tasks.

We make use of celery. It is a widely popular module for running complex distributed
task workflows on a variety of task and result backends. We are using the popular
rabbitmq, redis combination. Rabbitmq allows for a very scalable distribution of tasks in
clouds and clusters, while redis provides a more centralized, reliable temporary storage
for task stati and results.

The class :class:`UploadProcessing` allows to start, manage, and read the status of
a upload processing in a serializable form.

.. autoclass:: nomad.processing.UploadProcessing
"""

from typing import List, Any, Tuple, Dict
from celery import Celery, Task, chord, group
from celery.result import ResultBase, result_from_tuple
from celery.signals import after_setup_task_logger, after_setup_logger
from celery.utils.log import get_task_logger
from celery.canvas import Signature
import logging
import logstash
from datetime import datetime

import nomad.config as config
from nomad.files import Upload, UploadError
from nomad import files, utils
from nomad.parsing import parsers, parser_dict
from nomad.normalizing import normalizers
from nomad import search
import nomad.patch  # pylint: disable=unused-import

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


app = Celery('nomad.processing', backend=config.celery.backend_url, broker=config.celery.broker_url)
app.add_defaults(dict(
    accept_content=['json', 'pickle'],
    task_serializer=config.celery.serializer,
    result_serializer=config.celery.serializer,
))

ProcessingTaskResult = List[Tuple[str, Tuple[str, List[str]]]]
""" A list of parser/normalizer (tool, (status, errors)) tuples. """

ParseSpec = Tuple[str, str]
""" A tuple of (parser name, mainfile). """


class UploadProcessing():
    """
    Represents the processing of an uploaded file. It allows to start and manage a
    processing run, acts itself as an client accesible and serializable state of the
    processing.

    It is serializable (JSON, pickle). Iternaly stores
    :class:`~celery.results.AsyncResults` instance in serialized *tuple* form to
    keep connected to the results backend.

    Instances of this class represent the state of a processing and are handed
    from task to task through the processing workflow.

    Warning:
        You have to call :func:`forget` eventually to free all resources and the celery
        results backend.

        Anyhow, results will be deleted after 1 day, depending on `configuration
        <http://docs.celeryproject.org/en/latest/userguide/configuration.html#result-expires>`_.

    Arguments:
        upload_id: The id of the uploaded file in the object storage,
                   see also :mod:`nomad.files`.

    Attributes:
        parse_specs: List of (parser_name, mainfile) tuples.
        parse_results: Result of the parsers, currently bool indicating success.
        upload_hash: The hash of the uploaded file. E.g., used for archive/repo ids.
        status: Aggregated celery status for the whole process.
        task_name: Name of the currently running task.
        task_id: Id of the currently running task.
        cause: *None* or *Exception* that caused failure.
        result_tuple: Serialized form of the celery async_results tree for the processing.
    """
    def __init__(self, upload_id: str) -> None:
        self.upload_id = upload_id

        self.parse_specs: List[ParseSpec] = None
        self.processing_results: List[ProcessingTaskResult] = None
        self.upload_hash: str = None

        self.status: str = 'PENDING'
        self.task_name: str = None
        self.task_id: str = None
        self.cause: Exception = None
        self.result_tuple: Any = None

    @staticmethod
    def from_result_backend(upload_id, result_tuple):
        """ Loads the processing data from the results backend and returnes an updated instance. """
        processing = UploadProcessing(upload_id)
        processing.result_tuple = result_tuple
        return processing.updated()

    def start(self):
        """ Initiates the processing tasks via celery canvas. """
        assert not self._is_started, 'Cannot start a started or used processing.'

        # Keep the results of the last task is the workflow.
        # The last task is started by another task, therefore it
        # is not the end of the main task chain.
        finalize = close_upload.s()
        finalize_result = finalize.freeze()

        # start the main chain
        main_chain = open_upload.s(self) | distributed_parse.s(finalize)
        main_chain_result = main_chain.delay()

        # Create a singular result tree. This might not be the right way to do it.
        finalize_result.parent = main_chain_result

        # Keep the result as tuple that also includes all parents, i.e. the whole
        # serializable tree
        self.result_tuple = finalize_result.as_tuple()

    @property
    def _async_result(self) -> ResultBase:
        """
        The celery async_result in its regular usable, but not serializable form.

        We use the tuple form to allow serialization (i.e. storage). Keep in mind
        that the sheer `task_id` is not enough, because it does not contain
        the parent tasks, i.e. result tree.
        See `third comment <https://github.com/celery/celery/issues/1328>`_
        for details.
        """
        return result_from_tuple(self.result_tuple)

    @property
    def _is_started(self) -> bool:
        """ True, if the task is started. """
        return self.result_tuple is not None

    def _update(self, other: 'UploadProcessing') -> 'UploadProcessing':
        """ Updates all attributes from another instance. Returns itself. """
        self.parse_specs = other.parse_specs
        self.processing_results = other.processing_results
        self.upload_hash = other.upload_hash

        self.status = other.status
        self.task_name = other.task_name
        self.task_id = other.task_id
        self.cause = other.cause

        return self

    def updated(self) -> 'UploadProcessing':
        """ Consults the result backend and updates itself with the available results. """
        assert self._is_started, 'Run is not yet started.'

        async_result = self._async_result
        is_last_task = True
        while async_result is not None:
            if async_result.ready():
                status = async_result.status
                if status == 'SUCCESS' and not is_last_task:
                    status = 'PROGRESS'
                async_result.result.status = status
                return self._update(async_result.result)
            else:
                is_last_task = False
                async_result = async_result.parent
        return self

    def forget(self) -> None:
        """ Forget the results of a completed run; free all resources in the results backend. """
        assert self.ready(), 'Run is not completed.'

        async_result = self._async_result
        while async_result is not None:
            async_result.forget()
            async_result = async_result.parent

    def ready(self) -> bool:
        """ Returns: True if the task has been executed. """
        assert self._is_started, 'Run is not yet started.'

        return self._async_result.ready()

    def get(self, *args, **kwargs) -> 'UploadProcessing':
        """
        Blocks until the processing has finished. Forwards args, kwargs to
        *celery.result.get* for timeouts, etc.

        Returns: An upadted instance of itself with all the results.
        """
        assert self._is_started, 'Run is not yet started.'

        self._async_result.get(*args, **kwargs)
        return self.updated()

    def fail(self, e: Exception) -> 'UploadProcessing':
        """ Allows tasks to mark this processing as failed. All following task will do nothing. """
        self.cause = e
        self.status = 'FAILURE'
        return self

    def continue_with(self, task: Task) -> bool:
        """ Upadtes itself with information about the new current task. """
        assert self.status != 'SUCCESS', 'Cannot continue on completed workflow.'

        if self.status == 'FAILURE':
            return False
        else:
            self.status = 'STARTED'
            self.task_name = task.name
            self.task_id = task.request.id
            return True

    def _calc_processing(self, parse_spec: ParseSpec, results: ProcessingTaskResult):
        proc: Dict[str, Any] = {
            'parser': parse_spec[0],
            'mainfile': parse_spec[1]
        }

        if results is not None:
            pipeline = list()
            for result in results:
                stage: Dict[str, Any] = {
                    'stage': result[0],
                    'status': result[1][0]
                }
                if result[1][1] is not None:
                    stage['errors'] = result[1][1]
                pipeline.append(stage)

            proc['pipeline'] = pipeline

        return proc

    @property
    def calc_processings(self):
        if self.parse_specs is None:
            return None

        results = self.processing_results
        if results is None:
            results = [None for _ in self.parse_specs]
        return list(
            self._calc_processing(spec, result)
            for spec, result in zip(self.parse_specs, results))


@app.task(bind=True)
def open_upload(task: Task, processing: UploadProcessing) -> UploadProcessing:
    if not processing.continue_with(task):
        return processing

    try:
        upload = Upload(processing.upload_id)
        upload.open()
    except KeyError as e:
        logger.debug('Process request for non existing upload %s.' % processing.upload_id)
        return processing.fail(e)
    except UploadError as e:
        logger.debug('Could not open upload %s: %s' % (processing.upload_id, e))
        return processing.fail(e)

    logger.debug('Opened upload %s' % processing.upload_id)

    try:
        processing.upload_hash = upload.hash()
    except UploadError as e:
        logger.error('Could not create an upload hash %s: %s' % (processing.upload_id, e))
        return processing.fail(e)

    try:
        # TODO: deal with multiple possible parser specs
        processing.parse_specs = list()
        for filename in upload.filelist:
            for parser in parsers:
                if parser.is_mainfile(upload, filename):
                    parse_spec = (parser.name, upload.get_path(filename))
                    processing.parse_specs.append(parse_spec)
    except UploadError as e:
        logger.warning('Could find parse specs in open upload %s: %s' % (processing.upload_id, e))
        return processing.fail(e)

    return processing


@app.task(bind=True)
def close_upload(
        task, processing_results: List[ProcessingTaskResult], processing: UploadProcessing) \
        -> UploadProcessing:

    if not processing.continue_with(task):
        return processing

    processing.processing_results = processing_results

    try:
        upload = Upload(processing.upload_id)
    except KeyError as e:
        logger.warning('No upload %s' % processing.upload_id)
        return processing.fail(e)

    try:
        upload.close()
    except Exception as e:
        logger.error('Could not close upload %s: %s' % (processing.upload_id, e))
        return processing.fail(e)

    logger.debug('Closed upload %s' % processing.upload_id)

    return processing


@app.task(bind=True)
def distributed_parse(
        task: Task, processing: UploadProcessing, close_upload: Signature) -> UploadProcessing:

    if not processing.continue_with(task):
        chord([])(close_upload.clone(args=(processing,)))
        return processing

    parses = group(parse.s(processing, parse_spec) for parse_spec in processing.parse_specs)
    chord(parses)(close_upload.clone(args=(processing,)))
    return processing


@app.task()
def parse(processing: UploadProcessing, parse_spec: ParseSpec) -> ProcessingTaskResult:
    assert processing.upload_hash is not None

    upload_hash = processing.upload_hash
    parser, mainfile = parse_spec
    calc_hash = utils.hash(mainfile)
    results: ProcessingTaskResult = list()

    # parsing
    logger.debug('Start %s for %s/%s.' % (parser, upload_hash, mainfile))
    try:
        parser_backend = parser_dict[parser].run(mainfile)
        results.append((parser, parser_backend.status))
    except Exception as e:
        logger.warning(
            '%s stopped on %s/%s: %s' %
            (parser, upload_hash, mainfile, e), exc_info=e)
        results.append((parser, ('ParseFailed', [e.__str__()])))
        return results

    # normalization
    for normalizer in normalizers:
        logger.debug('Start %s for %s/%s.' % (normalizer, upload_hash, mainfile))
        try:
            normalizer(parser_backend).normalize()
            results.append((normalizer.__name__, parser_backend.status))
        except Exception as e:
            logger.warning(
                '%s stopped on %s/%s: %s' %
                (normalizer, upload_hash, mainfile, e), exc_info=e)
            results.append((normalizer, ('NormalizeFailed', [e.__str__()])))
            return results

    # update search
    try:
        search.Calc.add_from_backend(
            parser_backend,
            upload_hash=upload_hash,
            calc_hash=calc_hash,
            mainfile=mainfile,
            upload_time=datetime.now())
        results.append(('Indexer', ('IndexSuccess', [])))
    except Exception as e:
        logger.error(
            'Could not add %s/%s to search index: %s.' %
            (upload_hash, mainfile, e), exc_info=e)
        results.append(('Indexer', ('IndexFailed', [e.__str__()])))

    # calc data persistence
    archive_id = '%s/%s' % (upload_hash, calc_hash)
    try:
        with files.write_archive_json(archive_id) as out:
            parser_backend.write_json(out, pretty=True)
        results.append(('Storage', ('PersistenceSuccess', [])))
    except Exception as e:
        logger.error(
            'Could not write archive %s for paring %s with %s.' %
            (archive_id, mainfile, parser), exc_info=e)
        results.append(('Storage', ('PersistenceFailed', [e.__str__()])))

    logger.debug('Written results of %s for %s to %s.' % (parser, mainfile, archive_id))

    return results


@app.task()
def mul(x, y):
    return x * y
