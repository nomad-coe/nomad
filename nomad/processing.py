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

from typing import List, Any
from celery import Celery, Task, chord, group
from celery.result import AsyncResult, result_from_tuple
from celery.signals import after_setup_task_logger, after_setup_logger
from celery.utils.log import get_task_logger
from celery.canvas import Signature
import logging
import logstash
from datetime import datetime

import nomad.config as config
from nomad.files import Upload, UploadError
from nomad import files, utils, users
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


class CalcProcessing(dict):
    """
    Represents the processing of a singular calculation. It is used as argument and
    results for the task that processes an individual calculation from a mainfile.

    Arguments:
        upload_hash: The hash that identifies the upload in the archive.
        mainfile: The path to the mainfile in the upload.
        parser_name: The name of the parser to use/used.
        tmp_mainfile: The full path to the mainfile in the local fs.

    Attributes:
        pipeline: A list of sub processings for all parser, normalizers, indexing, storage.
        calc_hash: The mainfile hash that identifies the calc in the archive.
        archive_id: The id that identifies the archive via `upload_hash/calc_hash`.
    """
    def __init__(self, upload_hash, mainfile, parser_name, tmp_mainfile):
        self['upload_hash'] = upload_hash
        self['parser_name'] = parser_name
        self['mainfile'] = mainfile
        self['calc_hash'] = utils.hash(mainfile)
        self.tmp_mainfile = tmp_mainfile

    def append(self, task, status, errors=[]):
        if errors is None:
            errors = []
        stage = dict(task=task, status=status, errors=errors)
        self.setdefault('pipeline', []).append(stage)

    @property
    def parser_name(self):
        return self['parser_name']

    @property
    def mainfile(self):
        return self['mainfile']

    @property
    def pipeline(self):
        return self.get('pipeline', [])

    @property
    def upload_hash(self):
        return self.get('upload_hash', None)

    @property
    def calc_hash(self):
        return self.get('calc_hash', None)

    @property
    def archive_id(self):
        return '%s/%s' % (self.upload_hash, self.calc_hash)


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
        calc_processings: Information about identified calcs and their processing.
        upload_hash: The hash of the uploaded file. E.g., used for archive/repo ids.
        status: Aggregated celery status for the whole process.
        task_name: Name of the currently running task.
        task_id: Id of the currently running task.
        cause: *None* or *Exception* that caused failure.
        result_tuple: Serialized form of the celery async_results tree for the processing.
    """
    def __init__(self, upload_id: str) -> None:
        self.upload_id = upload_id

        self.calc_processing_task_ids: List[str] = []
        self.calc_processings: List[CalcProcessing] = None
        self.upload_hash: str = None

        self.status: str = 'PENDING'
        self.task_name: str = None
        self.task_id: str = None
        self.cause: Exception = None
        self.result_tuple: Any = None

        self._main = None

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
    def _async_result(self) -> AsyncResult:
        """
        The celery async_result in its regular usable, but not serializable form.

        We use the tuple form to allow serialization (i.e. storage). Keep in mind
        that the sheer `task_id` is not enough, because it does not contain
        the parent tasks, i.e. result tree.
        See `third comment <https://github.com/celery/celery/issues/1328>`_
        for details.
        """
        return result_from_tuple(self.result_tuple, app=app)

    @property
    def _is_started(self) -> bool:
        """ True, if the task is started. """
        return self.result_tuple is not None

    def _update(self, other: 'UploadProcessing') -> 'UploadProcessing':
        """ Updates all attributes from another instance. Returns itself. """
        self.calc_processing_task_ids = other.calc_processing_task_ids
        self.calc_processings = other.calc_processings
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
                self._update(async_result.result)
                break
            else:
                is_last_task = False
                async_result = async_result.parent

        self.calc_processings = []
        for calc_task_id in self.calc_processing_task_ids:
            calc_task_result = parse.AsyncResult(calc_task_id)
            if calc_task_result.ready() and calc_task_result.status == 'SUCCESS':
                self.calc_processings.append(calc_task_result.result)
            elif calc_task_result.state == 'PROGRESS':
                self.calc_processings.append(calc_task_result.info['processing'])

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
        processing.calc_processings = list()
        for filename in upload.filelist:
            for parser in parsers:
                if parser.is_mainfile(upload, filename):
                    calc_processing = CalcProcessing(
                        processing.upload_hash, filename, parser.name,
                        upload.get_path(filename))

                    processing.calc_processings.append(calc_processing)
    except UploadError as e:
        logger.warning('Could find parse specs in open upload %s: %s' % (processing.upload_id, e))
        return processing.fail(e)

    return processing


@app.task(bind=True)
def close_upload(
        task, calc_processings: List[CalcProcessing], processing: UploadProcessing) \
        -> UploadProcessing:

    if not processing.continue_with(task):
        return processing

    processing.calc_processings = calc_processings

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


def _report_progress(task, **kwargs):
    if not task.request.called_directly:
        task.update_state(state='PROGRESS', meta=kwargs)


@app.task(bind=True)
def distributed_parse(
        task: Task, processing: UploadProcessing, close_upload: Signature) -> UploadProcessing:
    if not processing.continue_with(task):
        chord([])(close_upload.clone(args=(processing,)))
        return processing

    # prepare the group of parallel calc processings
    parses = group(parse.s(calc_processing) for calc_processing in processing.calc_processings)
    # save the calc processing task ids to the overall processing
    processing.calc_processing_task_ids = list(child.task_id for child in parses.freeze().children)
    # initiate the chord that runs calc processings first, and close_upload afterwards
    chord(parses)(close_upload.clone(args=(processing,)))

    return processing


@app.task(bind=True)
def parse(self, processing: CalcProcessing) -> CalcProcessing:
    assert processing.upload_hash is not None

    upload_hash = processing.upload_hash
    parser, mainfile = processing.parser_name, processing.mainfile

    # parsing
    logger.debug('Start %s for %s/%s.' % (parser, upload_hash, mainfile))
    try:
        parser_backend = parser_dict[parser].run(processing.tmp_mainfile)
        processing.append(parser, *parser_backend.status)
        _report_progress(self, processing=processing)
    except Exception as e:
        logger.warning(
            '%s stopped on %s/%s: %s' %
            (parser, upload_hash, mainfile, e), exc_info=e)
        processing.append(parser, 'ParseFailed', [e.__str__()])
        return processing

    # normalization
    for normalizer in normalizers:
        normalizer_name = normalizer.__name__
        logger.debug('Start %s for %s/%s.' % (normalizer, upload_hash, mainfile))
        try:
            normalizer(parser_backend).normalize()
            processing.append(normalizer_name, *parser_backend.status)
            _report_progress(self, processing=processing)
        except Exception as e:
            logger.warning(
                '%s stopped on %s/%s: %s' %
                (normalizer, upload_hash, mainfile, e), exc_info=e)
            processing.append(normalizer_name, 'NormalizeFailed', [e.__str__()])
            return normalizer_name

    # update search
    try:
        search.Calc.add_from_backend(
            parser_backend,
            upload_hash=upload_hash,
            calc_hash=processing.calc_hash,
            mainfile=mainfile,
            upload_time=datetime.now())
        processing.append('Indexer', 'IndexSuccess')
        _report_progress(self, processing=processing)
    except Exception as e:
        logger.error(
            'Could not add %s/%s to search index: %s.' %
            (upload_hash, mainfile, e), exc_info=e)
        processing.append('Indexer', 'IndexFailed', [e.__str__()])

    # calc data persistence
    archive_id = processing.archive_id
    try:
        with files.write_archive_json(archive_id) as out:
            parser_backend.write_json(out, pretty=True)
        processing.append('Storage', 'PersistenceSuccess')
        _report_progress(self, processing=processing)
    except Exception as e:
        logger.error(
            'Could not write archive %s for paring %s with %s.' %
            (archive_id, mainfile, parser), exc_info=e)
        processing.append('Storage', 'PersistenceFailed', [e.__str__()])
        return processing

    logger.debug('Written results of %s for %s to %s.' % (parser, mainfile, archive_id))

    return processing


def handle_uploads(quit=False):
    """
    Listens for new uploads in files and initiates their processing.

    Arguments:
        quit: If true, will only handling one event and stop. Otherwise run forever.
    """
    @files.upload_put_handler
    def handle_upload_put(received_upload_id: str):
        logger = utils.get_logger(__name__, upload_id=received_upload_id)
        logger.debug('Initiate upload processing')
        try:
            upload = users.Upload.objects(id=received_upload_id).first()
            if upload is None:
                logger.error('Upload does not exist')
                raise Exception()

            if upload.upload_time is not None:
                logger.warn('Ignore upload notification, since file is already uploaded')
                raise StopIteration

            with logger.lnr_error('Save upload time'):
                upload.upload_time = datetime.now()
                upload.save()

            with logger.lnr_error('Start processing'):
                proc = UploadProcessing(received_upload_id)
                proc.start()
                upload.processing = proc.result_tuple
                upload.save()
        except Exception:
            pass

        if quit:
            raise StopIteration
        logger.debug('Initiated upload processing')

    logger.debug('Start upload put notification handler.')
    handle_upload_put(received_upload_id='provided by decorator')
