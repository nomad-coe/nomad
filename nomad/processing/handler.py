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

from datetime import datetime
from threading import Thread

from nomad import files, utils

from nomad.processing.tasks import extracting_task, cleanup_task, parse_all_task
from nomad.processing.state import UploadProc
from nomad.utils import get_logger, lnr


def start_processing(upload_id, proc: UploadProc=None) -> UploadProc:
    """ Starts the processing tasks via celery canvas. """

    if proc is not None:
        proc = UploadProc(**proc)
    else:
        proc = UploadProc(upload_id)

    # Keep the results of the last task is the workflow.
    # The last task is started by another task, therefore it
    # is not the end of the main task chain.
    finalize = cleanup_task.s()
    finalize_result = finalize.freeze()

    # start the main chain
    main_chain = extracting_task.s(proc) | parse_all_task.s(finalize)
    main_chain_result = main_chain.delay()

    # Create a singular result tree. This might not be the right way to do it.
    finalize_result.parent = main_chain_result

    # Keep the result as tuple that also includes all parents, i.e. the whole
    # serializable tree
    proc.celery_task_ids = finalize_result.as_tuple()
    proc.status = 'STARTED'

    return proc


def handle_uploads(quit=False):
    """
    Starts a daemon that will listen to files for new uploads. For each new
    upload it will initiate the processing and save the task in the upload user data,
    it will wait for processing to be completed and store the results in the upload
    user data.

    Arguments:
        quit: If true, will only handling one event and stop. Otherwise run forever.
    """

    @files.upload_put_handler
    def handle_upload_put(received_upload_id: str):
        from nomad.data import Upload
        logger = get_logger(__name__, upload_id=received_upload_id)
        logger.debug('Initiate upload processing')
        try:
            with lnr(logger, 'Could not load'):
                try:
                    upload = Upload.get(upload_id=received_upload_id)
                except KeyError as e:
                    logger.error('Upload does not exist')
                    raise e

            if upload.upload_time is not None:
                logger.warn('Ignore upload notification, since file is already uploaded')
                raise StopIteration

            with lnr(logger, 'Save upload time'):
                upload.upload_time = datetime.now()
                upload.save()

            with lnr(logger, 'Start processing'):
                proc = start_processing(received_upload_id, proc=upload.proc)
                assert proc.is_started
                upload.proc = proc
                upload.save()

        except Exception:
            logger.error('Exception while handling upload put notification.', exc_info=e)

        if quit:
            raise StopIteration

    utils.get_logger(__name__).debug('Start upload put notification handler.')
    handle_upload_put(received_upload_id='provided by decorator')


def handle_uploads_thread(quit=True):
    thread = Thread(target=lambda: handle_uploads(quit))
    thread.start()
    return thread
