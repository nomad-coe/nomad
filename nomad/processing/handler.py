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
from threading import Thread, Event

from nomad import files, utils

from nomad.processing.data import Upload
from nomad.utils import get_logger, lnr


def handle_uploads(ready=None, quit=False):
    """
    Starts a daemon that will listen to files for new uploads. For each new
    upload it will initiate the processing and save the task in the upload user data,
    it will wait for processing to be completed and store the results in the upload
    user data.

    Arguments:
        ready (Event): optional, will be set when thread is ready
        quit: If true, will only handling one event and stop. Otherwise run forever.
    """

    @files.upload_put_handler
    def handle_upload_put(received_upload_id: str):
        logger = get_logger(__name__, upload_id=received_upload_id)
        logger.debug('Initiate upload processing')
        try:
            with lnr(logger, 'Could not load'):
                try:
                    upload = Upload.get(received_upload_id)
                except KeyError as e:
                    logger.error('Upload does not exist')
                    raise e

            if upload.upload_time is not None:
                logger.warn('Ignore upload notification, since file is already uploaded')
                raise StopIteration

            with lnr(logger, 'Save upload time'):
                upload.upload_time = datetime.now()

            with lnr(logger, 'Start processing'):
                upload.process()

        except Exception as e:
            logger.error('Exception while handling upload put notification.', exc_info=e)

        if quit:
            raise StopIteration

    utils.get_logger(__name__).debug('Start upload put notification handler.')
    if ready is not None:
        ready.set()
    handle_upload_put(received_upload_id='provided by decorator')


def handle_uploads_thread(quit=True):
    """ Same as :func:`handle_uploads` but run in a separate thread. """
    ready = Event()
    thread = Thread(target=lambda: handle_uploads(ready=ready, quit=quit))
    thread.start()
    ready.wait()
    return thread
