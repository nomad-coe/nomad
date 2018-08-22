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
This file storage abstraction currently uses the object storage API
http://minio.io to manage and organize files. Object storage
organizes files in *buckets/ids*, with small amounts of *buckets* and virtually
unlimited numbers of *ids*. *Ids* can contain delimiters like `/` to mimic
filesystem structure. There is a 1024 utf-8 character limit on *id* length.

The file storage is organized in multiple buckets:

* *uploads*: used for uploaded user code input/output archives. Currently only .zip files
  are suported

Presigned URLs
--------------
Users (or GUI clients) can upload files directly to the object storage system. To avoid
authentication hassly, presigned URLs can be created that can be used directly to safely
*PUT* files.

.. autofunction:: nomad.files.get_presigned_upload_url
.. autofunction:: nomad.files.create_curl_upload_cmd

Uploads
-------
.. autoclass:: Upload

"""
from typing import Callable, List, Any, Generator, IO, TextIO, cast
import sys
import os
import os.path
from zipfile import ZipFile, BadZipFile
import shutil
from minio import Minio
import minio.error
import logging
import itertools
import hashlib
import base64
from contextlib import contextmanager
import gzip
import io
import json

import nomad.config as config

logger = logging.getLogger(__name__)

_client = None

if _client is None and 'sphinx' not in sys.modules:
    _client = Minio('%s:%s' % (config.minio.host, config.minio.port),
                    access_key=config.minio.accesskey,
                    secret_key=config.minio.secret,
                    secure=False)

    # ensure all neccessary buckets exist
    def ensure_bucket(name):
        try:
            _client.make_bucket(bucket_name=name)
            logger.info("Created uploads bucket with name %s." % name)
        except minio.error.BucketAlreadyOwnedByYou:
            pass

    ensure_bucket(config.files.uploads_bucket)
    ensure_bucket(config.files.archive_bucket)


def get_presigned_upload_url(upload_id: str) -> str:
    """
    Generates a presigned upload URL. Presigned URL allows users (and their client programs)
    to safely *PUT* a single file without further authorization or API to the *uploads* bucket
    using the given ``upload_id``. Example usages for presigned URLs include
    browser based uploads or simple *curl* commands (see also :func:`create_curl_upload_cmd`).

    Arguments:
        upload_id: The upload id for the uploaded file.

    Returns:
        The presigned URL string.
    """
    return _client.presigned_put_object(config.files.uploads_bucket, upload_id)


def create_curl_upload_cmd(presigned_url: str, file_dummy: str='<ZIPFILE>') -> str:
    """Creates a readymade curl command for uploading.

    Arguments:
        presigned_url: The presigned URL to base the command on.

    Kwargs:
        file_dummy: A placeholder for the file that the user/client has to replace.

    Returns:
        The curl shell command with correct method, url, headers, etc.
    """
    headers = 'Content-Type: application/octet-steam'
    return 'curl -X PUT "%s" -H "%s" -F file=@%s' % (presigned_url, headers, file_dummy)


def upload_put_handler(func: Callable[[str], None]) -> Callable[[], None]:
    def upload_notifications(events: List[Any]) -> Generator[str, None, None]:
        for event in events:
            for event_record in event['Records']:
                try:
                    event_name = event_record['eventName']
                    if event_name == 's3:ObjectCreated:Put':
                        upload_id = event_record['s3']['object']['key']
                        logger.debug('Received bucket upload event of for upload %s.' % upload_id)
                        yield upload_id
                        break  # only one per record, pls
                    else:
                        logger.debug('Unhanled bucket event %s.' % event_name)
                except KeyError:
                    logger.warning(
                        'Unhandled bucket event due to unexprected event format: %s' %
                        event_record)

    def wrapper(*args, **kwargs) -> None:
        logger.info('Start listening to uploads notifications.')

        _client.remove_all_bucket_notification(config.files.uploads_bucket)
        events = _client.listen_bucket_notification(
            config.files.uploads_bucket,
            events=['s3:ObjectCreated:*'])

        upload_ids = upload_notifications(events)
        for upload_id in upload_ids:
            try:
                func(upload_id)
            except StopIteration:
                # Using StopIteration to allow clients to stop handling of events.
                logging.debug(
                    'Handling of upload notifications was stopped via StopIteration.')
                return
            except Exception:
                pass

    return wrapper


class UploadError(Exception):
    def __init__(self, msg, cause):
        super().__init__(msg, cause)


class Upload():
    """
    Instances represent an uploaded file in the object storage. Class supports open/close,
    i.e. extract .zip files, and opening contained files. Some functions are only available
    for open (i.e. tmp. downloaded and extracted uploads) uploads.

    This class is also a context manager that opens and closes the upload respectively.

    Arguments:
        upload_id: The upload of this uploaded file.

    Attributes:
        upload_file: The path of the tmp version of this file for an open upload.
        upload_extract_dir: The path of the tmp directory with the extracted contents.
        filelist: A list of filenames relative to the .zipped upload root.
        metadata: The upload object storage metadata.
    """
    def __init__(self, upload_id: str) -> None:
        self.upload_id = upload_id
        self.upload_file: str = os.path.join(config.fs.tmp, 'uploads', upload_id)
        self.upload_extract_dir: str = os.path.join(config.fs.tmp, 'uploads_extracted', upload_id)
        self.filelist: List[str] = None

        try:
            self.metadata = _client.stat_object(config.files.uploads_bucket, upload_id).metadata
        except minio.error.NoSuchKey:
            raise KeyError(self.upload_id)

    # There is not good way to capsule decorators in a class:
    # https://medium.com/@vadimpushtaev/decorator-inside-python-class-1e74d23107f6
    class Decorators:
        @classmethod
        def handle_errors(cls, decorated):
            def wrapper(self, *args, **kwargs):
                try:
                    return decorated(self, *args, **kwargs)
                except Exception as e:
                    msg = 'Could not %s upload %s.' % (decorated.__name__, self.upload_id)
                    logger.error(msg, exc_info=e)
                    raise UploadError(msg, e)
            return wrapper

    @Decorators.handle_errors
    def hash(self) -> str:
        """ Calculates the first 28 bytes of a websafe base64 encoded SHA512 of the upload. """
        hash = hashlib.sha512()
        with open(self.upload_file, 'rb') as f:
            for data in iter(lambda: f.read(65536), b''):
                hash.update(data)

        return base64.b64encode(hash.digest(), altchars=b'-_')[0:28].decode('utf-8')

    @Decorators.handle_errors
    def open(self) -> None:
        """
        Opens the upload. This means the uploaed files gets tmp. downloaded and extracted.

        Raises:
            UploadError: If some IO went wrong.
            KeyError: If the upload does not exist.
        """
        os.makedirs(os.path.join(config.fs.tmp, 'uploads'), exist_ok=True)
        os.makedirs(os.path.join(config.fs.tmp, 'uploads_extracted'), exist_ok=True)

        try:
            _client.fget_object(config.files.uploads_bucket, self.upload_id, self.upload_file)
        except minio.error.NoSuchKey:
            raise KeyError(self.upload_id)

        zipFile = None
        try:
            zipFile = ZipFile(self.upload_file)
            zipFile.extractall(self.upload_extract_dir)
            self.filelist = [zipInfo.filename for zipInfo in zipFile.filelist]
        except BadZipFile as e:
            raise UploadError('Upload is not a zip file', e)
        finally:
            if zipFile is not None:
                zipFile.close()

    @Decorators.handle_errors
    def close(self) -> None:
        """
        Closes the upload. This means the tmp. files are deleted.

        Raises:
            UploadError: If some IO went wrong.
            KeyError: If the upload does not exist.
        """
        try:
            os.remove(self.upload_file)
            shutil.rmtree(self.upload_extract_dir)
        except FileNotFoundError:
            raise KeyError(self.upload_id)

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, exc_type, exc, exc_tb):
        self.close()

    @Decorators.handle_errors
    def open_file(self, filename: str, *args, **kwargs) -> IO[Any]:
        """ Opens a file within an open upload and returns a file like. """
        return open(self.get_path(filename), *args, **kwargs)

    def get_path(self, filename: str) -> str:
        """ Returns the tmp directory relative version of a filename. """
        return os.path.join(self.upload_extract_dir, filename)


@contextmanager
def write_archive_json(archive_id) -> Generator[TextIO, None, None]:
    """ Context manager that yiels a file-like to write the archive json. """
    binary_out = io.BytesIO()
    if config.files.compress_archive:
        gzip_wrapper = cast(TextIO, gzip.open(binary_out, 'wt'))
        out = gzip_wrapper
        metadata = {'Content-Encoding': 'gzip'}
    else:
        text_wrapper = io.TextIOWrapper(binary_out, encoding='utf-8')
        out = text_wrapper
        metadata = {}

    try:
        yield out
    finally:
        out.flush()
        binary_out.seek(0)
        length = len(binary_out.getvalue())

        _client.put_object(
            config.files.archive_bucket, archive_id, binary_out, length=length,
            content_type='application/json',
            metadata=metadata)

        out.close()
        binary_out.close()


def open_archive_json(archive_id) -> IO:
    """ Returns a file-like to read the archive json. """
    # The result already is a file-like and due to the Content-Encoding metadata is
    # will automatically be un-gzipped.
    try:
        return _client.get_object(config.files.archive_bucket, archive_id)
    except minio.error.NoSuchKey:
        raise KeyError()
