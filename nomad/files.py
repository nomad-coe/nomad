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
This module (and its main class :class:`Files`) represents an abstraction for NOMAD
file storage system.

Responsibilities: create, access files; create, receive, notify on, and access uploads.
"""
import os
from zipfile import ZipFile, BadZipFile
import shutil
from minio import Minio
import minio.error
import logging
import itertools

import nomad.config as config

LOGGER = logging.getLogger(__name__)

_client = Minio('%s:%s' % (config.minio.host, config.minio.port),
                access_key=config.minio.accesskey,
                secret_key=config.minio.secret,
                secure=False)

# ensure all neccessary buckets exist
try:
    _client.make_bucket(bucket_name=config.s3.uploads_bucket)
    LOGGER.info("Created uploads bucket with name %s." % config.s3.uploads_bucket)
except minio.error.BucketAlreadyOwnedByYou:
    LOGGER.debug("Uploads bucket with name %s already existed." % config.s3.uploads_bucket)


def get_presigned_upload_url(upload_id):
    return _client.presigned_put_object(config.s3.uploads_bucket, upload_id)


def create_curl_upload_cmd(presigned_url, file_dummy='<ZIPFILE>'):
    headers = 'Content-Type: application/octet-steam'
    return 'curl -X PUT "%s" -H "%s" -F file=@%s' % (presigned_url, headers, file_dummy)


def upload(upload_id):
    return Upload(upload_id)


def upload_put_handler(func):
    def upload_notifications(events):
        # The given events is a generator that will block and yield indefinetely.
        # Therefore, we have to use generator expressions and must not use list
        # comprehension. Same for chain vs chain.from_iterable.
        nested_event_records = (event['Records'] for event in events)
        event_records = itertools.chain.from_iterable(nested_event_records)

        for event_record in event_records:
            try:
                event_name = event_record['eventName']
                if event_name == 's3:ObjectCreated:Put':
                    LOGGER.debug('Received bucket upload event of type %s.' % event_name)
                    upload_id = event_record['s3']['object']['key']
                    yield upload_id
                else:
                    LOGGER.debug('Unhandled bucket event of type %s.' % event_name)
            except KeyError:
                LOGGER.warning(
                    'Unhandled bucket event due to unexprected event format: %s' %
                    event_record)

    def wrapper(*args, **kwargs):
        LOGGER.info('Start listening to uploads notifications.')

        events = _client.listen_bucket_notification(config.s3.uploads_bucket)

        upload_ids = upload_notifications(events)
        for upload_id in upload_ids:
            try:
                func(upload_id)
            except StopIteration:
                # Using StopIteration to allow clients to stop handling of events.
                logging.debug(
                    'Handling of upload notifications was stopped via StopIteration.')
                return
            except Exception as e:
                LOGGER.error(
                    'Unexpected exception in upload handler for upload:id:' %
                    upload_id, exc_info=e)

    return wrapper


class UploadError(Exception):
    IMPLEMENTATION_ERROR = 'implementation error'
    NOT_ZIP = 'upload is not a zip file'

    def __init__(self, msg, cause, code=IMPLEMENTATION_ERROR):
        super().__init__(msg, cause)
        self.code = code


class Upload():
    def __init__(self, upload_id):
        self.upload_id = upload_id
        self.upload_file = '%s/uploads/%s.zip' % (config.fs.tmp, upload_id)
        self.upload_extract_dir = '%s/uploads_extracted/%s' % (config.fs.tmp, upload_id)
        self.filelist = None

        try:
            _client.stat_object(config.s3.uploads_bucket, upload_id)
        except minio.error.NoSuchKey:
            raise KeyError(self.upload_id)

    # There is not good way to capsule decorators in a class:
    # https://medium.com/@vadimpushtaev/decorator-inside-python-class-1e74d23107f6
    class Decorators:
        @classmethod
        def log_upload_error(cls, decorated):
            def wrapper(self, *args, **kwargs):
                try:
                    return decorated(self, *args, **kwargs)
                except Exception as e:
                    msg = 'Could not %s upload %s.' % (decorated.__name__, self.upload_id)
                    LOGGER.error(msg, exc_info=e)
                    raise UploadError(msg, e)
            return wrapper

    @Decorators.log_upload_error
    def open(self):
        try:
            _client.fget_object(config.s3.uploads_bucket, self.upload_id, self.upload_file)
        except minio.error.NoSuchKey:
            raise KeyError(self.upload_id)

        zipFile = None
        try:
            zipFile = ZipFile(self.upload_file)
            zipFile.extractall(self.upload_extract_dir)
            self.filelist = [zipInfo.filename for zipInfo in zipFile.filelist]
        except BadZipFile as e:
            raise UploadError('Upload is not a zip file', e, UploadError.NOT_ZIP)
        finally:
            if zipFile is not None:
                zipFile.close()

    @Decorators.log_upload_error
    def close(self):
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

    @Decorators.log_upload_error
    def open_file(self, filename, *args, **kwargs):
        return open('%s/%s' % (self.upload_extract_dir, filename), *args, **kwargs)
