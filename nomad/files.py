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
from zipfile import ZipFile
import shutil
from minio import Minio
from minio.error import BucketAlreadyOwnedByYou
import logging
import itertools

import nomad.config as config

logger = logging.getLogger(__name__)

_client = Minio('%s:%s' % (config.minio.host, config.minio.port),
                access_key=config.minio.accesskey,
                secret_key=config.minio.secret,
                secure=False)

# ensure all neccessary buckets exist
try:
    _client.make_bucket(bucket_name=config.s3.uploads_bucket)
    logger.info("Created uploads bucket with name %s." % config.s3.uploads_bucket)
except BucketAlreadyOwnedByYou:
    logger.debug("Uploads bucket with name %s already existed." % config.s3.uploads_bucket)


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
                    logger.debug('Received bucket upload event of type %s.' % event_name)
                    upload_id = event_record['s3']['object']['key']
                    yield upload_id
                else:
                    logger.debug('Unhandled bucket event of type %s.' % event_name)
            except KeyError:
                logger.warning(
                    'Unhandled bucket event due to unexprected event format: %s' %
                    event_record)

    def wrapper(*args, **kwargs):
        logger.info('Start listening to uploads notifications.')

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
                logger.error(
                    'Unexpected exception in upload handler for upload:id:' %
                    upload_id, exc_info=e)

    return wrapper


class Upload():
    def __init__(self, upload_id):
        self.upload_id = upload_id
        self.upload_file = '%s/uploads/%s.zip' % (config.fs.tmp, upload_id)
        self.upload_extract_dir = '%s/uploads_extracted/%s' % (config.fs.tmp, upload_id)
        self.filelist = None

    def open(self):
        _client.fget_object(config.s3.uploads_bucket, self.upload_id, self.upload_file)
        zipFile = ZipFile(self.upload_file)
        zipFile.extractall(self.upload_extract_dir)
        self.filelist = [zipInfo.filename for zipInfo in zipFile.filelist]
        zipFile.close()

    def close(self):
        os.remove(self.upload_file)
        shutil.rmtree(self.upload_extract_dir)

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, exc_type, exc, exc_tb):
        self.close()

    def open_file(self, filename, *args, **kwargs):
        return open('%s/%s' % (self.upload_extract_dir, filename), *args, **kwargs)
