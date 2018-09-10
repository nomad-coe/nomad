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
This module comprises a set of persistent document classes that hold all user related
data. These are information about users, their uploads and datasets, the associated
calculations, and files


.. autoclass:: Calc
    :members:
.. autoclass:: Upload
    :members:
.. autoclass:: DataSet
.. autoclass:: User

"""

from typing import List, Any
import sys
from datetime import datetime
from elasticsearch.exceptions import NotFoundError
from mongoengine import \
    Document, EmailField, StringField, BooleanField, DateTimeField, \
    ListField, DictField, ReferenceField, IntField, connect
import mongoengine.errors
import logging

from nomad import config, files, utils
from nomad.repo import RepoCalc
from nomad.user import User, me
from nomad.processing.base import Proc, process, task, PENDING
from nomad.parsing import LocalBackend, parsers, parser_dict
from nomad.normalizing import normalizers
from nomad.utils import get_logger, lnr


class NotAllowedDuringProcessing(Exception): pass


class Calc(Proc):
    """
    Instances of this class represent calculations. This class manages the elastic
    search index entry, files, and archive for the respective calculation.

    It also contains the calculations processing and its state.

    The attribute list, does not include the various repository properties generated
    while parsing, including ``program_name``, ``program_version``, etc.

    Attributes:
        archive_id: the hash based archive id of the calc
        parser: the name of the parser used to process this calc
        upload_id: the id of the upload used to create this calculation
        mainfile: the mainfile (including path in upload) that was used to create this calc
        mainfile_tmp_path: path to the mainfile extracted for processing
    """
    archive_id = StringField(primary_key=True)
    upload_id = StringField()
    mainfile = StringField()
    parser = StringField()
    mainfile_tmp_path = StringField()

    meta: Any = {
        'indices': [
            'upload_id', 'mainfile', 'code', 'parser'
        ]
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._parser_backend = None
        self._upload = None

    @classmethod
    def get(cls, id):
        return cls.get_by_id(id, 'archive_id')

    def delete(self):
        """
        Delete this calculation and all associated data. This includes all files,
        the archive, and this search index entry.
        TODO is this needed? Or do we always delete hole uploads in bulk.
        """
        # delete the archive
        if self.archive_id is not None:
            files.delete_archive(self.archive_id)

        # delete the search index entry
        try:
            elastic_entry = RepoCalc.get(self.archive_id)
            if elastic_entry is not None:
                elastic_entry.delete()
        except NotFoundError:
            pass

        # delete this mongo document
        super().delete()

    def get_logger(self, **kwargs):
        upload_hash, calc_hash = self.archive_id.split('/')
        logger = super().get_logger()
        logger = logger.bind(
            upload_id=self.upload_id, mainfile=self.mainfile,
            upload_hash=upload_hash, calc_hash=calc_hash, **kwargs)
        return logger

    @property
    def json_dict(self):
        """ A json serializable dictionary representation. """
        data = {
            'archive_id': self.archive_id,
            'mainfile': self.mainfile,
            'upload_id': self.upload_id,
            'parser': self.parser
        }
        data.update(super().json_dict)
        return {key: value for key, value in data.items() if value is not None}

    @process
    def process(self):
        self._upload = Upload.get(self.upload_id)
        if self._upload is None:
            get_logger().error('calculation upload does not exist')

        try:
            self.parsing()
            self.normalizing()
            self.archiving()
        finally:
            self._upload.calc_proc_completed()

    @task
    def parsing(self):
        self._parser_backend = parser_dict[self.parser].run(self.mainfile_tmp_path)
        if self._parser_backend.status[0] != 'ParseSuccess':
            error = self._parser_backend.status[1]
            self.fail(error, level=logging.DEBUG)

    @task
    def normalizing(self):
        for normalizer in normalizers:
            normalizer_name = normalizer.__name__
            normalizer(self._parser_backend).normalize()
            if self._parser_backend.status[0] != 'ParseSuccess':
                error = self._parser_backend.status[1]
                self.fail(error, normalizer=normalizer_name, level=logging.WARNING)
                return
            self.get_logger().debug(
                'completed normalizer successfully', normalizer=normalizer_name)

    @task
    def archiving(self):
        upload_hash, calc_hash = self.archive_id.split('/')
        # persist to elastic search
        RepoCalc.create_from_backend(
            self._parser_backend,
            upload_hash=upload_hash,
            calc_hash=calc_hash,
            upload_id=self.upload_id,
            mainfile=self.mainfile,
            upload_time=self._upload.upload_time,
            staging=True,
            restricted=False,
            user_id=self._upload.user_id)

        # persist the archive
        with files.write_archive_json(self.archive_id) as out:
            self._parser_backend.write_json(out, pretty=True)


class Upload(Proc):
    """
    Represents uploads in the databases. Provides persistence access to the files storage,
    and processing state.

    Attributes:
        name: optional user provided upload name
        additional_metadata: optional user provided additional meta data
        upload_id: the upload id generated by the database
        in_staging: true if the upload is still in staging and can be edited by the uploader
        is_private: true if the upload and its derivitaves are only visible to the uploader
        presigned_url: the presigned url for file upload
        upload_time: the timestamp when the system realised the upload
        upload_hash: the hash of the uploaded file
        user_id: the id of the user that created this upload
    """
    id_field = 'upload_id'

    upload_id = StringField(primary_key=True)

    name = StringField(default=None)
    additional_metadata = DictField(default=None)

    in_staging = BooleanField(default=True)
    is_private = BooleanField(default=False)

    presigned_url = StringField()
    upload_command = StringField()
    upload_time = DateTimeField()
    upload_hash = StringField(default=None)

    processed_calcs = IntField(default=0)
    total_calcs = IntField(default=-1)

    user_id = StringField(required=True)

    meta: Any = {
        'indexes': [
            'upload_hash',
            'user_id'
        ]
    }

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._upload = None

    @classmethod
    def get(cls, id):
        return cls.get_by_id(id, 'upload_id')

    @classmethod
    def user_uploads(cls, user: User) -> List['Upload']:
        """ Returns all uploads for the given user. Currently returns all uploads. """
        return cls.objects()

    def get_logger(self, **kwargs):
        logger = super().get_logger()
        logger = logger.bind(upload_id=self.upload_id, **kwargs)
        return logger

    def delete(self):
        logger = self.get_logger(task='delete')

        if not (self.completed or self.is_stale or self.current_task == 'uploading'):
            raise NotAllowedDuringProcessing()

        with lnr(logger, 'delete upload file'):
            try:
                files.Upload(self.upload_id).delete()
            except KeyError:
                if self.current_task == 'uploading':
                    logger.debug(
                        'Upload exist, but file does not exist. '
                        'It was probably aborted and deleted.')
                else:
                    logger.debug('Upload exist, but uploaded file does not exist.')

        with lnr(logger, 'deleting calcs'):
            # delete archive files
            files.delete_archives(upload_hash=self.upload_hash)

            # delete repo entries
            RepoCalc.search().query('match', upload_id=self.upload_id).delete()

            # delete calc processings
            Calc.objects(upload_id=self.upload_id).delete()

        with lnr(logger, 'deleting upload'):
            super().delete()

    @classmethod
    def create(cls, **kwargs) -> 'Upload':
        """
        Creates a new upload for the given user, a user given name is optional.
        It will populate the record with a signed url and pending :class:`UploadProc`.
        The upload will be already saved to the database.
        """
        self = super().create(**kwargs)
        self.presigned_url = files.get_presigned_upload_url(self.upload_id)
        self.upload_command = files.create_curl_upload_cmd(self.presigned_url, 'your_file')
        self._continue_with('uploading')
        return self

    @property
    def is_stale(self) -> bool:
        if self.current_task == 'uploading' and self.upload_time is None:
            return (datetime.now() - self.create_time).days > 1
        else:
            return False

    @property
    def json_dict(self) -> dict:
        """ A json serializable dictionary representation. """
        data = {
            'name': self.name,
            'additional_metadata': self.additional_metadata,
            'upload_id': self.upload_id,
            'presigned_url': self.presigned_url,
            'upload_command': self.upload_command,
            'upload_time': self.upload_time.isoformat() if self.upload_time is not None else None,
            'is_stale': self.is_stale,
        }
        data.update(super().json_dict)
        return {key: value for key, value in data.items() if value is not None}

    @process
    def process(self):
        self.extracting()
        self.parse_all()

    @task
    def uploading(self):
        pass

    @task
    def extracting(self):
        logger = self.get_logger()
        try:
            self._upload = files.Upload(self.upload_id)
            self._upload.open()
            logger.debug('opened upload')
        except KeyError as e:
            self.fail('process request for non existing upload', level=logging.INFO)
            return

        try:
            self.upload_hash = self._upload.hash()
        except files.UploadError as e:
            self.fail('could not create upload hash', e)
            return

        if RepoCalc.upload_exists(self.upload_hash):
            self.fail('The same file was already uploaded and processed.', level=logging.INFO)
            return

    @task
    def parse_all(self):
        # TODO: deal with multiple possible parser specs
        self.total_calcs = 0
        for filename in self._upload.filelist:
            for parser in parsers:
                try:
                    if parser.is_mainfile(filename, lambda fn: self._upload.open_file(fn)):
                        tmp_mainfile = self._upload.get_path(filename)
                        calc = Calc.create(
                            archive_id='%s/%s' % (self.upload_hash, utils.hash(filename)),
                            mainfile=filename, parser=parser.name,
                            mainfile_tmp_path=tmp_mainfile,
                            upload_id=self.upload_id)

                        calc.process()
                        self.total_calcs += 1
                except Exception as e:
                    self.warning(
                        'exception while matching pot. mainfile',
                        mainfile=filename, exc_info=e)

        if self.total_calcs == 0:
            self.cleanup()

        # have to save the total_calcs information
        self.save()

    @task
    def cleanup(self):
        try:
            upload = files.Upload(self.upload_id)
        except KeyError as e:
            upload_proc.fail('Upload does not exist', exc_info=e)
            return

        upload.close()
        self.get_logger().debug('closed upload')

    def calc_proc_completed(self):
        processed_calcs, (total_calcs,) = self.incr_counter(
            'processed_calcs', other_fields=['total_calcs'])

        if processed_calcs == total_calcs:
            self.cleanup()

    @property
    def calcs(self):
        return Calc.objects(upload_id=self.upload_hash)
