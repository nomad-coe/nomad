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
"""

from typing import List, Any, ContextManager, Tuple, Generator
from mongoengine import StringField, DateTimeField, DictField
import logging
from structlog import wrap_logger
from contextlib import contextmanager

from nomad import utils, coe_repo, datamodel
from nomad.files import PathObject, ArchiveBasedStagingUploadFiles, ExtractError
from nomad.processing.base import Proc, Chord, process, task, PENDING, SUCCESS, FAILURE
from nomad.parsing import parsers, parser_dict
from nomad.normalizing import normalizers


class Calc(Proc, datamodel.Calc):
    """
    Instances of this class represent calculations. This class manages the elastic
    search index entry, files, and archive for the respective calculation.

    It also contains the calculations processing and its state.

    The attribute list, does not include the various repository properties generated
    while parsing, including ``program_name``, ``program_version``, etc.

    Attributes:
        calc_id: the calc_id of this calc
        parser: the name of the parser used to process this calc
        upload_id: the id of the upload used to create this calculation
        mainfile: the mainfile (including path in upload) that was used to create this calc
    """
    calc_id = StringField(primary_key=True)
    upload_id = StringField()
    mainfile = StringField()
    parser = StringField()

    meta: Any = {
        'indices': [
            'upload_id', 'mainfile', 'code', 'parser', 'tasks_status'
        ]
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._parser_backend = None
        self._upload: Upload = None
        self._upload_files: ArchiveBasedStagingUploadFiles = None
        self._calc_proc_logwriter = None
        self._calc_proc_logwriter_ctx: ContextManager = None

    @classmethod
    def get(cls, id):
        return cls.get_by_id(id, 'calc_id')

    @property
    def mainfile_file(self) -> PathObject:
        return self.upload_files.raw_file_object(self.mainfile)

    @property
    def upload(self) -> 'Upload':
        if not self._upload:
            self._upload = Upload.get(self.upload_id)
        return self._upload

    @property
    def upload_files(self) -> ArchiveBasedStagingUploadFiles:
        if not self._upload_files:
            self._upload_files = ArchiveBasedStagingUploadFiles(self.upload_id, is_authorized=lambda: True, local_path=self.upload.local_path)
        return self._upload_files

    def get_logger(self, **kwargs):
        logger = super().get_logger()
        logger = logger.bind(
            upload_id=self.upload_id, mainfile=self.mainfile, calc_id=self.calc_id, **kwargs)

        return logger

    def get_calc_logger(self, **kwargs):
        """
        Returns a wrapped logger that additionally saves all entries to the calculation
        processing log in the archive.
        """
        logger = self.get_logger(**kwargs)

        if self._calc_proc_logwriter is None:
            self._calc_proc_logwriter_ctx = self.upload_files.archive_log_file(self.calc_id, 'wt')
            self._calc_proc_logwriter = self._calc_proc_logwriter_ctx.__enter__()  # pylint: disable=E1101

        def save_to_calc_log(logger, method_name, event_dict):
            program = event_dict.get('normalizer', 'parser')
            event = event_dict.get('event', '')
            entry = '[%s] %s: %s' % (method_name, program, event)
            if len(entry) > 120:
                self._calc_proc_logwriter.write(entry[:120])
                self._calc_proc_logwriter.write('...')
            else:
                self._calc_proc_logwriter.write(entry)
            self._calc_proc_logwriter.write('\n')
            return event_dict

        return wrap_logger(logger, processors=[save_to_calc_log])

    @process
    def process_calc(self):
        logger = self.get_logger()
        if self.upload is None:
            logger.error('calculation upload does not exist')

        try:
            self.parsing()
            self.normalizing()
            self.archiving()
        finally:
            # close loghandler that was not closed due to failures
            try:
                if self._calc_proc_logwriter is not None:
                    self._calc_proc_logwriter.close()
                    self._calc_proc_logwriter = None
            except Exception as e:
                logger.error('could not close calculation proc log', exc_info=e)

            # inform parent proc about completion
            self.upload.completed_child()

    @task
    def parsing(self):
        context = dict(parser=self.parser, step=self.parser)
        logger = self.get_calc_logger(**context)
        parser = parser_dict[self.parser]

        with utils.timer(logger, 'parser executed', input_size=self.mainfile_file.size):
            self._parser_backend = parser.run(
                self.upload_files.raw_file_object(self.mainfile).os_path, logger=logger)

        self._parser_backend.openNonOverlappingSection('section_calculation_info')
        self._parser_backend.addValue('upload_id', self.upload_id)
        self._parser_backend.addValue('calc_id', self.calc_id)
        self._parser_backend.addValue('calc_hash', self.upload_files.calc_hash(self.mainfile))
        self._parser_backend.addValue('main_file', self.mainfile)
        self._parser_backend.addValue('parser_name', self.parser)

        if self._parser_backend.status[0] != 'ParseSuccess':
            logger.error(self._parser_backend.status[1])
            error = self._parser_backend.status[1]
            self._parser_backend.addValue('parse_status', 'ParseFailure')
            self.fail(error, level=logging.DEBUG, **context)
        else:
            self._parser_backend.addValue('parse_status', 'ParseSuccess')

        self._parser_backend.closeNonOverlappingSection('section_calculation_info')

        self._parser_backend.openNonOverlappingSection('section_repository_info')
        self._parser_backend.addValue('repository_archive_gid', '%s/%s' % (self.upload_id, self.calc_id))
        self._parser_backend.addValue(
            'repository_filepaths', self.upload_files.calc_files(self.mainfile))
        self._parser_backend.closeNonOverlappingSection('section_repository_info')

        self.add_processor_info(self.parser)

    @contextmanager
    def use_parser_backend(self, processor_name):
        self._parser_backend.reset_status()
        yield self._parser_backend
        self.add_processor_info(processor_name)

    def add_processor_info(self, processor_name: str) -> None:
        self._parser_backend.openContext('/section_calculation_info/0')
        self._parser_backend.openNonOverlappingSection('section_archive_processing_info')
        self._parser_backend.addValue('archive_processor_name', processor_name)

        if self._parser_backend.status[0] == 'ParseSuccess':
            warnings = getattr(self._parser_backend, '_warnings', [])
            if len(warnings) > 0:
                self._parser_backend.addValue('archive_processor_status', 'WithWarnings')
                self._parser_backend.addValue('archive_processor_warning_number', len(warnings))
                self._parser_backend.addArrayValues('archive_processor_warnings', [str(warning) for warning in warnings])
            else:
                self._parser_backend.addValue('archive_processor_status', 'Success')
        else:
            errors = self._parser_backend.tasks_status[1]
            self._parser_backend.addValue('archive_processor_error', str(errors))

        self._parser_backend.closeNonOverlappingSection('section_archive_processing_info')
        self._parser_backend.closeContext('/section_calculation_info/0')

    @task
    def normalizing(self):
        for normalizer in normalizers:
            normalizer_name = normalizer.__name__
            context = dict(normalizer=normalizer_name, step=normalizer_name)
            logger = self.get_calc_logger(**context)

            with utils.timer(
                    logger, 'normalizer executed', input_size=self.mainfile_file.size):
                with self.use_parser_backend(normalizer_name) as backend:
                    normalizer(backend).normalize(logger=logger)

            failed = self._parser_backend.status[0] != 'ParseSuccess'
            if failed:
                logger.error(self._parser_backend.status[1])
                error = self._parser_backend.status[1]
                self.fail(error, level=logging.WARNING, **context)
                break
            else:
                logger.debug(
                    'completed normalizer successfully', normalizer=normalizer_name)

    @task
    def archiving(self):
        logger = self.get_logger()

        # persist the repository metadata
        with utils.timer(logger, 'indexed', step='index'):
            self.upload_files.metadata.insert(self._parser_backend.metadata())

        # persist the archive
        with utils.timer(
                logger, 'archived', step='archive',
                input_size=self.mainfile_file.size) as log_data:
            with self.upload_files.archive_file(self.calc_id, 'wt') as out:
                self._parser_backend.write_json(out, pretty=True)

            log_data.update(archive_size=self.upload_files.archive_file_object(self.calc_id).size)

        # close loghandler
        if self._calc_proc_logwriter is not None:
            with utils.timer(
                    logger, 'archived log', step='archive_log',
                    input_size=self.mainfile_file.size) as log_data:
                self._calc_proc_logwriter_ctx.__exit__(None, None, None)  # pylint: disable=E1101
                self._calc_proc_logwriter = None

                log_data.update(log_size=self.upload_files.archive_log_file_object(self.calc_id).size)


class Upload(Chord, datamodel.Upload):
    """
    Represents uploads in the databases. Provides persistence access to the files storage,
    and processing state.

    Attributes:
        name: optional user provided upload name
        local_path: optional local path, e.g. for files that are already somewhere on the server
        metadata: optional user provided additional meta data
        upload_id: the upload id generated by the database
        upload_time: the timestamp when the system realised the upload
        user_id: the id of the user that created this upload
    """
    id_field = 'upload_id'

    upload_id = StringField(primary_key=True)

    name = StringField(default=None)
    local_path = StringField(default=None)
    metadata = DictField(default=None)
    upload_time = DateTimeField()
    user_id = StringField(required=True)

    meta: Any = {
        'indexes': [
            'user_id', 'tasks_status'
        ]
    }

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._upload_files: ArchiveBasedStagingUploadFiles = None

    @classmethod
    def get(cls, id):
        return cls.get_by_id(id, 'upload_id')

    @classmethod
    def user_uploads(cls, user: coe_repo.User) -> List['Upload']:
        """ Returns all uploads for the given user. Currently returns all uploads. """
        return cls.objects(user_id=str(user.user_id))

    @property
    def uploader(self):
        return coe_repo.User.from_user_id(self.user_id)

    def get_logger(self, **kwargs):
        logger = super().get_logger()
        logger = logger.bind(upload_id=self.upload_id, **kwargs)
        return logger

    @classmethod
    def create(cls, **kwargs) -> 'Upload':
        """
        Creates a new upload for the given user, a user given name is optional.
        It will populate the record with a signed url and pending :class:`UploadProc`.
        The upload will be already saved to the database.

        Arguments:
            user (coe_repo.User): The user that created the upload.
        """
        user: coe_repo.User = kwargs['user']
        del(kwargs['user'])
        if 'upload_id' not in kwargs:
            kwargs.update(upload_id=utils.create_uuid())
        kwargs.update(user_id=str(user.user_id))
        self = super().create(**kwargs)

        self._continue_with('uploading')

        return self

    def delete(self):
        """ Deletes this upload process state entry and its calcs. """
        Calc.objects(upload_id=self.upload_id).delete()
        super().delete()

    @process
    def delete_upload(self):
        """
        Deletes of the upload, including its processing state and
        staging files.
        """
        logger = self.get_logger()

        with utils.lnr(logger, 'staged upload delete failed'):
            with utils.timer(
                    logger, 'staged upload deleted', step='delete',
                    upload_size=self.upload_files.size):
                import time
                time.sleep(10)
                self.upload_files.delete()
                self.delete()

        return True  # do not save the process status on the delete upload

    @process
    def commit_upload(self):
        """
        Moves the upload out of staging to add it to the coe repository. It will
        pack the staging upload files in to public upload files, add entries to the
        coe repository db and remove this instance and its calculation from the
        processing state db.
        """
        logger = self.get_logger()

        with utils.lnr(logger, 'commit failed'):
            with utils.timer(
                    logger, 'upload added to repository', step='commit',
                    upload_size=self.upload_files.size):
                coe_repo.Upload.add(self, self.metadata)

            with utils.timer(
                    logger, 'staged upload files packed', step='commit',
                    upload_size=self.upload_files.size):
                self.upload_files.pack()

            with utils.timer(
                    logger, 'staged upload deleted', step='commit',
                    upload_size=self.upload_files.size):
                self.upload_files.delete()
                self.delete()

        return True  # do not save the process status on the delete upload

    @process
    def process_upload(self):
        self.extracting()
        self.parse_all()

    @task
    def uploading(self):
        pass

    @property
    def upload_files(self) -> ArchiveBasedStagingUploadFiles:
        if not self._upload_files:
            self._upload_files = ArchiveBasedStagingUploadFiles(self.upload_id, is_authorized=lambda: True, local_path=self.local_path)
        return self._upload_files

    @task
    def extracting(self):
        """
        Task performed before the actual parsing/normalizing. Extracting and bagging
        the uploaded files, computing all keys, create an *upload* entry in the NOMAD-coe
        repository db, etc.
        """
        # extract the uploaded file, this will also create a bagit bag.
        logger = self.get_logger()
        try:
            with utils.timer(
                    logger, 'upload extracted', step='extracting',
                    upload_size=self.upload_files.size):
                self.upload_files.extract()
        except KeyError:
            self.fail('processing requested for non existing upload', log_level=logging.ERROR)
            return
        except ExtractError:
            self.fail('bad .zip/.tar file', log_level=logging.INFO)
            return

    def match_mainfiles(self) -> Generator[Tuple[str, object], None, None]:
        """
        Generator function that matches all files in the upload to all parsers to
        determine the upload's mainfiles.

        Returns:
            Tuples of mainfile, filename, and parsers
        """
        for filename in self.upload_files.raw_file_manifest():
            for parser in parsers:
                try:
                    with self.upload_files.raw_file(filename) as mainfile_f:
                        if parser.is_mainfile(filename, lambda fn: mainfile_f):
                            yield filename, parser
                except Exception as e:
                    self.get_logger().error(
                        'exception while matching pot. mainfile',
                        mainfile=filename, exc_info=e)

    @task
    def parse_all(self):
        """
        Identified mainfile/parser combinations among the upload's files, creates
        respective :class:`Calc` instances, and triggers their processing.
        """
        logger = self.get_logger()

        # TODO: deal with multiple possible parser specs
        with utils.timer(
                logger, 'upload extracted', step='matching',
                upload_size=self.upload_files.size):
            total_calcs = 0
            for filename, parser in self.match_mainfiles():
                calc = Calc.create(
                    calc_id=self.upload_files.calc_id(filename),
                    mainfile=filename, parser=parser.name,
                    upload_id=self.upload_id)

                calc.process_calc()
                total_calcs += 1

        # have to save the total_calcs information for chord management
        self.spwaned_childred(total_calcs)

    def join(self):
        self.cleanup()

    @task
    def cleanup(self):
        # TODO issue #83
        with utils.timer(
                self.get_logger(), 'pack staging upload', step='cleaning',
                upload_size=self.upload_files.size):
            pass

    @property
    def processed_calcs(self):
        return Calc.objects(upload_id=self.upload_id, tasks_status__in=[SUCCESS, FAILURE]).count()

    @property
    def total_calcs(self):
        return Calc.objects(upload_id=self.upload_id).count()

    @property
    def failed_calcs(self):
        return Calc.objects(upload_id=self.upload_id, tasks_status=FAILURE).count()

    @property
    def pending_calcs(self):
        return Calc.objects(upload_id=self.upload_id, tasks_status=PENDING).count()

    def all_calcs(self, start, end, order_by='mainfile'):
        return Calc.objects(upload_id=self.upload_id)[start:end].order_by(order_by)

    @property
    def calcs(self):
        return Calc.objects(upload_id=self.upload_id)
