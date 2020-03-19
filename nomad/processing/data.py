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

'''
This module comprises a set of persistent document classes that hold all user related
data. These are information about users, their uploads and datasets, the associated
calculations, and files


.. autoclass:: Calc

.. autoclass:: Upload

'''

from typing import cast, List, Any, ContextManager, Tuple, Generator, Dict, cast, Iterable
from mongoengine import StringField, DateTimeField, DictField, BooleanField, IntField
import logging
from structlog import wrap_logger
from contextlib import contextmanager
import os.path
from datetime import datetime
from pymongo import UpdateOne
import hashlib
from structlog.processors import StackInfoRenderer, format_exc_info, TimeStamper
import json

from nomad import utils, config, infrastructure, search, datamodel
from nomad.files import PathObject, UploadFiles, ExtractError, ArchiveBasedStagingUploadFiles, PublicUploadFiles, StagingUploadFiles
from nomad.processing.base import Proc, process, task, PENDING, SUCCESS, FAILURE
from nomad.parsing import parser_dict, match_parser, Backend
from nomad.normalizing import normalizers


def _pack_log_event(logger, method_name, event_dict):
    try:
        log_data = dict(event_dict)
        log_data.update(**{
            key: value
            for key, value in getattr(logger, '_context', {}).items()
            if key not in ['service', 'release', 'upload_id', 'calc_id', 'mainfile', 'process_status']})
        log_data.update(logger=logger.name)

        return log_data
    except Exception:
        # raising an exception would cause an indefinite loop
        return event_dict


_log_processors = [
    StackInfoRenderer(),
    _pack_log_event,
    format_exc_info,
    TimeStamper(fmt="%Y-%m-%d %H:%M.%S", utc=False)]


class Calc(Proc):
    '''
    Instances of this class represent calculations. This class manages the elastic
    search index entry, files, and archive for the respective calculation.

    It also contains the calculations processing and its state.

    The attribute list, does not include the various metadata properties generated
    while parsing, including ``code_name``, ``code_version``, etc.

    Attributes:
        calc_id: the calc_id of this calc
        parser: the name of the parser used to process this calc
        upload_id: the id of the upload used to create this calculation
        mainfile: the mainfile (including path in upload) that was used to create this calc

        metadata: the metadata record wit calc and user metadata, see :class:`datamodel.EntryMetadata`
    '''
    calc_id = StringField(primary_key=True)
    upload_id = StringField()
    mainfile = StringField()
    parser = StringField()

    metadata = DictField()

    meta: Any = {
        'indexes': [
            'upload_id',
            ('upload_id', 'mainfile'),
            ('upload_id', 'parser'),
            ('upload_id', 'tasks_status'),
            ('upload_id', 'process_status'),
            ('upload_id', 'metadata.nomad_version'),
            'metadata.published',
            'metadata.datasets'
            'metadata.pid'
        ]
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._parser_backend: Backend = None
        self._upload: Upload = None
        self._upload_files: ArchiveBasedStagingUploadFiles = None
        self._calc_proc_logwriter = None
        self._calc_proc_logwriter_ctx: ContextManager = None

    @classmethod
    def from_entry_metadata(cls, entry_metadata):
        calc = Calc.create(
            calc_id=entry_metadata.calc_id,
            upload_id=entry_metadata.upload_id,
            mainfile=entry_metadata.mainfile,
            metadata=entry_metadata.m_to_dict(include_defaults=True))

        return calc

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
            self._upload.worker_hostname = self.worker_hostname
        return self._upload

    @property
    def upload_files(self) -> ArchiveBasedStagingUploadFiles:
        if not self._upload_files:
            self._upload_files = ArchiveBasedStagingUploadFiles(
                self.upload_id, is_authorized=lambda: True, upload_path=self.upload.upload_path)
        return self._upload_files

    def get_logger(self, **kwargs):
        '''
        Returns a wrapped logger that additionally saves all entries to the calculation
        processing log in the archive.
        '''
        logger = super().get_logger()
        logger = logger.bind(
            upload_id=self.upload_id, mainfile=self.mainfile, calc_id=self.calc_id, **kwargs)

        if self._calc_proc_logwriter_ctx is None:
            try:
                self._calc_proc_logwriter_ctx = self.upload_files.archive_log_file(self.calc_id, 'wt')
                self._calc_proc_logwriter = self._calc_proc_logwriter_ctx.__enter__()  # pylint: disable=E1101
            except KeyError:
                # cannot open log file
                pass

        if self._calc_proc_logwriter_ctx is None:
            return logger
        else:
            def save_to_calc_log(logger, method_name, event_dict):
                if self._calc_proc_logwriter is not None:
                    try:
                        dump_dict = dict(event_dict)
                        dump_dict.update(level=method_name.upper())
                        json.dump(dump_dict, self._calc_proc_logwriter, sort_keys=True)
                        self._calc_proc_logwriter.write('\n')

                    except Exception:
                        # Exceptions here will cause indefinite loop
                        pass

                return event_dict

            return wrap_logger(logger, processors=_log_processors + [save_to_calc_log])

    @process
    def re_process_calc(self):
        '''
        Processes a calculation again. This means there is already metadata and
        instead of creating it initially, we are just updating the existing
        records.
        '''
        parser = match_parser(self.upload_files.raw_file_object(self.mainfile).os_path, strict=False)

        if parser is None and not config.reprocess_unmatched:
            # Remove the logsfile and set a fake logwriter to avoid creating a log file,
            # because we will probably remove this calc and don't want to have ghost logfiles.
            if self._calc_proc_logwriter_ctx is not None:
                self._calc_proc_logwriter_ctx.__exit__(None, None, None)
                self.upload_files.archive_log_file_object(self.calc_id).delete()

            self._calc_proc_logwriter_ctx = open('/dev/null', 'wt')
            self._calc_proc_logwriter = self._calc_proc_logwriter_ctx.__enter__()
            self.get_logger().error(
                'no parser matches during re-process, will not re-process this calc')

            self.errors = ['no parser matches during re-process, will not re-process this calc']

            # mock the steps of actual processing
            self._continue_with('parsing')
            self._continue_with('normalizing')
            self._continue_with('archiving')
            self._complete()
            return

        logger = self.get_logger()
        if parser is None:
            self.get_logger().error('no parser matches during re-process, use the old parser')
            self.errors = ['no matching parser found during re-processing']
        if self.parser != parser.name:
            self.parser = parser.name
            logger.info(
                'different parser matches during re-process, use new parser',
                parser=parser.name)

        try:
            entry_metadata = datamodel.EntryMetadata.m_from_dict(self.metadata)
            entry_metadata.upload_id = self.upload_id
            entry_metadata.calc_id = self.calc_id
            entry_metadata.calc_hash = self.upload_files.calc_hash(self.mainfile)
            entry_metadata.mainfile = self.mainfile
            entry_metadata.nomad_version = config.version
            entry_metadata.nomad_commit = config.commit
            entry_metadata.last_processing = datetime.utcnow()
            entry_metadata.files = self.upload_files.calc_files(self.mainfile)
            self.metadata = entry_metadata.m_to_dict(include_defaults=True)

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

    @process
    def process_calc(self):
        '''
        Processes a new calculation that has no prior records in the mongo, elastic,
        or filesystem storage. It will create an initial set of (user) metadata.
        '''
        logger = self.get_logger()
        if self.upload is None:
            logger.error('calculation upload does not exist')

        try:
            # save preliminary minimum calc metadata in case processing fails
            # successful processing will replace it with the actual metadata
            calc_metadata = datamodel.EntryMetadata()
            calc_metadata.domain = parser_dict[self.parser].domain
            calc_metadata.upload_id = self.upload_id
            calc_metadata.calc_id = self.calc_id
            calc_metadata.calc_hash = self.upload_files.calc_hash(self.mainfile)
            calc_metadata.mainfile = self.mainfile
            calc_metadata.calc_hash = self.upload_files.calc_hash(self.mainfile)
            calc_metadata.nomad_version = config.version
            calc_metadata.nomad_commit = config.commit
            calc_metadata.last_processing = datetime.utcnow()
            calc_metadata.files = self.upload_files.calc_files(self.mainfile)
            calc_metadata.uploader = self.upload.user_id
            calc_metadata.upload_time = self.upload.upload_time
            calc_metadata.upload_name = self.upload.name
            self.metadata = calc_metadata.m_to_dict(include_defaults=True)  # TODO use embedded doc?

            if len(calc_metadata.files) >= config.auxfile_cutoff:
                self.warning(
                    'This calc has many aux files in its directory. '
                    'Have you placed many calculations in the same directory?')

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

    def fail(self, *errors, log_level=logging.ERROR, **kwargs):
        # in case of failure, index a minimum set of metadata and mark
        # processing failure
        try:
            entry_metadata = datamodel.EntryMetadata.m_from_dict(self.metadata)
            if self.parser is not None:
                parser = parser_dict[self.parser]
                if hasattr(parser, 'code_name'):
                    entry_metadata.code_name = parser.code_name

            entry_metadata.processed = False
            self.metadata = entry_metadata.m_to_dict(include_defaults=True)
            entry_metadata.a_elastic.index()
        except Exception as e:
            self.get_logger().error('could not index after processing failure', exc_info=e)

        super().fail(*errors, log_level=log_level, **kwargs)

    def on_process_complete(self, process_name):
        # the save might be necessary to correctly read the join condition from the db
        self.save()
        # in case of error, the process_name might be unknown
        if process_name == 'process_calc' or process_name == 're_process_calc' or process_name is None:
            self.upload.reload()
            self.upload.check_join()

    @task
    def parsing(self):
        ''' The *task* that encapsulates all parsing related actions. '''
        context = dict(parser=self.parser, step=self.parser)
        logger = self.get_logger(**context)
        parser = parser_dict[self.parser]
        self.metadata['parser_name'] = self.parser

        with utils.timer(logger, 'parser executed', input_size=self.mainfile_file.size):
            try:
                self._parser_backend = parser.run(
                    self.upload_files.raw_file_object(self.mainfile).os_path, logger=logger)
            except Exception as e:
                self.fail('parser failed with exception', exc_info=e, error=str(e), **context)
                return
            except SystemExit:
                self.fail('parser raised system exit', error='system exit', **context)
                return

        # add the non code specific calc metadata to the backend
        # all other quantities have been determined by parsers/normalizers
        self._parser_backend.openNonOverlappingSection('section_entry_info')
        self._parser_backend.addValue('upload_id', self.upload_id)
        self._parser_backend.addValue('calc_id', self.calc_id)
        self._parser_backend.addValue('calc_hash', self.metadata['calc_hash'])
        self._parser_backend.addValue('mainfile', self.mainfile)
        self._parser_backend.addValue('parser_name', self.parser)
        filepaths = self.metadata['files']
        self._parser_backend.addValue('number_of_files', len(filepaths))
        self._parser_backend.addValue('filepaths', filepaths)
        uploader = self.upload.uploader
        self._parser_backend.addValue(
            'entry_uploader_name', '%s, %s' % (uploader.first_name, uploader.last_name))
        self._parser_backend.addValue(
            'entry_uploader_id', str(uploader.user_id))
        self._parser_backend.addValue('entry_upload_time', int(self.upload.upload_time.timestamp()))
        self._parser_backend.closeNonOverlappingSection('section_entry_info')

        self.add_processor_info(self.parser)

        if self._parser_backend.status[0] != 'ParseSuccess':
            error = self._parser_backend.status[1]
            self.fail('parser failed', error=error, **context)

    @contextmanager
    def use_parser_backend(self, processor_name):
        self._parser_backend.reset_status()
        yield self._parser_backend
        self.add_processor_info(processor_name)

    def add_processor_info(self, processor_name: str) -> None:
        self._parser_backend.openContext('/section_entry_info/0')
        self._parser_backend.openNonOverlappingSection('section_archive_processing_info')
        self._parser_backend.addValue('archive_processor_name', processor_name)

        if self._parser_backend.status[0] == 'ParseSuccess':
            warnings = getattr(self._parser_backend, '_warnings', [])
            if len(warnings) > 0:
                self._parser_backend.addValue('archive_processor_status', 'WithWarnings')
                self._parser_backend.addValue('number_of_archive_processor_warnings', len(warnings))
                self._parser_backend.addArrayValues('archive_processor_warnings', [str(warning) for warning in warnings])
            else:
                self._parser_backend.addValue('archive_processor_status', 'Success')
        else:
            errors = self._parser_backend.status[1]
            self._parser_backend.addValue('archive_processor_error', str(errors))

        self._parser_backend.closeNonOverlappingSection('section_archive_processing_info')
        self._parser_backend.closeContext('/section_entry_info/0')

    @task
    def normalizing(self):
        ''' The *task* that encapsulates all normalizing related actions. '''
        for normalizer in normalizers:
            if normalizer.domain != parser_dict[self.parser].domain:
                continue

            normalizer_name = normalizer.__name__
            context = dict(normalizer=normalizer_name, step=normalizer_name)
            logger = self.get_logger(**context)

            with utils.timer(
                    logger, 'normalizer executed', input_size=self.mainfile_file.size):
                with self.use_parser_backend(normalizer_name) as backend:
                    try:
                        normalizer(backend).normalize(logger=logger)
                    except Exception as e:
                        self._parser_backend.finishedParsingSession('ParseFailure', [str(e)])
                        self.fail(
                            'normalizer failed with exception', exc_info=e, error=str(e), **context)
                        break
                    else:
                        if self._parser_backend.status[0] != 'ParseSuccess':
                            error = self._parser_backend.status[1]
                            self.fail('normalizer failed', error=error, **context)
                            break
                        else:
                            logger.debug(
                                'completed normalizer successfully', normalizer=normalizer_name)

    @task
    def archiving(self):
        ''' The *task* that encapsulates all archival related actions. '''
        logger = self.get_logger()

        entry_metadata = datamodel.EntryMetadata.m_from_dict(self.metadata)
        entry_metadata.apply_domain_metadata(self._parser_backend)
        entry_metadata.processed = True

        # persist the calc metadata
        with utils.timer(logger, 'saved calc metadata', step='metadata'):
            self.metadata = entry_metadata.m_to_dict(include_defaults=True)

        # index in search
        with utils.timer(logger, 'indexed', step='index'):
            entry_metadata.a_elastic.index()

        # persist the archive
        with utils.timer(
                logger, 'archived', step='archive',
                input_size=self.mainfile_file.size) as log_data:

            with self.upload_files.archive_file(self.calc_id, 'wt') as out:
                json.dump(self._parser_backend.resource.m_to_dict(
                    lambda section: section.m_def.name in datamodel.root_sections), out, indent=2)

            log_data.update(archive_size=self.upload_files.archive_file_object(self.calc_id).size)

        # close loghandler
        if self._calc_proc_logwriter is not None:
            with utils.timer(
                    logger, 'archived log', step='logs',
                    input_size=self.mainfile_file.size) as log_data:
                self._calc_proc_logwriter_ctx.__exit__(None, None, None)  # pylint: disable=E1101
                self._calc_proc_logwriter = None

                log_data.update(log_size=self.upload_files.archive_log_file_object(self.calc_id).size)

    def __str__(self):
        return 'calc %s calc_id=%s upload_id%s' % (super().__str__(), self.calc_id, self.upload_id)


class Upload(Proc):
    '''
    Represents uploads in the databases. Provides persistence access to the files storage,
    and processing state.

    Attributes:
        name: optional user provided upload name
        upload_path: the path were the uploaded files was stored
        temporary: True if the uploaded file should be removed after extraction
        upload_id: the upload id generated by the database
        upload_time: the timestamp when the system realised the upload
        user_id: the id of the user that created this upload
        published: Boolean that indicates the publish status
        publish_time: Date when the upload was initially published
        last_update: Date of the last publishing/re-processing
        joined: Boolean indicates if the running processing has joined (:func:`check_join`)
    '''
    id_field = 'upload_id'

    upload_id = StringField(primary_key=True)
    upload_path = StringField(default=None)
    temporary = BooleanField(default=False)
    embargo_length = IntField(default=36)

    name = StringField(default=None)
    upload_time = DateTimeField()
    user_id = StringField(required=True)
    published = BooleanField(default=False)
    publish_time = DateTimeField()
    last_update = DateTimeField()

    joined = BooleanField(default=False)

    meta: Any = {
        'indexes': [
            'user_id', 'tasks_status', 'process_status', 'published', 'upload_time'
        ]
    }

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._upload_files: ArchiveBasedStagingUploadFiles = None

    @property
    def metadata(self) -> dict:
        '''
        Getter, setter for user metadata. Metadata is pickled to and from the public
        bucket to allow sharing among all processes. Usually uploads do not have (much)
        user defined metadata, but users provide all metadata per upload as part of
        the publish process. This will change, when we introduce editing functionality
        and metadata will be provided through different means.
        '''
        try:
            upload_files = PublicUploadFiles(self.upload_id, is_authorized=lambda: True)
        except KeyError:
            return None
        return upload_files.user_metadata

    @metadata.setter
    def metadata(self, data: dict) -> None:
        upload_files = PublicUploadFiles(self.upload_id, is_authorized=lambda: True, create=True)
        upload_files.user_metadata = data

    @classmethod
    def get(cls, id: str, include_published: bool = True) -> 'Upload':
        return cls.get_by_id(id, 'upload_id')

    @classmethod
    def user_uploads(cls, user: datamodel.User, **kwargs) -> List['Upload']:
        ''' Returns all uploads for the given user. Kwargs are passed to mongo query. '''
        return cls.objects(user_id=str(user.user_id), **kwargs)

    @property
    def uploader(self):
        return datamodel.User.get(self.user_id)

    def get_logger(self, **kwargs):
        logger = super().get_logger()
        user = self.uploader
        user_name = '%s %s' % (user.first_name, user.last_name)
        # We are not using 'user_id' because logstash (?) will filter these entries ?!
        logger = logger.bind(
            upload_id=self.upload_id, upload_name=self.name, user_name=user_name,
            user=user.user_id, **kwargs)
        return logger

    @classmethod
    def create(cls, **kwargs) -> 'Upload':
        '''
        Creates a new upload for the given user, a user given name is optional.
        It will populate the record with a signed url and pending :class:`UploadProc`.
        The upload will be already saved to the database.

        Arguments:
            user: The user that created the upload.
        '''
        # use kwargs to keep compatibility with super method
        user: datamodel.User = kwargs['user']
        del(kwargs['user'])

        if 'upload_id' not in kwargs:
            kwargs.update(upload_id=utils.create_uuid())
        kwargs.update(user_id=user.user_id)
        self = super().create(**kwargs)

        self._continue_with('uploading')

        return self

    def delete(self):
        ''' Deletes this upload process state entry and its calcs. '''
        Calc.objects(upload_id=self.upload_id).delete()
        super().delete()

    def delete_upload_local(self):
        '''
        Deletes the upload, including its processing state and
        staging files. Local version without celery processing.
        '''
        logger = self.get_logger()

        with utils.lnr(logger, 'staged upload delete failed'):
            with utils.timer(
                    logger, 'upload deleted from index', step='index',
                    upload_size=self.upload_files.size):
                search.delete_upload(self.upload_id)

            with utils.timer(
                    logger, 'staged upload deleted', step='files',
                    upload_size=self.upload_files.size):
                self.upload_files.delete()

            self.delete()

    @process
    def delete_upload(self):
        '''
        Deletes of the upload, including its processing state and
        staging files. This starts the celery process of deleting the upload.
        '''
        self.delete_upload_local()

        return True  # do not save the process status on the delete upload

    @process
    def publish_upload(self):
        '''
        Moves the upload out of staging to the public area. It will
        pack the staging upload files in to public upload files.
        '''
        assert self.processed_calcs > 0

        logger = self.get_logger()
        logger.info('started to publish')

        with utils.lnr(logger, 'publish failed'):
            calcs = self.entries_metadata(self.metadata)

            with utils.timer(
                    logger, 'upload metadata updated', step='metadata',
                    upload_size=self.upload_files.size):

                def create_update(calc):
                    calc.published = True
                    calc.with_embargo = calc.with_embargo if calc.with_embargo is not None else False
                    return UpdateOne(
                        {'_id': calc.calc_id},
                        {'$set': {'metadata': calc.m_to_dict(include_defaults=True)}})

                Calc._get_collection().bulk_write([create_update(calc) for calc in calcs])

            if isinstance(self.upload_files, StagingUploadFiles):
                with utils.timer(
                        logger, 'staged upload files packed', step='pack',
                        upload_size=self.upload_files.size):
                    self.upload_files.pack(calcs)

            with utils.timer(
                    logger, 'index updated', step='index',
                    upload_size=self.upload_files.size):
                search.publish(calcs)

            if isinstance(self.upload_files, StagingUploadFiles):
                with utils.timer(
                        logger, 'staged upload deleted', step='delete staged',
                        upload_size=self.upload_files.size):
                    self.upload_files.delete()
                    self.published = True
                    self.publish_time = datetime.utcnow()
                    self.last_update = datetime.utcnow()
                    self.save()
            else:
                self.last_update = datetime.utcnow()
                self.save()

    @process
    def re_process_upload(self):
        '''
        A *process* that performs the re-processing of a earlier processed
        upload.

        Runs the distributed process of fully reparsing/renormalizing an existing and
        already published upload. Will renew the archive part of the upload and update
        mongo and elastic search entries.

        TODO this implementation does not do any re-matching. This will be more complex
        due to handling of new or missing matches.
        '''
        assert self.published

        logger = self.get_logger()
        logger.info('started to re-process')

        # mock the steps of actual processing
        self._continue_with('uploading')

        # extract the published raw files into a staging upload files instance
        self._continue_with('extracting')
        public_upload_files = cast(PublicUploadFiles, self.upload_files)
        staging_upload_files = public_upload_files.to_staging_upload_files(create=True)

        self._continue_with('parse_all')
        try:
            # check if a calc is already/still processing
            processing = Calc.objects(
                upload_id=self.upload_id,
                **Calc.process_running_mongoengine_query()).count()

            if processing > 0:
                logger.warn(
                    'processes are still/already running on calc, they will be resetted',
                    count=processing)

            # reset all calcs
            Calc._get_collection().update_many(
                dict(upload_id=self.upload_id),
                {'$set': Calc.reset_pymongo_update(worker_hostname=self.worker_hostname)})

            # process call calcs
            Calc.process_all(Calc.re_process_calc, dict(upload_id=self.upload_id), exclude=['metadata'])

            logger.info('completed to trigger re-process of all calcs')
        except Exception as e:
            # try to remove the staging copy in failure case
            logger.error('failed to trigger re-process of all calcs', exc_info=e)

            if staging_upload_files is not None and staging_upload_files.exists():
                staging_upload_files.delete()

            raise e

        # the packing and removing of the staging upload files, will be trigged by
        # the 'cleanup' task after processing all calcs

    @process
    def re_pack(self):
        ''' A *process* that repacks the raw and archive data based on the current embargo data. '''
        assert self.published

        # mock the steps of actual processing
        self._continue_with('uploading')
        self._continue_with('extracting')
        self._continue_with('parse_all')
        self._continue_with('cleanup')

        self.upload_files.re_pack(self.entries_metadata())
        self.joined = True
        self._complete()

    @process
    def process_upload(self):
        ''' A *process* that performs the initial upload processing. '''
        self.extracting()
        self.parse_all()

    @task
    def uploading(self):
        ''' A no-op *task* as a stand-in for receiving upload data. '''
        pass

    @property
    def upload_files(self) -> UploadFiles:
        upload_files_class = ArchiveBasedStagingUploadFiles if not self.published else PublicUploadFiles
        kwargs = dict(upload_path=self.upload_path) if not self.published else {}

        if not self._upload_files or not isinstance(self._upload_files, upload_files_class):
            self._upload_files = upload_files_class(
                self.upload_id, is_authorized=lambda: True, **kwargs)

        return self._upload_files

    @property
    def staging_upload_files(self) -> ArchiveBasedStagingUploadFiles:
        assert not self.published
        return cast(ArchiveBasedStagingUploadFiles, self.upload_files)

    @task
    def extracting(self):
        '''
        The *task* performed before the actual parsing/normalizing: extracting
        the uploaded files.
        '''
        # extract the uploaded file
        self._upload_files = ArchiveBasedStagingUploadFiles(
            upload_id=self.upload_id, is_authorized=lambda: True, create=True,
            upload_path=self.upload_path)

        logger = self.get_logger()
        try:
            with utils.timer(
                    logger, 'upload extracted', step='extracting',
                    upload_size=self.upload_files.size):
                self.upload_files.extract()

            if self.temporary:
                os.remove(self.upload_path)
                self.upload_path = None

        except KeyError:
            self.fail('processing requested for non existing upload', log_level=logging.ERROR)
            return
        except ExtractError:
            self.fail('bad .zip/.tar file', log_level=logging.INFO)
            return

    def _preprocess_files(self, path):
        '''
        Some files need preprocessing. Currently we need to add a stripped POTCAR version
        and always restrict/embargo the original.
        '''
        if os.path.basename(path).startswith('POTCAR'):
            # create checksum
            hash = hashlib.sha224()
            with open(self.staging_upload_files.raw_file_object(path).os_path, 'rb') as orig_f:
                for line in orig_f.readlines():
                    hash.update(line)

            checksum = hash.hexdigest()

            # created stripped POTCAR
            stripped_path = path + '.stripped'
            with open(self.staging_upload_files.raw_file_object(stripped_path).os_path, 'wt') as stripped_f:
                stripped_f.write('Stripped POTCAR file. Checksum of original file (sha224): %s\n' % checksum)
            os.system(
                '''
                    awk < %s >> %s '
                    BEGIN { dump=1 }
                    /End of Dataset/ { dump=1 }
                    dump==1 { print }
                    /END of PSCTR/ { dump=0 }'
                ''' % (
                    self.staging_upload_files.raw_file_object(path).os_path,
                    self.staging_upload_files.raw_file_object(stripped_path).os_path))

    def match_mainfiles(self) -> Generator[Tuple[str, object], None, None]:
        '''
        Generator function that matches all files in the upload to all parsers to
        determine the upload's mainfiles.

        Returns:
            Tuples of mainfile, filename, and parsers
        '''
        directories_with_match: Dict[str, str] = dict()
        upload_files = self.staging_upload_files
        for filename in upload_files.raw_file_manifest():
            self._preprocess_files(filename)
            try:
                parser = match_parser(upload_files.raw_file_object(filename).os_path)
                if parser is not None:
                    directory = os.path.dirname(filename)
                    if directory in directories_with_match:
                        # TODO this might give us the chance to store directory based relationship
                        # between calcs for the future?
                        pass
                    else:
                        directories_with_match[directory] = filename

                    yield filename, parser
            except Exception as e:
                self.get_logger().error(
                    'exception while matching pot. mainfile',
                    mainfile=filename, exc_info=e)

    @task
    def parse_all(self):
        '''
        The *task* used to identify mainfile/parser combinations among the upload's files, creates
        respective :class:`Calc` instances, and triggers their processing.
        '''
        logger = self.get_logger()

        with utils.timer(
                logger, 'upload extracted', step='matching',
                upload_size=self.upload_files.size):
            for filename, parser in self.match_mainfiles():
                calc = Calc.create(
                    calc_id=self.upload_files.calc_id(filename),
                    mainfile=filename, parser=parser.name,
                    worker_hostname=self.worker_hostname,
                    upload_id=self.upload_id)

                calc.process_calc()

    def on_process_complete(self, process_name):
        if process_name == 'process_upload' or process_name == 're_process_upload':
            self.check_join()

    def check_join(self):
        '''
        Performs an evaluation of the join condition and triggers the :func:`cleanup`
        task if necessary. The join condition allows to run the ``cleanup`` after
        all calculations have been processed. The upload processing stops after all
        calculation processings have been triggered (:func:`parse_all` or
        :func:`re_process_upload`). The cleanup task is then run within the last
        calculation process (the one that triggered the join by calling this method).
        '''
        total_calcs = self.total_calcs
        processed_calcs = self.processed_calcs

        self.get_logger().debug('check join', processed_calcs=processed_calcs, total_calcs=total_calcs)
        # check if process is not running anymore, i.e. not still spawining new processes to join
        # check the join condition, i.e. all calcs have been processed
        if not self.process_running and processed_calcs >= total_calcs:
            # this can easily be called multiple times, e.g. upload finished after all calcs finished
            modified_upload = self._get_collection().find_one_and_update(
                {'_id': self.upload_id, 'joined': {'$ne': True}},
                {'$set': {'joined': True}})
            if modified_upload is not None:
                self.get_logger().debug('join')
                self.cleanup()
            else:
                # the join was already done due to a prior call
                pass

    def reset(self):
        self.joined = False
        super().reset()

    @classmethod
    def reset_pymongo_update(cls, worker_hostname: str = None):
        update = super().reset_pymongo_update()
        update.update(joined=False)
        return update

    def _cleanup_after_processing(self):
        # send email about process finish
        user = self.uploader
        name = '%s %s' % (user.first_name, user.last_name)
        message = '\n'.join([
            'Dear %s,' % name,
            '',
            'your data %suploaded at %s has completed processing.' % (
                '"%s" ' % self.name if self.name else '', self.upload_time.isoformat()),  # pylint: disable=no-member
            'You can review your data on your upload page: %s' % config.gui_url(),
            '',
            'If you encouter any issues with your upload, please let us know and replay to this email.',
            '',
            'The nomad team'
        ])
        try:
            infrastructure.send_mail(
                name=name, email=user.email, message=message, subject='Processing completed')
        except Exception as e:
            # probably due to email configuration problems
            # don't fail or present this error to clients
            self.logger.error('could not send after processing email', exc_info=e)

    def _cleanup_after_re_processing(self):
        logger = self.get_logger()
        logger.info('started to repack re-processed upload')

        staging_upload_files = self.upload_files.to_staging_upload_files()

        with utils.timer(
                logger, 'reprocessed staged upload packed', step='delete staged',
                upload_size=self.upload_files.size):

            staging_upload_files.pack(self.entries_metadata(), skip_raw=True)

        with utils.timer(
                logger, 'reprocessed staged upload deleted', step='delete staged',
                upload_size=self.upload_files.size):

            staging_upload_files.delete()
            self.last_update = datetime.utcnow()
            self.save()

    @task
    def cleanup(self):
        '''
        The *task* that "cleans" the processing, i.e. removed obsolete files and performs
        pending archival operations. Depends on the type of processing.
        '''
        search.refresh()

        if self.current_process == 're_process_upload':
            self._cleanup_after_re_processing()
        else:
            self._cleanup_after_processing()

    def get_calc(self, calc_id) -> Calc:
        ''' Returns the upload calc with the given id or ``None``. '''
        return Calc.objects(upload_id=self.upload_id, calc_id=calc_id).first()

    @property
    def processed_calcs(self):
        '''
        The number of successfully or not successfully processed calculations. I.e.
        calculations that have finished processing.
        '''
        return Calc.objects(upload_id=self.upload_id, tasks_status__in=[SUCCESS, FAILURE]).count()

    @property
    def total_calcs(self):
        ''' The number of all calculations. '''
        return Calc.objects(upload_id=self.upload_id).count()

    @property
    def failed_calcs(self):
        ''' The number of calculations with failed processing. '''
        return Calc.objects(upload_id=self.upload_id, tasks_status=FAILURE).count()

    @property
    def pending_calcs(self) -> int:
        ''' The number of calculations with pending processing. '''
        return Calc.objects(upload_id=self.upload_id, tasks_status=PENDING).count()

    def all_calcs(self, start, end, order_by=None):
        '''
        Returns all calculations, paginated and ordered.

        Arguments:
            start: the start index of the requested page
            end: the end index of the requested page
            order_by: the property to order by
        '''
        query = Calc.objects(upload_id=self.upload_id)[start:end]
        return query.order_by(order_by) if order_by is not None else query

    @property
    def outdated_calcs(self):
        ''' All successfully processed and outdated calculations. '''
        return Calc.objects(
            upload_id=self.upload_id, tasks_status=SUCCESS,
            metadata__nomad_version__ne=config.version)

    @property
    def calcs(self):
        ''' All successfully processed calculations. '''
        return Calc.objects(upload_id=self.upload_id, tasks_status=SUCCESS)

    def entries_metadata(self, user_metadata: dict = None) -> Iterable[datamodel.EntryMetadata]:
        '''
        This is the :py:mod:`nomad.datamodel` transformation method to transform
        processing uploads into datamodel uploads. It will also implicitely transform
        all calculations of this upload.

        Arguments:
            user_metadata: A dict of user metadata that is applied to the resulting
                datamodel data and the respective calculations.
        '''
        # prepare user metadata per upload and per calc
        if user_metadata is not None:
            entries_metadata_dict: Dict[str, Any] = dict()
            upload_metadata: Dict[str, Any] = dict()

            upload_metadata.update(user_metadata)
            if 'calculations' in upload_metadata:
                del(upload_metadata['calculations'])

            for calc in user_metadata.get('calculations', []):  # pylint: disable=no-member
                entries_metadata_dict[calc['mainfile']] = calc

            user_upload_time = upload_metadata.get('upload_time', None)
            user_upload_name = upload_metadata.get('upload_name', None)

            def get_metadata(calc: Calc):
                entry_metadata = datamodel.EntryMetadata.m_from_dict(calc.metadata)
                entry_user_metadata = dict(upload_metadata)
                entry_user_metadata.pop('embargo_length', None)  # this is for uploads only
                entry_user_metadata.update(entries_metadata_dict.get(calc.mainfile, {}))
                entry_metadata.apply_user_metadata(entry_user_metadata)
                if entry_metadata.upload_time is None:
                    entry_metadata.upload_time = self.upload_time if user_upload_time is None else user_upload_time
                if entry_metadata.upload_name is None:
                    entry_metadata.upload_name = self.name if user_upload_name is None else user_upload_name

                return entry_metadata
        else:
            user_upload_time = None

            def get_metadata(calc: Calc):
                entry_metadata = datamodel.EntryMetadata.m_from_dict(calc.metadata)
                entry_metadata.upload_time = self.upload_time
                entry_metadata.upload_name = self.name

                return entry_metadata

        return [get_metadata(calc) for calc in Calc.objects(upload_id=self.upload_id)]

    def compress_and_set_metadata(self, metadata: Dict[str, Any]) -> None:
        '''
        Stores the given user metadata in the upload document. This is the metadata
        adhering to the API model (``UploadMetaData``). Most quantities can be stored
        for the upload and for each calculation. This method will try to move same values
        from the calculation to the upload to "compress" the data.
        '''
        self.embargo_length = min(metadata.get('embargo_length', 36), 36)

        compressed = {
            key: value for key, value in metadata.items() if key != 'calculations'}
        calculations: List[Dict[str, Any]] = []
        compressed['calculations'] = calculations

        for calc in metadata.get('calculations', []):
            compressed_calc: Dict[str, Any] = {}
            calculations.append(compressed_calc)
            for key, value in calc.items():
                if key in ['pid', 'mainfile', 'external_id']:
                    # these quantities are explicitly calc specific and have to stay with
                    # the calc
                    compressed_calc[key] = value
                else:
                    if key not in compressed:
                        compressed[key] = value
                    elif compressed[key].__repr__ != value.__repr__:
                        compressed_calc[key] = value
                    else:
                        compressed[key] = value

        self.metadata = compressed

    def __str__(self):
        return 'upload %s upload_id%s' % (super().__str__(), self.upload_id)
