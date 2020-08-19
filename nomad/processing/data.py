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

from typing import cast, List, Any, Tuple, Iterator, Dict, cast, Iterable
from mongoengine import StringField, DateTimeField, DictField, BooleanField, IntField
import logging
from structlog import wrap_logger
from contextlib import contextmanager
import os.path
from datetime import datetime
from pymongo import UpdateOne
import hashlib
from structlog.processors import StackInfoRenderer, format_exc_info, TimeStamper

from nomad import utils, config, infrastructure, search, datamodel
from nomad.files import PathObject, UploadFiles, ExtractError, ArchiveBasedStagingUploadFiles, PublicUploadFiles, StagingUploadFiles
from nomad.processing.base import Proc, process, task, PENDING, SUCCESS, FAILURE
from nomad.parsing import Backend
from nomad.parsing.parsers import parser_dict, match_parser
from nomad.normalizing import normalizers
from nomad.datamodel import EntryArchive
from nomad.archive import query_archive
from nomad.datamodel.encyclopedia import (
    EncyclopediaMetadata,
)
import phonopyparser.metainfo


section_metadata = datamodel.EntryArchive.section_metadata.name
section_workflow = datamodel.EntryArchive.section_workflow.name


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
            'parser',
            ('upload_id', 'mainfile'),
            ('upload_id', 'parser'),
            ('upload_id', 'tasks_status'),
            ('upload_id', 'process_status'),
            ('upload_id', 'metadata.nomad_version'),
            'metadata.processed',
            'metadata.last_processing',
            'metadata.published',
            'metadata.datasets',
            'metadata.pid'
        ]
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._parser_backend: Backend = None
        self._upload: Upload = None
        self._upload_files: ArchiveBasedStagingUploadFiles = None
        self._calc_proc_logs: List[Any] = None

        self._entry_metadata = None

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

    def apply_entry_metadata(self, entry_metadata: datamodel.EntryMetadata):
        self.metadata = entry_metadata.m_to_dict(
            include_defaults=True,
            categories=[datamodel.MongoMetadata])  # TODO use embedded doc?

    def create_metadata(self) -> datamodel.EntryMetadata:
        '''
        Returns a :class:`nomad.datamodel.EntryMetadata` with values from this
        processing object, not necessarely the user metadata nor the metadata from
        the archive.
        '''
        entry_metadata = datamodel.EntryMetadata()
        if self.parser is not None:
            entry_metadata.domain = parser_dict[self.parser].domain
        entry_metadata.upload_id = self.upload_id
        entry_metadata.calc_id = self.calc_id
        entry_metadata.mainfile = self.mainfile
        entry_metadata.nomad_version = config.meta.version
        entry_metadata.nomad_commit = config.meta.commit
        entry_metadata.uploader = self.upload.user_id
        entry_metadata.upload_time = self.upload.upload_time
        entry_metadata.upload_name = self.upload.name

        return entry_metadata

    def entry_metadata(self, upload_files: UploadFiles) -> datamodel.EntryMetadata:
        '''
        Returns a complete set of :class:`nomad.datamodel.EntryMetadata` including
        the user metadata and metadata from the archive.

        Arguments:
            upload_files:
                The :class:`nomad.files.UploadFiles` instance to read the archive from.
            cache:
                A boolean that indicates if the archive file should be left unclosed,
                e.g. if this method is called for many entries of the same upload.
        '''
        archive = upload_files.read_archive(self.calc_id)
        try:
            # instead of loading the whole archive, it should be enough to load the
            # parts that are referenced by section_metadata/EntryMetadata
            # TODO somehow it should determine which root setions too load from the metainfo
            # or configuration
            calc_archive = archive[self.calc_id]
            entry_archive_dict = {section_metadata: calc_archive[section_metadata].to_dict()}
            if section_workflow in calc_archive:
                entry_archive_dict[section_workflow] = calc_archive[section_workflow].to_dict()
            entry_metadata = datamodel.EntryArchive.m_from_dict(entry_archive_dict)[section_metadata]

        except KeyError:
            # Due hard processing failures, it might be possible that an entry might not
            # have an archive
            if self._entry_metadata is not None:
                entry_metadata = self._entry_metadata

            else:
                entry_metadata = self.create_metadata()

        entry_metadata.m_update_from_dict(self.metadata)

        return entry_metadata

    def user_metadata(self) -> datamodel.EntryMetadata:
        '''
        Returns a :class:`nomad.datamodel.EntryMetadata` with values from this
        processing object and the user metadata, not necessarely the metadata from
        the archive.
        '''
        entry_metadata = self.create_metadata()
        entry_metadata.m_update_from_dict(self.metadata)

        return entry_metadata

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

        if self._calc_proc_logs is None:
            self._calc_proc_logs = []

        def save_to_calc_log(logger, method_name, event_dict):
            try:
                # sanitize the event_dict, because all kinds of values might have been added
                dump_dict = {key: str(value) for key, value in event_dict.items()}
                dump_dict.update(level=method_name.upper())
                self._calc_proc_logs.append(dump_dict)

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
        logger = self.get_logger()

        if parser is None and not config.reprocess_unmatched:
            self.errors = ['no parser matches during re-process, will not re-process this calc']

            try:
                upload_files = PublicUploadFiles(self.upload_id, is_authorized=lambda: True)
                with upload_files.read_archive(self.calc_id) as archive:
                    self.upload_files.write_archive(self.calc_id, archive[self.calc_id].to_dict())

            except Exception as e:
                logger.error('could not copy archive for non matching, non reprocessed entry', exc_info=e)
                raise e

            # mock the steps of actual processing
            self._continue_with('parsing')
            self._continue_with('normalizing')
            self._continue_with('archiving')
            self._complete()
            return

        if parser is None:
            self.get_logger().warn('no parser matches during re-process, use the old parser')
            self.warnings = ['no matching parser found during re-processing']

        elif self.parser != parser.name:
            if parser_dict[self.parser].name == parser.name:
                # parser was just renamed
                self.parser = parser.name

            else:
                self.parser = parser.name
                logger.info(
                    'different parser matches during re-process, use new parser',
                    parser=parser.name)

        try:
            self._entry_metadata = self.user_metadata()
            self._entry_metadata.calc_hash = self.upload_files.calc_hash(self.mainfile)
            self._entry_metadata.nomad_version = config.meta.version
            self._entry_metadata.nomad_commit = config.meta.commit
            self._entry_metadata.last_processing = datetime.utcnow()
            self._entry_metadata.files = self.upload_files.calc_files(self.mainfile)

            self.parsing()
            self.normalizing()
            self.archiving()
        finally:
            # close loghandler that was not closed due to failures
            try:
                if self._parser_backend and self._parser_backend.resource:
                    self._parser_backend.resource.unload()
            except Exception as e:
                logger.error('could unload processing results', exc_info=e)

    def _setup_fallback_metadata(self):
        self._entry_metadata = self.create_metadata()
        self._entry_metadata.calc_hash = self.upload_files.calc_hash(self.mainfile)
        self._entry_metadata.last_processing = datetime.utcnow()
        self._entry_metadata.files = self.upload_files.calc_files(self.mainfile)

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
            self._setup_fallback_metadata()

            if len(self._entry_metadata.files) >= config.auxfile_cutoff:
                self.warning(
                    'This calc has many aux files in its directory. '
                    'Have you placed many calculations in the same directory?')

            self.parsing()
            self.normalizing()
            self.archiving()
        finally:
            # close loghandler that was not closed due to failures
            try:
                if self._parser_backend and self._parser_backend.resource:
                    self._parser_backend.resource.unload()
            except Exception as e:
                logger.error('could unload processing results', exc_info=e)

    def on_fail(self):
        # in case of failure, index a minimum set of metadata and mark
        # processing failure
        try:
            if self._entry_metadata is None:
                self._setup_fallback_metadata()

            self._entry_metadata.processed = False

            self.apply_entry_metadata(self._entry_metadata)
            if self._parser_backend and self._parser_backend.resource:
                backend = self._parser_backend
            else:
                backend = None
            self._entry_metadata.apply_domain_metadata(backend)

            self._entry_metadata.a_elastic.index()
        except Exception as e:
            self.get_logger().error(
                'could not index after processing failure', exc_info=e)

        try:
            self.write_archive(None)
        except Exception as e:
            self.get_logger().error(
                'could not write archive after processing failure', exc_info=e)

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
        self._entry_metadata.parser_name = self.parser

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

        if self._parser_backend.status[0] != 'ParseSuccess':
            error = self._parser_backend.status[1]
            self.fail('parser failed', error=error, **context)

    def process_phonon(self):
        """Function that is run for phonon calculation before cleanup.
        This task is run by the celery process that is calling the join for the
        upload.

        This function re-opens the Archive for this calculation to add method
        information from another referenced archive. Updates the method
        information in section_encyclopedia as well as the DFT domain metadata.
        """
        try:
            logger = self.get_logger(parser=self.parser, step=self.parser)

            # Open the archive of the phonon calculation.
            upload_files = StagingUploadFiles(self.upload_id, is_authorized=lambda: True)
            with upload_files.read_archive(self.calc_id) as archive:
                arch = query_archive(archive, {self.calc_id: self.calc_id})[self.calc_id]
                phonon_archive = EntryArchive.m_from_dict(arch)
            self._entry_metadata = phonon_archive.section_metadata
            self._calc_proc_logs = phonon_archive.processing_logs

            # Re-create a backend
            metainfo = phonopyparser.metainfo.m_env
            self._parser_backend = Backend(metainfo, logger=logger, domain="dft")
            self._parser_backend.entry_archive = phonon_archive

            # Read in the first referenced calculation. The reference is given as
            # an absolute path which needs to be converted into a path that is
            # relative to upload root.
            scc = self._parser_backend.entry_archive.section_run[0].section_single_configuration_calculation[0]
            relative_ref = scc.section_calculation_to_calculation_refs[0].calculation_to_calculation_external_url
            ref_id = upload_files.calc_id(relative_ref)
            with upload_files.read_archive(ref_id) as archive:
                arch = query_archive(archive, {ref_id: ref_id})[ref_id]
                ref_archive = EntryArchive.m_from_dict(arch)

            # Get encyclopedia method information directly from the referenced calculation.
            ref_enc_method = ref_archive.section_metadata.encyclopedia.method
            if ref_enc_method is None or len(ref_enc_method) == 0 or ref_enc_method.functional_type is None:
                raise ValueError("No method information available in referenced calculation.")
            self._parser_backend.entry_archive.section_metadata.encyclopedia.method = ref_enc_method

            # Overwrite old entry with new data. The metadata is updated with
            # new timestamp and method details taken from the referenced
            # archive.
            self._entry_metadata.last_processing = datetime.utcnow()
            self._entry_metadata.dft.xc_functional = ref_archive.section_metadata.dft.xc_functional
            self._entry_metadata.dft.basis_set = ref_archive.section_metadata.dft.basis_set
            self._entry_metadata.dft.update_group_hash()
            self._entry_metadata.encyclopedia.status = EncyclopediaMetadata.status.type.success
        except Exception as e:
            logger.error("Could not retrieve method information for phonon calculation.", exc_info=e)
            if self._entry_metadata is None:
                self._setup_fallback_metadata()
                self._entry_metadata.processed = False

            try:
                if self._entry_metadata.encyclopedia is None:
                    self._entry_metadata.encyclopedia = EncyclopediaMetadata()
                self._entry_metadata.encyclopedia.status = EncyclopediaMetadata.status.type.failure
            except Exception as e:
                logger.error("Could set encyclopedia status.", exc_info=e)

        finally:
            # persist the calc metadata
            with utils.timer(logger, 'saved calc metadata', step='metadata'):
                self.apply_entry_metadata(self._entry_metadata)

            # index in search
            with utils.timer(logger, 'indexed', step='index'):
                self._entry_metadata.a_elastic.index()

            # persist the archive
            with utils.timer(
                    logger, 'archived', step='archive',
                    input_size=self.mainfile_file.size) as log_data:

                archive_size = self.write_archive(self._parser_backend)
                log_data.update(archive_size=archive_size)

    @contextmanager
    def use_parser_backend(self, processor_name):
        self._parser_backend.reset_status()
        yield self._parser_backend

        if self._parser_backend.status[0] == 'ParseSuccess':
            warnings = getattr(self._parser_backend, '_warnings', [])

            if len(warnings) > 0:
                self.get_logger().warn(
                    'processor completed successful with warnings',
                    processor=processor_name, warnings=[str(warning) for warning in warnings])

            else:
                self.get_logger().info(
                    'processor completed successful',
                    processor=processor_name)

        else:
            errors = self._parser_backend.status[1]
            self.get_logger().error(
                'processor completed with failure',
                processor=processor_name, errors=str(errors))

    @task
    def normalizing(self):
        ''' The *task* that encapsulates all normalizing related actions. '''

        # allow normalizer to access and add data to the entry metadata
        self._parser_backend.entry_archive.m_add_sub_section(
            datamodel.EntryArchive.section_metadata, self._entry_metadata)

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
                        normalizer(backend.entry_archive).normalize(logger=logger)
                    except Exception as e:
                        self._parser_backend.finishedParsingSession('ParseFailure', [str(e)])
                        logger.error(
                            'normalizer failed with exception', exc_info=e,
                            error=str(e), **context)
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

        self._entry_metadata.apply_domain_metadata(self._parser_backend)
        self._entry_metadata.processed = True

        # persist the calc metadata
        with utils.timer(logger, 'saved calc metadata', step='metadata'):
            self.apply_entry_metadata(self._entry_metadata)

        # index in search
        with utils.timer(logger, 'indexed', step='index'):
            self._entry_metadata.a_elastic.index()

        # persist the archive
        with utils.timer(
                logger, 'archived', step='archive',
                input_size=self.mainfile_file.size) as log_data:

            archive_size = self.write_archive(self._parser_backend)
            log_data.update(archive_size=archive_size)

    def write_archive(self, backend: Backend):
        def filter_processing_logs(logs):
            if len(logs) > 100:
                return [
                    log for log in logs
                    if log.get('level') != 'DEBUG']
            return logs

        if self._calc_proc_logs is None:
            self._calc_proc_logs = []

        if backend is not None:
            entry_archive = backend.entry_archive.m_copy()
        else:
            entry_archive = datamodel.EntryArchive()

        if entry_archive.section_metadata is None:
            entry_archive.m_add_sub_section(datamodel.EntryArchive.section_metadata, self._entry_metadata)

        entry_archive.processing_logs = filter_processing_logs(self._calc_proc_logs)

        try:
            return self.upload_files.write_archive(self.calc_id, entry_archive.m_to_dict())
        except Exception as e:
            if backend is None:
                raise e

            # most likely failed due to domain data, try to write metadata and processing logs
            entry_archive = datamodel.EntryArchive()
            entry_archive.m_add_sub_section(datamodel.EntryArchive.section_metadata, self._entry_metadata)
            entry_archive.processing_logs = filter_processing_logs(self._calc_proc_logs)
            self.upload_files.write_archive(self.calc_id, entry_archive.m_to_dict())
            raise e

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

        with utils.lnr(logger, 'upload delete failed'):
            with utils.timer(
                    logger, 'upload deleted from index', step='index',
                    upload_size=self.upload_files.size):
                search.delete_upload(self.upload_id)

            with utils.timer(
                    logger, 'upload deleted', step='files',
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
            with self.entries_metadata(self.metadata) as calcs:

                with utils.timer(
                        logger, 'upload metadata updated', step='metadata',
                        upload_size=self.upload_files.size):

                    def create_update(calc):
                        calc.published = True
                        calc.with_embargo = calc.with_embargo if calc.with_embargo is not None else False
                        return UpdateOne(
                            {'_id': calc.calc_id},
                            {'$set': {'metadata': calc.m_to_dict(
                                include_defaults=True, categories=[datamodel.MongoMetadata])}})

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
        logger = self.get_logger()
        logger.info('started to re-process')

        # mock the steps of actual processing
        self._continue_with('uploading')

        # extract the published raw files into a staging upload files instance
        self._continue_with('extracting')

        if self.published:
            try:
                staging_upload_files = StagingUploadFiles(self.upload_id)
                # public files exist and there is a staging directory, it is probably old
                # and we delete it first
                staging_upload_files.delete()
                logger.warn('deleted old staging files')

            except KeyError as e:
                logger.info('reprocessing published files')
        else:
            logger.info('reprocessing staging files')

        staging_upload_files = self.upload_files.to_staging_upload_files(create=True)

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

            if self.published:
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

        self.upload_files.re_pack(self.user_metadata())
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

    def match_mainfiles(self) -> Iterator[Tuple[str, object]]:
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
            if modified_upload is None or modified_upload['joined'] is False:
                self.get_logger().info('join')

                # Before cleaning up, run an additional normalizer on phonon
                # calculations. TODO: This should be replaced by a more
                # extensive mechamism that supports more complex dependencies
                # between calculations.
                phonon_calculations = Calc.objects(upload_id=self.upload_id, parser="parsers/phonopy")
                for calc in phonon_calculations:
                    calc.process_phonon()

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
            'You can review your data on your upload page: %s' % config.gui_url(page='uploads'),
            '',
            'If you encounter any issues with your upload, please let us know and reply to this email.',
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
        if self.published:
            staging_upload_files = self.upload_files.to_staging_upload_files()
            logger.info('started to repack re-processed upload')

            with utils.timer(
                    logger, 'reprocessed staged upload packed', step='repack staged',
                    upload_size=self.upload_files.size):

                staging_upload_files.pack(self.user_metadata(), skip_raw=True)

            with utils.timer(
                    logger, 'reprocessed staged upload deleted', step='delete staged',
                    upload_size=self.upload_files.size):

                staging_upload_files.delete()
                self.last_update = datetime.utcnow()
                self.save()

        else:
            logger.info('no cleanup after re-processing unpublished upload')

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
            metadata__nomad_version__ne=config.meta.version)

    @property
    def calcs(self):
        ''' All successfully processed calculations. '''
        return Calc.objects(upload_id=self.upload_id, tasks_status=SUCCESS)

    @contextmanager
    def entries_metadata(
            self, user_metadata: dict = None) -> Iterator[Iterable[datamodel.EntryMetadata]]:
        '''
        This is the :py:mod:`nomad.datamodel` transformation method to transform
        processing upload's entries into list of :class:`nomad.datamodel.EntryMetadata` objects.

        Arguments:
            user_metadata: A dict of user metadata that is applied to the resulting
                datamodel data and the respective calculations.
        '''
        upload_files = self.upload_files

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
                entry_metadata = calc.entry_metadata(upload_files)
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
                entry_metadata = calc.entry_metadata(upload_files)
                entry_metadata.upload_time = self.upload_time
                entry_metadata.upload_name = self.name

                return entry_metadata

        try:
            # read all calc objects first to avoid missing curser errors
            yield [
                get_metadata(calc)
                for calc in list(Calc.objects(upload_id=self.upload_id))]

        finally:
            upload_files.close()

    def entry_ids(self) -> Iterable[str]:
        return [calc.calc_id for calc in Calc.objects(upload_id=self.upload_id)]

    def user_metadata(self) -> Iterable[datamodel.EntryMetadata]:
        return [calc.user_metadata() for calc in Calc.objects(upload_id=self.upload_id)]

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
