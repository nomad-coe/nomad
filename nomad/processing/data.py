#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

'''
This module comprises a set of persistent document classes that hold all user related
data. These are information about users, their uploads and datasets, the associated
calculations, and files


.. autoclass:: Calc

.. autoclass:: Upload

'''

from typing import cast, List, Any, Tuple, Iterator, Dict, cast, Iterable
from mongoengine import (
    StringField, DateTimeField, DictField, BooleanField, IntField, ListField)
import logging
from structlog import wrap_logger
from contextlib import contextmanager
import os.path
from datetime import datetime
from pymongo import UpdateOne
import hashlib
from structlog.processors import StackInfoRenderer, format_exc_info, TimeStamper
import yaml
import json
from functools import lru_cache
import urllib.parse
import requests

from nomad import utils, config, infrastructure, search, datamodel, metainfo, parsing
from nomad.files import (
    PathObject, UploadFiles, ExtractError, ArchiveBasedStagingUploadFiles,
    PublicUploadFiles, StagingUploadFiles)
from nomad.processing.base import Proc, process, task, PENDING, SUCCESS, FAILURE
from nomad.parsing.parsers import parser_dict, match_parser
from nomad.normalizing import normalizers
from nomad.datamodel import (
    EntryArchive, EditableUserMetadata, OasisMetadata, UserProvidableMetadata)
from nomad.archive import (
    write_partial_archive_to_mongo, delete_partial_archives_from_mongo)
from nomad.datamodel.encyclopedia import EncyclopediaMetadata


section_metadata = datamodel.EntryArchive.section_metadata.name
section_workflow = datamodel.EntryArchive.section_workflow.name


_editable_metadata: Dict[str, metainfo.Definition] = {}
_editable_metadata.update(**{
    quantity.name: quantity for quantity in UserProvidableMetadata.m_def.definitions})
_editable_metadata.update(**{
    quantity.name: quantity for quantity in EditableUserMetadata.m_def.definitions})

_oasis_metadata = {
    quantity.name: quantity for quantity in OasisMetadata.m_def.definitions}


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


def _normalize_oasis_upload_metadata(upload_id, metadata):
    # This is overwritten by the tests to do necessary id manipulations
    return upload_id, metadata


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
        'strict': False,
        'indexes': [
            'upload_id',
            'parser',
            ('upload_id', 'mainfile'),
            ('upload_id', 'parser'),
            ('upload_id', 'tasks_status'),
            ('upload_id', 'current_task'),
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
        self._parser_results: EntryArchive = None
        self._upload: Upload = None
        self._upload_files: ArchiveBasedStagingUploadFiles = None
        self._calc_proc_logs: List[Any] = None

        self._entry_metadata = None

    @classmethod
    def get(cls, id):
        return cls.get_by_id(id, 'calc_id')

    @property
    def entry_id(self) -> str:
        ''' Just an alias for calc_id. '''
        return self.calc_id

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
        processing object, not necessarily the user metadata nor the metadata from
        the archive.
        '''
        entry_metadata = datamodel.EntryMetadata()
        if self.parser is not None:
            parser = parser_dict[self.parser]
            if parser.domain:
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
        processing object and the user metadata, not necessarily the metadata from
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

                if method_name == 'error':
                    error = event_dict.get('event', None)
                    if error is not None:
                        self._entry_metadata.processing_errors.append(error)

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
        logger = self.get_logger()

        if config.reprocess_rematch:
            with utils.timer(logger, 'parser matching executed'):
                parser = match_parser(
                    self.upload_files.raw_file_object(self.mainfile).os_path, strict=False)

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

        else:
            parser = parser_dict.get(self.parser)

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
            self._entry_metadata.processing_errors = []
            self._entry_metadata.files = self.upload_files.calc_files(self.mainfile)

            self.parsing()
            self.normalizing()
            self.archiving()
        finally:
            # close loghandler that was not closed due to failures
            try:
                if self._parser_results and self._parser_results.m_resource:
                    self._parser_results.section_metadata = None
                    self._parser_results.m_resource.unload()
            except Exception as e:
                logger.error('could not unload processing results', exc_info=e)

    def _setup_fallback_metadata(self):
        self._entry_metadata = self.create_metadata()
        self._entry_metadata.calc_hash = self.upload_files.calc_hash(self.mainfile)
        self._entry_metadata.last_processing = datetime.utcnow()
        self._entry_metadata.processing_errors = []
        self._entry_metadata.files = self.upload_files.calc_files(self.mainfile)
        self._entry_metadata.parser_name = self.parser

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
                if self._parser_results and self._parser_results.m_resource:
                    self._parser_results.section_metadata = None
                    self._parser_results.m_resource.unload()
            except Exception as e:
                logger.error('could unload processing results', exc_info=e)

    def on_fail(self):
        # in case of failure, index a minimum set of metadata and mark
        # processing failure
        try:
            if self._entry_metadata is None:
                self._setup_fallback_metadata()

            self._entry_metadata.processed = False

            try:
                self.apply_entry_metadata(self._entry_metadata)
            except Exception as e:
                self.get_logger().error(
                    'could not apply entry metadata to entry', exc_info=e)

            try:
                self._entry_metadata.apply_domain_metadata(self._parser_results)
            except Exception as e:
                self.get_logger().error(
                    'could not apply domain metadata to entry', exc_info=e)

            self._entry_metadata.a_elastic.index()
        except Exception as e:
            self.get_logger().error(
                'could not index after processing failure', exc_info=e)

        try:
            self.write_archive(self._parser_results)
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
            if not config.process_reuse_parser:
                if isinstance(parser, parsing.FairdiParser):
                    try:
                        parser = parser.__class__()
                    except Exception as e:
                        self.fail(
                            'could not re-create parser instance',
                            exc_info=e, error=str(e), **context)
                        return
            try:
                self._parser_results = EntryArchive()
                # allow parsers to read/write metadata
                self._parser_results.m_add_sub_section(EntryArchive.section_metadata, self._entry_metadata)
                parser.parse(
                    self.upload_files.raw_file_object(self.mainfile).os_path,
                    self._parser_results, logger=logger)

            except Exception as e:
                self.fail('parser failed with exception', exc_info=e, error=str(e), **context)
                return
            except SystemExit:
                self.fail('parser raised system exit', error='system exit', **context)
                return

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
                arch = archive[self.calc_id]
                phonon_archive = EntryArchive.m_from_dict(arch.to_dict())
            self._entry_metadata = phonon_archive.section_metadata
            self._calc_proc_logs = phonon_archive.processing_logs

            # Re-create the parse results
            self._parser_results = phonon_archive

            # Read in the first referenced calculation. The reference is given as
            # an absolute path which needs to be converted into a path that is
            # relative to upload root.
            scc = self._parser_results.section_run[0].section_single_configuration_calculation[0]
            calculation_refs = scc.section_calculation_to_calculation_refs
            if calculation_refs is None:
                logger.error("No calculation_to_calculation references found")
                return

            relative_ref = scc.section_calculation_to_calculation_refs[0].calculation_to_calculation_external_url
            ref_id = upload_files.calc_id(relative_ref)
            with upload_files.read_archive(ref_id) as archive:
                arch = archive[ref_id]
                ref_archive = EntryArchive.m_from_dict(arch.to_dict())

            # Get encyclopedia method information directly from the referenced calculation.
            ref_enc_method = ref_archive.section_metadata.encyclopedia.method
            if ref_enc_method is None or len(ref_enc_method) == 0 or ref_enc_method.functional_type is None:
                logger.error("No method information available in referenced calculation.")
                return

            self._parser_results.section_metadata.encyclopedia.method = ref_enc_method

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
            with utils.timer(logger, 'calc metadata saved'):
                self.apply_entry_metadata(self._entry_metadata)

            # index in search
            with utils.timer(logger, 'calc metadata indexed'):
                self._entry_metadata.a_elastic.index()

            # persist the archive
            with utils.timer(
                    logger, 'calc archived',
                    input_size=self.mainfile_file.size) as log_data:

                archive_size = self.write_archive(self._parser_results)
                log_data.update(archive_size=archive_size)

    @task
    def normalizing(self):
        ''' The *task* that encapsulates all normalizing related actions. '''

        # allow normalizer to access and add data to the entry metadata
        if self._parser_results.section_metadata is None:
            self._parser_results.m_add_sub_section(
                datamodel.EntryArchive.section_metadata, self._entry_metadata)

        for normalizer in normalizers:
            if normalizer.domain != parser_dict[self.parser].domain:
                continue

            normalizer_name = normalizer.__name__
            context = dict(normalizer=normalizer_name, step=normalizer_name)
            logger = self.get_logger(**context)

            with utils.timer(logger, 'normalizer executed', input_size=self.mainfile_file.size):
                try:
                    normalizer(self._parser_results).normalize(logger=logger)
                    logger.info('normalizer completed successfull', **context)
                except Exception as e:
                    self.fail('normalizer failed with exception', exc_info=e, error=str(e), **context)

    def _read_metadata_from_file(self, logger):
        # metadata file name defined in nomad.config nomad_metadata.yaml/json
        # which can be placed in the directory containing the mainfile or somewhere up
        # highest priority is directory with mainfile
        metadata_file = config.metadata_file_name
        metadata_dir = os.path.dirname(self.mainfile_file.os_path)
        upload_raw_dir = self.upload_files._raw_dir.os_path

        metadata = {}
        metadata_part = None
        # apply the nomad files of the current directory and parent directories
        while True:
            metadata_part = self.upload.metadata_file_cached(
                os.path.join(metadata_dir, metadata_file))
            for key, val in metadata_part.items():
                if key in ['entries', 'oasis_datasets']:
                    continue
                metadata.setdefault(key, val)

            if metadata_dir == upload_raw_dir:
                break

            metadata_dir = os.path.dirname(metadata_dir)

        # Top-level nomad file can also contain an entries dict with entry
        # metadata per mainfile as key. This takes precedence of the other files.
        entries = metadata_part.get('entries', {})
        metadata_part = entries.get(self.mainfile, {})
        for key, val in metadata_part.items():
            metadata[key] = val

        if len(metadata) > 0:
            logger.info('Apply user metadata from nomad.yaml/json file')

        for key, val in metadata.items():
            if key in ['entries', 'skip_matching']:
                continue

            definition = _editable_metadata.get(key, None)
            if definition is None and self.upload.from_oasis:
                definition = _oasis_metadata.get(key, None)

            if key == 'uploader':
                if datamodel.User.get(self.upload.user_id).is_admin:
                    definition = datamodel.EntryMetadata.uploader

            if definition is None:
                logger.warn('Users cannot set metadata', quantity=key)
                continue

            try:
                self._entry_metadata.m_set(definition, val)
                if definition == datamodel.EntryMetadata.calc_id:
                    self.calc_id = val
            except Exception as e:
                logger.error(
                    'Could not apply user metadata from nomad.yaml/json file',
                    quantitiy=definition.name, exc_info=e)

    @task
    def archiving(self):
        ''' The *task* that encapsulates all archival related actions. '''
        logger = self.get_logger()

        self._entry_metadata.apply_domain_metadata(self._parser_results)
        self._entry_metadata.processed = True

        if self.upload.publish_directly:
            self._entry_metadata.published = False

        try:
            self._read_metadata_from_file(logger)
        except Exception as e:
            logger.error('could not process user metadata in nomad.yaml/json file', exc_info=e)

        # persist the calc metadata
        with utils.timer(logger, 'calc metadata saved'):
            self.apply_entry_metadata(self._entry_metadata)

        # index in search
        with utils.timer(logger, 'calc metadata indexed'):
            self._entry_metadata.a_elastic.index()

        # persist the archive
        with utils.timer(
                logger, 'calc archived',
                input_size=self.mainfile_file.size) as log_data:

            archive_size = self.write_archive(self._parser_results)
            log_data.update(archive_size=archive_size)

    def write_archive(self, archive: EntryArchive):
        # save the archive mongo entry
        try:
            if self._entry_metadata.processed:
                write_partial_archive_to_mongo(archive)
        except Exception as e:
            self.get_logger().error('could not write mongodb archive entry', exc_info=e)

        # add the processing logs to the archive
        def filter_processing_logs(logs):
            if len(logs) > 100:
                return [
                    log for log in logs
                    if log.get('level') != 'DEBUG']
            return logs

        if self._calc_proc_logs is None:
            self._calc_proc_logs = []

        if archive is not None:
            archive = archive.m_copy()
        else:
            archive = datamodel.EntryArchive()

        if archive.section_metadata is None:
            archive.m_add_sub_section(datamodel.EntryArchive.section_metadata, self._entry_metadata)

        archive.processing_logs = filter_processing_logs(self._calc_proc_logs)

        # save the archive msg-pack
        try:
            return self.upload_files.write_archive(self.calc_id, archive.m_to_dict())
        except Exception as e:
            # most likely failed due to domain data, try to write metadata and processing logs
            archive = datamodel.EntryArchive()
            archive.m_add_sub_section(datamodel.EntryArchive.section_metadata, self._entry_metadata)
            archive.processing_logs = filter_processing_logs(self._calc_proc_logs)
            self.upload_files.write_archive(self.calc_id, archive.m_to_dict())
            raise e

    def __str__(self):
        return 'calc %s calc_id=%s upload_id%s' % (super().__str__(), self.calc_id, self.upload_id)


class Upload(Proc):
    '''
    Represents uploads in the databases. Provides persistence access to the files storage,
    and processing state.

    Attributes:
        name: Optional user provided upload name.
        upload_path: The fs path were the uploaded files was stored during upload.
        temporary: True if the uploaded file should be removed after extraction.

        upload_id: The upload id generated by the database or the uploaded NOMAD deployment.
        upload_time: Datetime of the original upload independent of the NOMAD deployment
            it was first uploaded to.
        user_id: The id of the user that created this upload.
        published: Boolean that indicates that the upload is published on this NOMAD deployment.
        publish_time: Datetime when the upload was initially published on this NOMAD deployment.
        last_update: Datetime of the last modifying process run (publish, re-processing, upload).

        publish_directly: Boolean indicating that this upload should be published after initial processing.
        from_oasis: Boolean indicating that this upload is coming from another NOMAD deployment.
        oasis_id: The deployment id of the NOMAD that uploaded the upload.
        published_to: A list of deployment ids where this upload has been successfully uploaded to.

        joined: Boolean indicates if the running processing has joined (:func:`check_join`).
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

    publish_directly = BooleanField(default=False)
    from_oasis = BooleanField(default=False)
    oasis_deployment_id = StringField(default=None)
    published_to = ListField(StringField())

    joined = BooleanField(default=False)

    meta: Any = {
        'strict': False,
        'indexes': [
            'user_id', 'tasks_status', 'process_status', 'published', 'upload_time', 'create_time'
        ]
    }

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.publish_directly = self.publish_directly or self.from_oasis
        self._upload_files: ArchiveBasedStagingUploadFiles = None

    @lru_cache()
    def metadata_file_cached(self, path):
        for ext in config.metadata_file_extensions:
            full_path = '%s.%s' % (path, ext)
            if os.path.isfile(full_path):
                try:
                    with open(full_path) as f:
                        if full_path.endswith('.json'):
                            return json.load(f)
                        elif full_path.endswith('.yaml') or full_path.endswith('.yml'):
                            return yaml.load(f, Loader=getattr(yaml, 'FullLoader'))
                        else:
                            return {}
                except Exception as e:
                    self.get_logger().warn('could not parse nomad.yaml/json', path=path, exc_info=e)
                    # ignore the file contents if the file is not parsable
                    pass
        return {}

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
        logger = self.get_logger(upload_size=self.upload_files.size)

        with utils.lnr(logger, 'upload delete failed'):
            with utils.timer(logger, 'upload deleted from index'):
                search.delete_upload(self.upload_id)

            with utils.timer(logger, 'upload partial archives deleted'):
                calc_ids = [calc.calc_id for calc in Calc.objects(upload_id=self.upload_id)]
                delete_partial_archives_from_mongo(calc_ids)

            with utils.timer(logger, 'upload files deleted'):
                self.upload_files.delete()

            self.delete()

    @process
    def delete_upload(self):
        '''
        Deletes the upload, including its processing state and
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

        logger = self.get_logger(upload_size=self.upload_files.size)
        logger.info('started to publish')

        with utils.lnr(logger, 'publish failed'):
            with self.entries_metadata(self.metadata) as calcs:

                with utils.timer(logger, 'upload metadata updated'):
                    def create_update(calc):
                        calc.published = True
                        calc.with_embargo = calc.with_embargo if calc.with_embargo is not None else False
                        return UpdateOne(
                            {'_id': calc.calc_id},
                            {'$set': {'metadata': calc.m_to_dict(
                                include_defaults=True, categories=[datamodel.MongoMetadata])}})

                    Calc._get_collection().bulk_write([create_update(calc) for calc in calcs])

                if isinstance(self.upload_files, StagingUploadFiles):
                    with utils.timer(logger, 'staged upload files packed'):
                        self.upload_files.pack(calcs)

                with utils.timer(logger, 'index updated'):
                    search.publish(calcs)

                if isinstance(self.upload_files, StagingUploadFiles):
                    with utils.timer(logger, 'upload staging files deleted'):
                        self.upload_files.delete()
                        self.published = True
                        self.publish_time = datetime.utcnow()
                        self.last_update = datetime.utcnow()
                        self.save()
                else:
                    self.last_update = datetime.utcnow()
                    self.save()

    @process
    def publish_from_oasis(self):
        '''
        Uploads the already published upload to a different NOMAD deployment. This allows
        to push uploads from an OASIS to the central NOMAD.
        '''
        assert self.published, \
            'Only published uploads can be published to the central NOMAD.'
        assert config.oasis.central_nomad_deployment_id not in self.published_to, \
            'Upload is already published to the central NOMAD.'

        from nomad.cli.client.client import _create_client as create_client
        central_nomad_client = create_client(
            user=config.keycloak.username,
            password=config.keycloak.password,
            api_base_url=config.oasis.central_nomad_api_url,
            use_token=False)

        # compile oasis metadata for the upload
        upload_metadata = dict(upload_time=str(self.upload_time))
        upload_metadata_entries = {}
        upload_metadata_datasets = {}
        for calc in self.calcs:
            entry_metadata = dict(**{
                key: str(value) if isinstance(value, datetime) else value
                for key, value in calc.metadata.items()
                if key in _editable_metadata or key in _oasis_metadata})
            entry_metadata['calc_id'] = calc.calc_id
            if entry_metadata.get('with_embargo'):
                continue
            upload_metadata_entries[calc.mainfile] = entry_metadata
            if 'datasets' in entry_metadata:
                for dataset_id in entry_metadata['datasets']:
                    if dataset_id in upload_metadata_datasets:
                        continue

                    dataset = datamodel.Dataset.m_def.a_mongo.get(dataset_id=dataset_id)
                    upload_metadata_datasets[dataset_id] = dataset.m_to_dict()

        upload_metadata['entries'] = upload_metadata_entries
        upload_metadata['oasis_datasets'] = {
            dataset['name']: dataset for dataset in upload_metadata_datasets.values()}
        oasis_upload_id, upload_metadata = _normalize_oasis_upload_metadata(
            self.upload_id, upload_metadata)

        self.last_status_message = 'Compiled metadata to upload to the central NOMAD.'
        self.save()

        assert len(upload_metadata_entries) > 0, \
            'Only uploads with public contents can be published to the central NOMAD.'

        # add oasis metadata to the upload
        public_upload_files = cast(PublicUploadFiles, self.upload_files)
        public_upload_files.add_metadata_file(upload_metadata)
        file_to_upload = public_upload_files.public_raw_data_file

        self.last_status_message = 'Prepared the upload for uploading to central NOMAD.'
        self.save()

        # upload to central NOMAD
        oasis_admin_token = central_nomad_client.auth.get_auth().response().result.access_token
        upload_headers = dict(Authorization='Bearer %s' % oasis_admin_token)
        upload_parameters = dict(
            oasis_upload_id=oasis_upload_id,
            oasis_uploader_id=self.user_id,
            oasis_deployment_id=config.meta.deployment_id)
        upload_url = '%s/uploads/?%s' % (
            config.oasis.central_nomad_api_url,
            urllib.parse.urlencode(upload_parameters))

        with open(file_to_upload, 'rb') as f:
            response = requests.put(upload_url, headers=upload_headers, data=f)

        if response.status_code != 200:
            self.get_logger().error(
                'Could not upload to central NOMAD', status_code=response.status_code)
            self.last_status_message = 'Could not upload to central NOMAD.'
            return

        self.published_to.append(config.oasis.central_nomad_deployment_id)
        self.last_status_message = 'Successfully uploaded to central NOMAD.'

    @process
    def re_process_upload(self):
        '''
        A *process* that performs the re-processing of a earlier processed
        upload.

        Runs the distributed process of fully reparsing/re-normalizing an existing and
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
                with utils.timer(logger, 'upload staging files deleted'):
                    staging_upload_files.delete()
                    logger.warn('deleted old staging files')

            except KeyError as e:
                logger.info('reprocessing published files')
        else:
            logger.info('reprocessing staging files')

        with utils.timer(logger, 'upload extracted'):
            staging_upload_files = self.upload_files.to_staging_upload_files(create=True)

        self._continue_with('parse_all')
        try:
            if config.reprocess_match:
                with utils.timer(logger, 'calcs match on reprocess'):
                    with_embargo_values = Calc._get_collection().distinct(
                        'metadata.with_embargo', dict(upload_id=self.upload_id))
                    with_embargo = with_embargo_values != [False]

                    for filename, parser in self.match_mainfiles():
                        calc_id = staging_upload_files.calc_id(filename)
                        try:
                            Calc.get(calc_id)
                        except KeyError:
                            calc = Calc.create(
                                calc_id=calc_id,
                                mainfile=filename, parser=parser.name,
                                worker_hostname=self.worker_hostname,
                                create_time=self.upload_time,
                                upload_id=self.upload_id)

                            calc.metadata.update(
                                published=self.published,
                                with_embargo=with_embargo)

                            calc.save()

            with utils.timer(logger, 'calcs resetted'):
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

            with utils.timer(logger, 'calcs re-processing called'):
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
        self.upload_files.re_pack(self.user_metadata())

    @process
    def process_upload(self):
        ''' A *process* that performs the initial upload processing. '''
        self.extracting()

        if self.from_oasis:
            # we might need to add datasets from the oasis before processing and
            # adding the entries
            oasis_metadata_file = os.path.join(self.upload_files._raw_dir.os_path, config.metadata_file_name + '.json')
            with open(oasis_metadata_file, 'rt') as f:
                oasis_metadata = json.load(f)
            oasis_datasets = oasis_metadata.get('oasis_datasets', {})
            metadata_was_changed = False
            for oasis_dataset in oasis_datasets.values():
                try:
                    existing_dataset = datamodel.Dataset.m_def.a_mongo.get(
                        user_id=self.user_id, name=oasis_dataset['name'])
                except KeyError:
                    datamodel.Dataset(**oasis_dataset).a_mongo.save()
                else:
                    oasis_dataset_id = oasis_dataset['dataset_id']
                    if existing_dataset.dataset_id != oasis_dataset_id:
                        # A dataset for the same user with the same name was created
                        # in both deployments. We consider this to be the "same" dataset.
                        # These datasets have different ids and we need to migrate the provided
                        # dataset ids:
                        for entry in oasis_metadata['entries'].values():
                            entry_datasets = entry.get('datasets', [])
                            for index, dataset_id in enumerate(entry_datasets):
                                if dataset_id == oasis_dataset_id:
                                    entry_datasets[index] = existing_dataset.dataset_id
                                    metadata_was_changed = True

            if metadata_was_changed:
                with open(oasis_metadata_file, 'wt') as f:
                    json.dump(oasis_metadata, f)

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
    def staging_upload_files(self) -> StagingUploadFiles:
        return self.upload_files.to_staging_upload_files()

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
            with utils.timer(logger, 'upload extracted', upload_size=self.upload_files.size):
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
                    awk < '%s' >> '%s' '
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
        metadata = self.metadata_file_cached(
            os.path.join(self.upload_files.os_path, 'raw', config.metadata_file_name))
        skip_matching = metadata.get('skip_matching', False)
        entries_metadata = metadata.get('entries', {})

        directories_with_match: Dict[str, str] = dict()
        upload_files = self.upload_files.to_staging_upload_files()
        for filename in upload_files.raw_file_manifest():
            self._preprocess_files(filename)

            if skip_matching and filename not in entries_metadata:
                continue

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

        oasis_metadata = {}
        if self.from_oasis:
            oasis_metadata = self.metadata_file_cached(
                os.path.join(self.upload_files.os_path, 'raw', config.metadata_file_name)).get('entries', {})

        with utils.timer(logger, 'calcs processing called'):
            for filename, parser in self.match_mainfiles():
                oasis_entry_metadata = oasis_metadata.get(filename)
                if oasis_entry_metadata is not None:
                    calc_id = oasis_entry_metadata.get('calc_id')
                    if calc_id is None:
                        logger.warn('Oasis entry without id', mainfile=filename)
                        calc_id = self.upload_files.calc_id(filename)
                else:
                    calc_id = self.upload_files.calc_id(filename)

                calc = Calc.create(
                    calc_id=calc_id,
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
        calculation processing have been triggered (:func:`parse_all` or
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
                # extensive mechanism that supports more complex dependencies
                # between calculations.
                phonon_calculations = Calc.objects(upload_id=self.upload_id, parser="parsers/phonopy")
                for calc in phonon_calculations:
                    calc.process_phonon()

                self.cleanup()
            else:
                # the join was already done due to a prior call
                pass

    def reset(self, force=False):
        self.joined = False
        super().reset(force=force)

    @classmethod
    def reset_pymongo_update(cls, worker_hostname: str = None):
        update = super().reset_pymongo_update()
        update.update(joined=False)
        return update

    def _cleanup_after_processing(self):
        # send email about process finish
        if not self.publish_directly:
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

        if not self.publish_directly or self.processed_calcs == 0:
            return

        logger = self.get_logger(upload_size=self.upload_files.size)
        logger.info('started to publish upload directly')

        with utils.lnr(logger, 'publish failed'):
            metadata = self.metadata_file_cached(
                os.path.join(self.upload_files.os_path, 'raw', config.metadata_file_name))

            with self.entries_metadata(self.metadata) as calcs:
                with utils.timer(logger, 'upload metadata updated'):
                    def create_update(calc):
                        calc.published = True
                        calc.with_embargo = calc.with_embargo if calc.with_embargo is not None else False
                        return UpdateOne(
                            {'_id': calc.calc_id},
                            {'$set': {'metadata': calc.m_to_dict(
                                include_defaults=True, categories=[datamodel.MongoMetadata])}})

                    Calc._get_collection().bulk_write([create_update(calc) for calc in calcs])

                with utils.timer(logger, 'upload staging files packed'):
                    self.upload_files.pack(calcs)

                with utils.timer(logger, 'index updated'):
                    search.publish(calcs)

            with utils.timer(logger, 'upload staging files deleted'):
                self.upload_files.delete()

            if self.from_oasis:
                if metadata is not None:
                    self.upload_time = metadata.get('upload_time')

                if self.upload_time is None:
                    self.upload_time = datetime.utcnow()
                    logger.warn('oasis upload without upload time')

            self.publish_time = datetime.utcnow()
            self.published = True
            self.last_update = datetime.utcnow()
            self.save()

    def _cleanup_after_re_processing(self):
        logger = self.get_logger(upload_size=self.upload_files.size)
        if self.published:
            staging_upload_files = self.staging_upload_files
            logger.info('started to repack re-processed upload')

            with utils.timer(logger, 'staged upload files re-packed'):
                staging_upload_files.pack(self.user_metadata(), skip_raw=True)

            with utils.timer(logger, 'staged upload files deleted'):
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
        if not order_by:
            return query
        if type(order_by) == str:
            return query.order_by(order_by)
        assert type(order_by) == tuple, 'order_by must be a string or a tuple if set'
        return query.order_by(*order_by)

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
