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
This module contains functions to read data from NOMAD coe, external sources,
other/older nomad@FAIRDI instances to mass upload it to a new nomad@FAIRDI instance.

.. autoclass:: NomadCOEMigration
.. autoclass:: SourceCalc
"""

from typing import Generator, Tuple, List
import os.path
import json
import zipstream
import math
from mongoengine import Document, IntField, StringField, DictField
from passlib.hash import bcrypt

from nomad import utils, config
from nomad.coe_repo import User, Calc
from nomad.datamodel import CalcWithMetadata
from nomad.processing import FAILURE, SUCCESS


class SourceCalc(Document):
    """
    Mongo document used as a calculation, upload, and metadata db and index
    build from a given source db. Each :class:`SourceCacl` entry relates
    a pid, mainfile, upload "id" with each other for a corressponding calculation.
    It might alos contain the user metadata. The uploads are "id"ed via the
    specific path segment that identifies an upload on the CoE repo FS(s) without
    any prefixes (e.g. $EXTRACTED, /data/upload, etc.)
    """
    pid = IntField(primary_key=True)
    mainfile = StringField()
    upload = StringField()
    metadata = DictField()

    extracted_prefix = '$EXTRACTED/'
    sites = ['/data/nomad/extracted/', '/nomad/repository/extracted/']
    prefixes = [extracted_prefix] + sites

    meta = dict(indexes=['pid', 'upload'])

    _dataset_cache: dict = {}

    @staticmethod
    def index(source, drop: bool = False, with_metadata: bool = True, per_query: int = 100) \
            -> Generator[Tuple['SourceCalc', int], None, None]:
        """
        Creates a collection of :class:`SourceCalc` documents that represent source repo
        db entries.

        Arguments:
            source: The source db sql alchemy session
            drop: True to drop and create a new collection, update the existing otherwise,
                default is False.
            with_metadata: True to also grab all metadata and store it, default is True.
            per_query: The implementation tries to grab almost all data with a heavely joined
                query on the CoE snoflake/star shaped schema.
                The query cannot ask for the whole db at once: choose how many calculations
                should be read at a time to optimize for your application.

        Returns:
            yields tuples (:class:`SourceCalc`, #calcs_total)
        """
        if drop:
            SourceCalc.drop_collection()

        last_source_calc = SourceCalc.objects().order_by('-pid').first()
        start_pid = last_source_calc.pid if last_source_calc is not None else 0
        source_query = source.query(Calc)
        total = source_query.count()

        while True:
            calcs = source_query.filter(Calc.coe_calc_id > start_pid).order_by(Calc.coe_calc_id).limit(per_query)
            source_calcs = []
            for calc in calcs:
                if calc.calc_metadata is None or calc.calc_metadata.filenames is None:
                    yield None, total
                    continue  # dataset case

                filenames = json.loads(calc.calc_metadata.filenames.decode('utf-8'))
                filename = filenames[0]
                for prefix in SourceCalc.prefixes:
                    filename = filename.replace(prefix, '')
                segments = [file.strip('\\') for file in filename.split('/')]

                source_calc = SourceCalc(pid=calc.pid)
                source_calc.upload = segments[0]
                source_calc.mainfile = os.path.join(*segments[1:])
                if with_metadata:
                    source_calc.metadata = calc.to(CalcWithMetadata)
                source_calcs.append(source_calc)
                start_pid = source_calc.pid

                yield source_calc, total

            if len(source_calcs) == 0:
                break
            else:
                SourceCalc.objects.insert(source_calcs)


class ZipStreamFileAdaptor:
    def __init__(self, path_to_files: str) -> None:
        self.zip = zipstream.ZipFile()
        self.zip.write(path_to_files)
        self._current_chunk = bytes()
        self._current_chunk_index = 0

    def read(self, n: int) -> bytes:
        while (len(self._current_chunk) - self._current_chunk_index) < n:
            next_chunk = next(self.zip, None)
            left_over_chunk = self._current_chunk[self._current_chunk_index:]
            if next_chunk is None:
                self._current_chunk = bytes()
                self._current_chunk_index = 0
                return left_over_chunk

            self._current_chunk = left_over_chunk + next_chunk
            self._current_chunk_index = 0

        old_index = self._current_chunk_index
        self._current_chunk_index = self._current_chunk_index + n
        return self._current_chunk[old_index:self._current_chunk_index]


class NomadCOEMigration:
    """
    Drives a migration from the NOMAD coe repository db to nomad@FAIRDI. It is assumed
    that this class is never used on the worker or api service. It assumes the
    default coe repo connection as a connection to the source repository db.

    Attributes:
        source: SQLAlchemy session for the source NOMAD coe repository db.

    Arguments:
        sites: Directories that might contain uploads to migrate. Use to override defaults.
        pid_prefix: All PIDs for previously unknown calculations will get a PID higher
            than that. Use to override default.
    """

    default_sites = [
        '/nomad/repository/data/uploads',
        '/nomad/repository/data/extracted',
        '/data/nomad/uploaded/',
        '/data/nomad/extracted/']

    default_pid_prefix = int(1e7)

    archive_filename = 'archive.tar.gz'
    """ The standard name for tarred uploads in the CoE repository. """

    def __init__(
            self,
            sites: List[str] = default_sites,
            pid_prefix: int = default_pid_prefix) -> None:

        self.sites, self.pid_prefix = sites, pid_prefix
        self.logger = utils.get_logger(__name__)
        from nomad.infrastructure import repository_db
        self.source = repository_db

    def copy_users(self, target_db):
        """ Copy all users, keeping their ids, within a single transaction. """
        target_db.begin()
        for source_user in self.source.query(User).all():
            self.source.expunge(source_user)  # removes user from the source session
            target_db.merge(source_user)

        admin = target_db.query(User).filter_by(email='admin').first()
        if admin is None:
            admin = User(
                user_id=0, email='admin', firstname='admin', lastname='admin',
                password=bcrypt.encrypt(config.services.admin_password, ident='2y'))
            target_db.add(admin)
        target_db.commit()

    def migrate(self, *args):
        """
        Migrate the given uploads.

        It takes upload 'id's as args. Alternatively takes absolute paths to uploads.
        It tries to be as flexible as possible with those 'id's: looking at all
        configured sites, dealing with extracted and tarred/zipped uploads, dealing
        with paths to files and directories.

        Requires a build :func:`index` to look for existing data in the source db. This
        will be used to add user (and other, PID, ...) metadata and validate calculations.

        Uses PIDs of identified old calculations. Will create new PIDs for previously
        unknown uploads. New PIDs will be choosed from a `prefix++` range of ints.
        """

        migrated = 0
        upload_specs = args
        for upload_spec in upload_specs:
            # identify uploads
            if os.path.isabs(upload_spec):
                if os.path.exists(upload_spec):
                    upload_path = upload_spec
                else:
                    upload_path = None
            else:
                for site in self.sites:
                    potential_upload_path = os.path.join(site, upload_spec)
                    if os.path.exists(potential_upload_path):
                        upload_path = potential_upload_path
                        break

            if upload_path is None:
                self.logger.error('upload does not exist', upload_spec=upload_spec)
                continue

            if os.path.isfile(upload_path):
                upload_archive_f = open(upload_path, 'rb')
                upload_id = os.path.split(os.path.split(upload_path)[0])[1]
                upload_name = os.path.basename(upload_path)
            else:
                potential_upload_archive = os.path.join(upload_path, NomadCOEMigration.archive_filename)
                if os.path.isfile(potential_upload_archive):
                    upload_archive_f = open(potential_upload_archive, 'rb')
                    upload_id = os.path.split(os.path.split(potential_upload_archive)[0])[1]
                    upload_name = '%s.tar.gz' % upload_id
                else:
                    upload_id = os.path.split(upload_path)[1]
                    upload_archive_f = ZipStreamFileAdaptor(upload_path)
                    upload_name = '%s.zip' % upload_id

            # process and upload
            from nomad.client import create_client
            client = create_client()
            upload = client.uploads.upload(file=upload_archive_f, name=upload_name).response().result

            upload_logger = self.logger.bind(
                source_upload_id=upload_id,
                upload_id=upload.upload_id)

            # grab old metadata
            upload_metadata_calcs = list()
            metadata_dict = dict()
            upload_metadata = dict(calculations=upload_metadata_calcs)
            for source_calc in SourceCalc.objects(upload=upload_id):
                upload_calc_metadata = dict(
                    mainfile=source_calc.mainfile,
                    _pid=source_calc.pid)
                upload_calc_metadata.update(source_calc.user_metadata)  # TODO to upload calc metadata
                upload_metadata_calcs.append(upload_calc_metadata)
                source_metadata = CalcWithMetadata(**source_calc.metadata)
                source_metadata.__migrated = None
                metadata_dict[source_calc.mainfile] = source_metadata

            # wait for finished processing
            while upload.tasks_status not in [SUCCESS, FAILURE]:
                upload = client.uploads.get_upload(upload_id=upload.upload_id).response().result

            if upload.tasks_status == FAILURE:
                upload_logger.info('upload could not be processed', errors=upload.errors)

            # verify
            total_calcs = upload.calcs.pagination.total
            total_source_calcs = len(metadata_dict)
            unprocessed_calcs = 0
            migrated_calcs = 0
            calcs_with_diffs = 0
            new_calcs = 0
            missing_calcs = 0

            for page in range(0, math.ceil(total_calcs / 100)):
                upload = client.uploads.get_upload(
                    upload_id=upload.upload_id, per_page=100, page=page,
                    order_by='mainfile')

                for calc_proc in upload.calcs.results:
                    calc_logger = upload_logger.bind(
                        calc_id=calc_proc.calc_id,
                        mainfile=calc_proc.mainfile)

                    source_calc = metadata_dict.get(calc_proc.mainfile, None)
                    repo_calc = None
                    if calc_proc.tasks_status == SUCCESS:
                        repo_calc = client.uploads.get_repo(
                            upload_id=upload.upload_id,
                            calc_id=calc_proc.calc_id).response().result
                    else:
                        unprocessed_calcs += 1
                        calc_logger.info(
                            'could not process a calc%s.' %
                            ', that is in source' if source_calc is not None else '')

                        if source_calc is not None:
                            source_calc.__migrated = False
                        continue

                    if source_calc is None:
                        calc_logger.info('processed a calc that has no source')
                        new_calcs += 1
                        continue

                    # TODO add calc metadata to index
                    # TODO do the comparison
                    has_diff = False
                    if source_calc.mainfile != repo_calc['section_repository_info']['repository_filepaths'][0]:
                        has_diff = True
                        calc_logger.info('source target missmatch', quantity='mainfile')

                    if has_diff:
                        calcs_with_diffs += 1

                    source_calc.__migrated = True

            for source_calc in upload_metadata_calcs:
                if source_calc.__migrated is None:
                    missing_calcs += 1
                    upload_logger.info('no match for source calc', mainfile=source_calc.mainfile)
                elif source_calc.__migrated is False:
                    upload_logger.info('source calc not processed', mainfile=source_calc.mainfile)

            upload_metadata['calculations'] = [
                calc for calc in upload_metadata['calculations'] if calc.__migrated]

            # commit with metadata
            if total_calcs > unprocessed_calcs:
                upload = client.uploads.exec_upload_command(
                    upload_id=upload.upload_id, payload=dict(command='commit', metadata=upload_metadata)).response().result

                while upload.process_running:
                    client.uploads.get_upload(upload_id=upload.upload_id)

            # report
            upload_logger.info(
                'migrated upload',
                total_calcs=total_calcs,
                total_source_calcs=total_source_calcs,
                unprocessed_calcs=unprocessed_calcs,
                migrated_calcs=migrated_calcs,
                calcs_with_diffs=calcs_with_diffs,
                new_calcs=new_calcs,
                missing_calcs=new_calcs)
            migrated += 1

        return migrated

    def index(self, *args, **kwargs):
        """ see :func:`SourceCalc.index` """
        return SourceCalc.index(self.source, *args, **kwargs)
