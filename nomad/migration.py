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

from typing import Generator, Tuple, List, Iterable
import os.path
import zipstream
import zipfile
import math
from mongoengine import Document, IntField, StringField, DictField
from passlib.hash import bcrypt
from werkzeug.contrib.iterio import IterIO
import time
from bravado.exception import HTTPNotFound

from nomad import utils, config, infrastructure
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

    meta = dict(indexes=['upload'])

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
            yields tuples (:class:`SourceCalc`, #calcs_total[incl. datasets])
        """
        logger = utils.get_logger(__name__)
        if drop:
            SourceCalc.drop_collection()

        last_source_calc = SourceCalc.objects().order_by('-pid').first()
        start_pid = last_source_calc.pid if last_source_calc is not None else 0
        source_query = source.query(Calc)
        total = source_query.count() - SourceCalc.objects.count()

        while True:
            query_timer = utils.timer(logger, 'query source db')
            query_timer.__enter__()  # pylint: disable=E1101
            calcs: Iterable[Calc] = source_query \
                .filter(Calc.coe_calc_id > start_pid) \
                .order_by(Calc.coe_calc_id) \
                .limit(per_query)

            source_calcs = []
            for calc in calcs:
                query_timer.__exit__(None, None, None)  # pylint: disable=E1101
                try:
                    filenames = calc.files
                    if filenames is None or len(filenames) == 0:
                        continue  # dataset case

                    filename = filenames[0]
                    if len(filenames) == 1 and (filename.endswith('.tgz') or filename.endswith('.zip')):
                        continue  # also a dataset, some datasets have a downloadable archive

                    for prefix in SourceCalc.prefixes:
                        filename = filename.replace(prefix, '')
                    segments = [file.strip('\\') for file in filename.split('/')]

                    source_calc = SourceCalc(pid=calc.pid)
                    source_calc.upload = segments[0]
                    source_calc.mainfile = os.path.join(*segments[1:])
                    if with_metadata:
                        source_calc.metadata = calc.to_calc_with_metadata().__dict__
                    source_calcs.append(source_calc)
                    start_pid = source_calc.pid

                    yield source_calc, total
                except Exception as e:
                    logger.error('could not index', pid=calc.pid, exc_info=e)

            if len(source_calcs) == 0:
                break
            else:
                with utils.timer(logger, 'write index'):
                    SourceCalc.objects.insert(source_calcs)


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
        self._client = None
        self.source = infrastructure.repository_db

    @property
    def client(self):
        if self._client is None:
            from nomad.client import create_client
            self._client = create_client()

        return self._client

    def copy_users(self, target_db):
        """ Copy all users, keeping their ids, within a single transaction. """
        target_db.begin()
        for source_user in self.source.query(User).all():
            self.source.expunge(source_user)  # removes user from the source session
            target_db.merge(source_user)

        admin = target_db.query(User).filter_by(email='admin').first()
        if admin is None:
            admin = User(
                user_id=0, email='admin', first_name='admin', last_name='admin',
                password=bcrypt.encrypt(config.services.admin_password, ident='2y'))
            target_db.add(admin)
        target_db.commit()

    def set_new_pid_prefix(self, target_db, prefix=7000000):
        target_db.begin()
        target_db.execute('ALTER SEQUENCE calculations_calc_id_seq RESTART WITH %d' % prefix)
        target_db.commit()

    def _to_comparable_list(self, list):
        for item in list:
            if isinstance(item, dict):
                for key in item.keys():
                    if key.endswith('id'):
                        yield item.get(key)
            else:
                yield item

    def _validate(self, upload_id: str, calc_id: str, source_calc: CalcWithMetadata, logger) -> bool:
        """
        Validates the given processed calculation, assuming that the data in the given
        source_calc is correct.

        Returns:
            False, if the calculation differs from the source calc.
        """
        repo_calc = self.client.repo.get_repo_calc(
            upload_id=upload_id, calc_id=calc_id).response().result

        is_valid = True
        for key, target_value in repo_calc.items():
            if key in ['calc_id', 'upload_id', 'files', 'calc_hash']:
                continue

            source_value = getattr(source_calc, key, None)

            def report_mismatch():
                logger.info(
                    'source target missmatch', quantity=key,
                    source_value=source_value, target_value=target_value)

            if (source_value is None or target_value is None) and source_value != target_value:
                report_mismatch()
                is_valid = False
                continue

            if isinstance(target_value, list):
                source_list = list(self._to_comparable_list(source_value))
                target_list = list(self._to_comparable_list(target_value))
                if len(set(source_list).intersection(target_list)) != len(target_list):
                    report_mismatch()
                    is_valid = False
                continue

            if isinstance(source_value, str):
                source_value = source_value.lower()
                target_value = str(target_value).lower()

            if source_value != target_value:
                report_mismatch()
                is_valid = False

        return is_valid

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

        Returns: Yields a dictionary with status and statistics for each given upload.
        """

        upload_specs = args
        for upload_spec in upload_specs:
            # identify upload
            upload_path = None
            abs_upload_path = os.path.abspath(upload_spec)
            if os.path.exists(abs_upload_path):
                upload_path = upload_spec
            else:
                for site in self.sites:
                    potential_upload_path = os.path.join(site, upload_spec)
                    if os.path.exists(potential_upload_path):
                        upload_path = potential_upload_path
                        break

            if upload_path is None:
                error = 'upload does not exist'
                self.logger.error(error, upload_spec=upload_spec)
                yield dict(status=FAILURE, error=error)
                continue

            # prepare the upload by determining/creating an upload file, name, source upload id
            if os.path.isfile(upload_path):
                upload_archive_f = open(upload_path, 'rb')
                source_upload_id = os.path.split(os.path.split(upload_path)[0])[1]
                upload_name = os.path.basename(upload_path)
            else:
                potential_upload_archive = os.path.join(
                    upload_path, NomadCOEMigration.archive_filename)
                if os.path.isfile(potential_upload_archive):
                    upload_archive_f = open(potential_upload_archive, 'rb')
                    source_upload_id = os.path.split(os.path.split(potential_upload_archive)[0])[1]
                    upload_name = os.path.basename(potential_upload_archive)
                else:
                    source_upload_id = os.path.split(upload_path)[1]
                    zip_file = zipstream.ZipFile()
                    path_prefix = len(upload_path) + 1
                    for root, _, files in os.walk(upload_path):
                        for file in files:
                            zip_file.write(
                                os.path.join(root, file),
                                os.path.join(root[path_prefix:], file),
                                zipfile.ZIP_DEFLATED)
                    zip_file.write(upload_path)
                    upload_archive_f = IterIO(zip_file)
                    upload_name = '%s.zip' % source_upload_id

            # upload and process the upload file
            upload = self.client.uploads.upload(
                file=upload_archive_f, name=upload_name).response().result
            upload_archive_f.close()

            upload_logger = self.logger.bind(
                source_upload_id=source_upload_id, upload_id=upload.upload_id)

            # grab source metadata
            upload_metadata_calcs = list()
            metadata_dict = dict()
            upload_metadata = dict(calculations=upload_metadata_calcs)
            for source_calc in SourceCalc.objects(upload=source_upload_id):
                source_metadata = CalcWithMetadata(**source_calc.metadata)
                source_metadata.upload_id = upload.upload_id
                source_metadata.mainfile = source_calc.mainfile
                source_metadata.pid = source_calc.pid
                source_metadata.__migrated = False
                upload_metadata_calcs.append(source_metadata)
                metadata_dict[source_calc.mainfile] = source_metadata

            report = utils.POPO()
            report.total_source_calcs = len(metadata_dict)
            report.failed_calcs = 0
            report.migrated_calcs = 0
            report.calcs_with_diffs = 0
            report.new_calcs = 0
            report.missing_calcs = 0

            # wait for complete upload
            while upload.tasks_running:
                upload = self.client.uploads.get_upload(upload_id=upload.upload_id).response().result
                time.sleep(0.1)

            if upload.tasks_status == FAILURE:
                error = 'failed to process upload'
                report.missing_calcs = report.total_source_calcs
                report.total_calcs = 0
                upload_logger.error(error, process_errors=upload.errors)
                yield report
                continue
            else:
                report.total_calcs = upload.calcs.pagination.total

            # verify upload
            for page in range(1, math.ceil(report.total_calcs / 100) + 1):
                upload = self.client.uploads.get_upload(
                    upload_id=upload.upload_id, per_page=100, page=page,
                    order_by='mainfile').response().result

                for calc_proc in upload.calcs.results:
                    calc_logger = upload_logger.bind(
                        calc_id=calc_proc.calc_id,
                        mainfile=calc_proc.mainfile)

                    source_calc = metadata_dict.get(calc_proc.mainfile, None)
                    if calc_proc.tasks_status == SUCCESS:
                        if source_calc is None:
                            calc_logger.info('processed a calc that has no source')
                            report.new_calcs += 1
                            continue
                        else:
                            source_calc.__migrated = True
                            report.migrated_calcs += 1

                        if not self._validate(
                                upload.upload_id, calc_proc.calc_id, source_calc, calc_logger):
                            report.calcs_with_diffs += 1
                    else:
                        report.failed_calcs += 1
                        calc_logger.error(
                            'could not process a calc', process_errors=calc_proc.errors)
                        continue

            for source_calc in upload_metadata_calcs:
                if source_calc.__migrated is False:
                    report.missing_calcs += 1
                    upload_logger.info(
                        'no match or processed calc for source calc',
                        mainfile=source_calc.mainfile)

            # publish upload
            upload_metadata['calculations'] = [
                self._to_api_metadata(calc)
                for calc in upload_metadata['calculations'] if calc.__migrated]

            if report.total_calcs > report.failed_calcs:
                upload = self.client.uploads.exec_upload_operation(
                    upload_id=upload.upload_id,
                    payload=dict(operation='publish', metadata=upload_metadata)
                ).response().result

                while upload.process_running:
                    try:
                        upload = self.client.uploads.get_upload(
                            upload_id=upload.upload_id).response().result
                        time.sleep(0.1)
                    except HTTPNotFound:
                        # the proc upload will be deleted by the publish operation
                        break

            # report
            upload_logger.info('migrated upload', **report)
            yield report

    def _to_api_metadata(self, source: CalcWithMetadata) -> dict:
        """ Transforms to a dict that fullfils the API's uploade metadata model. """
        return dict(
            _upload_time=source.upload_time,
            _uploader=source.uploader['id'],
            _pid=source.pid,
            references=[ref['value'] for ref in source.references],
            datasets=[dict(
                id=ds['id'],
                _doi=ds.get('doi', {'value': None})['value'],
                _name=ds.get('name', None)) for ds in source.datasets],
            mainfile=source.mainfile,
            with_embargo=source.with_embargo,
            comment=source.comment,
            coauthors=list(user['id'] for user in source.coauthors),
            shared_with=list(user['id'] for user in source.shared_with)
        )

    def index(self, *args, **kwargs):
        """ see :func:`SourceCalc.index` """
        return SourceCalc.index(self.source, *args, **kwargs)
