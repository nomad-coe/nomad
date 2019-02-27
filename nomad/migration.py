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

from typing import Generator, Tuple, List, Iterable, Any, IO
import os.path
import zipstream
import zipfile
import math
from mongoengine import Document, IntField, StringField, DictField, ListField
import time
from bravado.exception import HTTPNotFound, HTTPBadRequest
import glob
import os
import runstats
import io

from nomad import utils, infrastructure, config
from nomad.coe_repo import User, Calc, LoginException
from nomad.datamodel import CalcWithMetadata
from nomad.processing import FAILURE, SUCCESS


default_pid_prefix = 7000000
""" The default pid prefix for new non migrated calculations """

max_package_size = 16 * 1024 * 1024 * 1024  # 16 GB
""" The maximum size of a package that will be used as an upload on nomad@FAIRDI """
use_stats_for_filestats_threshold = 1024


def iterable_to_stream(iterable, buffer_size=io.DEFAULT_BUFFER_SIZE):
    """
    Lets you use an iterable (e.g. a generator) that yields bytestrings as a read-only
    input stream.

    The stream implements Python 3's newer I/O API (available in Python 2's io module).
    For efficiency, the stream is buffered.
    """
    class IterStream(io.RawIOBase):
        def __init__(self):
            self.leftover = None
            self.iterator = iter(iterable)

        def readable(self):
            return True

        def readinto(self, b):
            requested_len = len(b)  # We're supposed to return at most this much
            while True:
                try:
                    chunk = next(self.iterator)
                except StopIteration:
                    if len(self.leftover) == 0:
                        return 0  # indicate EOF
                    chunk = self.leftover
                output, self.leftover = chunk[:requested_len], chunk[requested_len:]
                len_output = len(output)
                if len_output == 0:
                    continue  # do not prematurely indicate EOF
                b[:len_output] = output
                return len_output

    return io.BufferedReader(IterStream(), buffer_size=buffer_size)


class Package(Document):
    """
    A Package represents split origin NOMAD CoE uploads. We use packages as uploads
    in nomad@FAIRDI. Some of the uploads in nomad are very big (alfow lib) and need
    to be split down to yield practical (i.e. for mirrors) upload sizes. Therefore,
    uploads are split over multiple packages if one upload gets to large. A package
    always contains full directories of files to preserve *mainfile* *aux* file relations.
    """

    package_id = StringField(primary_key=True)
    """ A random UUID for the package. Could serve later is target upload id."""
    filenames = ListField(StringField(), default=[])
    """ The files in the package relative to the upload path """
    upload_path = StringField()
    """ The absolute path of the source upload """
    upload_id = StringField()
    """ The source upload_id. There might be multiple packages per upload (this is the point). """
    restricted = IntField()
    """ The restricted in month, 0 for unrestricted """
    size = IntField()
    """ The sum of all file sizes """

    meta = dict(indexes=['upload_id'])

    def open_package_upload_file(self) -> IO:
        """ Creates a streaming zip file from the files of this package. """
        zip_file = zipstream.ZipFile(compression=zipstream.ZIP_STORED, allowZip64=True)
        for filename in self.filenames:
            filepath = os.path.join(self.upload_path, filename)
            zip_file.write(filepath, filename)

        return iterable_to_stream(zip_file)  # type: ignore

    def create_package_upload_file(self) -> str:
        """  Creates a zip file for the package in tmp and returns its path. """
        upload_filepath = os.path.join(config.fs.nomad_tmp, '%s.zip' % self.package_id)
        if not os.path.isfile(upload_filepath):
            with zipfile.ZipFile(
                    upload_filepath, 'w',
                    compression=zipfile.ZIP_STORED, allowZip64=True) as zip_file:
                for filename in self.filenames:
                    filepath = os.path.join(self.upload_path, filename)
                    zip_file.write(filepath, filename)

        return upload_filepath

    @classmethod
    def index(cls, *upload_paths):
        """
        Creates Package objects for the given uploads in nomad. The given uploads are
        supposed to be path to the extracted upload directories. If the upload is already
        in the db, the upload is skipped entirely.
        """
        logger = utils.get_logger(__name__)

        for upload_path in upload_paths:
            try:
                stats = runstats.Statistics()
                upload_path = os.path.abspath(upload_path)
                upload_id = os.path.basename(upload_path)
                if cls.objects(upload_id=upload_id).first() is not None:
                    logger.info('upload already exists, skip', upload_id=upload_id)
                    continue

                restrict_files = glob.glob(os.path.join(upload_path, 'RESTRICTED_*'))
                month = 0
                for restrict_file in restrict_files:
                    restrict_file = os.path.basename(restrict_file)
                    try:
                        new_month = int(restrict_file[len('RESTRICTED_'):])
                        if new_month > month:
                            month = new_month
                    except Exception:
                        month = 36
                if month > 36:
                    month = 36
                restricted = month

                def create_package():
                    cls.timer = time.time()
                    package = Package(
                        package_id=utils.create_uuid(),
                        upload_id=upload_id,
                        upload_path=upload_path,
                        restricted=restricted,
                        size=0)
                    return package

                def save_package(package):
                    if len(package.filenames) == 0:
                        return

                    if package.size > max_package_size:
                        # a single directory seems to big for a package
                        logger.error(
                            'directory exceeds max package size', directory=upload_path, size=package.size)

                    package.save()
                    logger.info(
                        'created package',
                        size=package.size,
                        files=len(package.filenames),
                        package_id=package.package_id,
                        exec_time=time.time() - cls.timer,
                        upload_id=package.upload_id)

                package = create_package()

                for root, _, files in os.walk(upload_path):
                    directory_filenames: List[str] = []
                    directory_size = 0

                    if len(files) == 0:
                        continue

                    for file in files:
                        filepath = os.path.join(root, file)
                        filename = filepath[len(upload_path) + 1:]
                        directory_filenames.append(filename)
                        # getting file stats is pretty expensive with gpfs
                        # if an upload has more then 1000 files, its pretty likely that
                        # size patterns repeat ... goood enough
                        if len(stats) < use_stats_for_filestats_threshold:
                            filesize = os.path.getsize(filepath)
                            stats.push(filesize)
                        else:
                            filesize = stats.mean()
                        directory_size += filesize

                    if (package.size + directory_size) > max_package_size and package.size > 0:
                        save_package(package)
                        package = create_package()

                    for filename in directory_filenames:
                        package.filenames.append(filename)
                    package.size += directory_size

                    logger.debug('packaged directory', directory=root, size=directory_size)

                save_package(package)

                logger.info('completed upload', directory=upload_path, upload_id=upload_id)
            except Exception as e:
                logger.error(
                    'could create packae from upload',
                    upload_path=upload_path, upload_id=upload_id, exc_info=e)
                continue


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

    def copy_users(self):
        """ Copy all users. """
        for source_user in self.source.query(User).all():
            if source_user.user_id <= 2:
                # skip first two users to keep example users
                # they probably are either already the example users, or [root, Evgeny]
                continue

            create_user_payload = dict(
                user_id=source_user.user_id,
                email=source_user.email,
                first_name=source_user.first_name,
                last_name=source_user.last_name,
                password=source_user.password,
                created=source_user.created
            )

            try:
                create_user_payload.update(token=source_user.token)
            except LoginException:
                pass

            if source_user.affiliation is not None:
                create_user_payload.update(affiliation=dict(
                    name=source_user.affiliation.name,
                    address=source_user.affiliation.address))

            try:
                self.client.auth.create_user(payload=create_user_payload).response()
                self.logger.info('copied user', user_id=source_user.user_id)
            except HTTPBadRequest as e:
                self.logger.error('could not create user due to bad data', exc_info=e, user_id=source_user.user_id)

    def _to_comparable_list(self, list):
        for item in list:
            if isinstance(item, dict):
                for key in item.keys():
                    if key.endswith('id'):
                        yield item.get(key)
            else:
                yield item

    def _validate(self, repo_calc: dict, source_calc: CalcWithMetadata, logger) -> bool:
        """
        Validates the given processed calculation, assuming that the data in the given
        source_calc is correct.

        Returns:
            False, if the calculation differs from the source calc.
        """
        keys_to_validate = [
            'atoms', 'basis_set', 'xc_functional', 'system', 'crystal_system',
            'spacegroup', 'code_name', 'code_version']

        is_valid = True
        for key, target_value in repo_calc.items():
            if key not in keys_to_validate:
                continue

            source_value = getattr(source_calc, key, None)

            def check_mismatch() -> bool:
                # some exceptions
                if source_value == '3d' and target_value == 'bulk':
                    return True

                logger.info(
                    'source target missmatch', quantity=key,
                    source_value=source_value, target_value=target_value)
                return False

            if source_value is None and target_value is not None:
                continue

            if target_value is None and source_value is not None:
                is_valid &= check_mismatch()

            if isinstance(target_value, list):
                source_list = list(self._to_comparable_list(source_value))
                target_list = list(self._to_comparable_list(target_value))
                if len(set(source_list).intersection(target_list)) != len(target_list):
                    is_valid &= check_mismatch()
                continue

            if isinstance(source_value, str):
                source_value = source_value.lower()
                target_value = str(target_value).lower()

            if source_value != target_value:
                is_valid &= check_mismatch()

        return is_valid

    def _packages(
            self, source_upload_path: str,
            create: bool = False) -> Iterable[Tuple[Package, str]]:

        logger = self.logger.bind(source_upload_path=source_upload_path)

        if os.path.isfile(source_upload_path):
            # assume its a path to an archive files
            logger.error('currently no support for migrating archive files')
            return

        source_upload_id = os.path.basename(source_upload_path)
        logger = logger.bind(source_upload_id=source_upload_id)

        package_query = Package.objects(upload_id=source_upload_id)

        if package_query.count() == 0:
            if create:
                Package.index(source_upload_path)
                package_query = Package.objects(upload_id=source_upload_id)
                if package_query.count() == 0:
                    logger.error('no package exists, even after indexing')
                    return
            else:
                logger.error('no package exists for upload')
                return

        logger.debug('identified packages for source upload', n_packages=package_query.count())
        for package in package_query:
            yield package, source_upload_id,

    def migrate(
            self, upload_path, prefix: int = default_pid_prefix,
            create_packages: bool = False, local: bool = False):
        """
        Migrate the given uploads.

        It takes paths to extracted uploads as arguments.

        Requires :class:`Package` instances for the given upload paths. Those will
        be created, if they do not already exists. The packages determine the uploads
        for the target infrastructure.

        Requires a build :func:`index` to look for existing data in the source db. This
        will be used to add user (and other, PID, ...) metadata and validate calculations.

        Uses PIDs of identified old calculations. Will create new PIDs for previously
        unknown uploads. New PIDs will be choosed from a `prefix++` range of ints.

        Arguments:
            upload_path: A filepath to the upload directory.
            prefix: The PID prefix that should be used for new non migrated calcualtions.
            create_packages: If True, will create non existing packages.
                Will skip with errors otherwise.
            local: Instead of streaming an upload, create a local file and use local_path on the upload.

        Returns: Yields a dictionary with status and statistics for each given upload.
        """
        from nomad.client import stream_upload_with_client

        if prefix is not None:
            self.logger.info('set pid prefix', pid_prefix=prefix)
            self.client.admin.exec_pidprefix_command(payload=dict(prefix=prefix)).response()

        for package, source_upload_id in self._packages(upload_path, create=create_packages):
            package_id = package.package_id

            logger = self.logger.bind(
                package_id=package_id, source_upload_id=source_upload_id, upload_path=upload_path)

            # upload and process the upload file
            with utils.timer(logger, 'upload completed'):
                try:
                    if local:
                        upload_filepath = package.create_package_upload_file()
                        self.logger.debug('created package upload file')
                        upload = self.client.uploads.upload(
                            name=package_id, local_path=upload_filepath).response().result
                    else:
                        upload_f = package.open_package_upload_file()
                        self.logger.debug('opened package upload file')
                        upload = stream_upload_with_client(self.client, upload_f, name=package_id)
                except Exception as e:
                    self.logger.error('could not upload package', exc_info=e)
                    continue

            logger.debug('package file uploaded')

            logger = logger.bind(
                source_upload_id=source_upload_id, upload_id=upload.upload_id)

            # grab source metadata
            upload_metadata_calcs: List[Any] = list()
            metadata_dict = dict()
            upload_metadata = dict(calculations=upload_metadata_calcs)
            for source_calc in SourceCalc.objects(upload=source_upload_id):
                source_metadata = CalcWithMetadata(**source_calc.metadata)
                source_metadata.upload_id = upload.upload_id
                source_metadata.mainfile = source_calc.mainfile
                source_metadata.pid = source_calc.pid
                source_metadata.__migrated = False  # type: ignore
                upload_metadata_calcs.append(source_metadata)
                metadata_dict[source_calc.mainfile] = source_metadata
            logger.debug('loaded source metadata', calcs=len(upload_metadata_calcs))

            report = utils.POPO()
            report.total_source_calcs = len(metadata_dict)
            report.failed_calcs = 0
            report.migrated_calcs = 0
            report.calcs_with_diffs = 0
            report.new_calcs = 0
            report.missing_calcs = 0

            # wait for complete upload
            with utils.timer(logger, 'upload processing completed'):
                while upload.tasks_running:
                    sleep_time = 0.1
                    upload = self.client.uploads.get_upload(upload_id=upload.upload_id).response().result
                    time.sleep(sleep_time)
                    if sleep_time < 60:
                        sleep_time *= 2

            if upload.tasks_status == FAILURE:
                error = 'failed to process upload'
                report.missing_calcs = report.total_source_calcs
                report.total_calcs = 0
                logger.error(error, process_errors=upload.errors)
                yield report
                continue
            else:
                report.total_calcs = upload.calcs.pagination.total

            timer = utils.timer(logger, 'validation completed')
            timer.__enter__()  # type: ignore, pylint: disable=no-member
            # verify upload: check for processing errors
            per_page = 200
            for page in range(1, math.ceil(report.total_calcs / per_page) + 1):
                upload = self.client.uploads.get_upload(
                    upload_id=upload.upload_id, per_page=per_page, page=page,
                    order_by='mainfile').response().result

                for calc_proc in upload.calcs.results:
                    calc_logger = logger.bind(
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

                    else:
                        report.failed_calcs += 1
                        calc_logger.error(
                            'could not process a calc', process_errors=calc_proc.errors)
                        continue

            # very upload: calc data
            for page in range(1, math.ceil(report.total_calcs / per_page) + 1):
                search = self.client.repo.search(
                    page=page, per_page=per_page, upload_id=upload.upload_id,
                    order_by='mainfile').response().result
                for calc in search.results:
                    source_calc = metadata_dict.get(calc_proc.mainfile, None)
                    calc_logger = logger.bind(calc_id=calc['calc_id'], mainfile=calc['mainfile'])
                    if not self._validate(calc, source_calc, calc_logger):
                        report.calcs_with_diffs += 1

            for source_calc in upload_metadata_calcs:
                if source_calc.__migrated is False:
                    report.missing_calcs += 1
                    logger.info(
                        'no match or processed calc for source calc',
                        mainfile=source_calc.mainfile)
            timer.__exit__(None, None, None)  # type: ignore, pylint: disable=no-member

            # publish upload
            with utils.timer(logger, 'upload published'):
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
            logger.info('migrated upload', **report)
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
            coauthors=list(int(user['id']) for user in source.coauthors),
            shared_with=list(int(user['id']) for user in source.shared_with)
        )

    def index(self, *args, **kwargs):
        """ see :func:`SourceCalc.index` """
        return SourceCalc.index(self.source, *args, **kwargs)

    def package(self, *args, **kwargs):
        """ see :func:`Package.add` """
        return Package.index(*args, **kwargs)
