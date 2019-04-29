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

from typing import Generator, Tuple, List, Iterable, Any, Dict
import multiprocessing
import multiprocessing.pool
import time
import os
import os.path
import zipfile
import tarfile
import math
from mongoengine import Document, IntField, StringField, DictField
import datetime
from bravado.exception import HTTPNotFound, HTTPBadRequest, HTTPGatewayTimeout
import os
import runstats
import io
import threading
from contextlib import contextmanager
import shutil
import json
import random

from nomad import utils, infrastructure, files, config
from nomad.coe_repo import User, Calc, LoginException
from nomad.datamodel import CalcWithMetadata
from nomad.processing import FAILURE


default_pid_prefix = 7000000
""" The default pid prefix for new non migrated calculations """

max_package_size = 32 * 1024 * 1024 * 1024  # 32 GB
""" The maximum size of a package that will be used as an upload on nomad@FAIRDI """
use_stats_for_filestats_threshold = 1024

default_comment = 'entry with unknown provernance'
default_uploader = dict(id=1)


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


Directory = Tuple[List[str], str, int]


def create_package_zip(
        upload_id: str, upload_path: str, package_id: str, package_path: str, compress: bool,
        package_filepaths: Iterable[str]) -> None:
    logger = utils.get_logger(
        __name__, source_upload_id=upload_id,
        source_upload_path=upload_path, package_id=package_id)

    package_zip = zipfile.ZipFile(
        package_path, 'w',
        compression=zipfile.ZIP_DEFLATED if compress else zipfile.ZIP_STORED)

    try:
        for filepath in package_filepaths:
            package_zip.write(filepath, filepath[len(upload_path):])
    except Exception as e:
        logger.error('could not write file to zip', filepath=filepath, exc_info=e)
    finally:
        package_zip.close()

    logger.info('package zip created')


class Package(Document):
    """
    A Package represents split origin NOMAD CoE uploads. We use packages as uploads
    in nomad@FAIRDI. Some of the uploads in nomad are very big (alfow lib) and need
    to be split down to yield practical (i.e. for mirrors) upload sizes. Therefore,
    uploads are split over multiple packages if one upload gets to large. A package
    always contains full directories of files to preserve *mainfile* *aux* file relations.
    Package have a package entry in mongo and a .zip file with the raw data.
    """

    package_id = StringField(primary_key=True)
    """ A random UUID for the package. Could serve later its target upload id."""
    package_path = StringField(required=True)
    """ The path to the package .zip file. """
    upload_path = StringField(required=True)
    """ The absolute path of the source upload """
    upload_id = StringField(required=True)
    """ The source upload_id. There might be multiple packages per upload (this is the point). """
    restricted = IntField(default=0)
    """ The restricted in month, 0 for unrestricted. """
    size = IntField()
    """ The sum of all file sizes. """
    files = IntField()
    """ The number of files. """
    packages = IntField(default=-1)
    """ The number of packages in the same upload. """

    migration_version = IntField(default=-1)
    """ The version of the last successful migration of this package """
    report = DictField()
    """ The report of the last successful migration of this package """

    migration_failure = StringField()
    """ String that describe the cause for last failed migration attempt """

    meta = dict(indexes=['upload_id', 'migration_version'])

    @classmethod
    def aggregate_reports(cls, migration_version: str = None) -> 'Report':
        """
        Returns an aggregated report over all package reports of the given migration version,
        or all packages.
        """
        if migration_version is not None:
            query = cls.objects(migration_version__gte=migration_version, report__exists=True)
        else:
            query = cls.objects(report__exists=True)

        report = Report()
        for package in query:
            report.add(Report(package.report))

        return report

    @classmethod
    def get_packages(
            cls, upload_path: str, target_dir: str, create: bool = False,
            compress: bool = False, parallel: int = 1) -> Iterable['Package']:
        """
        Will get packages for the given upload_path. Creates the package zip files and
        package index entries if ``create`` is True. But, either will only be created if
        it does not already exist. Yields the Package objects.
        """
        upload_id = os.path.basename(upload_path)
        logger = utils.get_logger(__name__, source_upload_path=upload_path, source_upload_id=upload_id)

        if not os.path.isdir(upload_path):
            logger.error('upload path is not a directory')
            return []

        upload_directory = files.DirectoryObject(target_dir, upload_id, create=True, prefix=True)
        restricted = 0

        # The packages number is written after all packages of an upload have been created.
        # this should allow to abort mid upload packaging and continue later by removing
        # all started packages first.
        is_packaged = cls.objects(upload_id=upload_id, packages__ne=-1).count() != 0

        async_results: List[multiprocessing.pool.AsyncResult] = []
        pool = multiprocessing.Pool(parallel)
        pool.__enter__()

        if not is_packaged:
            if not create:
                return None

            cls.objects(upload_id=upload_id).delete()

            def create_package_entry():
                package_id = utils.create_uuid()
                return Package(
                    package_id=package_id,
                    upload_path=upload_path,
                    upload_id=upload_id,
                    package_path=upload_directory.join_file('%s.zip' % package_id).os_path)

            def close_package(package_size: int, package_filepaths: List[str]):
                package_entry_to_close = package_entry

                def save_package_entry(*args) -> None:
                    package_entry_to_close.size = package_size
                    package_entry_to_close.files = len(package_filepaths)
                    package_entry_to_close.save()

                    logger.debug(
                        'created package',
                        package_id=package_entry.package_id, size=package_size)

                def handle_package_error(*args) -> None:
                    logger.error(
                        'could not create package zip due to unexpected exception',
                        exc_info=args[0])

                while len(async_results) > parallel:
                    async_results[:] = [
                        async_result for async_result in async_results
                        if not async_result.ready()]
                    time.sleep(0.1)

                async_result = pool.apply_async(
                    create_package_zip,
                    args=(
                        upload_id, upload_path, package_entry.package_id,
                        package_entry.package_path, compress, package_filepaths),
                    callback=save_package_entry, error_callback=handle_package_error)

                async_results.append(async_result)

            package_entry = create_package_entry()
            package_size = 0
            package_filepaths = []
            with cls.upload_iterator(upload_path) as directory:
                for filepaths, parent_directory, size in directory:
                    for filepath in filepaths:
                        basepath = os.path.basename(filepath)
                        if basepath.startswith('RESTRICTED'):
                            restricted = 36
                            try:
                                restricted = min(36, int(basepath[len('RESTRICTED_'):]))
                            except Exception:
                                pass

                        package_filepaths.append(os.path.join(parent_directory, filepath))

                    if size > max_package_size:
                        logger.warn(
                            'directory exceeds max package size',
                            package_id=package_entry.package_id, size=package_size)

                    package_size += size
                    if package_size > max_package_size:
                        close_package(package_size, package_filepaths)
                        package_size, package_filepaths = 0, []
                        package_entry = create_package_entry()

                if len(package_filepaths) > 0:
                    close_package(package_size, package_filepaths)

            # wait for all zip processes to complete
            while not all(async_result.ready() for async_result in async_results):
                time.sleep(0.1)

            pool.__exit__(None, None, None)

            package_query = cls.objects(upload_id=upload_id)
            package_query.update(restricted=restricted, packages=package_query.count())
            logger.debug(
                'packaged upload', source_upload_id=upload_id, source_upload_path=upload_path,
                restricted=restricted)

            return package_query
        else:
            return cls.objects(upload_id=upload_id)

    @classmethod
    @contextmanager
    def upload_iterator(cls, upload_path: str) -> Generator[Generator[Directory, None, None], None, None]:
        """
        A contextmanager that opens the given upload and provides a generator for
        directories. Directories are tuple of an iterable of upload relative filepaths
        and the directory size.
        """
        potential_archive_path = os.path.join(upload_path, 'archive.tar.gz')
        if os.path.isfile(potential_archive_path):
            with cls.extracted_archive(potential_archive_path) as extracted_archive:
                yield cls.iterate_upload_directory(extracted_archive)
        else:
            yield cls.iterate_upload_directory(upload_path)

    @classmethod
    def iterate_upload_directory(cls, upload_path) -> Generator[Directory, None, None]:
        """
        Interprets the given upload path as a directory. Files path are given as upload
        path relative paths.
        """
        stats = runstats.Statistics()
        for root, _, files in os.walk(upload_path):
            directory_filepaths: List[str] = []
            directory_size = 0

            if len(files) == 0:
                continue

            if len(files) < 20 and any(file.endswith('.tar.gz') for file in files):
                # TODO the OQMD case, files are managed as bunch of .tar.gz files
                for file in files:
                    archive_path = os.path.join(root, file)
                    prefix = os.path.dirname(archive_path)[len(upload_path) + 1:]
                    with cls.extracted_archive(archive_path) as extracted_archive:
                        for paths, _, size in cls.iterate_upload_directory(extracted_archive):
                            yield [os.path.join(prefix, path) for path in paths], upload_path, size
                continue

            for file in files:
                filepath = os.path.join(root, file)
                filename = filepath[len(upload_path) + 1:]
                directory_filepaths.append(filename)
                # getting file stats is pretty expensive with gpfs
                # if an upload has more then 1000 files, its pretty likely that
                # size patterns repeat ... goood enough
                if len(stats) < use_stats_for_filestats_threshold:
                    try:
                        filesize = os.path.getsize(filepath)
                    except Exception:
                        # if there are individual files that cannot be accessed, we fully ignore them
                        # they are most likely just broken links
                        pass

                    stats.push(filesize)
                else:
                    filesize = stats.mean()
                directory_size += filesize

            yield directory_filepaths, upload_path, directory_size

    @classmethod
    @contextmanager
    def extracted_archive(cls, archive_path: str) -> Generator[str, None, None]:
        """
        Temporarily extracts the given archive and returns directory with the extracted
        data.
        """
        tmp_directory = os.path.join(config.fs.local_tmp, utils.create_uuid())
        os.mkdir(tmp_directory)

        with tarfile.TarFile.open(archive_path) as tar_file:
            tar_file.extractall(tmp_directory)
        # try to fix permissions, do not care if command fails
        os.system('chmod -Rf 0755 %s/*' % tmp_directory)

        yield tmp_directory

        shutil.rmtree(tmp_directory)


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

    migration_version = IntField(default=-1)

    extracted_prefix = '$EXTRACTED/'
    sites = ['/data/nomad/extracted/', '/nomad/repository/extracted/']
    prefixes = [extracted_prefix] + sites

    meta = dict(indexes=['upload', 'mainfile', 'migration_version'])

    _dataset_cache: dict = {}

    @staticmethod
    def missing(use_cache=False):
        """
        Produces data about non migrated calcs
        """
        tmp_data_path = '/tmp/nomad_migration_missing.json'
        if os.path.exists(tmp_data_path) and use_cache:
            with open(tmp_data_path, 'rt') as f:
                data = utils.POPO(**json.load(f))
        else:
            data = utils.POPO(step=0)

        try:
            # get source_uploads that have non migrated calcs
            if data.step < 1 or not use_cache:
                import re
                data.source_uploads = SourceCalc._get_collection() \
                    .find({'migration_version': {'$lt': 0}, 'mainfile': {'$not': re.compile(r'^aflowlib_data.*')}}) \
                    .distinct('upload')
                data.step = 1

            if data.step < 2 or not use_cache:
                source_uploads = []
                data.packages = utils.POPO()
                data.uploads_with_no_package = []
                for source_upload in data.source_uploads:
                    package = Package.objects(upload_id=source_upload).first()
                    if package is None:
                        data.uploads_with_no_package.append(source_upload)
                    else:
                        calcs = SourceCalc.objects(upload=source_upload).count()
                        packages = Package.objects(upload_id=source_upload).count()
                        source_uploads.append(dict(
                            id=source_upload, package_count=packages,
                            packages=package.packages, calcs=calcs,
                            path=package.upload_path))
                        source_uploads = sorted(source_uploads, key=lambda k: k['calcs'], reverse=True)
                data.source_uploads = source_uploads
                data.step = 2
        finally:
            with open(tmp_data_path, 'wt') as f:
                json.dump(data, f)

        return data

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

                    # this is taken from metadata.location and has inconsistent directory prefix,
                    # but is more accurate than taking the first file as mainfile, which
                    # also is sometimes not the actual mainfile.
                    if calc.mainfile is not None:
                        calc_mainfile = os.path.basename(calc.mainfile)
                        if calc_mainfile != os.path.basename(source_calc.mainfile):
                            source_calc.mainfile = os.path.join(
                                os.path.dirname(source_calc.mainfile), calc_mainfile)

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


NO_PROCESSED_CALCS = 0
FAILED_PROCESSING = 1
FAILED_PUBLISH = 2


class NomadCOEMigration:
    """
    Drives a migration from the NOMAD coe repository db to nomad@FAIRDI.

    Arguments:
        migration_version: The migration version. Only packages/calculations with
            no migration version or a lower migration version are migrated.
        package_directory: The directory that packages are/get stored in.
        compress_packages: True to use compression on package creation.
        threads: Number of threads to run migration in parallel.
        quiet: Prints stats if not quiet
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
            migration_version: int = 0,
            package_directory: str = None,
            compress_packages: bool = False,
            threads: int = 1, quiet: bool = False) -> None:
        self.logger = utils.get_logger(__name__, migration_version=migration_version)

        self.migration_version = migration_version
        self.package_directory = package_directory if package_directory is not None else config.fs.migration_packages
        self.compress_packages = compress_packages
        self._client = None
        self._threads = threads
        self._quiet = quiet

        self.source = infrastructure.repository_db

    @property
    def client(self):
        if self._client is None:
            from nomad.client import create_client
            self._client = create_client()

        return self._client

    def report(self):
        """ returns an aggregated report over all prior migrated packages """
        return Package.aggregate_reports(migration_version=self.migration_version)

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

    expected_differences = {
        '0d': 'molecule / cluster',
        '3d': 'bulk',
        '2d': '2d / surface',
        '+u': 'gga'
    }

    def validate(self, repo_calc: dict, source_calc: CalcWithMetadata, logger) -> bool:
        """
        Validates the given processed calculation, assuming that the data in the given
        source_calc is correct.

        Returns:
            False, if the calculation differs from the source calc.
        """
        keys_to_validate = [
            'atoms', 'basis_set', 'xc_functional', 'system', 'crystal_system',
            'spacegroup', 'code_name', 'code_version']

        def to_comparable_list(list):
            for item in list:
                if isinstance(item, dict):
                    for key in item.keys():
                        if key.endswith('id'):
                            yield item.get(key)
                else:
                    yield item

        is_valid = True
        for key, target_value in repo_calc.items():
            if key not in keys_to_validate:
                continue

            source_value = getattr(source_calc, key, None)

            def check_mismatch() -> bool:
                # some exceptions
                if isinstance(source_value, str) and \
                        source_value in NomadCOEMigration.expected_differences and \
                        target_value == NomadCOEMigration.expected_differences.get(source_value):
                    return True

                logger.info(
                    'source target missmatch', quantity=key,
                    source_value=source_value, target_value=target_value,
                    value_diff='%s->%s' % (str(source_value), str(target_value)))
                return False

            if source_value is None and target_value is not None:
                continue

            if target_value is None and source_value is not None:
                is_valid &= check_mismatch()

            if isinstance(target_value, list):
                source_list = list(to_comparable_list(source_value))
                target_list = list(to_comparable_list(target_value))
                if len(source_list) != len(target_list):
                    is_valid &= check_mismatch()
                elif any(a != b for a, b in zip(source_list, target_list)):
                    is_valid &= check_mismatch()
                continue

            if isinstance(source_value, str):
                source_value = source_value.lower()
                target_value = str(target_value).lower()

            if source_value != target_value:
                is_valid &= check_mismatch()

        return is_valid

    def surrogate_metadata(self, source: CalcWithMetadata):
        """
        Compute metadata from the given metadata that can be used for new calcs of the
        same upload.
        """
        return CalcWithMetadata(
            uploader=source.uploader,
            with_embargo=source.with_embargo,
            upload_time=source.upload_time,
            coauthors=source.coauthors,
            shared_with=source.shared_with,
            comment=source.comment,
            references=source.references,
            datasets=source.datasets)

    def set_pid_prefix(self, prefix: int = default_pid_prefix):
        """
        Sets the repo db pid counter to the given values. Allows to create new calcs
        without interfering with migration calcs with already existing PIDs.
        """
        self.logger.info('set pid prefix', pid_prefix=prefix)
        self.client.admin.exec_pidprefix_command(payload=dict(prefix=prefix)).response()

    def call_api(self, operation: str, *args, **kwargs) -> Any:
        """
        Calls nomad via the bravado client. It deals with a very busy nomad and catches,
        backsoff, and retries on gateway timouts. It also circumvents bravados/jsonschemas
        thread safety issues using a global lock on client usage.

        Arguments:
            operation: Comma separated string of api, endpoint, operation,
                e.g. 'uploads.get_upload'.
        """
        op_path = operation.split('.')
        op = self.client
        for op_path_segment in op_path:
            op = getattr(op, op_path_segment)

        sleep = utils.SleepTimeBackoff()
        while True:
            try:
                NomadCOEMigration._client_lock.acquire(blocking=True)
                return op(*args, **kwargs).response().result
            except HTTPGatewayTimeout:
                sleep()
            except Exception as e:
                raise e
            finally:
                NomadCOEMigration._client_lock.release()

    def migrate(
            self, *args, delete_failed: str = '',
            create_packages: bool = False, only_republish: bool = False,
            wait: int = 0) -> utils.POPO:
        """
        Migrate the given uploads.

        It takes paths to extracted uploads as arguments.

        Requires :class:`Package` instances for the given upload paths. Those will
        be created, if they do not already exists. The packages determine the uploads
        for the target infrastructure.

        Requires a build :func:`index` to look for existing data in the source db. This
        will be used to add user (and other, PID, ...) metadata and validate calculations.

        Uses PIDs of identified old calculations. Will create new PIDs for previously
        unknown uploads. See :func:`set_pid_prefix` on how to avoid conflicts.

        Arguments:
            upload_path: A filepath to the upload directory.
            delete_failed: String from ``N``, ``U``, ``P`` to determine that uploads with
                no processed calcs (N), failed upload processing (U), or failed publish
                operation (P) should be deleted after the migration attempt.
            create_packages: If True, it will attempt to create upload packages if they
                do not exists.
            only_republish: If the package exists and is published, it will be republished.
                Nothing else. Useful to reindex/recreate coe repo, etc.
            offset: Will add a random sleep before migrating each package between 0 and
                ``wait`` seconds.

        Returns: Dictionary with statistics on the migration.
        """

        cv = threading.Condition()
        overall_report = Report()
        threads = []

        def print_report():
            if not self._quiet:
                print(overall_report)

        def migrate_package(package: Package):
            logger = self.logger.bind(
                package_id=package.package_id, source_upload_id=package.upload_id)

            if package.migration_version is not None and package.migration_version >= self.migration_version:
                if only_republish:
                    self.republish_package(package)
                else:
                    self.logger.info(
                        'package already migrated, skip it',
                        package_id=package.package_id, source_upload_id=package.upload_id)

                package_report = package.report
                overall_report.skipped_packages += 1
            else:
                try:
                    if wait > 0:
                        self.logger.info('wait for a random amount of time')
                        time.sleep(random.randint(0, wait))

                    package_report = self.migrate_package(package, delete_failed=delete_failed)

                except Exception as e:
                    package_report = Report()
                    package_report.failed_packages = 1
                    event = 'unexpected exception while migrating packages'
                    package.migration_failure = event + ': ' + str(e)
                    logger.error(event, exc_info=e)
                finally:
                    package.report = package_report
                    package.migration_version = self.migration_version
                    package.save()

            with cv:
                try:
                    overall_report.add(package_report)

                    migrated_all_packages = all(
                        p.migration_version == self.migration_version
                        for p in Package.objects(upload_id=package.upload_id))

                    if migrated_all_packages:
                        missing_calcs = SourceCalc.objects(
                            upload=package.upload_id, migration_version__ne=self.migration_version).count()
                        total_source_calcs = SourceCalc.objects(upload=package.upload_id).count()

                        overall_report.missing_calcs += missing_calcs
                        overall_report.total_source_calcs += total_source_calcs

                        logger.info('migrated upload')

                    print_report()
                except Exception as e:
                    logger.error('unexpected exception while migrating packages', exc_info=e)

                self._threads += 1
                cv.notify()

        for arg in args:
            packages = Package.get_packages(
                arg, self.package_directory, create=create_packages,
                compress=self.compress_packages)

            if packages is None:
                self.logger.error('there are no packages for upload', upload_source_path=arg)
                continue

            for package in packages:
                with cv:
                    cv.wait_for(lambda: self._threads > 0)
                    self._threads -= 1
                    thread = threading.Thread(target=lambda: migrate_package(package))
                    threads.append(thread)
                    thread.start()

        for thread in threads:
            thread.join()

        return overall_report

    _client_lock = threading.Lock()

    def republish_package(self, package: Package) -> None:

        source_upload_id = package.upload_id
        package_id = package.package_id

        logger = self.logger.bind(package_id=package_id, source_upload_id=source_upload_id)

        uploads = self.call_api('uploads.get_uploads', name=package_id)
        if len(uploads) > 1:
            self.logger.warning('upload name is not unique')
        if len(uploads) == 0:
            self.logger.info('upload does not exist')
            return

        for upload in uploads:
            if not upload.published:
                self.logger.info('upload is not published, therefore cannot re-publish')
                continue

            upload = self.call_api(
                'uploads.exec_upload_operation', upload_id=upload.upload_id,
                payload=dict(operation='publish'))

            sleep = utils.SleepTimeBackoff()
            while upload.process_running:
                upload = self.call_api('uploads.get_upload', upload_id=upload.upload_id)
                sleep()

            if upload.tasks_status == FAILURE:
                event = 'could not re publish upload'
                logger.error(event, process_errors=upload.errors)
            else:
                logger.info('republished upload')

    def migrate_package(self, package: Package, delete_failed: str = '') -> 'Report':
        """ Migrates the given package. For other params see :func:`migrate`. """

        source_upload_id = package.upload_id
        package_id = package.package_id

        logger = self.logger.bind(package_id=package_id, source_upload_id=source_upload_id)
        logger.debug('start to process package')

        report = Report()
        report.total_packages += 1

        # check if the package is already uploaded
        upload = None
        try:
            uploads = self.call_api('uploads.get_uploads', name=package_id)
            if len(uploads) > 1:
                event = 'duplicate upload name'
                package.migration_failure(event)
                report.failed_packages += 1
                return report
            elif len(uploads) == 1:
                upload = uploads[0]

        except Exception as e:
            event = 'could not verify if upload already exists'
            logger.error(event, exc_info=e)
            package.migration_failure(event)
            report.failed_packages += 1
            return report

        # upload and process the upload file
        if upload is None:
            with utils.timer(logger, 'upload completed'):
                try:
                    upload = self.call_api(
                        'uploads.upload', name=package_id, local_path=package.package_path)
                except Exception as e:
                    event = 'could not upload package'
                    logger.error(event, exc_info=e)
                    package.migration_failure = event + ': ' + str(e)
                    report.failed_packages += 1
                    return report
        else:
            self.logger.info('package was already uploaded')
            # get more details than the get_uploads call provided
            upload = self.call_api('uploads.get_upload', upload_id=upload.upload_id)

        logger = logger.bind(
            source_upload_id=source_upload_id, upload_id=upload.upload_id)

        def delete_upload(reason: int):
            delete = \
                (reason == NO_PROCESSED_CALCS and 'N' in delete_failed) or \
                (reason == FAILED_PROCESSING and 'U' in delete_failed) or \
                (reason == FAILED_PUBLISH and 'P' in delete_failed)

            upload_to_delete = upload

            if delete:
                upload_to_delete = self.call_api(
                    'uploads.delete_upload', upload_id=upload_to_delete.upload_id)

                sleep = utils.SleepTimeBackoff()
                while upload_to_delete.process_running:
                    try:
                        upload_to_delete = self.call_api(
                            'uploads.get_upload', upload_id=upload_to_delete.upload_id)
                        sleep()
                    except HTTPNotFound:
                        # the proc upload will be deleted by the delete operation
                        break
                logger.info('deleted upload after migration failure')
            else:
                logger.warning(
                    'will keep upload after migration failure for debugging',
                    reason=reason, delete_failed=delete_failed)

        # grab source calcs, while waiting for upload
        source_calcs = dict()
        surrogate_source_calc_with_metadata = None
        with utils.timer(logger, 'loaded source metadata'):
            with zipfile.ZipFile(package.package_path) as zf:
                for filenames_chunk in utils.chunks(zf.namelist(), 1000):
                    for source_calc in SourceCalc.objects(
                            upload=source_upload_id, mainfile__in=filenames_chunk):

                        source_calc_with_metadata = CalcWithMetadata(**source_calc.metadata)
                        source_calc_with_metadata.pid = source_calc.pid
                        source_calc_with_metadata.mainfile = source_calc.mainfile
                        source_calcs[source_calc.mainfile] = (source_calc, source_calc_with_metadata)

                        # establish a surrogate for new calcs
                        if surrogate_source_calc_with_metadata is None:
                            surrogate_source_calc_with_metadata = \
                                self.surrogate_metadata(source_calc_with_metadata)

        # try to find a surrogate outside the package, if necessary
        if surrogate_source_calc_with_metadata is None:
            source_calc = SourceCalc.objects(upload=source_upload_id).first()
            if source_calc is not None:
                source_calc_with_metadata = CalcWithMetadata(**source_calc.metadata)
                surrogate_source_calc_with_metadata = \
                    self.surrogate_metadata(source_calc_with_metadata)

        # wait for complete upload
        with utils.timer(logger, 'upload processing completed'):
            sleep = utils.SleepTimeBackoff()
            while upload.tasks_running:
                upload = self.call_api('uploads.get_upload', upload_id=upload.upload_id)
                sleep()

        if upload.tasks_status == FAILURE:
            event = 'failed to process upload'
            logger.error(event, process_errors=upload.errors)
            package.migration_failure = event + ': ' + str(upload.errors)
            report.failed_packages += 1
            delete_upload(FAILED_PROCESSING)
            return report
        else:
            report.total_calcs += upload.calcs.pagination.total

        calc_mainfiles = []
        upload_total_calcs = upload.calcs.pagination.total

        # check for processing errors
        with utils.timer(logger, 'checked upload processing'):
            per_page = 10000
            for page in range(1, math.ceil(upload_total_calcs / per_page) + 1):
                upload = self.call_api(
                    'uploads.get_upload', upload_id=upload.upload_id, per_page=per_page,
                    page=page)

                for calc_proc in upload.calcs.results:
                    calc_logger = logger.bind(
                        calc_id=calc_proc.calc_id,
                        mainfile=calc_proc.mainfile)

                    calc_mainfiles.append(calc_proc.mainfile)

                    if calc_proc.tasks_status == FAILURE:
                        report.failed_calcs += 1
                        calc_logger.info(
                            'could not process a calc', process_errors=calc_proc.errors)
                        continue

        # verify upload against source
        calcs_in_search = 0
        with utils.timer(logger, 'verified upload against source calcs'):
            scroll_id = 'first'
            while scroll_id is not None:
                scroll_args: Dict[str, Any] = dict(scroll=True)
                if scroll_id != 'first':
                    scroll_args['scroll_id'] = scroll_id

                search = self.call_api('repo.search', upload_id=upload.upload_id, **scroll_args)

                scroll_id = search.scroll_id

                for calc in search.results:
                    calcs_in_search += 1
                    source_calc, source_calc_with_metadata = source_calcs.get(
                        calc['mainfile'], (None, None))

                    if source_calc is not None:
                        report.migrated_calcs += 1

                        calc_logger = logger.bind(calc_id=calc['calc_id'], mainfile=calc['mainfile'])
                        if calc.get('processed', False):
                            try:
                                if not self.validate(
                                        calc, source_calc_with_metadata, calc_logger):
                                    report.calcs_with_diffs += 1
                            except Exception as e:
                                calc_logger.warning(
                                    'unexpected exception during validation', exc_info=e)
                                report.calcs_with_diffs += 1
                    else:
                        calc_logger.info('processed a calc that has no source')
                        report.new_calcs += 1
                        # guessing the metadata from other calcs in upload/package
                        if surrogate_source_calc_with_metadata is not None:
                            new_calc_with_metadata = CalcWithMetadata(**surrogate_source_calc_with_metadata.to_dict())
                            new_calc_with_metadata.mainfile = calc['mainfile']
                        else:
                            calc_logger.warning('could not determine any metadata for new calc')
                            create_time_epoch = os.path.getctime(package.upload_path)
                            new_calc_with_metadata = CalcWithMetadata(
                                upload_time=datetime.datetime.fromtimestamp(create_time_epoch),
                                with_embargo=package.restricted > 0,
                                comment=default_comment,
                                uploader=default_uploader,
                                mainfile=calc['mainfile'])
                            surrogate_source_calc_with_metadata = new_calc_with_metadata

                        source_calcs[calc['mainfile']] = (None, new_calc_with_metadata)

            if len(calc_mainfiles) != calcs_in_search:
                logger.error('missmatch between processed calcs and calcs found with search')

        # publish upload
        if len(calc_mainfiles) > 0:
            with utils.timer(logger, 'upload published'):
                upload_metadata = dict(with_embargo=(package.restricted > 0))
                upload_metadata['calculations'] = [
                    self._to_api_metadata(source_calc_with_metadata)
                    for _, source_calc_with_metadata in source_calcs.values()]

                upload = self.call_api(
                    'uploads.exec_upload_operation', upload_id=upload.upload_id,
                    payload=dict(operation='publish', metadata=upload_metadata))

                sleep = utils.SleepTimeBackoff()
                while upload.process_running:
                    upload = self.call_api('uploads.get_upload', upload_id=upload.upload_id)
                    sleep()

                if upload.tasks_status == FAILURE:
                    event = 'could not publish upload'
                    logger.error(event, process_errors=upload.errors)
                    report.failed_calcs = report.total_calcs
                    report.migrated_calcs = 0
                    report.calcs_with_diffs = 0
                    report.new_calcs = 0
                    report.failed_packages += 1
                    package.migration_failure = event + ': ' + str(upload.errors)

                    delete_upload(FAILED_PUBLISH)
                    SourceCalc.objects(upload=source_upload_id, mainfile__in=calc_mainfiles) \
                        .update(migration_version=-1)
                else:
                    SourceCalc.objects(upload=source_upload_id, mainfile__in=calc_mainfiles) \
                        .update(migration_version=self.migration_version)
        else:
            delete_upload(NO_PROCESSED_CALCS)
            logger.info('no successful calcs, skip publish')

        logger.info('migrated package', **report)
        return report

    def _to_api_metadata(self, calc_with_metadata: CalcWithMetadata) -> dict:
        """ Transforms to a dict that fullfils the API's uploade metadata model. """

        return dict(
            _upload_time=calc_with_metadata.upload_time,
            _uploader=calc_with_metadata.uploader['id'],
            _pid=calc_with_metadata.pid,
            references=[ref['value'] for ref in calc_with_metadata.references],
            datasets=[dict(
                id=ds['id'],
                _doi=ds.get('doi', {'value': None})['value'],
                _name=ds.get('name', None)) for ds in calc_with_metadata.datasets],
            mainfile=calc_with_metadata.mainfile,
            with_embargo=calc_with_metadata.with_embargo,
            comment=calc_with_metadata.comment,
            coauthors=list(int(user['id']) for user in calc_with_metadata.coauthors),
            shared_with=list(int(user['id']) for user in calc_with_metadata.shared_with)
        )

    def source_calc_index(self, *args, **kwargs):
        """ see :func:`SourceCalc.index` """
        return SourceCalc.index(self.source, *args, **kwargs)

    def package_index(self, upload_path, **kwargs) -> None:
        """
        Creates Package objects and respective package zip files for the given uploads.
        The given uploads are supposed to be path to the extracted upload directories.
        If the upload is already in the index, it is not recreated.
        """
        logger = utils.get_logger(__name__)

        try:
            for package_entry in Package.get_packages(
                    upload_path, self.package_directory, create=True,
                    compress=self.compress_packages, **kwargs):

                logger.info(
                    'package in index',
                    source_upload_path=upload_path,
                    source_upload_id=package_entry.upload_id,
                    package_id=package_entry.package_id)
        except Exception as e:
            logger.error(
                'could not create package from upload',
                upload_path=upload_path, exc_info=e)


class Report(utils.POPO):
    def __init__(self, *args, **kwargs):
        self.total_packages = 0
        self.failed_packages = 0
        self.skipped_packages = 0
        self.total_calcs = 0  # the calcs that have been found by the target
        self.total_source_calcs = 0  # the calcs in the source index
        self.failed_calcs = 0  # calcs that have been migrated with failed processing
        self.migrated_calcs = 0   # the calcs from the source, successfully added to the target
        self.calcs_with_diffs = 0  # the calcs from the source, successfully added to the target with different metadata
        self.new_calcs = 0  # the calcs successfully added to the target that were not found in the source
        self.missing_calcs = 0  # the calcs in the source, that could not be added to the target due to failure or not founding the calc

        super().__init__(*args, **kwargs)

    def add(self, other: 'Report') -> None:
        for key, value in other.items():
            self[key] = self.get(key, 0) + value

    def __str__(self):
        return (
            'packages: {:,}, skipped: {:,}, source calcs: {:,}, migrated: {:,}, '
            'failed: {:,}, missing: {:,}, new: {:,}'.format(
                self.total_packages, self.skipped_packages,
                self.total_source_calcs, self.migrated_calcs,
                self.failed_calcs, self.missing_calcs,
                self.new_calcs))
