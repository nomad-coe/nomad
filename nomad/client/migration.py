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

import click
import time
import datetime
import os
import os.path
import re
import shutil
import multiprocessing
import queue
import json

from nomad import config, infrastructure
from nomad.migration import NomadCOEMigration, SourceCalc, Package

from .main import cli


def _Migration(**kwargs) -> NomadCOEMigration:
    return NomadCOEMigration(**kwargs)


def _setup():
    pass


@cli.group(help='Migrate data from NOMAD CoE to nomad@FAIRDI')
@click.option('-h', '--host', default=config.migration_source_db.host, help='The migration repository source db host, default is "%s".' % config.migration_source_db.host)
@click.option('-p', '--port', default=config.migration_source_db.port, help='The migration repository source db port, default is %d.' % config.migration_source_db.port)
@click.option('-u', '--user', default=config.migration_source_db.user, help='The migration repository source db user, default is %s.' % config.migration_source_db.user)
@click.option('-w', '--password', default=config.migration_source_db.password, help='The migration repository source db password.')
@click.option('-db', '--dbname', default=config.migration_source_db.dbname, help='The migration repository source db name, default is %s.' % config.migration_source_db.dbname)
@click.option('--migration-version', default=0, type=int, help='The version number, only packages with lower or no number will be migrated.')
@click.option('--package-directory', default=config.fs.migration_packages, help='The directory used as bucket for upload packages, default is %s.' % config.fs.migration_packages)
@click.option('--compress-packages', is_flag=True, help='Turn on compression for creating migration packages')
def migration(
        host, port, user, password, dbname, migration_version, package_directory, compress_packages):
    global _setup

    def _setup():
        infrastructure.setup_logging()
        infrastructure.setup_repository_db(
            readony=True, host=host, port=port, user=user, password=password, dbname=dbname)
        infrastructure.setup_mongo()

    global _Migration

    def _Migration(**kwargs):
        return NomadCOEMigration(
            migration_version=migration_version, package_directory=package_directory,
            compress_packages=compress_packages, **kwargs)


@migration.command(help='Create/update the coe repository db migration index')
@click.option('--drop', help='Drop the existing index, otherwise it will only add new data.', is_flag=True)
@click.option('--with-metadata', help='Extract metadata for each calc and add it to the index.', is_flag=True)
@click.option('--per-query', default=100, help='We index many objects with one query. Default is 100.')
def index(drop, with_metadata, per_query):
    _setup()
    start = time.time()
    indexed_total = 0
    indexed_calcs = 0
    for calc, total in _Migration().source_calc_index(drop=drop, with_metadata=with_metadata, per_query=int(per_query)):
        indexed_total += 1
        indexed_calcs += 1 if calc is not None else 0
        eta = total * ((time.time() - start) / indexed_total)
        print(
            'indexed: %8d, calcs: %8d, total: %8d, ETA: %s\r' %
            (indexed_total, indexed_calcs, total, datetime.timedelta(seconds=eta)), end='')
    print('done')


@migration.command(help='Reset migration version to start a new migration.')
@click.option('--delete-packages', is_flag=True, help='Also remove all packages.')
def reset(delete_packages: bool):
    infrastructure.setup_logging()
    infrastructure.setup_mongo()

    SourceCalc.objects(migration_version__ne=-1).update(migration_version=-1)
    if delete_packages:
        for subdir in os.listdir(config.fs.migration_packages):
            shutil.rmtree(os.path.join(config.fs.migration_packages, subdir))
        Package.objects().delete()
    else:
        Package.objects(migration_version__ne=-1).update(migration_version=-1)


def determine_upload_paths(paths, pattern=None):
    if pattern is not None:
        assert len(paths) == 1, "Can only apply pattern on a single directory."
        path = paths[0]
        if pattern == "ALL":
            paths = [os.path.join(path, directory) for directory in os.listdir(path)]
        else:
            paths = []
            compiled_pattern = re.compile(pattern)
            directories = os.listdir(path)
            directories.sort()
            for sub_directory in directories:
                if re.fullmatch(compiled_pattern, sub_directory):
                    paths.append(os.path.join(path, sub_directory))

    return paths


@migration.command(help='Add an upload folder to the package index.')
@click.argument('upload-paths', nargs=-1)
@click.option('--pattern', default=None, type=str, help='Interpret the paths as directory and migrate those subdirectory that match the given regexp')
@click.option('--parallel', default=1, type=int, help='Use the given amount of parallel processes to process uploads. Default is 1.')
@click.option('--parallel-zip', default=1, type=int, help='Use the given amount of parallel processes to pack packages. Default is 1.')
def package(upload_paths, pattern, parallel, parallel_zip):
    upload_paths = determine_upload_paths(upload_paths, pattern)
    upload_path_queue = multiprocessing.Queue(len(upload_paths))

    print('Package %d uploads with %d/%d processes.' % (len(upload_paths), parallel, parallel_zip))

    for upload_path in upload_paths:
        upload_path_queue.put(upload_path)

    def package_paths():
        infrastructure.setup_logging()
        infrastructure.setup_mongo()

        migration = _Migration()

        try:
            while True:
                upload_path = upload_path_queue.get()
                migration.package_index(upload_path, parallel=parallel_zip)
        except queue.Empty:
            pass

    processes = []
    for _ in range(0, parallel):
        process = multiprocessing.Process(target=package_paths)
        process.start()
        processes.append(process)

    for process in processes:
        process.join()
    upload_path_queue.close()


@migration.command(help='Get an report over all migrated packages.')
def report():
    infrastructure.setup_logging()
    infrastructure.setup_mongo()

    report = _Migration().report()
    print(report)


@migration.command(help='Copy users from source into empty target db')
def copy_users(**kwargs):
    _setup()
    _Migration().copy_users()


@migration.command(help='Set the repo db PID calc counter.')
@click.argument('prefix', nargs=1, type=int, default=7000000)
def pid_prefix(prefix: int):
    infrastructure.setup_logging()
    _Migration().set_pid_prefix(prefix=prefix)


@migration.command(help='Upload the given upload locations. Uses the existing index to provide user metadata')
@click.argument('upload-paths', nargs=-1)
@click.option('--pattern', default=None, type=str, help='Interpret the paths as directory and migrate those subdirectory that match the given regexp')
@click.option('--delete-failed', default='', type=str, help='String from N, U, P to determine if empty (N), failed (U), or failed to publish (P) uploads should be deleted or kept for debugging.')
@click.option('--parallel', default=1, type=int, help='Use the given amount of parallel processes. Default is 1.')
@click.option('--create-packages', is_flag=True, help='Indicate that packages should be created, if they do not already exist.')
def upload(
        upload_paths: list, pattern: str, parallel: int, delete_failed: str,
        create_packages: bool):

    infrastructure.setup_logging()
    infrastructure.setup_mongo()

    _Migration(threads=parallel).migrate(
        *determine_upload_paths(upload_paths, pattern), delete_failed=delete_failed,
        create_packages=create_packages)


@migration.command(help='Get an report about not migrated calcs.')
def missing():
    infrastructure.setup_logging()
    infrastructure.setup_mongo()

    report = SourceCalc.missing()
    print(json.dumps(report, indent=2))
