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

from nomad import config, infrastructure
from nomad.migration import NomadCOEMigration

from .main import cli


_migration: NomadCOEMigration = None


def _setup():
    pass


@cli.group(help='Migrate data from NOMAD CoE to nomad@FAIRDI')
@click.option('-h', '--host', default=config.migration_source_db.host, help='The migration repository source db host, default is "%s".' % config.migration_source_db.host)
@click.option('-p', '--port', default=config.migration_source_db.port, help='The migration repository source db port, default is %d.' % config.migration_source_db.port)
@click.option('-u', '--user', default=config.migration_source_db.user, help='The migration repository source db user, default is %s.' % config.migration_source_db.user)
@click.option('-w', '--password', default=config.migration_source_db.password, help='The migration repository source db password.')
@click.option('-db', '--dbname', default=config.migration_source_db.dbname, help='The migration repository source db name, default is %s.' % config.migration_source_db.dbname)
def migration(host, port, user, password, dbname):
    global _setup

    def _setup():
        infrastructure.setup_logging()
        infrastructure.setup_repository_db(
            readony=True, host=host, port=port, user=user, password=password, dbname=dbname)
        infrastructure.setup_mongo()

        global _migration
        _migration = NomadCOEMigration()


@migration.command(help='Create/update the coe repository db migration index')
@click.option('--drop', help='Drop the existing index, otherwise it will only add new data.', is_flag=True)
@click.option('--with-metadata', help='Extract metadata for each calc and add it to the index.', is_flag=True)
@click.option('--per-query', default=100, help='We index many objects with one query. Default is 100.')
def index(drop, with_metadata, per_query):
    _setup()
    start = time.time()
    indexed_total = 0
    indexed_calcs = 0
    for calc, total in _migration.index(drop=drop, with_metadata=with_metadata, per_query=int(per_query)):
        indexed_total += 1
        indexed_calcs += 1 if calc is not None else 0
        eta = total * ((time.time() - start) / indexed_total)
        print(
            'indexed: %8d, calcs: %8d, total: %8d, ETA: %s\r' %
            (indexed_total, indexed_calcs, total, datetime.timedelta(seconds=eta)), end='')
    print('done')


@migration.command(help='Add an upload folder to the package index.')
@click.argument('upload-paths', nargs=-1)
def package(upload_paths):
    infrastructure.setup_logging()
    infrastructure.setup_mongo()

    migration = NomadCOEMigration()
    migration.package(*upload_paths)


@migration.command(help='Copy users from source into empty target db')
def copy_users(**kwargs):
    _setup()
    _migration.copy_users()


@migration.command(help='Set the repo db PID calc counter.')
@click.argument('prefix', nargs=1, type=int, default=7000000)
def pid_prefix(prefix: int):
    infrastructure.setup_logging()
    migration = NomadCOEMigration()
    migration.set_pid_prefix(prefix=prefix)


@migration.command(help='Upload the given upload locations. Uses the existing index to provide user metadata')
@click.argument('paths', nargs=-1)
@click.option('--create-packages', help='Allow migration to create package entries on the fly.', is_flag=True)
@click.option('--local', help='Create local upload files.', is_flag=True)
@click.option('--delete-local', help='Delete created local upload files after upload.', is_flag=True)
@click.option('--parallel', default=1, type=int, help='Use the given amount of parallel processes. Default is 1.')
@click.option('--migration-version', default=0, type=int, help='The version number, only packages with lower or no number will be migrated.')
def upload(
        paths: list, create_packages, local: bool, delete_local: bool, parallel: int,
        migration_version: int):

    infrastructure.setup_logging()
    infrastructure.setup_mongo()

    migration = NomadCOEMigration(migration_version=migration_version, threads=parallel)
    migration.migrate(
        *paths, local=local, delete_local=delete_local, create_packages=create_packages)
