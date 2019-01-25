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

import pytest
import os
from bravado.client import SwaggerClient

from nomad import infrastructure, coe_repo

from nomad.migration import NomadCOEMigration, SourceCalc
from nomad.infrastructure import repository_db_connection

from .bravado_flaks import FlaskTestHttpClient
from tests.conftest import create_repository_db
from tests.test_api import client as flask_client, create_auth_headers  # noqa pylint: disable=unused-import
from tests.test_client import client as bravado_client  # noqa pylint: disable=unused-import

test_source_db_name = 'test_nomad_fair_migration_source'
test_target_db_name = 'test_nomad_fair_migration_target'


@pytest.fixture(scope='module')
def source_repo(monkeysession, repository_db):
    try:
        with repository_db_connection(dbname='postgres', with_trans=False) as con:
            with con.cursor() as cursor:
                cursor.execute("CREATE DATABASE %s ;" % test_source_db_name)
    except Exception:
        pass

    with repository_db_connection(dbname=test_source_db_name, with_trans=False) as con:
        with con.cursor() as cur:
            cur.execute(
                'DROP SCHEMA IF EXISTS public CASCADE;'
                'CREATE SCHEMA IF NOT EXISTS public;'
                'GRANT ALL ON SCHEMA public TO postgres;'
                'GRANT ALL ON SCHEMA public TO public;')
            sql_file = os.path.join(os.path.dirname(infrastructure.__file__), 'empty_repository_db.sql')
            cur.execute(open(sql_file, 'r').read())
            cur.execute(
                'TRUNCATE TABLE public.users CASCADE;'
                "INSERT INTO public.users VALUES (1, 'one', 'one', 'one', 'one', NULL, NULL, NULL);"
                "INSERT INTO public.users VALUES (2, 'two', 'two', 'two', 'two', NULL, NULL, NULL);"
                "INSERT INTO public.calculations VALUES (NULL, NULL, NULL, NULL, 0, false, 1, NULL); "
                "INSERT INTO public.calculations VALUES (NULL, NULL, NULL, NULL, 0, false, 2, NULL); "
                "INSERT INTO public.codefamilies VALUES (1, 'test_code'); "
                "INSERT INTO public.codeversions VALUES (1, 1, 'test_version'); "
                "INSERT INTO public.metadata VALUES (1, NULL, NULL, NULL, NULL, 'formula1', '2019-01-01 12:00:00', NULL, decode('[\"$EXTRACTED/upload/test/mainfile.json\"]', 'escape'), 1, NULL); "
                "INSERT INTO public.metadata VALUES (1, NULL, NULL, NULL, NULL, 'formula2', '2015-01-01 13:00:00', NULL, decode('[\"$EXTRACTED/upload/test/mainfile.json\"]', 'escape'), 2, NULL); "
                "INSERT INTO public.spacegroups VALUES (1, 255); "
                "INSERT INTO public.spacegroups VALUES (2, 255); "
                "INSERT INTO public.user_metadata VALUES (1, 0, 'label1'); "
                "INSERT INTO public.user_metadata VALUES (2, 1, 'label2'); "
                "INSERT INTO public.ownerships VALUES (1, 1); "
                "INSERT INTO public.ownerships VALUES (2, 2); ")

    with create_repository_db(monkeysession, exists=True, readonly=True, dbname=test_source_db_name) as db:
        yield db


@pytest.fixture(scope='function')
def target_repo(repository_db):
    with create_repository_db(readonly=False, exists=False, dbname=test_target_db_name) as db:
        db.execute('TRUNCATE users CASCADE;')
        yield db
        db.execute('TRUNCATE uploads CASCADE;')


@pytest.fixture(scope='function')
def migration(source_repo, target_repo):
    migration = NomadCOEMigration(sites=['tests/data/migration'])
    yield migration


def test_copy_users(migration, target_repo):
    migration.copy_users(target_repo)
    assert target_repo.query(coe_repo.User).count() == 3
    assert target_repo.query(coe_repo.User).filter_by(user_id=1).first().email == 'one'
    assert target_repo.query(coe_repo.User).filter_by(user_id=2).first().email == 'two'


def perform_index(migration, has_indexed, with_metadata, **kwargs):
    has_source_calc = False
    for source_calc, total in SourceCalc.index(migration.source, with_metadata=with_metadata, **kwargs):
        assert source_calc.pid is not None
        assert source_calc.mainfile == 'test/mainfile.json'
        assert source_calc.upload == 'upload'
        has_source_calc = True
        assert total == 2

    assert has_source_calc == has_indexed

    test_calc = SourceCalc.objects(mainfile='test/mainfile.json', upload='upload').first()
    assert test_calc is not None

    if with_metadata:
        assert test_calc.metadata['uploader'] == 1
        assert test_calc.metadata['comment'] == 'label1'


@pytest.mark.parametrize('with_metadata', [False, True])
def test_create_index(migration, mockmongo, with_metadata: bool):
    perform_index(migration, has_indexed=True, drop=True, with_metadata=with_metadata)


@pytest.mark.parametrize('with_metadata', [True, False])
def test_update_index(migration, mockmongo, with_metadata: bool):
    perform_index(migration, has_indexed=True, drop=True, with_metadata=with_metadata)
    perform_index(migration, has_indexed=False, drop=False, with_metadata=with_metadata)


def test_migrate(migration, target_repo, flask_client, worker, monkeysession):
    """
    Parameters to test
    - missing upload, extracted, archive, broken archive
    - upload process failure
    - upload with no parsable files
    - calculations with process errors
    - matching, non matching calculations
    - to few calculations
    - to many caclualtions
    - not in the index

    All with two calcs, two users (for coauthors)
    """
    # source repo is the infrastructure repo
    migration.index(drop=True, with_metadata=True)
    # source repo is the infrastructure repo
    migration.copy_users(target_repo)

    # target repo is the infrastructure repo
    def create_client():
        admin = target_repo.query(coe_repo.User).filter_by(email='admin').first()
        http_client = FlaskTestHttpClient(flask_client, headers=create_auth_headers(admin))
        return SwaggerClient.from_url('/swagger.json', http_client=http_client)

    old_repo = infrastructure.repository_db
    monkeysession.setattr('nomad.infrastructure.repository_db', target_repo)
    monkeysession.setattr('nomad.client.create_client', create_client)
    migrated = 0
    try:
        migrated = migration.migrate('upload1')
    except Exception as e:
        raise e
    finally:
        monkeysession.setattr('nomad.infrastructure.repository_db', old_repo)

    assert migrated == 1
