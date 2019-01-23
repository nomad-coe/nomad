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

from nomad import infrastructure, coe_repo

from nomad.migration import NomadCOEMigration, SourceCalc
from nomad.infrastructure import repository_db_connection

from tests.conftest import create_repository_db

test_source_db_name = 'test_nomad_fair_migration_source'
test_target_db_name = 'test_nomad_fair_migration_target'


class TestNomadCOEMigration:
    @pytest.fixture(scope='module')
    def source_repo(self, monkeysession, repository_db):
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
                    "INSERT INTO public.metadata VALUES (NULL, NULL, NULL, NULL, NULL, 'formula', '2015-05-27 11:56:37', NULL, decode('[\"$EXTRACTED/upload1/test/out.xml\"]', 'escape'), 1, NULL); "
                    "INSERT INTO public.user_metadata VALUES (1, 0, 'label1'); "
                    "INSERT INTO public.ownerships VALUES (1, 1); ")

        with create_repository_db(monkeysession, exists=True, readonly=True, dbname=test_source_db_name) as db:
            yield db

    @pytest.fixture(scope='function')
    def target_repo(self, repository_db):
        with create_repository_db(readonly=False, exists=False, dbname=test_target_db_name) as db:
            db.execute('TRUNCATE users CASCADE;')
            yield db
            db.execute('TRUNCATE uploads CASCADE;')

    @pytest.fixture(scope='function')
    def migration(self, source_repo, target_repo):
        migration = NomadCOEMigration()
        yield migration

    def test_copy_users(self, migration, target_repo):
        migration.copy_users(target_repo)
        assert target_repo.query(coe_repo.User).count() == 2
        assert target_repo.query(coe_repo.User).filter_by(user_id=1).first().email == 'one'
        assert target_repo.query(coe_repo.User).filter_by(user_id=2).first().email == 'two'

    def test_index(self, migration, mockmongo):
        SourceCalc.index(migration.source)
        test_calc = SourceCalc.objects(mainfile='test/out.xml', upload='upload1').first()

        assert test_calc is not None
        assert test_calc.metadata['_uploader'] == 1
        assert test_calc.metadata['comment'] == 'label1'
