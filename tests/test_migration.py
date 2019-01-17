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

from nomad.migration import NomadCOEMigration
from nomad.infrastructure import repository_db_connection


test_source_db_name = 'migration_source'


class TestNomadCOEMigration:

    @pytest.fixture(scope='session')
    def source_repo(self):
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
                    "INSERT INTO public.users VALUES (2, 'two', 'two', 'two', 'two', NULL, NULL, NULL);")

        yield dict(dbname=test_source_db_name)

        with repository_db_connection(dbname='postgres', with_trans=False) as con:
            with con.cursor() as cur:
                cur.execute(
                    "UPDATE pg_database SET datallowconn = 'false' WHERE datname = '%s'; "
                    "ALTER DATABASE %s CONNECTION LIMIT 1; "
                    "SELECT pg_terminate_backend(pid) FROM pg_stat_activity WHERE datname = '%s'; " %
                    (test_source_db_name, test_source_db_name, test_source_db_name))
                cur.execute("DROP DATABASE %s;" % test_source_db_name)

    @pytest.fixture(scope='function')
    def target_repo(self, clean_repository_db):
        # delete the default test users
        clean_repository_db.query(coe_repo.user.Session).delete()
        clean_repository_db.query(coe_repo.User).delete()
        return clean_repository_db

    @pytest.fixture(scope='function')
    def migration(self, source_repo, target_repo):
        migration = NomadCOEMigration(**source_repo)
        yield migration

        migration.source_repo_db.expunge_all()
        migration.source_repo_db.invalidate()
        migration.source_repo_db.close_all()
        migration.source_connection.close()
        migration.source_connection.engine.dispose()

    def test_copy_users(self, migration, clean_repository_db):
        migration.copy_users()
        assert clean_repository_db.query(coe_repo.User).count() == 2
        assert clean_repository_db.query(coe_repo.User).filter_by(user_id=1).first().email == 'one'
        assert clean_repository_db.query(coe_repo.User).filter_by(user_id=2).first().email == 'two'
