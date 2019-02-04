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
import os.path
from bravado.client import SwaggerClient
import json

from nomad import infrastructure, coe_repo, datamodel

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
    """
    Fixture for an example migration source db with:
    - two user
    - two calculations (1 per user)
    - one calculation with all metadata (dataset, ref, comment, coauther, sharewith)
    """
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

            schema_sql_file, example_data_sql_file = (
                os.path.join(os.path.dirname(infrastructure.__file__), 'empty_repository_db.sql'),
                os.path.join('tests', 'data', 'migration', 'example_source_db.sql'))

            for sql_file in [schema_sql_file, example_data_sql_file]:
                with open(sql_file, 'r') as f:
                    cur.execute(f.read())

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
        assert source_calc.mainfile in ['1/template.json', '2/template.json']
        assert source_calc.upload == 'upload'
        has_source_calc = True
        assert total == 3  # 2 calcs + 1 dataset

    assert has_source_calc == has_indexed

    test_calc = SourceCalc.objects(mainfile='1/template.json', upload='upload').first()
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


@pytest.fixture(scope='function')
def migrate_infra(migration, target_repo, flask_client, worker, monkeysession):
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
    indexed = list(migration.index(drop=True, with_metadata=True))
    assert len(indexed) == 2
    # source repo is the infrastructure repo
    migration.copy_users(target_repo)
    migration.set_new_pid_prefix(target_repo)

    # target repo is the infrastructure repo
    def create_client():
        admin = target_repo.query(coe_repo.User).filter_by(email='admin').first()
        http_client = FlaskTestHttpClient(flask_client, headers=create_auth_headers(admin))
        return SwaggerClient.from_url('/swagger.json', http_client=http_client)

    old_repo = infrastructure.repository_db
    monkeysession.setattr('nomad.infrastructure.repository_db', target_repo)
    monkeysession.setattr('nomad.client.create_client', create_client)

    yield migration

    monkeysession.setattr('nomad.infrastructure.repository_db', old_repo)


mirgation_test_specs = [
    ('baseline', dict(migrated=2, source=2)),
    ('archive', dict(migrated=2, source=2)),
    ('new_upload', dict(new=2)),
    ('new_calc', dict(migrated=2, source=2, new=1)),
    ('missing_calc', dict(migrated=1, source=2, missing=1)),
    ('missmatch', dict(migrated=2, source=2, diffs=1)),
    ('failed_calc', dict(migrated=1, source=2, diffs=0, missing=1, failed=1, errors=1)),
    ('failed_upload', dict(migrated=0, source=2, missing=2, errors=1))
]


@pytest.mark.filterwarnings("ignore:SAWarning")
@pytest.mark.parametrize('test, assertions', mirgation_test_specs)
@pytest.mark.timeout(30)
def test_migrate(migrate_infra, test, assertions, caplog):
    uploads_path = os.path.join('tests', 'data', 'migration', test)
    reports = list(migrate_infra.migrate(
        *[os.path.join(uploads_path, dir) for dir in os.listdir(uploads_path)]))

    assert len(reports) == 1
    report = reports[0]
    assert report['total_calcs'] == assertions.get('migrated', 0) + assertions.get('new', 0) + assertions.get('failed', 0)

    # assert if new, diffing, migrated calcs where detected correctly
    assert report['total_source_calcs'] == assertions.get('source', 0)
    assert report['migrated_calcs'] == assertions.get('migrated', 0)
    assert report['calcs_with_diffs'] == assertions.get('diffs', 0)
    assert report['new_calcs'] == assertions.get('new', 0)
    assert report['missing_calcs'] == assertions.get('missing', 0)

    # assert if migrated calcs have correct user metadata
    repo_db = infrastructure.repository_db
    if assertions.get('migrated', 0) > 0:
        calc_1 = repo_db.query(coe_repo.Calc).get(1)
        assert calc_1 is not None
        metadata = calc_1.to(datamodel.CalcWithMetadata)
        assert metadata.pid <= 2
        assert metadata.uploader == 1
        assert metadata.upload_time.isoformat() == '2019-01-01T12:00:00+00:00'
        assert len(metadata.datasets) == 1
        assert metadata.datasets[0]['id'] == 3
        assert metadata.datasets[0]['name'] == 'test_dataset'
        assert metadata.datasets[0]['dois'][0] == 'internal_ref'
        assert metadata.comment == 'label1'
        assert len(metadata.coauthors) == 1
        assert metadata.coauthors[0] == 2
        assert len(metadata.references) == 1
        assert metadata.references[0] == 'external_ref'

    if assertions.get('migrated', 0) > 1:
        calc_2 = repo_db.query(coe_repo.Calc).get(2)
        assert calc_1 is not None
        metadata = calc_2.to(datamodel.CalcWithMetadata)
        assert len(metadata.shared_with) == 1
        assert metadata.shared_with[0] == 1

    # assert pid prefix of new calcs
    if assertions.get('new', 0) > 0:
        assert repo_db.query(coe_repo.Calc).get(7000000) is not None

    errors = 0
    for record in caplog.get_records(when='call'):
        if record.levelname in ['ERROR', 'CRITICAL']:
            record_data = json.loads(record.getMessage())
            if 'source_upload_id' in record_data:
                errors += 1

    assert errors == assertions.get('errors', 0)
