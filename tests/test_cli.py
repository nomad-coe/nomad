
#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import pytest
import click.testing
import json
import datetime
import time

from nomad import processing as proc, files
from nomad.search import search
from nomad.cli import cli
from nomad.cli.cli import POPO
from nomad.processing import Upload, Entry, ProcessStatus
from nomad.utils.exampledata import ExampleData

# TODO there is much more to test


def invoke_cli(*args, **kwargs):
    return click.testing.CliRunner().invoke(*args, obj=POPO(), **kwargs)


@pytest.mark.usefixtures('reset_config', 'nomad_logging')
class TestCli:
    def test_help(self, example_mainfile):

        start = time.time()
        result = invoke_cli(
            cli, ['--help'], catch_exceptions=False)
        assert result.exit_code == 0
        assert time.time() - start < 1


@pytest.mark.usefixtures('reset_config', 'nomad_logging')
class TestParse:
    def test_parser(self, example_mainfile):
        _, mainfile_path = example_mainfile
        result = invoke_cli(
            cli, ['parse', mainfile_path], catch_exceptions=False)
        assert result.exit_code == 0


@pytest.mark.usefixtures('reset_config', 'no_warn', 'mongo_infra', 'elastic_infra', 'raw_files_infra')
class TestAdmin:
    def test_reset(self, reset_infra):
        result = invoke_cli(
            cli, ['admin', 'reset', '--i-am-really-sure'], catch_exceptions=False)
        assert result.exit_code == 0

    def test_reset_not_sure(self):
        result = invoke_cli(
            cli, ['admin', 'reset'], catch_exceptions=False)
        assert result.exit_code == 1

    # TODO this has somekind of raise condition in it and the test fails every other time
    # on the CI/CD
    # def test_clean(self, published):
    #     upload_id = published.upload_id

    #     Upload.objects(upload_id=upload_id).delete()
    #     assert published.upload_files.exists()
    #     assert Entry.objects(upload_id=upload_id).first() is not None
    #     search.refresh()
    #     assert search.SearchRequest().search_parameter('upload_id', upload_id).execute()['total'] > 0
    #     # TODO test new index pair
    #     # assert es_search(owner=None, query=dict(upload_id=upload_id)).pagination.total == 0
    #     # assert es_search(owner=None, query=dict(upload_id=upload_id)).pagination.total == 0

    #     result = invoke_cli(
    #         cli, ['admin', 'clean', '--force', '--skip-es'], catch_exceptions=False)

    #     assert result.exit_code == 0
    #     assert not published.upload_files.exists()
    #     assert Entry.objects(upload_id=upload_id).first() is None
    #     search.refresh()
    #     assert search.SearchRequest().search_parameter('upload_id', upload_id).execute()['total'] > 0
    #     # TODO test new index pair
    #     # assert es_search(owner=None, query=dict(upload_id=upload_id)).pagination.total == 0

    @pytest.mark.parametrize('publish_time,dry,lifted', [
        (datetime.datetime.now(), False, False),
        (datetime.datetime(year=2012, month=1, day=1), True, False),
        (datetime.datetime(year=2012, month=1, day=1), False, True)])
    def test_lift_embargo(self, published, publish_time, dry, lifted):
        upload_id = published.upload_id
        published.publish_time = publish_time
        published.save()
        entry = Entry.objects(upload_id=upload_id).first()

        assert published.upload_files.exists()
        assert published.with_embargo

        assert search(owner='public', query=dict(upload_id=upload_id)).pagination.total == 0

        result = invoke_cli(
            cli, ['admin', 'lift-embargo'] + (['--dry'] if dry else []),
            catch_exceptions=False)

        assert result.exit_code == 0
        published.block_until_complete()
        assert not published.with_embargo == lifted
        assert (search(owner='public', query=dict(upload_id=upload_id)).pagination.total > 0) == lifted
        if lifted:
            with files.UploadFiles.get(upload_id=upload_id).read_archive(entry_id=entry.entry_id) as archive:
                assert entry.entry_id in archive

    def test_delete_entry(self, published):
        upload_id = published.upload_id
        entry = Entry.objects(upload_id=upload_id).first()

        result = invoke_cli(
            cli, ['admin', 'entries', 'rm', entry.entry_id], catch_exceptions=False)

        assert result.exit_code == 0
        assert 'deleting' in result.stdout
        assert Upload.objects(upload_id=upload_id).first() is not None
        assert Entry.objects(entry_id=entry.entry_id).first() is None


def transform_for_index_test(entry):
    entry.comment = 'specific'
    return entry


@pytest.mark.usefixtures('reset_config', 'no_warn')
class TestAdminUploads:

    def test_query_mongo(self, published):
        upload_id = published.upload_id

        query = dict(upload_id=upload_id)
        result = invoke_cli(
            cli, ['admin', 'uploads', '--entries-mongo-query', json.dumps(query), 'ls'],
            catch_exceptions=False)

        assert result.exit_code == 0
        assert '1 uploads selected' in result.stdout

    def test_ls(self, published):
        upload_id = published.upload_id

        result = invoke_cli(
            cli, ['admin', 'uploads', 'ls', upload_id], catch_exceptions=False)

        assert result.exit_code == 0
        assert '1 uploads selected' in result.stdout

    def test_ls_query(self, published):
        upload_id = published.upload_id

        result = invoke_cli(
            cli, ['admin', 'uploads', '--entries-es-query', f'{{"upload_id":"{upload_id}"}}', 'ls'], catch_exceptions=False)
        assert result.exit_code == 0
        assert '1 uploads selected' in result.stdout

    def test_rm(self, published):
        upload_id = published.upload_id

        result = invoke_cli(
            cli, ['admin', 'uploads', 'rm', upload_id], catch_exceptions=False)

        assert result.exit_code == 0
        assert 'deleting' in result.stdout
        assert Upload.objects(upload_id=upload_id).first() is None
        assert Entry.objects(upload_id=upload_id).first() is None

    def test_index(self, published):
        upload_id = published.upload_id
        entry = Entry.objects(upload_id=upload_id).first()
        entry.comment = 'specific'
        entry.save()

        assert search(owner='all', query=dict(comment='specific')).pagination.total == 0

        result = invoke_cli(
            cli, ['admin', 'uploads', 'index', upload_id], catch_exceptions=False)
        assert result.exit_code == 0
        assert 'index' in result.stdout

        assert search(owner='all', query=dict(comment='specific')).pagination.total == 1

    def test_index_with_transform(self, published):
        upload_id = published.upload_id
        assert search(owner='all', query=dict(comment='specific')).pagination.total == 0

        result = invoke_cli(
            cli, [
                'admin', 'uploads', 'index',
                '--transformer', 'tests.test_cli.transform_for_index_test',
                upload_id],
            catch_exceptions=False)
        assert result.exit_code == 0
        assert 'index' in result.stdout

        assert search(owner='all', query=dict(comment='specific')).pagination.total == 1

    def test_re_process(self, published, monkeypatch):
        monkeypatch.setattr('nomad.config.meta.version', 'test_version')
        upload_id = published.upload_id
        entry = Entry.objects(upload_id=upload_id).first()
        assert entry.nomad_version != 'test_version'

        result = invoke_cli(
            cli, ['admin', 'uploads', 'process', '--parallel', '2', upload_id], catch_exceptions=False)

        assert result.exit_code == 0
        assert 'processing' in result.stdout
        entry.reload()
        assert entry.nomad_version == 'test_version'

    def test_re_pack(self, published, monkeypatch):
        upload_id = published.upload_id
        entry = Entry.objects(upload_id=upload_id).first()
        assert published.with_embargo
        published.embargo_length = 0
        published.save()

        result = invoke_cli(
            cli, ['admin', 'uploads', 're-pack', upload_id], catch_exceptions=False)

        assert result.exit_code == 0
        assert 're-pack' in result.stdout
        entry.reload()
        upload_files = files.PublicUploadFiles(upload_id)
        for path_info in upload_files.raw_directory_list(recursive=True, files_only=True):
            with upload_files.raw_file(path_info.path) as f:
                f.read()
        for entry in Entry.objects(upload_id=upload_id):
            with upload_files.read_archive(entry.entry_id) as archive:
                assert entry.entry_id in archive

        published.reload()
        assert published.process_status == ProcessStatus.SUCCESS

    def test_chown(self, published: Upload, test_user, other_test_user):
        upload_id = published.upload_id
        assert published.main_author == test_user.user_id
        with published.entries_metadata() as entries_metadata:
            for entry_metadata in entries_metadata:
                assert entry_metadata.main_author.user_id == test_user.user_id

        result = invoke_cli(
            cli, ['admin', 'uploads', 'chown', other_test_user.username, upload_id], catch_exceptions=False)

        assert result.exit_code == 0
        assert 'changing' in result.stdout

        published.block_until_complete()

        assert published.main_author == other_test_user.user_id
        with published.entries_metadata() as entries_metadata:
            for entry_metadata in entries_metadata:
                assert entry_metadata.main_author.user_id == other_test_user.user_id

    @pytest.mark.parametrize('with_entries,success,failure', [
        (True, False, False),
        (False, False, False),
        (True, True, False),
        (False, False, True)])
    def test_reset(self, non_empty_processed, with_entries, success, failure):
        upload_id = non_empty_processed.upload_id

        upload = Upload.objects(upload_id=upload_id).first()
        entry = Entry.objects(upload_id=upload_id).first()
        assert upload.process_status == ProcessStatus.SUCCESS
        assert entry.process_status == ProcessStatus.SUCCESS

        args = ['admin', 'uploads', 'reset']
        if with_entries: args.append('--with-entries')
        if success: args.append('--success')
        if failure: args.append('--failure')
        args.append(upload_id)
        result = invoke_cli(cli, args, catch_exceptions=False)

        assert result.exit_code == 0
        assert 'reset' in result.stdout
        upload = Upload.objects(upload_id=upload_id).first()
        entry = Entry.objects(upload_id=upload_id).first()

        expected_state = ProcessStatus.READY
        if success: expected_state = ProcessStatus.SUCCESS
        if failure: expected_state = ProcessStatus.FAILURE
        assert upload.process_status == expected_state
        if not with_entries:
            assert entry.process_status == ProcessStatus.SUCCESS
        else:
            assert entry.process_status == expected_state

    @pytest.mark.parametrize('indexed', [True, False])
    def test_integrity_entry_index(self, test_user, mongo, elastic, indexed):
        data = ExampleData(main_author=test_user)
        data.create_upload(upload_id='test_upload')
        data.create_entry(upload_id='test_upload')
        data.save(with_es=indexed, with_files=False)

        result = invoke_cli(cli, 'admin uploads integrity entry-index', catch_exceptions=True)

        assert result.exit_code == 0
        assert ('test_upload' in result.output) != indexed


@pytest.mark.usefixtures('reset_config')
class TestClient:

    def test_upload(self, non_empty_example_upload, admin_user, proc_infra, client_with_api_v1):
        result = invoke_cli(
            cli,
            [
                'client', '-u', admin_user.username, '--token-via-api',
                'upload', '--upload-name', 'test_upload', '--local-path',
                non_empty_example_upload],
            catch_exceptions=False)

        assert result.exit_code == 0, result.output
        assert '1/0/1' in result.output
        assert proc.Upload.objects(upload_name='test_upload').first() is not None

    def test_local(self, published_wo_user_metadata, client_with_api_v1):
        result = invoke_cli(
            cli,
            ['client', 'local', published_wo_user_metadata.successful_entries[0].entry_id],
            catch_exceptions=True)

        assert result.exit_code == 0, result.output

    @pytest.mark.skip('Disabled. Tested code is temporaely commented.')
    def test_statistics(self):

        result = invoke_cli(
            cli, ['client', 'statistics-table'], catch_exceptions=True)

        assert result.exit_code == 0, result.output
        assert 'Calculations, e.g. total energies' in result.output
        assert 'Geometries' in result.output
        assert 'Bulk crystals' in result.output
        assert '2D / Surfaces' in result.output
        assert 'Atoms / Molecules' in result.output
        assert 'DOS' in result.output
        assert 'Band structures' in result.output
