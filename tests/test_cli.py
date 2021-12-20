
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

from typing import cast
import pytest
import click.testing
import json
import datetime
import time
import zipfile
import os.path
import os
import re

from nomad import search, processing as proc, files, config
from nomad.cli import cli
from nomad.cli.cli import POPO
from nomad.processing import Upload, Calc
from nomad.processing.base import SUCCESS

from tests.app.flask.test_app import BlueprintClient
from tests.app.flask.conftest import (  # pylint: disable=unused-import
    test_user_bravado_client, client, session_client, admin_user_bravado_client)  # pylint: disable=unused-import
from tests.app.conftest import test_user_auth, admin_user_auth  # pylint: disable=unused-import
from tests.processing.test_data import run_processing

# TODO there is much more to test


@pytest.mark.usefixtures('reset_config', 'nomad_logging')
class TestCli:
    def test_help(self, example_mainfile):

        start = time.time()
        result = click.testing.CliRunner().invoke(
            cli, ['--help'], catch_exceptions=False)
        assert result.exit_code == 0
        assert time.time() - start < 1


@pytest.mark.usefixtures('reset_config', 'nomad_logging')
class TestParse:
    def test_parser(self, example_mainfile):
        _, mainfile_path = example_mainfile
        result = click.testing.CliRunner().invoke(
            cli, ['parse', mainfile_path], catch_exceptions=False)
        assert result.exit_code == 0


@pytest.mark.usefixtures('reset_config', 'no_warn', 'mongo_infra', 'elastic_infra', 'raw_files_infra')
class TestAdmin:
    def test_reset(self, reset_infra):
        result = click.testing.CliRunner().invoke(
            cli, ['admin', 'reset', '--i-am-really-sure'], catch_exceptions=False)
        assert result.exit_code == 0

    def test_reset_not_sure(self):
        result = click.testing.CliRunner().invoke(
            cli, ['admin', 'reset'], catch_exceptions=False)
        assert result.exit_code == 1

    # def test_clean(self, published):
    #     upload_id = published.upload_id

    #     Upload.objects(upload_id=upload_id).delete()
    #     assert published.upload_files.exists()
    #     assert Calc.objects(upload_id=upload_id).first() is not None
    #     assert search.SearchRequest().search_parameter('upload_id', upload_id).execute()['total'] > 0

    #     result = click.testing.CliRunner().invoke(
    #         cli, ['admin', 'clean', '--force', '--skip-es'], catch_exceptions=False)

    #     assert result.exit_code == 0
    #     assert not published.upload_files.exists()
    #     assert Calc.objects(upload_id=upload_id).first() is None
    #     assert search.SearchRequest().search_parameter('upload_id', upload_id).execute()['total'] > 0

    @pytest.mark.parametrize('upload_time,dry,lifted', [
        (datetime.datetime.now(), False, False),
        (datetime.datetime(year=2012, month=1, day=1), True, False),
        (datetime.datetime(year=2012, month=1, day=1), False, True)])
    def test_lift_embargo(self, published, upload_time, dry, lifted):
        upload_id = published.upload_id
        published.upload_time = upload_time
        published.save()
        calc = Calc.objects(upload_id=upload_id).first()

        assert published.upload_files.exists()
        assert calc.metadata['with_embargo']
        assert search.SearchRequest().owner('public').search_parameter('upload_id', upload_id).execute()['total'] == 0
        with pytest.raises(Exception):
            with files.UploadFiles.get(upload_id=upload_id).read_archive(calc_id=calc.calc_id):
                pass

        result = click.testing.CliRunner().invoke(
            cli, ['admin', 'lift-embargo'] + (['--dry'] if dry else []),
            catch_exceptions=False)

        assert result.exit_code == 0
        assert not Calc.objects(upload_id=upload_id).first().metadata['with_embargo'] == lifted
        assert (search.SearchRequest().owner('public').search_parameter('upload_id', upload_id).execute()['total'] > 0) == lifted
        if lifted:
            with files.UploadFiles.get(upload_id=upload_id).read_archive(calc_id=calc.calc_id) as archive:
                assert calc.calc_id in archive

    def test_delete_entry(self, published):
        upload_id = published.upload_id
        calc = Calc.objects(upload_id=upload_id).first()

        result = click.testing.CliRunner().invoke(
            cli, ['admin', 'entries', 'rm', calc.calc_id], catch_exceptions=False)

        assert result.exit_code == 0
        assert 'deleting' in result.stdout
        assert Upload.objects(upload_id=upload_id).first() is not None
        assert Calc.objects(calc_id=calc.calc_id).first() is None


def transform_for_index_test(calc):
    calc.comment = 'specific'
    return calc


@pytest.mark.usefixtures('reset_config')
class TestAdminUploads:

    @pytest.mark.parametrize('codes, count', [
        (['VASP'], 1),
        (['doesNotExist'], 0),
        (['VASP', 'doesNotExist'], 1)])
    def test_uploads_code(self, published, codes, count, no_warn):
        codes_args = []
        for code in codes:
            codes_args.append('--code')
            codes_args.append(code)
        result = click.testing.CliRunner().invoke(
            cli, ['admin', 'uploads'] + codes_args + ['ls'], catch_exceptions=False)

        assert result.exit_code == 0
        assert '%d uploads selected' % count in result.stdout

    def test_query_mongo(self, published, no_warn):
        upload_id = published.upload_id

        query = dict(upload_id=upload_id)
        result = click.testing.CliRunner().invoke(
            cli, ['admin', 'uploads', '--query-mongo', 'ls', json.dumps(query)],
            catch_exceptions=False)

        assert result.exit_code == 0
        assert '1 uploads selected' in result.stdout

    def test_ls(self, published, no_warn):
        upload_id = published.upload_id

        result = click.testing.CliRunner().invoke(
            cli, ['admin', 'uploads', 'ls', upload_id], catch_exceptions=False)

        assert result.exit_code == 0
        assert '1 uploads selected' in result.stdout

    def test_ls_query(self, published, no_warn):
        upload_id = published.upload_id

        result = click.testing.CliRunner().invoke(
            cli, ['admin', 'uploads', 'ls', '{"match":{"upload_id":"%s"}}' % upload_id], catch_exceptions=False)

        assert result.exit_code == 0
        assert '1 uploads selected' in result.stdout

    def test_rm(self, published, no_warn):
        upload_id = published.upload_id

        result = click.testing.CliRunner().invoke(
            cli, ['admin', 'uploads', 'rm', upload_id], catch_exceptions=False)

        assert result.exit_code == 0
        assert 'deleting' in result.stdout
        assert Upload.objects(upload_id=upload_id).first() is None
        assert Calc.objects(upload_id=upload_id).first() is None

    def test_index(self, published, no_warn):
        upload_id = published.upload_id
        calc = Calc.objects(upload_id=upload_id).first()
        calc.metadata['comment'] = 'specific'
        calc.save()

        assert search.SearchRequest().search_parameters(comment='specific').execute()['total'] == 0

        result = click.testing.CliRunner().invoke(
            cli, ['admin', 'uploads', 'index', upload_id], catch_exceptions=False)
        assert result.exit_code == 0
        assert 'index' in result.stdout

        assert search.SearchRequest().search_parameters(comment='specific').execute()['total'] == 1

    def test_index_with_transform(self, published, no_warn):
        upload_id = published.upload_id
        assert search.SearchRequest().search_parameters(comment='specific').execute()['total'] == 0

        result = click.testing.CliRunner().invoke(
            cli, [
                'admin', 'uploads', 'index',
                '--transformer', 'tests.test_cli.transform_for_index_test',
                upload_id],
            catch_exceptions=False)
        assert result.exit_code == 0
        assert 'index' in result.stdout

        assert search.SearchRequest().search_parameters(comment='specific').execute()['total'] == 1

    def test_re_process(self, published, monkeypatch, no_warn):
        monkeypatch.setattr('nomad.config.meta.version', 'test_version')
        upload_id = published.upload_id
        calc = Calc.objects(upload_id=upload_id).first()
        assert calc.metadata['nomad_version'] != 'test_version'

        result = click.testing.CliRunner().invoke(
            cli, ['admin', 'uploads', 're-process', '--parallel', '2', upload_id], catch_exceptions=False)

        assert result.exit_code == 0
        assert 're-processing' in result.stdout
        calc.reload()
        assert calc.metadata['nomad_version'] == 'test_version'

    def test_re_pack(self, published, monkeypatch, no_warn):
        upload_id = published.upload_id
        calc = Calc.objects(upload_id=upload_id).first()
        assert calc.metadata['with_embargo']
        calc.metadata['with_embargo'] = False
        calc.save()

        result = click.testing.CliRunner().invoke(
            cli, ['admin', 'uploads', 're-pack', '--parallel', '2', upload_id], catch_exceptions=False)

        assert result.exit_code == 0
        assert 're-pack' in result.stdout
        calc.reload()
        upload_files = files.PublicUploadFiles(upload_id)
        for raw_file in upload_files.raw_file_manifest():
            with upload_files.raw_file(raw_file) as f:
                f.read()
        for calc in Calc.objects(upload_id=upload_id):
            with upload_files.read_archive(calc.calc_id) as archive:
                assert calc.calc_id in archive

        published.reload()
        assert published.tasks_status == SUCCESS

    def test_chown(self, published, test_user, other_test_user, no_warn):
        upload_id = published.upload_id
        calc = Calc.objects(upload_id=upload_id).first()
        assert calc.metadata['uploader'] == other_test_user.user_id

        result = click.testing.CliRunner().invoke(
            cli, ['admin', 'uploads', 'chown', test_user.username, upload_id], catch_exceptions=False)

        assert result.exit_code == 0
        assert 'changing' in result.stdout

        upload = Upload.objects(upload_id=upload_id).first()
        calc.reload()

        assert upload.user_id == test_user.user_id
        assert calc.metadata['uploader'] == test_user.user_id

    def test_edit(self, published, no_warn):
        upload_id = published.upload_id

        def assert_calcs(publish, with_embargo):
            calcs = Calc.objects(upload_id=upload_id)
            for calc in calcs:
                assert calc.metadata['published'] == publish
                assert calc.metadata['with_embargo'] == with_embargo

            for calc in search.search(owner=None, query=dict(upload_id=upload_id)).data:
                assert calc['published'] == publish
                assert calc['with_embargo'] == with_embargo

        assert_calcs(True, True)

        def perform_test(publish, with_embargo):
            if publish:
                params = ['--publish', 'with-embargo' if with_embargo else 'no-embargo']
            else:
                assert not with_embargo
                params = ['--unpublish']

            result = click.testing.CliRunner().invoke(
                cli, ['admin', 'uploads', 'edit'] + params, catch_exceptions=False)

            assert result.exit_code == 0
            assert 'editing' in result.stdout
            assert_calcs(publish, with_embargo)

        perform_test(False, False)
        perform_test(True, False)
        perform_test(True, True)

    @pytest.mark.parametrize('with_calcs,success,failure', [
        (True, False, False),
        (False, False, False),
        (True, True, False),
        (False, False, True)])
    def test_reset(self, non_empty_processed, with_calcs, success, failure, no_warn):
        upload_id = non_empty_processed.upload_id

        upload = Upload.objects(upload_id=upload_id).first()
        calc = Calc.objects(upload_id=upload_id).first()
        assert upload.tasks_status == proc.SUCCESS
        assert calc.tasks_status == proc.SUCCESS

        args = ['admin', 'uploads', 'reset']
        if with_calcs: args.append('--with-calcs')
        if success: args.append('--success')
        if failure: args.append('--failure')
        args.append(upload_id)
        result = click.testing.CliRunner().invoke(cli, args, catch_exceptions=False)

        assert result.exit_code == 0
        assert 'reset' in result.stdout
        upload = Upload.objects(upload_id=upload_id).first()
        calc = Calc.objects(upload_id=upload_id).first()

        expected_state = proc.PENDING
        if success: expected_state = proc.SUCCESS
        if failure: expected_state = proc.FAILURE
        assert upload.tasks_status == expected_state
        if not with_calcs:
            assert calc.tasks_status == proc.SUCCESS
        else:
            assert calc.tasks_status == expected_state

    @pytest.mark.parametrize('paths,entries', [
        pytest.param(['0/POTCAR'], ['public'], id='single'),
        pytest.param(['0/POTCAR.gz', '0/POTCAR.xz', '0/POTCAR'], ['public'], id='multiple'),
        pytest.param(['0/POTCAR', 'POTCAR.xz'], ['public'], id='root'),
        pytest.param(['0/POTCAR'], ['embargo'], id='embargo')
    ])
    def test_quarantine_raw_files(self, paths, entries, test_user, proc_infra, no_warn):
        upload_path = os.path.join(config.fs.tmp, 'upload.zip')
        nomad_json = {}
        for i, entry in enumerate(entries):
            if entry == 'embargo':
                nomad_json.setdefault('entries', {})[f'{i}/archive.json'] = dict(with_embargo=True)
        with zipfile.ZipFile(upload_path, 'w') as zf:
            for path in paths:
                with zf.open(path, 'w') as f:
                    f.write(b'content')

            for i, _ in enumerate(entries):
                zf.write('tests/data/parsers/archive.json', f'{i}/archive.json')

            with zf.open('nomad.json', 'w') as f:
                f.write(json.dumps(nomad_json, indent=2).encode())

        # process upload
        upload = run_processing(('test_upload_id', upload_path), test_user)
        upload.publish_upload()
        try:
            upload.block_until_complete(interval=.01)
        except Exception:
            pass

        public_upload_files = cast(files.PublicUploadFiles, upload.upload_files)

        # run command
        exec_results = click.testing.CliRunner().invoke(
            cli, ['admin', 'uploads', 'quarantine-raw-files', '--', 'test_upload_id'],
            catch_exceptions=False)

        # assert results
        assert re.search(r'Moving .* to quarantine', exec_results.output) is not None, exec_results.output
        assert re.search(r'could not move files', exec_results.output) is None, exec_results.output

        # assert files
        quarantined_zip = public_upload_files._raw_file_object('quarantined')
        with zipfile.ZipFile(quarantined_zip.os_path, 'r') as zf:
            for path in paths:
                assert path in zf.namelist()

    @pytest.mark.parametrize('entries,with_empty,results', [
        pytest.param(['embargo', 'embargo'], True, [r'removed empty zip.*public'], id='embargo'),
        pytest.param(['public', 'public'], True, [r'removed empty zip.*restricted'], id='public'),
        pytest.param(['embargo'], False, [r'non empty public file'], id='non-empty-public'),
        pytest.param(['embargo', 'embargo'], None, [r'was already removed.*public'], id='embargo-wo-empty'),
        pytest.param(['public', 'public'], None, [r'was already removed.*restricted'], id='public-wo-empty'),
        pytest.param(['public', 'embargo'], None, [r'inconsistent upload'], id='mixed'),
        pytest.param(['public', 'public', 'other'], True, [r'quarantined'], id='mixed-other'),
        pytest.param(['public', 'public', 'other-large'], True, [r'inconsistent pack'], id='mixed-other-large')
    ])
    def test_prepare_migration(self, entries, results, with_empty, proc_infra, test_user):
        # create upload
        upload_path = os.path.join(config.fs.tmp, 'upload.zip')
        nomad_json = dict(entries={})
        entries_json = nomad_json['entries']
        with zipfile.ZipFile(upload_path, 'w') as zf:
            for i, entry in enumerate(entries):
                if entry in ['other', 'other-large']:
                    continue

                mainfile = f'{i}/archive.json'
                zf.write('tests/data/parsers/archive.json', mainfile)
                entries_json[mainfile] = {
                    'with_embargo': entry == 'embargo'
                }

            with zf.open('nomad.json', 'w') as f:
                f.write(json.dumps(nomad_json).encode())

        # process upload
        upload = run_processing(('test_upload_id', upload_path), test_user)
        upload.publish_upload()
        try:
            upload.block_until_complete(interval=.01)
        except Exception:
            pass

        public_upload_files = cast(files.PublicUploadFiles, upload.upload_files)

        # Add other files
        if any(entry in ['other', 'other-large'] for entry in entries):
            restricted_zip = public_upload_files._raw_file_object('restricted')
            with zipfile.ZipFile(restricted_zip.os_path, 'w') as zf:
                for i, entry in enumerate(entries):
                    if entry not in ['other', 'other-large']:
                        continue
                    with zf.open(f'{i}/not_nomad.txt', 'w') as f:
                        if entry == 'other-large':
                            f.write(b'I am not a nomad mainfile.' * 256)
                        else:
                            f.write(b'I am not a nomad mainfile.')

        # create empty restricted or publish ... procesing isn't doing that anymore, but
        # did in the past
        if with_empty is not None:
            for access in ['public', 'restricted']:
                zip_path = public_upload_files._raw_file_object(access)
                if zip_path.exists():
                    continue
                with zipfile.ZipFile(zip_path.os_path, 'w') as zf:
                    if not with_empty:
                        with zf.open('file.txt', 'w') as f:
                            f.write(b'Content')
                    pass

        # run command
        exec_results = click.testing.CliRunner().invoke(
            cli, ['admin', 'uploads', 'prepare-migration', '--', 'test_upload_id'],
            catch_exceptions=False)

        # assert results
        for result in results:
            assert re.search(result, exec_results.output) is not None, exec_results.output

        # assert files
        has_quarantined = re.search(r'quarantined', exec_results.output) is not None
        assert public_upload_files._raw_file_object('quarantined').exists() == has_quarantined


@pytest.mark.usefixtures('reset_config')
class TestClient:

    def test_upload(self, test_user_bravado_client, non_empty_example_upload, proc_infra):
        result = click.testing.CliRunner().invoke(
            cli,
            ['client', 'upload', '--offline', '--name', 'test_upload', non_empty_example_upload],
            catch_exceptions=False)

        assert result.exit_code == 0
        assert '1/0/1' in result.output
        assert proc.Upload.objects(name='test_upload').first() is not None

    def test_local(self, client, published, admin_user_bravado_client, monkeypatch):
        def requests_get(url, stream, headers):
            assert stream
            rv = client.get(url[url.index('/api/raw'):], headers=headers)
            assert rv.status_code == 200
            return POPO(iter_content=lambda *args, **kwargs: [bytes(rv.data)])

        monkeypatch.setattr('requests.get', requests_get)
        result = click.testing.CliRunner().invoke(
            cli,
            ['client', 'local', '%s/%s' % (published.upload_id, list(published.calcs)[0].calc_id)],
            catch_exceptions=False)

        assert result.exit_code == 0

    def test_mirror_dry(self, published, admin_user_bravado_client, monkeypatch):
        monkeypatch.setattr('nomad.cli.client.mirror.__in_test', True)

        result = click.testing.CliRunner().invoke(
            cli, ['client', 'mirror', '--dry'], catch_exceptions=False)

        assert result.exit_code == 0
        assert published.upload_id in result.output
        assert published.upload_files.os_path in result.output

    @pytest.mark.parametrize('move, link', [(True, False), (False, True), (False, False)])
    def test_mirror(self, published, admin_user_bravado_client, monkeypatch, move, link):
        ref_search_results = search.flat(
            search.SearchRequest().search_parameters(
                upload_id=published.upload_id).execute_paginated()['results'][0])

        monkeypatch.setattr('nomad.cli.client.mirror.__in_test', True)

        if move:
            result = click.testing.CliRunner().invoke(
                cli, ['client', 'mirror', '--move'], catch_exceptions=False)
        elif link:
            result = click.testing.CliRunner().invoke(
                cli, ['client', 'mirror', '--link'], catch_exceptions=False)
        else:
            result = click.testing.CliRunner().invoke(
                cli, ['client', 'mirror', '--source-mapping', '.volumes/test_fs:.volumes/test_fs'], catch_exceptions=False)

        assert result.exit_code == 0
        assert published.upload_id in result.output
        assert published.upload_files.os_path in result.output
        assert proc.Upload.objects(upload_id=published.upload_id).count() == 1
        assert proc.Calc.objects(upload_id=published.upload_id).count() == 1
        new_search = search.SearchRequest().search_parameters(upload_id=published.upload_id).execute_paginated()
        calcs_in_search = new_search['pagination']['total']
        assert calcs_in_search == 1

        new_search_results = search.flat(new_search['results'][0])
        for key in new_search_results.keys():
            if key not in ['upload_time', 'last_processing', 'dft.labels.label', 'owners', 'authors', 'uploader', 'coauthors', 'shared_with']:
                # There is a sub second change due to date conversions (?).
                # Labels have arbitrary order.
                assert json.dumps(new_search_results[key]) == json.dumps(ref_search_results[key])

        published.upload_files.exists
        proc.Upload.objects(upload_id=published.upload_id).first().upload_files.exists

    def test_mirror_staging(self, non_empty_processed, admin_user_bravado_client, monkeypatch):
        ref_search_results = search.flat(
            search.SearchRequest().search_parameters(
                upload_id=non_empty_processed.upload_id).execute_paginated()['results'][0])

        monkeypatch.setattr('nomad.cli.client.mirror.__in_test', True)

        result = click.testing.CliRunner().invoke(
            cli, ['client', 'mirror', '--staging', '--link'], catch_exceptions=False)

        assert result.exit_code == 0
        assert non_empty_processed.upload_id in result.output
        assert non_empty_processed.upload_files.os_path in result.output
        assert proc.Upload.objects(upload_id=non_empty_processed.upload_id).count() == 1
        assert proc.Calc.objects(upload_id=non_empty_processed.upload_id).count() == 1
        new_search = search.SearchRequest().search_parameters(upload_id=non_empty_processed.upload_id).execute_paginated()
        calcs_in_search = new_search['pagination']['total']
        assert calcs_in_search == 1

        new_search_results = search.flat(new_search['results'][0])
        for key in new_search_results.keys():
            if key not in ['upload_time', 'last_processing']:  # There is a sub second change due to date conversions (?)
                assert json.dumps(new_search_results[key]) == json.dumps(ref_search_results[key])

        non_empty_processed.upload_files.exists
        proc.Upload.objects(upload_id=non_empty_processed.upload_id).first().upload_files.exists

    def test_mirror_files_only(self, published, admin_user_bravado_client, monkeypatch):
        monkeypatch.setattr('nomad.cli.client.mirror.__in_test', True)

        result = click.testing.CliRunner().invoke(
            cli, ['client', 'mirror', '--files-only'], catch_exceptions=False)

        assert result.exit_code == 0, result.output
        assert published.upload_id in result.output
        assert published.upload_files.os_path in result.output

        published.upload_files.exists

    def test_mirror_datasets(self, client, published_wo_user_metadata, test_user_auth, admin_user_bravado_client, monkeypatch):
        # use the API to create dataset and DOI
        api = BlueprintClient(client, '/api')
        rv = api.post(
            '/repo/edit', headers=test_user_auth, content_type='application/json',
            data=json.dumps({
                'actions': {
                    'datasets': [{
                        'value': 'test_dataset'
                    }]
                }
            }))
        assert rv.status_code == 200

        rv = api.post('/datasets/test_dataset', headers=test_user_auth)
        assert rv.status_code == 200

        # perform the mirror
        monkeypatch.setattr('nomad.cli.client.mirror.__in_test', True)

        result = click.testing.CliRunner().invoke(
            cli, ['client', 'mirror'], catch_exceptions=False)

        assert result.exit_code == 0, result.output
        assert published_wo_user_metadata.upload_id in result.output
        assert published_wo_user_metadata.upload_files.os_path in result.output

        published_wo_user_metadata.upload_files.exists

    def test_statistics(self, client, proc_infra, admin_user_bravado_client):

        result = click.testing.CliRunner().invoke(
            cli, ['client', 'statistics-table'], catch_exceptions=True)

        assert result.exit_code == 0, result.output
        assert 'Calculations, e.g. total energies' in result.output
        assert 'Geometries' in result.output
        assert 'Bulk crystals' in result.output
        assert '2D / Surfaces' in result.output
        assert 'Atoms / Molecules' in result.output
        assert 'DOS' in result.output
        assert 'Band structures' in result.output
