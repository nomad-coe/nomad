
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
import click.testing
import json

from nomad import utils, search, processing as proc
from nomad.cli import cli
from nomad.processing import Upload, Calc

# TODO there is much more to test


@pytest.mark.usefixtures('reset_config', 'no_warn')
class TestAdmin:
    def test_clean(self, published):
        upload_id = published.upload_id

        Upload.objects(upload_id=upload_id).delete()
        assert published.upload_files.exists()
        assert Calc.objects(upload_id=upload_id).first() is not None
        assert search.SearchRequest().search_parameter('upload_id', upload_id).execute()['total'] > 0

        result = click.testing.CliRunner().invoke(
            cli, ['admin', 'clean', '--force', '--skip-es'], catch_exceptions=False, obj=utils.POPO())

        assert result.exit_code == 0
        assert not published.upload_files.exists()
        assert Calc.objects(upload_id=upload_id).first() is None
        assert search.SearchRequest().search_parameter('upload_id', upload_id).execute()['total'] > 0

    def test_index(self, published):
        upload_id = published.upload_id
        calc = Calc.objects(upload_id=upload_id).first()
        calc.metadata['comment'] = 'specific'
        calc.save()

        assert search.SearchRequest().search_parameter('comment', 'specific').execute()['total'] == 0

        result = click.testing.CliRunner().invoke(
            cli, ['admin', 'index', '--threads', '2'], catch_exceptions=False, obj=utils.POPO())
        assert result.exit_code == 0
        assert 'index' in result.stdout

        assert search.SearchRequest().search_parameter('comment', 'specific').execute()['total'] == 1


@pytest.mark.usefixtures('reset_config', 'no_warn')
class TestAdminUploads:

    @pytest.mark.parametrize('codes, count', [
        (['VASP'], 1),
        (['doesNotExist'], 0),
        (['VASP', 'doesNotExist'], 1)])
    def test_uploads_code(self, published, codes, count):
        codes_args = []
        for code in codes:
            codes_args.append('--code')
            codes_args.append(code)
        result = click.testing.CliRunner().invoke(
            cli, ['admin', 'uploads'] + codes_args + ['ls'], catch_exceptions=False, obj=utils.POPO())

        assert result.exit_code == 0
        assert '%d uploads selected' % count in result.stdout

    def test_ls(self, published):
        upload_id = published.upload_id

        result = click.testing.CliRunner().invoke(
            cli, ['admin', 'uploads', 'ls', upload_id], catch_exceptions=False, obj=utils.POPO())

        assert result.exit_code == 0
        assert '1 uploads selected' in result.stdout

    def test_rm(self, published):
        upload_id = published.upload_id

        result = click.testing.CliRunner().invoke(
            cli, ['admin', 'uploads', 'rm', upload_id], catch_exceptions=False, obj=utils.POPO())

        assert result.exit_code == 0
        assert 'deleting' in result.stdout
        assert Upload.objects(upload_id=upload_id).first() is None
        assert Calc.objects(upload_id=upload_id).first() is None

    def test_index(self, published):
        upload_id = published.upload_id
        calc = Calc.objects(upload_id=upload_id).first()
        calc.metadata['comment'] = 'specific'
        calc.save()

        assert search.SearchRequest().search_parameters(comment='specific').execute()['total'] == 0

        result = click.testing.CliRunner().invoke(
            cli, ['admin', 'uploads', 'index', upload_id], catch_exceptions=False, obj=utils.POPO())
        assert result.exit_code == 0
        assert 'index' in result.stdout

        assert search.SearchRequest().search_parameters(comment='specific').execute()['total'] == 1

    def test_re_process(self, published, monkeypatch):
        monkeypatch.setattr('nomad.config.version', 'test_version')
        upload_id = published.upload_id
        calc = Calc.objects(upload_id=upload_id).first()
        assert calc.metadata['nomad_version'] != 'test_version'

        result = click.testing.CliRunner().invoke(
            cli, ['admin', 'uploads', 're-process', '--parallel', '2', upload_id], catch_exceptions=False, obj=utils.POPO())

        assert result.exit_code == 0
        assert 're-processing' in result.stdout
        calc.reload()
        assert calc.metadata['nomad_version'] == 'test_version'


@pytest.mark.usefixtures('reset_config')
class TestClient:

    def test_upload(self, test_user_bravado_client, non_empty_example_upload, proc_infra):
        result = click.testing.CliRunner().invoke(
            cli,
            ['client', 'upload', '--offline', '--name', 'test_upload', non_empty_example_upload],
            catch_exceptions=False, obj=utils.POPO())

        assert result.exit_code == 0
        assert '1/0/1' in result.output
        assert proc.Upload.objects(name='test_upload').first() is not None

    def test_local(self, client, published, admin_user_bravado_client, monkeypatch):
        def requests_get(url, stream, headers):
            assert stream
            rv = client.get(url[url.index('/raw'):], headers=headers)
            assert rv.status_code == 200
            return utils.POPO(iter_content=lambda *args, **kwargs: [bytes(rv.data)])

        monkeypatch.setattr('requests.get', requests_get)
        result = click.testing.CliRunner().invoke(
            cli,
            ['client', 'local', '%s/%s' % (published.upload_id, list(published.calcs)[0].calc_id)],
            catch_exceptions=False, obj=utils.POPO())

        assert result.exit_code == 0

    def test_mirror_dry(self, published, admin_user_bravado_client, monkeypatch):
        monkeypatch.setattr('nomad.cli.client.mirror.__in_test', True)

        result = click.testing.CliRunner().invoke(
            cli, ['client', 'mirror', '--dry'], catch_exceptions=False, obj=utils.POPO())

        assert result.exit_code == 0
        assert published.upload_id in result.output
        assert published.upload_files.os_path in result.output

    @pytest.mark.parametrize('move', [True, False])
    def test_mirror(self, published, admin_user_bravado_client, monkeypatch, move):
        ref_search_results = search.SearchRequest().search_parameters(upload_id=published.upload_id).execute_paginated()['results'][0]

        monkeypatch.setattr('nomad.cli.client.mirror.__in_test', True)

        if move:
            result = click.testing.CliRunner().invoke(
                cli, ['client', 'mirror', '--move'], catch_exceptions=False, obj=utils.POPO())
        else:
            result = click.testing.CliRunner().invoke(
                cli, ['client', 'mirror', '--source-mapping', '.volumes/test_fs:.volumes/test_fs'], catch_exceptions=False, obj=utils.POPO())

        assert result.exit_code == 0
        assert published.upload_id in result.output
        assert published.upload_files.os_path in result.output
        assert proc.Upload.objects(upload_id=published.upload_id).count() == 1
        assert proc.Calc.objects(upload_id=published.upload_id).count() == 1
        new_search = search.SearchRequest().search_parameters(upload_id=published.upload_id).execute_paginated()
        calcs_in_search = new_search['pagination']['total']
        assert calcs_in_search == 1

        new_search_results = new_search['results'][0]
        for key in new_search_results.keys():
            if key not in ['upload_time']:  # There is a sub second change due to date conversions (?)
                assert json.dumps(new_search_results[key]) == json.dumps(ref_search_results[key])

        published.upload_files.exists

    def test_mirror_files_only(self, published, admin_user_bravado_client, monkeypatch):
        monkeypatch.setattr('nomad.cli.client.mirror.__in_test', True)

        result = click.testing.CliRunner().invoke(
            cli, ['client', 'mirror', '--files-only'], catch_exceptions=False, obj=utils.POPO())

        assert result.exit_code == 0, result.output
        assert published.upload_id in result.output
        assert published.upload_files.os_path in result.output

        published.upload_files.exists
