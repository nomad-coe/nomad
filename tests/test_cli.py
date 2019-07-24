
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

import click.testing

from nomad import utils, search
from nomad.cli import cli
from nomad.processing import Upload, Calc


# TODO there is much more to test

class TestAdmin:
    def test_clean(self, published, no_warn, reset_config):
        upload_id = published.upload_id

        Upload.objects(upload_id=upload_id).delete()
        assert published.upload_files.exists()
        assert Calc.objects(upload_id=upload_id).first() is not None
        assert search.entry_search(search_parameters=dict(upload_id=upload_id))['pagination']['total'] > 0

        result = click.testing.CliRunner().invoke(
            cli, ['admin', 'clean', '--force', '--skip-es'], catch_exceptions=False, obj=utils.POPO())

        assert result.exit_code == 0
        assert not published.upload_files.exists()
        assert Calc.objects(upload_id=upload_id).first() is None
        assert search.entry_search(search_parameters=dict(upload_id=upload_id))['pagination']['total'] > 0


class TestAdminUploads:

    def test_ls(self, published, no_warn, reset_config):
        upload_id = published.upload_id

        result = click.testing.CliRunner().invoke(
            cli, ['admin', 'uploads', 'ls', upload_id], catch_exceptions=False, obj=utils.POPO())

        assert result.exit_code == 0
        assert '1 uploads selected' in result.stdout

    def test_rm(self, published, no_warn, reset_config):
        upload_id = published.upload_id

        result = click.testing.CliRunner().invoke(
            cli, ['admin', 'uploads', 'rm', upload_id], catch_exceptions=False, obj=utils.POPO())

        assert result.exit_code == 0
        assert 'deleting' in result.stdout
        assert Upload.objects(upload_id=upload_id).first() is None
        assert Calc.objects(upload_id=upload_id).first() is None

    def test_re_process(self, published, no_warn, monkeypatch, reset_config):
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
