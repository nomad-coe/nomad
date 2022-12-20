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
import os
import yaml

from nomad import config

from .utils import assert_log


@pytest.fixture
def with_config():
    old_values = config.fs.public, config.fs.archive_version_suffix
    yield config
    config.fs.public, config.fs.archive_version_suffix = old_values


def test_apply(with_config, caplog):
    config._apply('fs_public', 'test_value')
    assert config.fs.public == 'test_value'

    config._apply('fs_archive_version_suffix', 'test_value')
    assert config.fs.archive_version_suffix == 'test_value'

    config._apply('does_not_exist', 'test_value')
    assert_log(caplog, 'ERROR', 'does_not_exist does not exist')

    config._apply('fs_does_not_exist', 'test_value')
    assert_log(caplog, 'ERROR', 'fs_does_not_exist does not exist')

    config._apply('services_max_entry_download', 'not_a_number')
    assert_log(caplog, 'ERROR', 'cannot set')

    config._apply('nounderscore', 'test_value')
    assert_log(caplog, 'ERROR', 'nounderscore does not exist')


def test_env(with_config, monkeypatch):
    monkeypatch.setattr('os.environ', dict(NOMAD_FS_PUBLIC='test_value'))
    os.environ['NOMAD_FS_PUBLIC'] = 'test_value'
    config._apply_env_variables()
    assert config.fs.public == 'test_value'


def test_nomad_yaml(raw_files, with_config, monkeypatch, caplog):
    config_data = {
        'fs': {
            'public': 'test_value',
            'archive_version_suffix': 'test_value',
            'does_not_exist': 'test_value'
        },
        'does_not_exist': 'test_value',
        'services': {
            'max_entry_download': 'not_a_number'
        }
    }

    test_nomad_yaml = os.path.join(config.fs.tmp, 'nomad_test.yaml')
    monkeypatch.setattr('os.environ', dict(NOMAD_CONFIG=test_nomad_yaml))
    with open(test_nomad_yaml, 'w') as file:
        yaml.dump(config_data, file)

    config.load_config()

    os.remove(test_nomad_yaml)

    assert config.fs.public == 'test_value'
    assert config.fs.archive_version_suffix == 'test_value'
    assert_log(caplog, 'ERROR', 'does_not_exist does not exist')
    assert_log(caplog, 'ERROR', 'fs_does_not_exist does not exist')
    assert_log(caplog, 'ERROR', 'cannot set')
