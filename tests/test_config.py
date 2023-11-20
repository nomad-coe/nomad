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
from pydantic import parse_obj_as

from nomad import config
from nomad.utils import flatten_dict

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
    assert_log(caplog, 'ERROR', 'config key does not exist: does_not_exist')

    config._apply('fs_does_not_exist', 'test_value')
    assert_log(caplog, 'ERROR', 'config key does not exist: fs_does_not_exist')

    config._apply('services_max_entry_download', 'not_a_number')
    assert_log(caplog, 'ERROR', 'cannot set')

    config._apply('nounderscore', 'test_value')
    assert_log(caplog, 'ERROR', 'config key does not exist: nounderscore')


def test_env(with_config, monkeypatch):
    monkeypatch.setattr('os.environ', dict(NOMAD_FS_PUBLIC='test_value'))
    os.environ['NOMAD_FS_PUBLIC'] = 'test_value'
    config._apply_env_variables()
    assert config.fs.public == 'test_value'


def load_config(config_dict, monkeypatch):
    '''Loads the given dictionary into the current config.
    '''
    test_nomad_yaml = os.path.join(config.fs.tmp, 'nomad_test.yaml')
    monkeypatch.setattr('os.environ', dict(NOMAD_CONFIG=test_nomad_yaml))
    with open(test_nomad_yaml, 'w') as file:
        yaml.dump(config_dict, file)
    config.load_config()
    os.remove(test_nomad_yaml)


@pytest.mark.parametrize('config_dict', [
    pytest.param({'services': {'max_entry_download': 123}}, id='set integer'),
    pytest.param({'services': {'api_host': 'http://myhost'}}, id='set string'),
    pytest.param({'ui': {'theme': {'title': 'mytitle'}}}, id='set nested'),
])
def test_config(raw_files, with_config, monkeypatch, caplog, config_dict):
    load_config(config_dict, monkeypatch)
    for key, value in config_dict.items():
        flat_real = flatten_dict(getattr(config, key).dict())
        flat_expected = flatten_dict(value)
        for key2, value2 in flat_expected.items():
            assert flat_real[key2] == value2


@pytest.mark.parametrize('config_dict, error', [
    pytest.param({'does_not_exist': 'test_value'}, 'config key does not exist: does_not_exist', id='undefined field root'),
    pytest.param({'fs': {'does_not_exist': 'test_value'}}, 'config key does not exist: fs_does_not_exist', id='undefined field child'),
    pytest.param(
        {'services': {'max_entry_download': 'not_a_number'}},
        "cannot set config setting services_max_entry_download=not_a_number: 1 validation error for ParsingModel[int]",
        id='incompatible-value'
    ),
])
def test_config_error(raw_files, with_config, monkeypatch, caplog, config_dict, error):
    load_config(config_dict, monkeypatch)
    assert_log(caplog, 'ERROR', error)


def test_parser_plugins():
    from nomad import config
    from nomad.config import Parser
    parsers = [
        plugin for plugin in config.plugins.options.values()
        if isinstance(plugin, Parser)]
    assert len(parsers) == 69


def test_plugin_polymorphism():
    plugins_yaml = parse_obj_as(config.Plugins, yaml.safe_load('''
        options:
            schema:
                plugin_type: schema
                name: test
                python_package: nomad.datamodel.metainfo.simulation
    '''))

    plugins_config = config.Plugins(options=dict(parser=config.Parser(
        name='parsers/abinit',
        python_package='electronicparsers.abinit',
        parser_class_name='electronicparsers.abinit.parser.AbinitParser'
    )))

    from nomad.config import _merge
    _merge(plugins_config.options, plugins_yaml.options)
    plugins = plugins_config

    assert isinstance(plugins.options['schema'], config.Schema)
    assert isinstance(plugins.options['parser'], config.Parser)
