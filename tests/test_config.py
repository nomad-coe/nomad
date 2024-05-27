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
import yaml
import os

from nomad.utils import flatten_dict
from nomad.config import load_config
from nomad.config.models.plugins import Schema, Parser

from pydantic import ValidationError

from .utils import assert_log


def load_test_config(conf_yaml, conf_env, mockopen=None, monkeypatch=None):
    if conf_env:
        monkeypatch.setattr('os.environ', conf_env)
    config_file = os.environ.get('NOMAD_CONFIG', 'nomad.yaml')
    if conf_yaml:
        mockopen.write(config_file, yaml.dump(conf_yaml))
        old = os.path.exists
        monkeypatch.setattr(
            'os.path.exists', lambda x: True if x == config_file else old(x)
        )
    return load_config()


def get_config_env(config):
    if not config:
        return {}
    return {
        f'NOMAD_{key.upper()}': value
        for key, value in flatten_dict(config, '_').items()
    }


def load_format(config, format):
    conf_yaml = config if format == 'yaml' else None
    conf_env = get_config_env(config) if format == 'env' else None
    return conf_yaml, conf_env


def assert_config(config, config_expected):
    flattened = flatten_dict(config_expected)
    for key, val in flattened.items():
        root = config
        for part in key.split('.'):
            if isinstance(root, dict):
                root = root[part]
            else:
                root = getattr(root, part)
        assert root == val


def test_config_file_change(mockopen, monkeypatch):
    """Tests that changing the config file path works."""
    conf_yaml = {'fs': {'public': 'test'}}
    conf_env = {'NOMAD_CONFIG': 'test.yaml'}
    config = load_test_config(conf_yaml, conf_env, mockopen, monkeypatch)
    assert_config(config, conf_yaml)


@pytest.mark.parametrize(
    'config_dict',
    [
        pytest.param({'fs': {'public': 'test'}}, id='nested string'),
        pytest.param({'north': {'hub_ip': '1.2.3.4'}}, id='underscore in field name'),
    ],
)
@pytest.mark.parametrize('format', ['yaml', 'env'])
def test_config_success(config_dict, format, mockopen, monkeypatch):
    """Tests that config variables are correctly loaded."""
    conf_yaml, conf_env = load_format(config_dict, format)
    config_obj = load_test_config(conf_yaml, conf_env, mockopen, monkeypatch)
    assert_config(config_obj, config_dict)


@pytest.mark.parametrize(
    'config_dict, warning, formats_with_warning',
    [
        pytest.param(
            {'fs': {'does': 'not exist'}},
            'The following unsupported keys were found in your configuration, e.g. nomad.yaml: "does".',
            ['yaml', 'env'],
            id='non-existing nested field',
        ),
        pytest.param(
            {'does': 'not exist'},
            'The following unsupported keys were found in your configuration, e.g. nomad.yaml: "does".',
            ['yaml'],
            id='non-existing top-level field',
        ),
    ],
)
@pytest.mark.parametrize('format', ['yaml', 'env'])
def test_config_warning(
    config_dict, format, warning, formats_with_warning, caplog, mockopen, monkeypatch
):
    """Tests that extra fields create a warning message."""
    conf_yaml, conf_env = load_format(config_dict, format)
    load_test_config(conf_yaml, conf_env, mockopen, monkeypatch)
    assert_log(caplog, 'WARNING', warning, negate=format not in formats_with_warning)


@pytest.mark.parametrize(
    'config_dict, error',
    [
        pytest.param(
            {'celery': {'timeout': 'not_a_number'}},
            (
                r'1 validation error for Config\n'
                r'celery -> timeout\n'
                r'  value is not a valid integer \(type=type_error\.integer\)'
            ),
            id='invalid type',
        ),
    ],
)
@pytest.mark.parametrize('format', ['yaml', 'env'])
def test_config_error(config_dict, format, error, mockopen, monkeypatch):
    """Tests that validation errors raise exceptions."""
    conf_yaml, conf_env = load_format(config_dict, format)
    with pytest.raises(ValidationError, match=error):
        load_test_config(conf_yaml, conf_env, mockopen, monkeypatch)


@pytest.mark.parametrize(
    'conf_yaml, conf_env, value',
    [
        pytest.param(None, None, '.volumes/fs/public', id='default'),
        pytest.param(
            {'fs': {'public': 'yaml'}}, {}, 'yaml', id='yaml overrides default'
        ),
        pytest.param(
            None, {'NOMAD_FS_PUBLIC': 'env'}, 'env', id='env overrides default'
        ),
        pytest.param(
            {'fs': {'public': 'yaml'}},
            {'NOMAD_FS_PUBLIC': 'env'},
            'env',
            id='env overrides yaml and default',
        ),
    ],
)
def test_config_priority(conf_yaml, conf_env, value, mockopen, monkeypatch):
    """Tests that the priority between model defaults, yaml and environment
    variables is correctly handled."""
    config = load_test_config(conf_yaml, conf_env, mockopen, monkeypatch)
    assert config.fs.public == value


@pytest.mark.parametrize(
    'conf_yaml, conf_env, conf_expected',
    [
        pytest.param(
            {
                'plugins': {
                    'options': {
                        'normalizers/simulation/dos': {
                            'name': 'yaml',
                        }
                    }
                }
            },
            {'plugins': {'include': ['normalizers/simulation/dos']}},
            {
                'plugins': {
                    'entry_points': {
                        'include': ['normalizers/simulation/dos'],
                        'options': {
                            'normalizers/simulation/dos': {
                                'name': 'yaml',
                                'python_package': 'dosnormalizer',
                                'description': 'This is the normalizer for DOS in NOMAD.\n',
                                'plugin_documentation_url': None,
                                'plugin_source_code_url': None,
                                'normalizer_class_name': 'dosnormalizer.DosNormalizer',
                                'plugin_type': 'normalizer',
                            }
                        },
                    }
                }
            },
            id='dictionary: merges',
        ),
        pytest.param(
            {'plugins': {'include': ['a']}},
            {'plugins': {'include': ['b']}},
            {'plugins': {'entry_points': {'include': ['b']}}},
            id='list: overrides',
        ),
        pytest.param(
            {'celery': {'timeout': 100}},
            {'celery': {'timeout': 200}},
            {'celery': {'timeout': 200}},
            id='scalar: overrides',
        ),
    ],
)
def test_config_merge(conf_yaml, conf_env, conf_expected, mockopen, monkeypatch):
    """Tests that configs are correctly merged: dictionaries should be merged,
    everything else overridden."""
    config = load_test_config(
        conf_yaml, get_config_env(conf_env), mockopen, monkeypatch
    )
    config.load_plugins()
    assert_config(config, conf_expected)


@pytest.mark.parametrize(
    'conf_yaml, conf_expected',
    [
        pytest.param(
            {
                'plugins': {
                    'include': ['a'],
                    'exclude': ['a'],
                },
            },
            {'plugins': {'entry_points': {'include': ['a'], 'exclude': ['a']}}},
            id='only old values',
        ),
        pytest.param(
            {
                'plugins': {
                    'entry_points': {'include': ['b'], 'exclude': ['b']},
                },
            },
            {'plugins': {'entry_points': {'include': ['b'], 'exclude': ['b']}}},
            id='only new values',
        ),
        pytest.param(
            {
                'plugins': {
                    'include': ['a'],
                    'exclude': ['a'],
                    'entry_points': {'include': ['b'], 'exclude': ['b']},
                },
            },
            {'plugins': {'entry_points': {'include': ['a'], 'exclude': ['a']}}},
            id='old include and exclude have precedence: non-empty lists',
        ),
        pytest.param(
            {
                'plugins': {
                    'include': [],
                    'exclude': [],
                    'entry_points': {'include': ['b'], 'exclude': ['b']},
                },
            },
            {'plugins': {'entry_points': {'include': [], 'exclude': []}}},
            id='old include and exclude have precedence: empty lists',
        ),
        pytest.param(
            {
                'plugins': {
                    'include': None,
                    'exclude': None,
                    'entry_points': {'include': ['b'], 'exclude': ['b']},
                },
            },
            {'plugins': {'entry_points': {'include': None, 'exclude': None}}},
            id='old include and exclude have precedence: None',
        ),
        pytest.param(
            {
                'plugins': {
                    'options': {'parsers/vasp': {'mainfile_name_re': 'a'}},
                    'entry_points': {
                        'options': {'parsers/vasp': {'mainfile_name_re': 'b'}}
                    },
                },
            },
            {
                'plugins': {
                    'entry_points': {
                        'options': {
                            'parsers/vasp': {
                                'mainfile_name_re': 'a',
                                'python_package': 'electronicparsers.vasp',
                            }
                        }
                    }
                }
            },
            id='old, new and default options are merged with old config having precendence over new values.',
        ),
    ],
)
def test_plugin_entry_points(conf_yaml, conf_expected, mockopen, monkeypatch):
    """Tests that any conflicts between old and new plugin configs are resolved
    correctly."""
    config = load_test_config(conf_yaml, None, mockopen, monkeypatch)
    config.load_plugins()
    assert_config(config, conf_expected)


def test_parser_plugins():
    config = load_config()
    config.load_plugins()
    parsers = [
        entry_point
        for entry_point in config.plugins.entry_points.options.values()
        if isinstance(entry_point, Parser)
    ]
    assert len(parsers) == 71


def test_plugin_polymorphism(mockopen, monkeypatch):
    plugins = {
        'plugins': {
            'options': {
                'schema': {
                    'plugin_type': 'schema',
                    'name': 'test',
                    'python_package': 'runschema',
                },
                'parser': {
                    'plugin_type': 'parser',
                    'name': 'parsers/abinit',
                    'python_package': 'electronicparsers.abinit',
                    'parser_class_name': 'electronicparsers.abinit.parser.AbinitParser',
                },
            }
        }
    }
    config = load_test_config(plugins, None, mockopen, monkeypatch)
    config.load_plugins()
    assert isinstance(config.plugins.entry_points.options['schema'], Schema)
    assert isinstance(config.plugins.entry_points.options['parser'], Parser)
