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

import os
import re

import pytest
import yaml
from pydantic import ValidationError

from nomad.config import load_config
from nomad.config.models.plugins import ParserEntryPoint, SchemaPackageEntryPoint
from nomad.utils import flatten_dict

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
            'The following unsupported keys were found in the nomad configuration '
            '(nomad.yaml, defaults.yaml or environment variables): FS: "does"',
            ['yaml', 'env'],
            id='non-existing nested field',
        ),
        pytest.param(
            {'does': 'not exist'},
            'The following unsupported keys were found in the nomad configuration '
            '(nomad.yaml, defaults.yaml or environment variables): Config: "does"',
            ['yaml'],
            id='non-existing top-level field',
        ),
    ],
)
@pytest.mark.parametrize('format', ['yaml', 'env'])
def test_config_warning(
    config_dict,
    format,
    warning,
    formats_with_warning,
    log_output,
    mockopen,
    monkeypatch,
):
    """Tests that extra fields create a warning message."""
    conf_yaml, conf_env = load_format(config_dict, format)
    load_test_config(conf_yaml, conf_env, mockopen, monkeypatch)
    assert_log(
        log_output, 'WARNING', warning, negate=format not in formats_with_warning
    )


@pytest.mark.parametrize(
    'config_dict, error',
    [
        pytest.param(
            {'services': {'api_timeout': 'not_a_number'}},
            (
                '1 validation error for Config\nservices.api_timeout\n  '
                'Input should be a valid integer, unable to parse string as an '
                "integer [type=int_parsing, input_value='not_a_number', input_type=str]"
            ),
            id='invalid type',
        ),
    ],
)
@pytest.mark.parametrize('format', ['yaml', 'env'])
def test_config_error(config_dict, format, error, mockopen, monkeypatch):
    """Tests that validation errors raise exceptions."""
    conf_yaml, conf_env = load_format(config_dict, format)
    with pytest.raises(ValidationError, match=re.escape(error)):
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
                # TODO: the NORTH part is only needed to reset the tool
                'north': {'tools': {'include': []}},
                'plugins': {
                    'options': {
                        'dosnormalizer:dos_normalizer_entry_point': {
                            'name': 'yaml',
                        }
                    }
                },
            },
            {'plugins': {'include': ['dosnormalizer:dos_normalizer_entry_point']}},
            {
                'plugins': {
                    'entry_points': {
                        'include': ['dosnormalizer:dos_normalizer_entry_point'],
                        'options': {
                            'dosnormalizer:dos_normalizer_entry_point': {
                                'name': 'yaml',
                                'plugin_package': 'dosnormalizer',
                                'description': 'Normalizer for the DOS data.',
                                'entry_point_type': 'normalizer',
                            }
                        },
                    }
                }
            },
            id='dictionary: merges',
        ),
        pytest.param(
            {'north': {'tools': {'include': []}}, 'plugins': {'include': ['a']}},
            {'plugins': {'include': ['b']}},
            {'plugins': {'entry_points': {'include': ['b']}}},
            id='list: overrides',
        ),
        pytest.param(
            {'services': {'api_timeout': 100}},
            {'services': {'api_timeout': 200}},
            {'services': {'api_timeout': 200}},
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
                'north': {'tools': {'include': []}},
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
                'north': {'tools': {'include': []}},
                'plugins': {
                    'entry_points': {'include': ['b'], 'exclude': ['b']},
                },
            },
            {'plugins': {'entry_points': {'include': ['b'], 'exclude': ['b']}}},
            id='only new values',
        ),
        pytest.param(
            {
                'north': {'tools': {'include': []}},
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
                'north': {'tools': {'include': []}},
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
                'north': {'tools': {'include': []}},
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
                'north': {'tools': {'include': []}},
                'plugins': {
                    'options': {
                        'electronicparsers:vasp_parser_entry_point': {
                            'mainfile_name_re': 'a'
                        }
                    },
                    'entry_points': {
                        'options': {
                            'electronicparsers:vasp_parser_entry_point': {
                                'mainfile_name_re': 'b'
                            }
                        }
                    },
                },
            },
            {
                'north': {'tools': {'include': []}},
                'plugins': {
                    'entry_points': {
                        'options': {
                            'electronicparsers:vasp_parser_entry_point': {
                                'mainfile_name_re': 'a',
                                'plugin_package': 'electronicparsers',
                            }
                        }
                    }
                },
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
        if isinstance(entry_point, ParserEntryPoint)
    ]
    assert len(parsers) == 68


def test_plugin_polymorphism(mockopen, monkeypatch):
    plugins = {
        'plugins': {
            'options': {
                'schema': {
                    'entry_point_type': 'schema_package',
                    'name': 'test',
                    'plugin_package': 'runschema',
                },
                'parser': {
                    'entry_point_type': 'parser',
                    'name': 'parsers/abinit',
                    'plugin_package': 'electronicparsers',
                },
            }
        }
    }
    config = load_test_config(plugins, None, mockopen, monkeypatch)
    config.load_plugins()
    assert isinstance(
        config.plugins.entry_points.options['schema'], SchemaPackageEntryPoint
    )
    assert isinstance(config.plugins.entry_points.options['parser'], ParserEntryPoint)


@pytest.mark.parametrize(
    'conf_yaml, conf_expected',
    [
        pytest.param(
            None,
            {
                'uploads': {
                    'pagination': {
                        'page_size': 10,
                        'order_by': 'upload_create_time',
                        'order': 'desc',
                    },
                    'entries': {
                        'pagination': {
                            'page_size': 5,
                            'order_by': 'process_status',
                            'order': 'asc',
                        }
                    },
                }
            },
            id='default pagination values',
        ),
        pytest.param(
            {
                'uploads': {
                    'pagination': {
                        'page_size': 12,
                        'order_by': 'upload_id',
                        'order': 'asc',
                    },
                    'entries': {
                        'pagination': {
                            'page_size': 6,
                            'order_by': 'entry_id',
                            'order': 'desc',
                        }
                    },
                }
            },
            {
                'uploads': {
                    'pagination': {
                        'page_size': 12,
                        'order_by': 'upload_id',
                        'order': 'asc',
                    },
                    'entries': {
                        'pagination': {
                            'page_size': 6,
                            'order_by': 'entry_id',
                            'order': 'desc',
                        }
                    },
                }
            },
            id='all pagination values changed',
        ),
        pytest.param(
            {'projects': {'pagination': {'page_size': 12}}},
            {'uploads': {'pagination': {'page_size': 12}}},
            id='using alias projects',
        ),
    ],
)
def test_pagination(conf_yaml, conf_expected, mockopen, monkeypatch):
    config = load_test_config(conf_yaml, None, mockopen, monkeypatch)
    assert_config(config, conf_expected)


@pytest.mark.parametrize(
    'entry_point_id, custom_id_url_safe, expected_id_url_safe, should_raise',
    [
        pytest.param(
            'nomad_parser_vasp.parsers:VASPRunParser',
            None,
            'nomad_parser_vasp.parsers-VASPRunParser',
            False,
            id='auto-generate: typical entry point id',
        ),
        pytest.param(
            'simple_id',
            None,
            'simple_id',
            False,
            id='auto-generate: simple id unchanged',
        ),
        pytest.param(
            'some.entry:point',
            'my_custom_id',
            'my_custom_id',
            False,
            id='custom: valid with underscores',
        ),
        pytest.param(
            'some.entry:point',
            'my-custom-id',
            'my-custom-id',
            False,
            id='custom: valid with hyphens',
        ),
        pytest.param(
            'some.entry:point',
            'test%20string',
            'test%20string',
            False,
            id='custom: valid with percent-encoding',
        ),
        pytest.param(
            'some.entry:point',
            '_private_id',
            '_private_id',
            False,
            id='custom: valid starting with underscore',
        ),
        pytest.param(
            'some.entry:point',
            'id123',
            'id123',
            False,
            id='custom: valid with numbers',
        ),
        pytest.param(
            'some.entry:point',
            'valid.id',
            'valid.id',
            False,
            id='custom: valid with dots',
        ),
        pytest.param(
            'some.entry:point',
            'invalid id',
            None,
            True,
            id='custom: invalid with spaces',
        ),
        pytest.param(
            'some.entry:point',
            'invalid:id',
            None,
            True,
            id='custom: invalid with colons',
        ),
    ],
)
def test_id_url_safe(
    entry_point_id,
    custom_id_url_safe,
    expected_id_url_safe,
    should_raise,
    mockopen,
    monkeypatch,
):
    """Tests URL-safe identifier generation, validation, and assignment."""
    config_dict = {
        'plugins': {
            'options': {
                entry_point_id: {
                    'entry_point_type': 'parser',
                    'id_url_safe': custom_id_url_safe,
                }
            },
        }
    }
    conf_yaml, conf_env = load_format(config_dict, 'yaml')
    config = load_test_config(conf_yaml, conf_env, mockopen, monkeypatch)

    if should_raise:
        with pytest.raises(ValueError):
            config.load_plugins()
    else:
        config.load_plugins()
        assert (
            config.plugins.entry_points.options[entry_point_id].id_url_safe
            == expected_id_url_safe
        )


@pytest.mark.parametrize(
    'entry_points, collides',
    [
        pytest.param(
            {
                'options': {
                    'pkg.mod:Class': {
                        'id_url_safe': 'A',
                        'entry_point_type': 'schema_package',
                    },
                    'pkg-mod_Class': {
                        'id_url_safe': 'A',
                        'entry_point_type': 'schema_package',
                    },
                }
            },
            True,
            id='same custom id_url_safe and entry point type',
        ),
        pytest.param(
            {
                'options': {
                    'pkg.mod:Class': {
                        'id_url_safe': 'A',
                        'entry_point_type': 'parser',
                    },
                    'pkg-mod_Class': {
                        'id_url_safe': 'A',
                        'entry_point_type': 'schema_package',
                    },
                }
            },
            False,
            id='same custom id_url_safe, different entry point types',
        ),
        pytest.param(
            {
                'exclude': ['pkg.mod:Class'],
                'options': {
                    'pkg.mod:Class': {
                        'id_url_safe': 'A',
                        'entry_point_type': 'schema_package',
                    },
                    'pkg-mod_Class': {
                        'id_url_safe': 'A',
                        'entry_point_type': 'schema_package',
                    },
                },
            },
            False,
            id='no clash if inactive',
        ),
    ],
)
def test_id_url_safe_collision(entry_points, collides, mockopen, monkeypatch):
    """Tests that URL-safe identifier collisions are detected."""
    config_dict = {
        'plugins': {
            'entry_points': entry_points,
        }
    }

    conf_yaml, conf_env = load_format(config_dict, 'yaml')
    config = load_test_config(conf_yaml, conf_env, mockopen, monkeypatch)

    if collides:
        with pytest.raises(ValueError):
            config.load_plugins()
    else:
        config.load_plugins()


@pytest.mark.parametrize(
    'conf_yaml, conf_expected',
    [
        pytest.param(
            {'keycloak': {'server_url': 'http://example.com/auth'}},
            {'keycloak': {'server_url': 'http://example.com/auth'}},
            id='keycloak-no-slash',
        ),
        pytest.param(
            {'keycloak': {'server_url': 'http://example.com/auth/'}},
            {'keycloak': {'server_url': 'http://example.com/auth'}},
            id='keycloak-with-slash',
        ),
    ],
)
def test_normalized_url(conf_yaml, conf_expected, mockopen, monkeypatch):
    config = load_test_config(conf_yaml, None, mockopen, monkeypatch)
    assert_config(config, conf_expected)


# Tests for `Auth`


@pytest.mark.parametrize(
    ('conf_yaml', 'conf_expected'),
    [
        pytest.param(
            {'auth': {'authorized_users': None}},
            {'auth': {'authorized_users': None}},
            id='none',
        ),
        pytest.param(
            {'auth': {'authorized_users': []}},
            {'auth': {'authorized_users': []}},
            id='empty',
        ),
        pytest.param(
            {'auth': {'authorized_users': ['alice', 'bob@example.com']}},
            {'auth': {'authorized_users': ['alice', 'bob@example.com']}},
            id='already-normalized',
        ),
        pytest.param(
            {
                'auth': {
                    'authorized_users': [' Alice ', 'BOB@example.com ', '  CHARLIE  ']
                }
            },
            {'auth': {'authorized_users': ['alice', 'bob@example.com', 'charlie']}},
            id='strip-and-lower',
        ),
        pytest.param(
            {'auth': {'authorized_users': ['Alice', 'alice ', ' ALICE']}},
            {'auth': {'authorized_users': ['alice', 'alice', 'alice']}},
            id='duplicates-not-deduplicated',
        ),
    ],
)
def test_authorized_users(conf_yaml, conf_expected, mockopen, monkeypatch):
    config = load_test_config(conf_yaml, None, mockopen, monkeypatch)
    assert_config(config, conf_expected)
