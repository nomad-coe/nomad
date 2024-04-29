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

"""
This module describes all configurable parameters for the nomad python code. The
configuration is used for all executed python code including API, worker, CLI, and other
scripts. To use the configuration in your own scripts or new modules, simply import
this module.

All parameters are structured into objects for two reasons. First, to have
categories. Second, to allow runtime manipulation that is not effected
by python import logic. The categories are choosen along infrastructure components:
``mongo``, ``elastic``, etc.

This module also provides utilities to read the configuration from environment variables
and .yaml files. This is done automatically on import. The precedence is env over .yaml
over defaults.
"""

import os
import sys
import yaml
import logging
import os.path
from typing import Dict, Any
from nomad.config.models.config import Config

# use std python logger, since logging is not configured while loading configuration
logger = logging.getLogger(__name__)


def _load_config_yaml() -> Dict[str, Any]:
    """
    Loads the configuration from a YAML file.
    """
    config_file = os.environ.get('NOMAD_CONFIG', 'nomad.yaml')
    config_data: Dict[str, Any] = {}

    if os.path.exists(config_file):
        with open(config_file, 'r') as stream:
            try:
                config_data = yaml.load(stream, Loader=yaml.SafeLoader)
            except yaml.YAMLError as e:
                logger.error(f'cannot read nomad config: {e}')

    return config_data


def _load_config_env() -> Dict[str, Any]:
    """
    Loads the configuration from environment variables.

    TODO: The current syntax that uses underscores to separate different fields
    can lead to ambiguities. The models also cannot be used to intelligently
    decide the target field, because the model may contain discriminated unions,
    in which case the model __fields__ is undefined. This function simply splits
    the path at the first underscore and uses the first part as target section
    and the second part as field name.
    """

    def add_deep(data, path, value):
        parts = path.split('_', 1)
        root = data
        for i, part in enumerate(parts):
            if i == len(parts) - 1:
                root[part] = value
            else:
                new = root.get(part, {})
                root[part] = new
                root = new

    config_data: Dict[str, Any] = {}
    prefix = 'NOMAD_'
    for key, value in os.environ.items():
        if key.startswith(prefix) and key != 'NOMAD_CONFIG':
            add_deep(config_data, key[len(prefix) :].lower(), value)

    return config_data


def _merge(*args) -> Dict[str, Any]:
    """
    Recursively merge the given dictionaries one by one.

    When two dictionaries have conflicting keys, the values are overwritten,
    except when dealing with dictionaries, in which case they are combined.
    """

    def merge_dicts(dict1, dict2):
        merged = dict(dict1)  # Make a shallow copy of dict1

        for key, value in dict2.items():
            if key in merged:
                if isinstance(merged[key], dict) and isinstance(value, dict):
                    merged[key] = merge_dicts(merged[key], value)
                else:
                    merged[key] = value
            else:
                merged[key] = value

        return merged

    root = args[0]
    for config in args[1:]:
        root = merge_dicts(root, config)
    return root


_plugins = None


def load_config() -> Config:
    """Custom config loader. Used instead of Pydantic BaseSettings because of
    custom merging logic and custom loading of environment variables.
    """
    with open(os.path.join(os.path.dirname(__file__), 'defaults.yaml'), 'r') as stream:
        config_default = yaml.load(stream, Loader=yaml.SafeLoader)
    config_yaml = _load_config_yaml()
    config_env = _load_config_env()
    config_final = _merge(config_default, config_yaml, config_env)

    # The plugin config is stored for later when it is lazy-loaded.
    global _plugins
    _plugins = config_final['plugins']
    del config_final['plugins']

    return Config.parse_obj(config_final)


config = load_config()

# Expose config fields under this module for backwards compatibility
_module = sys.modules[__name__]
_fields = Config.__fields__
for field_name in _fields.keys():
    setattr(_module, field_name, getattr(config, field_name))
