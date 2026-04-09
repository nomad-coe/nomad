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
by python import logic. The categories are chosen along infrastructure components:
``mongo``, ``elastic``, etc.

This module also provides utilities to read the configuration from environment variables
and .yaml files. This is done automatically on import. The precedence is env over .yaml
over defaults.
"""

import json
import logging
import os
import sys
from typing import Any

import yaml

from nomad.config.models.config import Config

# use std python logger, since logging is not configured while loading configuration
logger = logging.getLogger(__name__)


def _load_config_yaml(files: list[str] | None = None) -> dict[str, Any]:
    """
    Loads the configuration from one or more YAML files. Files are merged in order,
    with later files overwriting values from earlier ones.
    """
    if not files:
        config_file = os.environ.get('NOMAD_CONFIG', 'nomad.yaml')
        if os.path.exists(config_file):
            files_to_load = [config_file]
        else:
            files_to_load = []
    else:
        files_to_load = files

    final_config_data: dict[str, Any] = {}
    for config_file in files_to_load:
        with open(config_file) as stream:
            try:
                config_data = yaml.load(stream, Loader=yaml.SafeLoader)
                if config_data:
                    final_config_data = _merge(final_config_data, config_data)
            except yaml.YAMLError as e:
                logger.error(f'Cannot read nomad config file {config_file}: {e}')

    return final_config_data


def _load_config_env() -> dict[str, Any]:
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

    config_data: dict[str, Any] = {}
    prefix = 'NOMAD_'
    for key, value in os.environ.items():
        if key == 'NOMAD_CONFIG' or not key.startswith(prefix):
            continue

        key = key[len(prefix) :].lower()
        # Some environment variables starting with NOMAD_ are unavoidable
        # in docker/kubernetes environments. We should ignore them here,
        # before they cause a warning later when the config is validated.
        if all([not key.startswith(field) for field in Config.model_fields.keys()]):
            continue

        try:
            value = json.loads(value)
        except json.decoder.JSONDecodeError:
            pass
        add_deep(config_data, key, value)

    return config_data


def _merge(*args) -> dict[str, Any]:
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
        if config:
            root = merge_dicts(root, config)
    return root


_plugins = None


def load_config(files: list[str] | None = None) -> Config:
    """Custom config loader. Used instead of Pydantic BaseSettings because of
    custom merging logic and custom loading of environment variables.
    """
    with open(os.path.join(os.path.dirname(__file__), 'defaults.yaml')) as stream:
        config_default = yaml.load(stream, Loader=yaml.SafeLoader)
    config_yaml = _load_config_yaml(files)
    config_env = _load_config_env()
    config_final = _merge(config_default, config_yaml, config_env)

    # The plugin config is stored for later when it is lazy-loaded.
    global _plugins
    _plugins = config_final['plugins']
    del config_final['plugins']

    validated = Config.model_validate(config_final)
    validated.archive.initialize()

    return validated


def load_and_set_config(files: list[str] | None = None) -> Config:
    """
    Loads the configuration from the specified files and updates the global 'config'
    object and module-level attributes.

    This function is necessary for runtime configuration changes, e.g., via CLI flags,
    as it ensures all parts of the application see the updated configuration.
    """
    new_config = load_config(files=files)
    globals()['config'] = new_config

    _module = sys.modules[__name__]
    _fields = Config.model_fields
    for field_name in _fields.keys():
        setattr(_module, field_name, getattr(new_config, field_name))

    return new_config


config = load_config()

# Expose config fields under this module for backwards compatibility
_module = sys.modules[__name__]
_fields = Config.model_fields
for field_name in _fields.keys():
    setattr(_module, field_name, getattr(config, field_name))
