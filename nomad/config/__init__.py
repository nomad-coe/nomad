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

'''
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
'''

import logging
import os
import os.path
import yaml
import warnings
from typing import List, Any, Union

from nomad.config.models import (
    NomadSettings, Services, Meta, Oasis, NORTH, RabbitMQ, Celery, FS, Elastic, Keycloak, Mongo, Logstash, Tests, Mail, Normalize, Resources, Client, DataCite, GitLab, Process, Reprocess, RFC3161Timestamp, BundleExport, BundleImport, Archive, UI
)


warnings.filterwarnings('ignore', message='numpy.dtype size changed')
warnings.filterwarnings('ignore', message='numpy.ufunc size changed')
warnings.filterwarnings('ignore', category=DeprecationWarning)


def api_url(ssl: bool = True, api: str = 'api', api_host: str = None, api_port: int = None):
    '''
    Returns the url of the current running nomad API. This is for server-side use.
    This is not the NOMAD url to use as a client, use `nomad.config.client.url` instead.
    '''
    if api_port is None:
        api_port = services.api_port
    if api_host is None:
        api_host = services.api_host
    protocol = 'https' if services.https and ssl else 'http'
    host_and_port = api_host
    if api_port not in [80, 443]:
        host_and_port += ':' + str(api_port)
    base_path = services.api_base_path.strip('/')
    return f'{protocol}://{host_and_port}/{base_path}/{api}'


def gui_url(page: str = None):
    base = api_url(True)[:-3]
    if base.endswith('/'):
        base = base[:-1]

    if page is not None:
        return '%s/gui/%s' % (base, page)

    return '%s/gui' % base


def rabbitmq_url():
    return 'pyamqp://%s:%s@%s//' % (rabbitmq.user, rabbitmq.password, rabbitmq.host)


def north_url(ssl: bool = True):
    return api_url(ssl=ssl, api='north', api_host=north.hub_host, api_port=north.hub_port)


def hub_url():
    return f'http://{north.hub_host}:{north.hub_port}{services.api_base_path}/north/hub'


services = Services()
meta = Meta(deployment_url=api_url())
oasis = Oasis()
north = NORTH()
rabbitmq = RabbitMQ()
celery = Celery()
fs = FS()
elastic = Elastic()
keycloak = Keycloak()
mongo = Mongo()
logstash = Logstash()
tests = Tests()
mail = Mail()
normalize = Normalize()
resources = Resources()
client = Client()
datacite = DataCite()
gitlab = GitLab()
process = Process()
reprocess = Reprocess()
rfc3161_timestamp = RFC3161Timestamp()
bundle_export = BundleExport()
bundle_import = BundleImport()
archive = Archive()
ui = UI()


def normalize_loglevel(value, default_level=logging.INFO):
    plain_value = value
    if plain_value is None:
        return default_level
    else:
        try:
            return int(plain_value)
        except ValueError:
            return getattr(logging, plain_value)


_transformations = {
    'console_log_level': normalize_loglevel,
    'logstash_level': normalize_loglevel
}


# use std python logger, since logging is not configured while loading configuration
logger = logging.getLogger(__name__)


def _check_config():
    '''Used to check that the current configuration is valid. Should only be
    called once after the final config is loaded.

    Raises:
        AssertionError: if there is a contradiction or invalid values in the
            config file settings.
    '''
    # TODO more if this should be translated into pydantic validations.

    # The AFLOW symmetry information is checked once on import
    proto_symmetry_tolerance = normalize.prototype_symmetry_tolerance
    symmetry_tolerance = normalize.symmetry_tolerance
    if proto_symmetry_tolerance != symmetry_tolerance:
        raise AssertionError(
            "The AFLOW prototype information is outdated due to changed tolerance "
            "for symmetry detection. Please update the AFLOW prototype information "
            "by running the CLI command 'nomad admin ops prototype-update "
            "--matches-only'."
        )

    if normalize.springer_db_path and not os.path.exists(normalize.springer_db_path):
        normalize.springer_db_path = None

    if keycloak.public_server_url is None:
        keycloak.public_server_url = keycloak.server_url

    def get_external_path(path):
        if fs.external_working_directory and not os.path.isabs(path):
            return os.path.join(fs.external_working_directory, path)
        return path

    if fs.staging_external is None:
        fs.staging_external = get_external_path(fs.staging)

    if fs.public_external is None:
        fs.public_external = get_external_path(fs.public)

    if fs.north_home_external is None:
        fs.north_home_external = get_external_path(fs.north_home)

    ui.north.enabled = north.enabled


def _merge(a: Union[dict, NomadSettings], b: dict, path: List[str] = None) -> Union[dict, NomadSettings]:
    '''
    Recursively merges b into a. Will add new key-value pairs, and will
    overwrite existing key-value pairs. Notice that this mutates the original
    dictionary/model a and if you want to return a copy, you will want to first
    (deep)copy the original value.
    '''
    def has(target, key):
        return key in target if isinstance(target, dict) else hasattr(target, key)

    def set(target, key, value):
        if isinstance(target, dict):
            target[key] = value
        else:
            setattr(target, key, value)

    def get(target, key):
        return target[key] if isinstance(target, dict) else getattr(target, key)

    if path is None: path = []
    for key in b:
        value = b[key]
        if has(a, key):
            child = get(a, key)
            if isinstance(value, dict) and isinstance(child, (NomadSettings, dict)):
                _merge(child, value, path + [str(key)])
            else:
                set(a, key, value)
        else:
            set(a, key, value)
    return a


def _apply(key, value, raise_error: bool = True) -> None:
    '''
    Changes the config according to given key and value. The first part of a key
    (with ``_`` as a separator) is interpreted as a group of settings. E.g. ``fs_staging``
    leading to ``config.fs.staging``.
    '''
    full_key = key
    try:
        group_key, config_key = full_key.split('_', 1)
    except Exception:
        if raise_error:
            logger.error(f'config key does not exist: {full_key}')
        return

    current = globals()

    current_value: Any = None
    if group_key not in current:
        if key not in current:
            if raise_error:
                logger.error(f'config key does not exist: {full_key}')
            return
    else:
        current = current[group_key]
        if not isinstance(current, NomadSettings):
            if raise_error:
                logger.error(f'config key does not exist: {full_key}')
            return

        try:
            current_value = getattr(current, config_key)
        except AttributeError:
            if raise_error:
                logger.error(f'config key does not exist: {full_key}')
            return

        key = config_key

    try:
        if isinstance(value, dict):
            value = _merge(current_value, value)
        elif current_value is not None and not isinstance(value, type(current_value)):
            value = _transformations.get(full_key, type(current_value))(value)
        setattr(current, key, value)
        logger.info(f'set config setting {full_key}={value}')
    except Exception as e:
        logger.error(f'cannot set config setting {full_key}={value}: {e}')


def _apply_env_variables():
    kwargs = {
        key[len('NOMAD_'):].lower(): value
        for key, value in os.environ.items()
        if key.startswith('NOMAD_') and key != 'NOMAD_CONFIG'}

    for key, value in kwargs.items():
        _apply(key, value, raise_error=False)


def _apply_nomad_yaml():
    config_file = os.environ.get('NOMAD_CONFIG', 'nomad.yaml')

    if not os.path.exists(config_file):
        return

    with open(config_file, 'r') as stream:
        try:
            config_data = yaml.load(stream, Loader=getattr(yaml, 'FullLoader'))
        except yaml.YAMLError as e:
            logger.error(f'cannot read nomad config: {e}')
            return

    if not config_data:
        return

    for key, value in config_data.items():
        if isinstance(value, dict):
            group_key = key
            for key, value in value.items():
                _apply(f'{group_key}_{key}', value)
        else:
            _apply(key, value)


def load_config():
    '''
    Loads the configuration from nomad.yaml and environment.
    '''
    _apply_nomad_yaml()
    _apply_env_variables()
    _check_config()


load_config()
