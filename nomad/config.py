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

.. autoclass:: nomad.config.NomadConfig
.. autofunction:: nomad.config.apply
.. autofunction:: nomad.config.load_config
"""

import logging
import os
import os.path
import yaml
import warnings

from nomad import gitinfo


warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")


class NomadConfig(dict):
    """
    A class for configuration categories. It is a dict subclass that uses attributes as
    key/value pairs.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def __getattr__(self, name):
        if name in self:
            return self[name]
        else:
            raise AttributeError("No such attribute: " + name)

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        if name in self:
            del self[name]
        else:
            raise AttributeError("No such attribute: " + name)


CELERY_WORKER_ROUTING = 'worker'
CELERY_QUEUE_ROUTING = 'queue'

rabbitmq = NomadConfig(
    host='localhost',
    user='rabbitmq',
    password='rabbitmq'
)


def rabbitmq_url():
    return 'pyamqp://%s:%s@%s//' % (rabbitmq.user, rabbitmq.password, rabbitmq.host)


celery = NomadConfig(
    max_memory=64e6,  # 64 GB
    timeout=1800,  # 1/2 h
    acks_late=True,
    routing=CELERY_QUEUE_ROUTING,
    priorities={
        'Upload.process_upload': 5,
        'Upload.delete_upload': 9,
        'Upload.publish_upload': 10
    }
)

fs = NomadConfig(
    tmp='.volumes/fs/tmp',
    staging='.volumes/fs/staging',
    public='.volumes/fs/public',
    migration_packages='.volumes/fs/migration_packages',
    local_tmp='/tmp',
    prefix_size=2,
    working_directory=os.getcwd()
)

elastic = NomadConfig(
    host='localhost',
    port=9200,
    index_name='nomad_fairdi_calcs'
)

keycloak = NomadConfig(
    server_url='https://repository.nomad-coe.eu/fairdi/keycloak/auth/',
    realm_name='fairdi_nomad_test',
    username='admin',
    password='password',
    client_id='nomad_api_dev',
    client_secret='**********'
)

mongo = NomadConfig(
    host='localhost',
    port=27017,
    db_name='nomad_fairdi'
)

logstash = NomadConfig(
    enabled=True,
    host='localhost',
    tcp_port='5000',
    level=logging.DEBUG
)

services = NomadConfig(
    api_host='localhost',
    api_port=8000,
    api_base_path='/fairdi/nomad/latest',
    api_secret='defaultApiSecret',
    api_chaos=0,
    admin_user_id='00000000-0000-0000-0000-000000000000',
    not_processed_value='not processed',
    unavailable_value='unavailable',
    https=False,
    upload_limit=10,
    force_raw_file_decoding=False
)

tests = NomadConfig(
    default_timeout=30
)


def api_url(ssl: bool = True):
    return '%s://%s/%s/api' % (
        'https' if services.https and ssl else 'http',
        services.api_host.strip('/'),
        services.api_base_path.strip('/'))


def gui_url():
    base = api_url(True)[:-3]
    if base.endswith('/'):
        base = base[:-1]
    return '%s/gui' % base


mail = NomadConfig(
    enabled=False,
    with_login=False,
    host='',
    port=8995,
    user='',
    password='',
    from_address='webmaster@nomad-coe.eu',
    cc_address='webmaster@nomad-coe.eu'
)

normalize = NomadConfig(
    system_classification_with_clusters_threshold=50
)

client = NomadConfig(
    user='leonard.hofstadter@nomad-fairdi.tests.de',
    password='password',
    url='http://localhost:8000/fairdi/nomad/latest/api'
)

datacite = NomadConfig(
    mds_host='https://mds.datacite.org',
    enabled=False,
    prefix='10.17172',
    user='*',
    password='*'
)

version = '0.7.2'
commit = gitinfo.commit
release = 'devel'
domain = 'DFT'
service = 'unknown nomad service'
auxfile_cutoff = 100
parser_matching_size = 9128
console_log_level = logging.WARNING
max_upload_size = 32 * (1024 ** 3)
raw_file_strip_cutoff = 1000


def normalize_loglevel(value, default_level=logging.INFO):
    plain_value = value
    if plain_value is None:
        return default_level
    else:
        try:
            return int(plain_value)
        except ValueError:
            return getattr(logging, plain_value)


transformations = {
    'console_log_level': normalize_loglevel,
    'logstash_level': normalize_loglevel
}


# use std python logger, since logging is not configured while loading configuration
logger = logging.getLogger(__name__)


def apply(key, value) -> None:
    """
    Changes the config according to given key and value. The keys are interpreted as paths
    to config values with ``_`` as a separator. E.g. ``fs_staging`` leading to
    ``config.fs.staging``
    """
    path = list(reversed(key.split('_')))
    child_segment = None
    current_value = None
    child_config = globals()
    child_key = None

    try:
        while len(path) > 0:
            if child_segment is None:
                child_segment = path.pop()
            else:
                child_segment += '_' + path.pop()

            if child_segment in child_config:
                current_value = child_config[child_segment]

            if current_value is None:
                if len(path) == 0:
                    raise KeyError

                continue
            if isinstance(current_value, NomadConfig):
                child_config = current_value
                current_value = None
                child_segment = None
            else:
                if len(path) > 0:
                    raise KeyError()

                child_key = child_segment
                break

        if child_key is None or current_value is None:
            raise KeyError()
    except KeyError:
        return

    if not isinstance(value, type(current_value)):
        try:
            value = transformations.get(key, type(current_value))(value)
        except Exception as e:
            logger.error(
                'config key %s value %s has wrong type: %s' % (key, str(value), str(e)))

    child_config[child_key] = value


def load_config(config_file: str = os.environ.get('NOMAD_CONFIG', 'nomad.yaml')) -> None:
    """
    Loads the configuration from the ``config_file`` and environment.

    Arguments:
        config_file: Override the configfile, default is file stored in env variable
            NOMAD_CONFIG or ``nomad.yaml``.
    """
    # load yaml and override defaults
    if os.path.exists(config_file):
        with open(config_file, 'r') as stream:
            try:
                config_data = yaml.load(stream, Loader=getattr(yaml, 'FullLoader'))
            except yaml.YAMLError as e:
                logger.error('cannot read nomad config', exc_info=e)

        def adapt(config, new_config, child_key=None):
            for key, value in new_config.items():
                if key in config:
                    if child_key is None:
                        qualified_key = key
                    else:
                        qualified_key = '%s_%s' % (child_key, key)

                    current_value = config[key]
                    if isinstance(value, dict) and isinstance(current_value, NomadConfig):
                        adapt(current_value, value, qualified_key)
                    else:
                        if not isinstance(value, type(current_value)):
                            try:
                                value = transformations.get(qualified_key, type(current_value))(value)
                            except Exception as e:
                                logger.error(
                                    'config key %s value %s has wrong type: %s' % (key, str(value), str(e)))
                        else:
                            config[key] = value
                            logger.info('override config key %s with value %s' % (key, str(value)))
                else:
                    logger.error('config key %s does not exist' % key)

        adapt(globals(), config_data)

    # load env and override yaml and defaults
    kwargs = {
        key[len('NOMAD_'):].lower(): value
        for key, value in os.environ.items()
        if key.startswith('NOMAD_')
    }

    for key, value in kwargs.items():
        apply(key, value)


load_config()
