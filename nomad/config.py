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
This module is used to store all configuration values. It makes use of
*namedtuples* to create key sensitive configuration objects.
"""

import os
import logging
from collections import namedtuple
import warnings

warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")


FilesConfig = namedtuple(
    'FilesConfig', ['uploads_bucket', 'raw_bucket', 'archive_bucket', 'staging_bucket', 'public_bucket'])
""" API independent configuration for the object storage. """

CeleryConfig = namedtuple('Celery', ['broker_url'])
""" Used to configure the RabbitMQ for celery. """

FSConfig = namedtuple('FSConfig', ['tmp', 'objects'])
""" Used to configure file stystem access. """

RepositoryDBConfig = namedtuple('RepositoryDBConfig', ['host', 'port', 'dbname', 'user', 'password'])
""" Used to configure access to NOMAD-coe repository db. """

ElasticConfig = namedtuple('ElasticConfig', ['host', 'port', 'index_name', 'coe_repo_calcs_index_name'])
""" Used to configure elastic search. """

MongoConfig = namedtuple('MongoConfig', ['host', 'port', 'db_name'])
""" Used to configure mongo db. """

LogstashConfig = namedtuple('LogstashConfig', ['enabled', 'host', 'tcp_port', 'level'])
""" Used to configure and enable/disable the ELK based centralized logging. """

NomadServicesConfig = namedtuple('NomadServicesConfig', ['api_host', 'api_port', 'api_base_path', 'api_secret', 'admin_password', 'upload_url', 'disable_reset'])
""" Used to configure nomad services: worker, handler, api """

MailConfig = namedtuple('MailConfig', ['host', 'port', 'user', 'password', 'from_address'])
""" Used to configure how nomad can send email """

files = FilesConfig(
    uploads_bucket='uploads',
    raw_bucket=os.environ.get('NOMAD_FILES_RAW_BUCKET', 'raw'),
    archive_bucket='archive',
    staging_bucket='staging',
    public_bucket='public'
)

rabbit_host = os.environ.get('NOMAD_RABBITMQ_HOST', 'localhost')
rabbit_port = os.environ.get('NOMAD_RABBITMQ_PORT', None)
rabbit_user = 'rabbitmq'
rabbit_password = 'rabbitmq'
rabbit_url = 'pyamqp://%s:%s@%s//' % (rabbit_user, rabbit_password, rabbit_host)


def get_loglevel_from_env(key, default_level=logging.INFO):
    plain_value = os.environ.get(key, None)
    if plain_value is None:
        return default_level
    else:
        try:
            return int(plain_value)
        except ValueError:
            return getattr(logging, plain_value, default_level)


celery = CeleryConfig(
    broker_url=rabbit_url
)

fs = FSConfig(
    tmp=os.environ.get('NOMAD_FILES_TMP_DIR', '.volumes/fs/tmp'),
    objects=os.environ.get('NOMAD_FILES_OBJECTS_DIR', '.volumes/fs/objects')
)
elastic = ElasticConfig(
    host=os.environ.get('NOMAD_ELASTIC_HOST', 'localhost'),
    port=int(os.environ.get('NOMAD_ELASTIC_PORT', 9200)),
    index_name=os.environ.get('NOMAD_ELASTIC_INDEX_NAME', 'calcs'),
    coe_repo_calcs_index_name='coe_repo_calcs'
)
repository_db = RepositoryDBConfig(
    host=os.environ.get('NOMAD_COE_REPO_DB_HOST', 'localhost'),
    port=int(os.environ.get('NOMAD_COE_REPO_DB_PORT', 5432)),
    dbname=os.environ.get('NOMAD_COE_REPO_DB_NAME', 'nomad'),
    user=os.environ.get('NOMAD_COE_REPO_DB_USER', 'postgres'),
    password=os.environ.get('NOMAD_COE_REPO_PASSWORD', 'nomad')
)
mongo = MongoConfig(
    host=os.environ.get('NOMAD_MONGO_HOST', 'localhost'),
    port=int(os.environ.get('NOMAD_MONGO_PORT', 27017)),
    db_name=os.environ.get('NOMAD_MONGO_DB_NAME', 'users')
)
logstash = LogstashConfig(
    enabled=True,
    host=os.environ.get('NOMAD_LOGSTASH_HOST', 'localhost'),
    tcp_port=int(os.environ.get('NOMAD_LOGSTASH_TCPPORT', '5000')),
    level=get_loglevel_from_env('NOMAD_LOGSTASH_LEVEL', default_level=logging.DEBUG)
)
services = NomadServicesConfig(
    api_host=os.environ.get('NOMAD_API_HOST', 'localhost'),
    api_port=int(os.environ.get('NOMAD_API_PORT', 8000)),
    api_base_path=os.environ.get('NOMAD_API_BASE_PATH', '/nomad/api'),
    api_secret=os.environ.get('NOMAD_API_SECRET', 'defaultApiSecret'),
    admin_password=os.environ.get('NOMAD_API_ADMIN_PASSWORD', 'password'),
    upload_url=os.environ.get('NOMAD_UPLOAD_URL', 'http://localhost/nomad/uploads'),
    disable_reset=os.environ.get('NOMAD_API_DISABLE_RESET', 'false') == 'true'
)
migration_source_db = RepositoryDBConfig(
    host=os.environ.get('NOMAD_MIGRATION_SOURCE_DB_HOST', 'db-repository.nomad.esc'),
    port=int(os.environ.get('NOMAD_MIGRATION_SOURCE_DB_PORT', 5432)),
    dbname=os.environ.get('NOMAD_MIGRATION_SOURCE_DB_NAME', 'nomad_prod'),
    user=os.environ.get('NOMAD_MIGRATION_SOURCE_USER', 'nomadlab'),
    password=os.environ.get('NOMAD_MIGRATION_SOURCE_PASSWORD', '*')
)
mail = MailConfig(
    host=os.environ.get('NOMAD_SMTP_HOST', ''),  # empty or None host disables email
    port=int(os.environ.get('NOMAD_SMTP_PORT', 8995)),
    user=os.environ.get('NOMAD_SMTP_USER', None),
    password=os.environ.get('NOMAD_SMTP_PASSWORD', None),
    from_address=os.environ.get('NOMAD_MAIL_FROM', 'webmaster@nomad-coe.eu')
)

console_log_level = get_loglevel_from_env('NOMAD_CONSOLE_LOGLEVEL', default_level=logging.INFO)
service = os.environ.get('NOMAD_SERVICE', 'unknown nomad service')
release = os.environ.get('NOMAD_RELEASE', 'devel')
