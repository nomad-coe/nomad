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

from collections import namedtuple

S3Config = namedtuple('S3', ['uploads_bucket', 'repository_bucket', 'archive_bucket'])
""" API independent configuration for the object storage. """

RabitMQConfig = namedtuple('RabbitMQ', ['host', 'port', 'user', 'password'])
""" Used to configure the RabbitMQ used by celery as a task backend. """

MinioConfig = namedtuple('Minio', ['host', 'port', 'accesskey', 'secret'])
""" Used to configure the minio object storage API. """

FSConfig = namedtuple('FSConfig', ['tmp'])
""" Used to configure file stystem access. """

LogstashConfig = namedtuple('LogstashConfig', ['enabled', 'host', 'tcp_port'])
""" Used to configure and enable/disable the ELK based centralized logging. """

s3 = S3Config(
    uploads_bucket='uploads',
    repository_bucket='repository',
    archive_bucket='archive'
)
rabbitmq = RabitMQConfig(
    host='localhost',
    port=None,
    user='rabbitmq',
    password='rabbitmq'
)
minio = MinioConfig(
    host='localhost',
    port=9007,
    accesskey='AKIAIOSFODNN7EXAMPLE',
    secret='wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY'
)
fs = FSConfig(
    tmp='.volumes/fs'
)
logstash = LogstashConfig(
    enabled=False,
    host='localhost',
    tcp_port=5000
)
