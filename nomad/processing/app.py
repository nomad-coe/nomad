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

from celery import Celery
from celery.signals import after_setup_task_logger, after_setup_logger
import logging

from nomad import config, utils
import nomad.patch  # pylint: disable=unused-import

if config.logstash.enabled:
    def initialize_logstash(logger=None, loglevel=logging.DEBUG, **kwargs):
        utils.add_logstash_handler(logger)
        return logger

    after_setup_task_logger.connect(initialize_logstash)
    after_setup_logger.connect(initialize_logstash)

app = Celery('nomad.processing', backend=config.celery.backend_url, broker=config.celery.broker_url)
app.add_defaults(dict(
    accept_content=['json', 'pickle'],
    task_serializer=config.celery.serializer,
    result_serializer=config.celery.serializer,
))
