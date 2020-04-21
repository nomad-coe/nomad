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
This module contains the Celery configuration.
"""
from celery import Celery
from nomad import config

# Celery is configured to use redis as a results backend. Although the results
# are not forwarded within the processing pipeline, celery requires the results
# backend to be configured in order to use chained tasks.
app = Celery(
    'nomad.processing',
    backend=config.redis_url(),
    broker=config.rabbitmq_url(),
)
app.conf.update(accept_content=["myjson"])
app.conf.update(task_serializer="myjson")
app.conf.update(result_serializer="myjson")
app.conf.update(worker_hijack_root_logger=False)
app.conf.update(worker_max_memory_per_child=config.celery.max_memory)
if config.celery.routing == config.CELERY_WORKER_ROUTING:
    app.conf.update(worker_direct=True)
app.conf.task_queue_max_priority = 10
