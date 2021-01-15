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

import click
import asyncio
from concurrent import futures as concurrent_futures

from nomad import config
from .admin import admin


@admin.group(help='Run a nomad service locally (outside docker).')
def run():
    pass


@run.command(help='Run the nomad development worker.')
def worker():
    run_worker()


@run.command(help='Run the nomad development app with all apis.')
@click.option('--debug', help='Does run flask in debug.', is_flag=True)
@click.option('--with-chaos', type=int, default=0, help='Enable a certain percentage of chaos.')
def app(debug: bool, with_chaos: int):
    config.services.api_chaos = with_chaos
    run_app(debug=debug)


def run_app(**kwargs):
    config.meta.service = 'app'
    from uvicorn import Server, Config

    uv_config = Config(
        'nomad.app.main:app', host='127.0.0.1',
        port=config.services.api_port, log_level='info')
    server = Server(config=uv_config)
    server.run()


def run_worker():
    config.meta.service = 'worker'
    from nomad import processing
    processing.app.worker_main(['worker', '--loglevel=INFO', '-Q', 'celery,uploads,calcs'])


@run.command(help='Run both app and worker.')
def appworker():
    executor = concurrent_futures.ProcessPoolExecutor(2)
    loop = asyncio.get_event_loop()
    loop.run_in_executor(executor, run_app)
    loop.run_in_executor(executor, run_worker)
