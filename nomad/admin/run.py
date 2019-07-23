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

import click
import asyncio
from concurrent.futures import ProcessPoolExecutor

from nomad import config
from .cli import cli


@cli.group(help='Run a nomad service locally (outside docker).')
def run():
    pass


@run.command(help='Run the nomad development worker.')
def worker():
    run_worker()


@run.command(help='Run the nomad development api.')
@click.option('--debug', help='Does run flask in debug.', is_flag=True)
def api(debug: bool):
    run_api(debug=debug)


def run_api(**kwargs):
    config.service = 'api'
    from nomad import infrastructure
    from nomad.api.__main__ import run_dev_server
    infrastructure.setup()
    run_dev_server(port=8000, **kwargs)


def run_worker():
    config.service = 'worker'
    from nomad import processing
    processing.app.worker_main(['worker', '--loglevel=INFO', '-Q', 'celery,uploads,calcs'])


@run.command(help='Run both api and worker.')
def apiworker():
    executor = ProcessPoolExecutor(2)
    loop = asyncio.get_event_loop()
    loop.run_in_executor(executor, run_api)
    loop.run_in_executor(executor, run_worker)
