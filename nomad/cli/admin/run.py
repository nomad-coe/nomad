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
from nomad import config
from .admin import admin
from nomad import utils


@admin.group(help='Run a nomad service locally (outside docker).')
def run():
    pass


@run.command(help='Run the jupyter hub.')
def hub():
    run_hub()


@run.command(help='Run the nomad development worker.')
def worker():
    run_worker()


@run.command(help='Run the nomad development app with all apis.')
@click.option('--debug', help='Does run app in debug.', is_flag=True)
def app(debug: bool):
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
    processing.app.worker_main(['worker', '--loglevel=INFO', '-Q', 'celery'])


def run_hub():
    from jupyterhub.app import main
    import sys
    import os
    import subprocess

    if 'JUPYTERHUB_CRYPT_KEY' not in os.environ:
        crypt_key = config.north.jupyterhub_crypt_key
        if crypt_key is None:
            crypt_key = subprocess.check_output('openssl rand -hex 32'.split(' ')).decode().strip('\n')
        os.environ['JUPYTERHUB_CRYPT_KEY'] = crypt_key

    config.meta.service = 'hub'
    sys.exit(main(argv=['-f', 'nomad/jupyterhub_config.py', '--Application.log_level=INFO']))


def task_app():
    logger = utils.get_logger('app')
    try:
        run_app()
    except Exception as error:
        logger.exception(error)


def task_worker():
    logger = utils.get_logger('worker')
    try:
        run_worker()
    except Exception as error:
        logger.exception(error)


@run.command(help='Run both app and worker.')
def appworker():
    from concurrent import futures as concurrent_futures
    import asyncio

    executor = concurrent_futures.ProcessPoolExecutor(2)
    loop = asyncio.get_event_loop()
    loop.run_in_executor(executor, task_app)
    loop.run_in_executor(executor, task_worker)
