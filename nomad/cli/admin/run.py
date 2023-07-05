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

from nomad import config, utils

from .admin import admin


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
@click.option('--with-gui', help='The app will configure the gui for production and service it.', is_flag=True)
@click.option('--host', type=str, help='Passed to uvicorn host parameter.')
@click.option('--port', type=int, help='Passed to uvicorn host parameter.')
@click.option('--log-config', type=str, help='Passed to uvicorn log-config parameter.')
@click.option('--workers', type=int, help='Passed to uvicorn workers parameter.')
def app(with_gui: bool, **kwargs):
    run_app(with_gui=with_gui, **kwargs)


@run.command(help='Run server that collects and submits logs to the central Nomad instance.')
def logtransfer():
    config.meta.service = 'logtransfer'

    from nomad.logtransfer import start_logtransfer_service
    start_logtransfer_service()


def run_app(with_gui: bool = False, **kwargs):
    config.meta.service = 'app'

    if with_gui:
        import os.path
        import glob
        import shutil

        gui_folder = os.path.abspath(
            os.path.join(os.path.dirname(__file__), '../../app/static/gui'))
        run_gui_folder = os.path.join(gui_folder, '../.gui_configured')

        # copy
        shutil.rmtree(run_gui_folder, ignore_errors=True)
        shutil.copytree(gui_folder, run_gui_folder)

        # replace base path in all GUI files
        source_file_globs = [
            '**/*.json',
            '**/*.html',
            '**/*.js',
            '**/*.js.map',
            '**/*.css']
        for source_file_glob in source_file_globs:
            source_files = glob.glob(os.path.join(run_gui_folder, source_file_glob), recursive=True)
            for source_file in source_files:
                with open(source_file, 'rt') as f:
                    file_data = f.read()
                file_data = file_data.replace('/fairdi/nomad/latest', config.services.api_base_path)
                with open(source_file, 'wt') as f:
                    f.write(file_data)

        # App and gui are served from the same server, same port. Replace the base urls with
        # relative paths
        config.ui.app_base = f'{config.services.api_base_path.rstrip("/")}'
        config.ui.north_base = f'{config.services.api_base_path.rstrip("/")}/north'

    from uvicorn import Server, Config
    from nomad.utils import get_logger

    uv_config = Config(
        'nomad.app.main:app',
        log_level='info',
        **{k: v for k, v in kwargs.items() if v is not None})

    server = Server(config=uv_config)
    get_logger(__name__).info('created uvicorn server', data=uv_config.__dict__)
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
