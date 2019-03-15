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

import os.path
import os
import sys
import click
import asyncio
from concurrent.futures import ProcessPoolExecutor
from celery.task.control import revoke

from nomad import config, infrastructure, processing, utils

from .main import cli


@cli.group(help='Processing related functions')
def proc():
    pass

@proc.command(help='List processing tasks')
def ls():
    infrastructure.setup_logging()
    infrastructure.setup_mongo()

    def ls(query):
        for proc in query:
            print(proc)


    ls(processing.Calc.objects(process_status=processing.PROCESS_RUNNING))
    ls(processing.Upload.objects(process_status=processing.PROCESS_RUNNING))

@proc.command(help='Stop all running processing')
@click.option('--calcs', is_flag=True, help='Only stop calculation processing')
@click.option('--kill', is_flag=True, help='Use the kill signal and force task failure')
def stop_all(calcs: bool, kill: bool):
    infrastructure.setup_logging()
    infrastructure.setup_mongo()

    logger = utils.get_logger(__name__)

    def stop_all(query):
        for proc in query:
            logger_kwargs = dict(upload_id=proc.upload_id)
            if isinstance(proc, processing.Calc):
                logger_kwargs.update(calc_id=proc.calc_id)

            logger.info(
                'send terminate celery task', celery_task_id=proc.celery_task_id,
                kill=kill, **logger_kwargs)

            kwargs = {}
            if kill:
                kwargs.update(signal='SIGKILL')
            revoke(proc.celery_task_id, terminate=True, **kwargs)
            if kill:
                logger.info(
                    'fail proc', celery_task_id=proc.celery_task_id, kill=kill,
                    **logger_kwargs)

                proc.fail('process terminate via nomad cli')

    stop_all(processing.Calc.objects(process_status=processing.PROCESS_RUNNING))
    if not calcs:
        stop_all(processing.Upload.objects(process_status=processing.PROCESS_RUNNING))


@cli.command(help='Attempts to reset the nomad.')
def reset():
    from .main import create_client
    create_client().admin.exec_reset_command().response()


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


@cli.command(help='Runs tests and linting. Useful before commit code.')
@click.option('--skip-tests', help='Do not test, just do code checks.', is_flag=True)
def qa(skip_tests: bool):
    os.chdir(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
    ret_code = 0
    if not skip_tests:
        click.echo('Run tests ...')
        ret_code += os.system('python -m pytest -svx tests')
    click.echo('Run code style checks ...')
    ret_code += os.system('python -m pycodestyle --ignore=E501,E701 nomad tests')
    click.echo('Run linter ...')
    ret_code += os.system('python -m pylint --load-plugins=pylint_mongoengine nomad tests')
    click.echo('Run static type checks ...')
    ret_code += os.system('python -m mypy --ignore-missing-imports --follow-imports=silent --no-strict-optional nomad tests')

    sys.exit(ret_code)
