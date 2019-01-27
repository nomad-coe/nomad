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

from nomad import config

from .main import cli


@cli.command(help='Attempts to reset the nomad.')
def reset():
    from .main import create_client
    create_client().admin.exec_reset_command().response()


@cli.group(help='Run a nomad service locally (outside docker).')
def run():
    pass


@run.command(help='Run the nomad development worker.')
def worker():
    config.service = 'nomad_worker'
    from nomad import processing
    processing.app.worker_main(['worker', '--loglevel=INFO'])


@run.command(help='Run the nomad development api.')
@click.option('--debug', help='Does run flask in debug.', is_flag=True)
def api(debug: bool):
    config.service = 'nomad_api'
    from nomad import infrastructure
    from nomad.api.__main__ import run_dev_server
    infrastructure.setup()
    run_dev_server(debug=debug, port=8000)


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
