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

import sys
import requests
import click
import logging
from bravado.requests_client import RequestsClient
from bravado.client import SwaggerClient

from nomad import config


api_base = 'http://%s:%d/%s' % (config.services.api_host, config.services.api_port, config.services.api_base_path)
user = 'leonard.hofstadter@nomad-fairdi.tests.de'
pw = 'password'


def create_client():
    return _create_client()


def _create_client(
        host: str = config.services.api_host,
        port: int = config.services.api_port,
        base_path: str = config.services.api_base_path,
        user: str = user, password: str = pw):
    """ A factory method to create the client. """

    http_client = RequestsClient()
    if user is not None:
        http_client.set_basic_auth(host, user, pw)

    client = SwaggerClient.from_url(
        'http://%s:%d%s/swagger.json' % (host, port, base_path),
        http_client=http_client)

    return client


def handle_common_errors(func):
    def wrapper(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except requests.exceptions.ConnectionError:
            click.echo(
                '\nCould not connect to nomad at %s. '
                'Check connection and host/port options.' % api_base)
            sys.exit(0)
    return wrapper


@click.group()
@click.option('-h', '--host', default=config.services.api_host, help='The host nomad runs on, default is "%s".' % config.services.api_host)
@click.option('-p', '--port', default=config.services.api_port, help='the port nomad runs with, default is %d.' % config.services.api_port)
@click.option('-u', '--user', default=None, help='the user name to login, default no login.')
@click.option('-w', '--password', default=None, help='the password use to login.')
@click.option('-v', '--verbose', help='sets log level to debug', is_flag=True)
def cli(host: str, port: int, verbose: bool, user: str, password: str):
    if verbose:
        config.console_log_level = logging.DEBUG
    else:
        config.console_log_level = logging.WARNING

    global api_base
    api_base = 'http://%s:%d/nomad/api' % (host, port)

    global create_client

    def create_client():  # pylint: disable=W0612
        if user is not None:
            return _create_client(host=host, port=port, user=user, password=password)
        else:
            return _create_client(host=host, port=port)
