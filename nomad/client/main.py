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

import os
import sys
import requests
import click
import logging
from bravado.requests_client import RequestsClient
from bravado.client import SwaggerClient
from urllib.parse import urlparse

from nomad import config, utils, infrastructure


def create_client():
    return _create_client()


def _create_client(*args, **kwargs):
    return __create_client(*args, **kwargs)


def __create_client(user: str = config.client.user, password: str = config.client.password):
    """ A factory method to create the client. """
    host = urlparse(config.client.url).netloc.split(':')[0]
    http_client = RequestsClient()
    if user is not None:
        http_client.set_basic_auth(host, user, password)

    client = SwaggerClient.from_url(
        '%s/swagger.json' % config.client.url,
        http_client=http_client)

    utils.get_logger(__name__).info('created bravado client', user=user)

    return client


def handle_common_errors(func):
    def wrapper(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except requests.exceptions.ConnectionError:
            click.echo(
                '\nCould not connect to nomad at %s. '
                'Check connection and url.' % config.client.url)
            sys.exit(0)
    return wrapper


@click.group()
@click.option('-n', '--url', default=config.client.url, help='The URL where nomad is running, default is "%s".' % config.client.url)
@click.option('-u', '--user', default=None, help='the user name to login, default is "%s" login.' % config.client.user)
@click.option('-w', '--password', default=config.client.password, help='the password used to login.')
@click.option('-v', '--verbose', help='sets log level to info', is_flag=True)
@click.option('--debug', help='sets log level to debug', is_flag=True)
def cli(url: str, verbose: bool, debug: bool, user: str, password: str):
    if debug:
        config.console_log_level = logging.DEBUG
    elif verbose:
        config.console_log_level = logging.INFO
    else:
        config.console_log_level = logging.WARNING

    config.service = os.environ.get('NOMAD_SERVICE', 'client')
    infrastructure.setup_logging()

    logger = utils.get_logger(__name__)

    logger.info('Used nomad is %s' % url)
    logger.info('Used user is %s' % user)

    config.client.url = url

    global _create_client

    def _create_client(*args, **kwargs):  # pylint: disable=W0612
        if user is not None:
            logger.info('create client', user=user)
            return __create_client(user=user, password=password)
        else:
            logger.info('create anonymous client')
            return __create_client()
