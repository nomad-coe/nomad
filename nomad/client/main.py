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


_default_url = 'http://%s:%d/%s' % (config.services.api_host, config.services.api_port, config.services.api_base_path.strip('/'))
_nomad_url = os.environ.get('NOMAD_URL', _default_url)
_user = os.environ.get('NOMAD_USER', 'leonard.hofstadter@nomad-fairdi.tests.de')
_pw = os.environ.get('NOMAD_PASSWORD', 'password')


def get_nomad_url():
    return _nomad_url


def create_client():
    return _create_client()


def _create_client(*args, **kwargs):
    return __create_client(*args, **kwargs)


def __create_client(user: str = _user, password: str = _pw):
    """ A factory method to create the client. """
    host = urlparse(_nomad_url).netloc.split(':')[0]
    http_client = RequestsClient()
    if user is not None:
        http_client.set_basic_auth(host, user, password)

    client = SwaggerClient.from_url(
        '%s/swagger.json' % _nomad_url,
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
                'Check connection and url.' % _nomad_url)
            sys.exit(0)
    return wrapper


@click.group()
@click.option('-n', '--url', default=_nomad_url, help='The URL where nomad is running "%s".' % _nomad_url)
@click.option('-u', '--user', default=None, help='the user name to login, default no login.')
@click.option('-w', '--password', default=_pw, help='the password use to login.')
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

    global _nomad_url
    _nomad_url = url

    global _create_client

    def _create_client(*args, **kwargs):  # pylint: disable=W0612
        if user is not None:
            logger.info('create client', user=user)
            return __create_client(user=user, password=password)
        else:
            logger.info('create anonymous client')
            return __create_client()
