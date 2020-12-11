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

import sys
import requests
import click
from bravado import requests_client as bravado_requests_client
from bravado import client as bravado_client
from urllib import parse as urllib_parse

from nomad import config as nomad_config
from nomad import utils
from nomad import client as nomad_client

from nomad.cli.cli import cli


def create_client():
    return _create_client()


def _create_client(*args, **kwargs):
    return __create_client(*args, **kwargs)


def __create_client(
        user: str = nomad_config.client.user,
        password: str = nomad_config.client.password,
        api_base_url: str = nomad_config.client.url,
        ssl_verify: bool = True, use_token: bool = True):
    ''' A factory method to create the client. '''
    if not ssl_verify:
        import warnings
        warnings.filterwarnings("ignore")
    http_client = bravado_requests_client.RequestsClient(ssl_verify=ssl_verify)

    client = bravado_client.SwaggerClient.from_url(
        '%s/swagger.json' % api_base_url, http_client=http_client)

    utils.get_logger(__name__).info('created bravado client', user=user)

    if user is not None:
        host = urllib_parse.urlparse(api_base_url).netloc
        if use_token:
            http_client.authenticator = nomad_client.KeycloakAuthenticator(
                host=host,
                user=user,
                password=password,
                server_url=nomad_config.keycloak.server_url,
                realm_name=nomad_config.keycloak.realm_name,
                client_id=nomad_config.keycloak.client_id)
        else:
            http_client.set_basic_auth(
                host=host,
                username=user,
                password=password)

        utils.get_logger(__name__).info('set bravado client authentication', user=user)

    return client


def handle_common_errors(func):
    def wrapper(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except requests.exceptions.ConnectionError:
            click.echo(
                '\nCould not connect to nomad at %s. '
                'Check connection and url.' % nomad_config.client.url)
            sys.exit(0)
    return wrapper


@cli.group(help='Commands that use the nomad API to do useful things')
@click.option('-n', '--url', default=nomad_config.client.url, help='The URL where nomad is running, default is "%s".' % nomad_config.client.url)
@click.option('-u', '--user', default=None, help='the user name to login, default is "%s" login.' % nomad_config.client.user)
@click.option('-w', '--password', default=nomad_config.client.password, help='the password used to login.')
@click.option('--no-ssl-verify', help='disables SSL verificaton when talking to nomad.', is_flag=True)
@click.option('--no-token', is_flag=True, help='replaces token with basic auth, e.g. to work with v0.6.x or older API versions')
@click.pass_context
def client(ctx, url: str, user: str, password: str, no_ssl_verify: bool, no_token: bool):
    logger = utils.get_logger(__name__)

    logger.info('Used nomad is %s' % url)
    logger.info('Used user is %s' % user)

    nomad_config.client.url = url

    ctx.obj.user = user

    global _create_client

    def _create_client(*args, **kwargs):  # pylint: disable=W0612
        if user is not None:
            logger.info('create client', user=user)
            return __create_client(
                user=user, password=password, ssl_verify=not no_ssl_verify,
                use_token=not no_token)
        else:
            logger.info('create anonymous client')
            return __create_client(ssl_verify=not no_ssl_verify)
