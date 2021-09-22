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

from bravado import requests_client as bravado_requests_client
from bravado import client as bravado_client
from urllib import parse as urllib_parse

from nomad import config as nomad_config
from nomad import utils
from nomad import client as nomad_client


def _create_client(*args, **kwargs):
    return __create_client(*args, **kwargs)


def __create_client(
        user: str = None, password: str = None, api_base_url: str = None,
        ssl_verify: bool = True, use_token: bool = True):
    ''' A factory method to create the client. '''
    # Deferred assigning of default values (instead of using Pythons default arguments),
    # because config might have been changed after import.
    if user is None:
        user = nomad_config.client.user
    if password is None:
        password = nomad_config.client.password
    if api_base_url is None:
        api_base_url = nomad_config.client.url

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
