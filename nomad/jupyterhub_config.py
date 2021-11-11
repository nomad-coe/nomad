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

'''
Configuration for JuypterHUB. We use JuypterHUB to run remote tools, like jupyter, on
NOMAD servers.
'''

from dockerspawner.dockerspawner import DockerSpawner
from jupyterhub.handlers.base import BaseHandler
from oauthenticator.generic import GenericOAuthenticator
from traitlets import Unicode
import os
import os.path
import requests
import pathlib
import json

from nomad import config

# TODO The AnonymousLogin kinda works, but it won't logout or remove old containers,
# if the user is still logged in non anonymously.


class AnonymousLoginHandler(BaseHandler):
    async def get(self):
        # TODO somehow set a fake user to avoid login
        # - read cookie and use cookie, or create, set, and use cookie
        await self.login_user(data='anonymous')
        self.redirect(self.get_next_url())


class NomadAuthenticator(GenericOAuthenticator):
    '''
    A custom OAuth authenticator. It can be used with NOMAD's keycloak. It imeplement
    a `pre_spawn_start` hook that passes user upload mounts to the spawner.
    '''

    nomad_api_url = Unicode(
        os.environ.get('NOMAD_CLIENT_URL', 'https://nomad-lab.eu/prod/rae/api'),
        config=True,
        help='Nomad Api Url. This is used to get all staging upload ids for the loggedin user.',
    )

    def get_handlers(self, app):
        ''' Add an additional handler for anonymous logins. '''
        handlers = super().get_handlers(app)
        return handlers + [(r'/anonymous_login', AnonymousLoginHandler)]

    async def authenticate(self, handler, data=None):
        if data == 'anonymous':
            return 'anonymous'

        return await super().authenticate(handler, data=data)

    async def pre_spawn_start(self, user, spawner):
        '''
        Uses the user credentials to request all staging uploads and pass the
        respective path as volume host mounts to the spawner.
        '''

        # This is guacamole specific
        # linuxserver/webtop guacamole-lite based guacamole client use SUBFOLDER to
        # confiugure base path
        if not spawner.environment:
            spawner.environment = {}
        spawner.environment['SUBFOLDER'] = f'{config.north.hub_base_path}/user/{user.name}/'

        if user.name == 'anonymous':
            return

        auth_state = await user.get_auth_state()
        if not auth_state:
            self.log.warn('Authentication state is not configured!')
            return

        access_token = auth_state['access_token']
        api_headers = {'Authorization': f'Bearer {access_token}'}

        try:
            uploads_response = requests.get(
                f'{self.nomad_api_url}/v1/uploads?is_published=false&per_page=100',
                headers=api_headers)
        except Exception as e:
            self.log.error('Cannot access Nomad API: %s', e)
            return

        if uploads_response.status_code != 200:
            self.log.error('Cannot get user uploads: %s', uploads_response.text)
            return

        user_dir = os.path.abspath(f'{config.north.users_fs}/{user.name}')
        if not os.path.exists(user_dir):
            pathlib.Path(user_dir).mkdir(parents=True, exist_ok=True)

        shared_dir = os.path.abspath(config.north.shared_fs)
        if not os.path.exists(shared_dir):
            pathlib.Path(user_dir).mkdir(parents=True, exist_ok=True)

        volumes = {
            user_dir: f'/home/jovyan/work',
            shared_dir: f'/home/jovyan/shared'
        }

        for upload in uploads_response.json()['data']:
            if 'upload_files_server_path' in upload:
                upload_id = upload['upload_id']
                upload_server_path = upload['upload_files_server_path']
                volumes[f'{upload_server_path}/raw'] = f'/home/jovyan/uploads/{upload_id}'

        self.log.debug('Configure spawner with nomad volumes: %s', volumes)

        spawner.volumes = volumes


c = get_config()  # type: ignore  # pylint: disable=undefined-variable

nomad_keycloak = f'{config.keycloak.server_url.rstrip("/")}/realms/{config.keycloak.realm_name}'

c.JupyterHub.allow_named_servers = True

c.JupyterHub.port = config.north.hub_port
c.JupyterHub.base_url = config.north.hub_base_path
c.JupyterHub.authenticator_class = NomadAuthenticator
c.Authenticator.enable_auth_state = True
c.GenericOAuthenticator.login_service = 'keycloak'
c.GenericOAuthenticator.client_id = 'nomad_public'
c.GenericOAuthenticator.authorize_url = f'{nomad_keycloak}/protocol/openid-connect/auth'
c.GenericOAuthenticator.token_url = f'{nomad_keycloak}/protocol/openid-connect/token'
c.GenericOAuthenticator.userdata_url = f'{nomad_keycloak}/protocol/openid-connect/userinfo'
c.GenericOAuthenticator.userdata_params = {'state': 'state'}
c.GenericOAuthenticator.username_key = 'preferred_username'
c.GenericOAuthenticator.userdata_method = 'GET'
c.GenericOAuthenticator.scope = ['openid', 'profile']
c.NomadAuthenticator.nomad_api_url = config.north.nomad_api_url

c.Authenticator.auto_login = True


class NomadDockerSpawner(DockerSpawner):

    async def start(self, image=None, extra_create_kwargs=None, extra_host_config=None):
        tool = self.container_name.split('-')[-1]
        if tool in tools:
            tools[tool](self)
        else:
            configure_default(self)

        return await super().start(image, extra_create_kwargs, extra_host_config)


# launch with docker
c.JupyterHub.spawner_class = NomadDockerSpawner
# c.JupyterHub.spawner_class = 'dockerspawner.DockerSpawner'
c.DockerSpawner.image = 'jupyter/base-notebook'
c.DockerSpawner.remove = True
if config.north.docker_network:
    c.DockerSpawner.network_name = config.north.docker_network
c.JupyterHub.hub_ip_connect = config.north.hub_ip_connect
if config.north.hub_ip:
    c.JupyterHub.hub_ip = config.north.hub_ip


def configure_default(spawner: DockerSpawner):
    spawner.image = 'jupyter/base-notebook'


def create_configure_from_tool_json(tool_json):
    def configure(spawner: DockerSpawner):
        spawner.image = tool_json['image']

    return configure


tools_json_path = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    '../dependencies/nomad-remote-tools-hub/tools.json')

with open(tools_json_path, 'rt') as f:
    tools_json = json.load(f)

tools = {
    key: create_configure_from_tool_json(value)
    for key, value in tools_json.items()
}
