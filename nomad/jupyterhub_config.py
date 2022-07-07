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
import os
import os.path
import requests
import json

from nomad import config, infrastructure

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

    def get_handlers(self, app):
        ''' Add an additional handler for anonymous logins. '''
        handlers = super().get_handlers(app)
        return handlers + [(r'/anonymous_login', AnonymousLoginHandler)]

    async def authenticate(self, handler, data=None):
        if data == 'anonymous':
            return 'anonymous'

        try:
            return await super().authenticate(handler, data=data)
        except Exception as e:
            # The default authenticate has failed, e.g. due to missing credentials.
            # Check if there is bearer authorization and try to identify the user with
            # this. Otherwise, propagate the exception.
            authorization = handler.request.headers.get('Authorization', None)
            if not (authorization and authorization.lower().startswith('bearer ')):
                raise e

        # Use the access token to authenticate the user
        access_token = authorization[7:]
        payload = infrastructure.keycloak.decode_access_token(access_token)

        authenticated = {
            'name': payload['preferred_username'],
            'auth_state': {
                'access_token': access_token,
                'refresh_token': None,  # It is unclear how we can get the refresh
                                        # token. The only way is through a proper oauth
                                        # flow. But we can't do that if the hub is only
                                        # used via API and not the hub pages.
                                        # For now the assumption is that
                                        # - we implicitly "refresh" the access token with
                                        #   each API call that causes this to run
                                        # - we assume that jhub we automatically go through
                                        #   a oauth flow, once the access token is expired
                                        #   and it does not find a stored refresh token.
                'oauth_user': payload,
                'scope': []}}
        # This will save the users authstate
        await handler.auth_to_user(authenticated)

        return authenticated

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
        spawner.environment['SUBFOLDER'] = f'{config.services.api_base_path.rstrip("/")}/north/user/{user.name}/'
        spawner.environment['JUPYTERHUB_CLIENT_API_URL'] = f'{config.north_url()}/hub/api'

        if user.name == 'anonymous':
            return

        try:
            auth_state = await user.get_auth_state()
            if not auth_state:
                self.log.warn('Authentication state is not configured!')
                return
            access_token = auth_state['access_token']
            api_headers = {'Authorization': f'Bearer {access_token}'}

            uploads_response = requests.get(
                f'{config.api_url().rstrip("/")}/v1/uploads?is_published=false&per_page=100',
                headers=api_headers)
        except Exception as e:
            self.log.error('Cannot access Nomad API: %s', e)
            return

        if uploads_response.status_code != 200:
            self.log.error('Cannot get user uploads: %s', uploads_response.text)
            return

        volumes = {}

        def add_volume(host_path, mount_path):
            host_path = os.path.abspath(host_path)
            volumes[host_path] = mount_path

        add_volume(os.path.join(config.north.users_fs, user.name), f'/prefix/work')
        add_volume(os.path.join(config.north.shared_fs), f'/prefix/shared')

        for upload in uploads_response.json()['data']:
            if 'upload_files_server_path' in upload:
                upload_id = upload['upload_id']
                upload_server_path = upload['upload_files_server_path']
                add_volume(f'{upload_server_path}/raw', f'/prefix/uploads/{upload_id}')

        self.log.debug('Configure spawner with nomad volumes: %s', volumes)

        spawner.volumes = volumes
        spawner.nomad_username = user.name


c = get_config()  # type: ignore  # pylint: disable=undefined-variable

# Allow named single-user servers per user (Default: False)
c.JupyterHub.allow_named_servers = True

# If named servers are enabled, default name of server to spawn or open, e.g. by
# user-redirect. (Default: '')
c.JupyterHub.default_server_name = 'jupyter'

# TODO: This is temporary.  Place everything behind a single origin aka nginx proxy
c.JupyterHub.tornado_settings = {
    'headers': {
        'Access-Control-Allow-Origin': '*',
        'Access-Control-Allow-Headers': '*',
        'Access-Control-Allow-Methods': '*'
    },
}

# The public facing URL of the whole JupyterHub application. This is the address on which
# the proxy will bind. (Default: 'http://:8000')
c.JupyterHub.bind_url = f'http://:9000/{config.services.api_base_path.strip("/")}/north'

nomad_public_keycloak = f'{config.keycloak.public_server_url.rstrip("/")}/realms/{config.keycloak.realm_name}'
nomad_keycloak = f'{config.keycloak.server_url.rstrip("/")}/realms/{config.keycloak.realm_name}'

c.JupyterHub.authenticator_class = NomadAuthenticator

c.GenericOAuthenticator.login_service = 'keycloak'
c.GenericOAuthenticator.client_id = 'nomad_public'
c.GenericOAuthenticator.authorize_url = f'{nomad_public_keycloak}/protocol/openid-connect/auth'
c.GenericOAuthenticator.token_url = f'{nomad_keycloak}/protocol/openid-connect/token'
c.GenericOAuthenticator.userdata_url = f'{nomad_keycloak}/protocol/openid-connect/userinfo'
c.GenericOAuthenticator.userdata_params = {'state': 'state'}
c.GenericOAuthenticator.username_key = 'preferred_username'
# c.GenericOAuthenticator.userdata_method = 'GET'
c.GenericOAuthenticator.scope = ['openid', 'profile']

c.Authenticator.auto_login = True
c.Authenticator.enable_auth_state = True


class NomadDockerSpawner(DockerSpawner):

    async def start(self, image=None, extra_create_kwargs=None, extra_host_config=None):
        self.log.debug(f'Configuring spawner for container {self.container_name}')
        tool = self.container_name.split('-')[-1]
        if tool in tools:
            tools[tool](self)
            self.log.debug(f'Configured spawner for {tool}')
        else:
            self.log.error(f'{tool} is not a tool, raise an error')
            raise NotImplementedError('You cannot launch non tool containers.')

        return await super().start(image, extra_create_kwargs, extra_host_config)

    def _docker(self, method, *args, **kwargs):
        if config.north.windows:
            tries = 0
            max_tries = 1
            if method == 'port':
                max_tries = 3
            while tries < max_tries:
                result = super()._docker(method, *args, **kwargs)
                if result is not None:
                    break
                import time
                time.sleep(3)
                tries += 1

            return result
        else:
            return super()._docker(method, *args, **kwargs)


# launch with docker
c.JupyterHub.spawner_class = NomadDockerSpawner
c.DockerSpawner.image = 'jupyter/base-notebook'
c.DockerSpawner.remove = True

# Prefix for container names. See name_template for full container name for a particular
# user's server. (Default: 'jupyter')
c.DockerSpawner.prefix = 'noamd_oasis_north'

if config.north.docker_network:
    c.DockerSpawner.network_name = config.north.docker_network

if config.north.hub_ip:
    c.JupyterHub.hub_ip = config.north.hub_ip

if config.north.hub_connect_ip:
    c.JupyterHub.hub_connect_ip = config.north.hub_connect_ip

if config.north.hub_connect_url:
    c.DockerSpawner.hub_connect_url = config.north.hub_connect_url


def create_configure_from_tool_json(tool_json):
    def configure(spawner: DockerSpawner):
        spawner.image = tool_json['image']
        if 'cmd' in tools_json:
            spawner.cmd = tools_json['cmd']

        for key, value in spawner.volumes.items():
            if value.startswith('/prefix') and 'mount_path' in tool_json:
                spawner.volumes[key] = value.replace('/prefix', tool_json['mount_path'])

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
