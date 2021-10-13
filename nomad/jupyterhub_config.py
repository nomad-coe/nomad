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

from oauthenticator.generic import GenericOAuthenticator
from traitlets import Unicode
import os
import requests

from nomad import config


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

    async def pre_spawn_start(self, user, spawner):
        '''
        Uses the user credentials to request all staging uploads and pass the
        respective path as volume host mounts to the spawner.
        '''
        auth_state = await user.get_auth_state()
        if not auth_state:
            self.log.warn('Authentication state is not configured!')
            return

        access_token = auth_state['access_token']
        api_headers = {'Authorization': f'Bearer {access_token}'}

        uploads_response = requests.get(
            f'{self.nomad_api_url}/v1/uploads?is_published=false&per_page=100',
            headers=api_headers)
        if not uploads_response.status_code == 200:
            self.log.error('Cannot get user uploads: %s', uploads_response.text)
            return

        volumes = {}
        for upload in uploads_response.json()['data']:
            if 'upload_files_server_path' in upload:
                upload_id = upload['upload_id']
                upload_server_path = upload['upload_files_server_path']
                volumes[f'{upload_server_path}/raw'] = f'/home/jovyan/uploads/{upload_id}'

        self.log.debug('Configure spawner with nomad volumes: %s', volumes)

        spawner.volumes = volumes


c = get_config()  # type: ignore  # pylint: disable=undefined-variable

nomad_keycloak = f'{config.keycloak.server_url.rstrip("/")}/realms/{config.keycloak.realm_name}'

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
c.NomadAuthenticator.nomad_api_url = config.api_url(ssl=False)

c.Authenticator.auto_login = True

# launch with docker
c.JupyterHub.spawner_class = 'dockerspawner.DockerSpawner'
c.DockerSpawner.image = 'jupyter/base-notebook'
c.DockerSpawner.remove = True
if config.north.docker_network:
    c.DockerSpawner.network_name = config.north.docker_network
c.DockerSpawner.hub_ip_connect = config.north.docker_host_ip
