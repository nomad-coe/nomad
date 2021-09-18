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

import requests
from keycloak import KeycloakOpenID
import time
from bravado import requests_client as bravado_requests_client

from nomad import config


class APIError(Exception): pass


def url(path):
    ''' Returns the full NOMAD API url for the given api path. '''
    return f'{config.client.url}/v1/{path}'


# This class is somewhat questionable, because there might be very similar functionality
# already in requests. But it is somewhere hidden in OAuth flow implementations.
class Auth(requests.auth.AuthBase):
    '''
    A request Auth class that can be used to authenticate in request callcs like this:

    .. code::

        requests.get('https://nomad-lab.eu/prod/rae/api/v1/entries', auth=Auth())

    This Auth will given username/password to acquire a token, it will keep and manage
    the token, and apply it to the headers requests.

    Arguments:
        user: Optional user name or email, default is take from ``config.client.user``
        password: Optional password, default is taken from ``config.client.password``
    '''
    def __init__(self, user: str = config.client.user, password: str = config.client.password):
        self.user = user
        self._password = password

        self._token = None

    def __call__(self, request):
        if self._token is None:
            response = requests.get(
                url('auth/token'),
                params=dict(username=self.user, password=self._password))

            if response.status_code != 200:
                response_json = response.json()
                raise APIError(
                    f'Could not acquire authentication token: '
                    f'{response_json.get("description") or response_json.get("detail") or "unknown reason"} '
                    f'({response_json.get("code", response.status_code)})')

            self._token = response.json()['access_token']

        # TODO check if token is still valid and refresh

        request.headers['Authorization'] = f'Bearer {self._token}'
        return request


class KeycloakAuthenticator(bravado_requests_client.Authenticator):
    def __init__(self, host, user, password, **kwargs):
        super().__init__(host=host)
        self.user = user
        self.password = password
        self.token = None
        self.__oidc = KeycloakOpenID(**kwargs)

    def apply(self, request=None):
        if self.token is None:
            self.token = self.__oidc.token(username=self.user, password=self.password)
            self.token['time'] = time.time()
        elif self.token['expires_in'] < int(time.time()) - self.token['time'] + 10:
            try:
                self.token = self.__oidc.refresh_token(self.token['refresh_token'])
                self.token['time'] = time.time()
            except Exception:
                self.token = self.__oidc.token(username=self.user, password=self.password)
                self.token['time'] = time.time()

        if request:
            request.headers.setdefault('Authorization', 'Bearer %s' % self.token['access_token'])
            return request
        else:
            return dict(Authorization='Bearer %s' % self.token['access_token'])
