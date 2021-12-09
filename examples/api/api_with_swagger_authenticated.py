"""
This is a brief example on how to authenticate with the public nomad@FAIRDI API.
"""

from bravado.requests_client import RequestsClient, Authenticator
from bravado.client import SwaggerClient
from urllib.parse import urlparse
from keycloak import KeycloakOpenID
from time import time

nomad_url = 'http://nomad-lab.eu/prod/rae/api'
user = 'yourusername'
password = 'yourpassword'


# an authenticator for NOMAD's keycloak user management
class KeycloakAuthenticator(Authenticator):
    def __init__(self, user, password):
        super().__init__(host=urlparse(nomad_url).netloc.split(':')[0])
        self.user = user
        self.password = password
        self.token = None
        self.__oidc = KeycloakOpenID(
            server_url='https://nomad-lab.eu/fairdi/keycloak/auth/',
            realm_name='fairdi_nomad_prod',
            client_id='nomad_public')

    def apply(self, request):
        if self.token is None:
            self.token = self.__oidc.token(username=self.user, password=self.password)
            self.token['time'] = time()
        elif self.token['expires_in'] < int(time()) - self.token['time'] + 10:
            try:
                self.token = self.__oidc.refresh_token(self.token['refresh_token'])
                self.token['time'] = time()
            except Exception:
                self.token = self.__oidc.token(username=self.user, password=self.password)
                self.token['time'] = time()

        request.headers.setdefault('Authorization', 'Bearer %s' % self.token['access_token'])

        return request


# create the bravado client
http_client = RequestsClient()
http_client.authenticator = KeycloakAuthenticator(user=user, password=password)
client = SwaggerClient.from_url('%s/swagger.json' % nomad_url, http_client=http_client)

# simple search request to print number of user entries
print(client.repo.search(owner='user').response().result.pagination.total)
