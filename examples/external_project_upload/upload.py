"""
This is a brief example demonstrating the public nomad@FAIRDI API for doing operations
that might be necessary to integrate external project data.
"""

from bravado.requests_client import RequestsClient, Authenticator
from bravado.client import SwaggerClient
from keycloak import KeycloakOpenID
from urllib.parse import urlparse
import time
import os.path
import sys

nomad_url = 'http://nomad-lab.eu/prod/rae/api'
user = 'youruser'
password = 'yourpassword'

upload_file = os.path.join(os.path.dirname(__file__), 'example.zip')


# an authenticator for NOMAD's keycloak user management
class KeycloakAuthenticator(Authenticator):
    def __init__(self, user, password):
        super().__init__(host=urlparse(nomad_url).netloc)
        self.user = user
        self.password = password
        self.token = None
        self.__oidc = KeycloakOpenID(
            server_url='https://nomad-lab.eu/fairdi/keycloak/auth/',
            realm_name='fairdi_nomad_test',
            client_id='nomad_public')

    def apply(self, request):
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

        request.headers.setdefault('Authorization', 'Bearer %s' % self.token['access_token'])

        return request


# create the bravado client
http_client = RequestsClient()
http_client.authenticator = KeycloakAuthenticator(user=user, password=password)
client = SwaggerClient.from_url('%s/swagger.json' % nomad_url, http_client=http_client)

# upload data
print('uploading  a file with "external_id/AcAg/vasp.xml" inside ...')
with open(upload_file, 'rb') as f:
    upload = client.uploads.upload(file=f, publish_directly=True).response().result

print('processing ...')
while upload.process_running:
    upload = client.uploads.get_upload(upload_id=upload.upload_id).response().result
    time.sleep(5)
    print('processed: %d, failures: %d' % (upload.processed_calcs, upload.failed_calcs))

# check if processing was a success
if upload.process_status != 'SUCCESS':
    print('something went wrong')
    print('errors: %s' % str(upload.errors))
    # try to delete the unsuccessful upload
    client.uploads.delete_upload(upload_id=upload.upload_id).response().result
    sys.exit(1)
