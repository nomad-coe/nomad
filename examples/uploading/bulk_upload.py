'''
The scenario is that you have a lot of <name>-*.zip files for upload. This script
will add a nomad.json to the uploads, upload them 1 at a time, watch the processing, repeat.
'''

from typing import Dict, Any
from bravado.requests_client import RequestsClient, Authenticator
from bravado.client import SwaggerClient
from keycloak import KeycloakOpenID
from urllib.parse import urlparse
import time
import os.path
import sys
import zipfile

from nomad.client import KeycloakAuthenticator
from nomad import config

nomad_url = 'http://nomad-lab.eu/prod/rae/api'
user = 'youruser'
password = 'yourpassword'
uploader_id = None


# create the bravado client
http_client = RequestsClient()
http_client.authenticator = KeycloakAuthenticator(
    host=urlparse(nomad_url).netloc,
    user=user,
    password=password,
    server_url=config.keycloak.server_url,
    realm_name=config.keycloak.realm_name,
    client_id=config.keycloak.client_id)
client = SwaggerClient.from_url('%s/swagger.json' % nomad_url, http_client=http_client)


def upload(
        path: str, local_path: bool = False, metadata_path: str = None,
        publish_directly: bool = False, uploader_id: str = None):
    '''
    Arguments:
        path: The file path to the upload file.
        local: If this file should be uploaded with &local_path
        metadata_path: Optional nomad.(yaml|json) metadata file
        publish_directly: If the upload should be published directly
    '''

    assert os.path.isfile(path), f'The {path} is not a file'

    # add metadata
    if metadata_path is not None:
        assert os.path.isfile(metadata_path), f'The {metadata_path} is not a file'
        assert os.path.basename(metadata_path).endswith('.json'), f'The {metadata_path} is not a nomad metadata file'
        assert path.endswith('.zip'), 'Adding nomad metadata is only supported for .zip files'

        with zipfile.ZipFile(path, 'a') as zip:
            zip.write(metadata_path, 'nomad.json')

    # upload
    print(f'uploading {path}')
    kwargs: Dict[str, Any] = {}
    if publish_directly:
        kwargs['publish_directly'] = True

    if uploader_id is not None:
        kwargs['uploader_id'] = uploader_id

    if local_path:
        upload = client.uploads.upload(local_path=path, **kwargs).response().result
    else:
        with open(path, 'rb') as f:
            upload = client.uploads.upload(file=f, **kwargs).response().result

    print(f'processing {path}')
    while upload.tasks_running:
        upload = client.uploads.get_upload(upload_id=upload.upload_id).response().result
        time.sleep(3)
        print('processed: %d, failures: %d' % (upload.processed_calcs, upload.failed_calcs))

    # check if processing was a success
    if upload.tasks_status != 'SUCCESS':
        print('something went wrong')
        print('errors: %s' % str(upload.errors))
        # try to delete the unsuccessful upload
        client.uploads.delete_upload(upload_id=upload.upload_id).response().result

        return False

    return True


if __name__ == '__main__':
    metadata_path = None
    if sys.argv[1].endswith('json'):
        metadata_path = sys.argv[1]
        paths = sys.argv[2:]
    else:
        paths = sys.argv[1:]

    for path in paths:
        upload(
            path, metadata_path=metadata_path, local_path=True, publish_directly=True,
            uploader_id=uploader_id)
