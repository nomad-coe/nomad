"""
This is a brief example demonstrating the public nomad@FAIRDI API for doing operations
that might be necessary to integrate external project data.
"""

import requests
import os.path
import time

from nomad.config import config
from nomad.client import Auth

nomad_url = config.client.url
user = 'youruser'
password = 'yourpassword'

upload_file = os.path.join(os.path.dirname(__file__), 'example.zip')

# create an auth object
auth = Auth(user=user, password=password)

# upload data
print('uploading  a file with "external_id/AcAg/vasp.xml" inside ...')
with open(upload_file, 'rb') as f:
    response = requests.post(
        f'{nomad_url}/v1/uploads', auth=auth, data=f, params=dict(publish_directly=True)
    )
    assert response.status_code == 200
    upload = response.json()['data']

print('processing ...')
while upload['process_running']:
    response = requests.get(f'{nomad_url}/v1/uploads/{upload["upload_id"]}', auth=auth)
    assert response.status_code == 200
    upload = response.json()['data']
    time.sleep(5)
    print(
        'processed: %d, failures: %d'
        % (upload['processed_entries_count'], upload['failed_entries_count'])
    )

# check if processing was a success
if upload['process_status'] != 'SUCCESS':
    print('something went wrong')
    print('errors: %s' % str(upload.errors))
    # try to delete the unsuccessful upload
    requests.delete(f'{nomad_url}/v1/uploads/{upload["upload_id"]}', auth=auth)
