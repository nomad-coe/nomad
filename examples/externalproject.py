"""
This is a brief example demonstrating the public nomad@FAIRDI API for doing operations
that might be necessary to integrate external project data.
"""

from bravado.requests_client import RequestsClient
from bravado.client import SwaggerClient
from bravado.exception import HTTPNotFound
from urllib.parse import urlparse
import time
import os.path
import sys

nomad_url = 'http://enc-staging-nomad.esc.rzg.mpg.de/fairdi/nomad/v0.3.0/api'
user = 'leonard.hofstadter@nomad-fairdi.tests.de'
password = 'password'

upload_file = 'external_project_example.zip'

# create the bravado client
host = urlparse(nomad_url).netloc.split(':')[0]
http_client = RequestsClient()
http_client.set_basic_auth(host, user, password)
client = SwaggerClient.from_url('%s/swagger.json' % nomad_url, http_client=http_client)

# upload data
with open(upload_file, 'rb') as f:
    upload = client.uploads.upload(file=f).response().result
while upload.tasks_running:
    upload = client.uploads.get_upload(upload_id=upload.upload_id).response().result
    time.sleep(5)
    print('processed: %d, failures: %d' % (upload.processed_calcs, upload.failed_calcs))

#   check if processing was a success
if upload.tasks_status != 'SUCCESS':
    print('something went wrong')
    print('errors: %s' % str(upload.errors))
    # delete the unsuccessful upload
    client.uploads.delete_upload(upload_id=upload.upload_id).response().result
    sys.exit(1)

# publish data
client.uploads.exec_upload_operation(upload_id=upload.upload_id, payload={
    'operation': 'publish',
    'metadata': {
        'comment': 'Data from a cool external project',
        'references': ['http://external.project.eu'],
        # 'coauthors': ['sheldon.cooper@ucla.edu'],  this does not yet work with emails
        # 'external_id': 'external_id'  this does also not work, but we could implement something like this
    }
})
while upload.process_running:
    try:
        upload = client.uploads.get_upload(upload_id=upload.upload_id).response().result
        time.sleep(1)
    except HTTPNotFound:
        # upload gets deleted from the upload staging area once published
        break

# search for data
result = client.repo.search(paths='external_id').response().result
if result.pagination.total == 0:
    print('not found')
    sys.exit(1)
elif result.pagination.total > 1:
    print('my ids are not specific enough, bummer ... or did I uploaded stuff multiple times?')
# The results key holds an array with the current page data
calc = result.results[0]

# download data
#   via api
client.raw.get(upload_id=calc['upload_id'], path=calc['mainfile']).response()
#   via download
#   just the 'mainfile'
url = '%s/raw/%s/%s' % (nomad_url, calc['upload_id'], calc['mainfile'])
#   all files
url = '%s/raw/%s/%s/*' % (nomad_url, calc['upload_id'], os.path.dirname(calc['mainfile']))
