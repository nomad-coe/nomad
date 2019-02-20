"""
This is a brief example demonstrating the public nomad@FAIRDI API for doing operations
that might be necessary to integrate external project data.
"""

# This does not assume many specific python packages. Only the bravado
# library that allows to use swagger-based ReST APIs is required.
# It can be install via `pip install bravado`
from bravado.requests_client import RequestsClient
from bravado.client import SwaggerClient
from bravado.exception import HTTPNotFound
from urllib.parse import urlparse
import time
import os.path
import sys

nomad_url = 'http://enc-staging-nomad.esc.rzg.mpg/fairdi/nomad/v0.3.0/api'
user = 'leonard.hofstadter@nomad-fairdi.tests.de'
password = 'password'

# lets assume we have a test file from our external project
# with (among others) `/external_id/BrSiTi/vasp.xml.gz`
upload_file = 'externa_project_example.tgz'

# create the bravado client
host = urlparse(nomad_url).netloc.split(':')[0]
http_client = RequestsClient()
http_client.set_basic_auth(host, user, password)
client = SwaggerClient.from_url('%s/swagger.json' % nomad_url, http_client=http_client)

# upload data
#   create the upload by uploaded a .zip file
upload = client.uploads.upload(file=upload_file).response().result
#   constantly polling the upload to get updates on the processing
while upload.processing_running:
    upload = client.uploads.get_upload(upload_id=upload.upload_id).response().result
    time.sleep(5)
    print('processed: %d, failures: %d' % (upload.processed_calcs, upload.failed_calcs))

#   check if processing was a success
if upload.tasks_status != 'SUCCESS':
    print('something went wrong')
    print('errors: %s' % str(upload.errors))
    # delete the unsuccessful upload
    client.uploads.delete_upload(upload_id=upload.upload_id)
    sys.exit(1)

# publish data
#   In the upload staging area the data is only visible to you. It has to be published
#   to get into the public nomad.
#   Therefore, the later search and download steps will work without publishing, but only if the
#   client is authenticated with your user account. You should do that for testing stuff out,
#   because there is no user-based deleting of published data.
#   The publish step also allows you to provide additional metadata, see below.
client.uploads.exec_upload_operation(upload_id=upload.upload_id, payload={
    'operation': 'publish',
    'metadata': {
        'comment': 'Data from materials project',
        'references': ['http://materials-project.gov'],
        # 'coauthors': ['person@lbl.gov', '...'],  this does not yet work with emails
        # 'external_id': 'a/mp/id'  this does also not work, but we could implement something like this
    }
})
while upload.processing_running:
    try:
        upload = client.uploads.get_upload(upload_id=upload.upload_id).response().result
        time.sleep(1)
    except HTTPNotFound:
        # upload gets deleted from the upload staging area once published
        break

# search for data
#   The `paths` searchkey is a text search (works like google), where tokens are separated
#   by '/'. You can basically search for (multiple) parts within a file path. Note
#   that the following will also match `/prefix/tag1/something/else/tag2/vasp.xml`. It will
#   not match `/tag2/tag3/vaps.xml` and also not `/tag1/tags2.xml`.
result = client.repo.get_calcs(paths='tag1/tag2').response().result
#   The results are paginated. The pagination key holds an object with total, page, per_page
#   kind of information
if result.pagination.total == 0:
    print('not found')
    sys.exit(1)
elif result.pagination.total > 1:
    print('my ids are not specific enough, bummer ...')
    sys.exit(1)
else:
    # The results key holds an array with the current page data
    calc = result.results[0]

# download data
#   via api
client.raw.get(upload_id=calc.upload_id, path=calc.mainfile).response()
#   via download
#   just the 'mainfile'
url = '%s/raw/%s/%s' % (nomad_url, calc.upload_id, calc.mainfile)
#   all files
url = '%s/raw/%s/%s/*' % (nomad_url, calc.upload_id, os.path.dirname(calc.mainfile))
