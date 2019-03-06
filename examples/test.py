"""
This is a brief example demonstrating the public nomad@FAIRDI API for doing operations
that might be necessary to integrate external project data.
"""

from bravado.requests_client import RequestsClient
from bravado.client import SwaggerClient
import math
from urllib.parse import urlparse
from concurrent.futures import ThreadPoolExecutor

# nomad_url = 'http://enc-staging-nomad.esc.rzg.mpg.de/fairdi/nomad/migration/api'
nomad_url = 'http://localhost:8000/nomad/api/'
user = 'admin'
password = 'password'

upload_file = 'external_project_example.zip'

# create the bravado client
host = urlparse(nomad_url).netloc.split(':')[0]
http_client = RequestsClient()
http_client.set_basic_auth(host, user, password)
client = SwaggerClient.from_url('%s/swagger.json' % nomad_url, http_client=http_client)

uploads = [upload.upload_id for upload in client.uploads.get_uploads().response().result]

executor = ThreadPoolExecutor(max_workers=10)


def run(upload_id):
    upload = client.uploads.get_upload(upload_id=upload_id).response().result
    upload_total_calcs = upload.calcs.pagination.total
    per_page = 200
    for page in range(1, math.ceil(upload_total_calcs / per_page) + 1):
        search = client.repo.search(
            page=page, per_page=per_page, order_by='mainfile',
            upload_id=upload_id).response().result

        print(search.pagination.page)


for upload in uploads:
    executor.submit(lambda: run(upload))
