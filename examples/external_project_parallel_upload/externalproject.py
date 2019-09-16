"""
This example shows how to read files from many sources (here .tar.gz files),
chunk the data into even sized uploads and upload/process them in parallel.
"""

from bravado.requests_client import RequestsClient
from bravado.client import SwaggerClient
from urllib.parse import urlparse
import time
import os.path
import sys

# config
nomad_url = 'http://labdev-nomad.esc.rzg.mpg.de/fairdi/nomad/testing/api'
user = 'leonard.hofstadter@nomad-fairdi.tests.de'
password = 'password'
approx_upload_size = 1 * 1024  # 32 * 1024^3
parallel_uploads = 6

# create the bravado client
host = urlparse(nomad_url).netloc.split(':')[0]
http_client = RequestsClient()
http_client.set_basic_auth(host, user, password)
client = SwaggerClient.from_url('%s/swagger.json' % nomad_url, http_client=http_client)

def source_generator():
    """
    Yields all data sources.
    """
    yield os.path.join(os.path.dirname(__file__), 'example-1.tar.gz')
    yield os.path.join(os.path.dirname(__file__), 'example-2.tar.gz')
    yield os.path.join(os.path.dirname(__file__), 'example-3.tar.gz')

def source_file_generator(source):
    """
    Yields [filepath, file] tuples from :func:`source_generator`
    """
    pass

def calc_generator(source_file):
    pass

def upload_generator(zip_streams):
    """
    Yields nomad uploads that are already uploaded, but are still processing.
    """
    size = 0
      streamed_files: Set[str] = set()

    def generator():
        """ Stream a zip file with all files using zipstream. """
        def iterator():
            """
            Replace the directory based iter of zipstream with an iter over all given
            files.
            """
            for zipped_filename, upload_filename, upload_files in files:
                if zipped_filename in streamed_files:
                    continue
                streamed_files.add(zipped_filename)

                # Write a file to the zipstream.
                try:
                    with upload_files.raw_file(upload_filename, 'rb') as f:
                        def iter_content():
                            while True:
                                data = f.read(1024 * 64)
                                if not data:
                                    break
                                yield data

                        yield dict(arcname=zipped_filename, iterable=iter_content())
                except KeyError:
                    # files that are not found, will not be returned
                    pass
                except Restricted:
                    # due to the streaming nature, we cannot raise 401 here
                    # we just leave it out in the download
                    pass

        compression = zipfile.ZIP_DEFLATED if compress else zipfile.ZIP_STORED
        zip_stream = zipstream.ZipFile(mode='w', compression=compression, allowZip64=True)
        zip_stream.paths_to_write = iterator()

        for chunk in zip_stream:
            yield chunk

    response = Response(stream_with_context(generator()), mimetype='application/zip')
    response.headers['Content-Disposition'] = 'attachment; filename={}'.format(zipfile_name)
    return response


def upload_and_process():
    """
    Uses the chain of generators to upload data sequentially, awaits processing in
    parallel.
    """

upload_file = os.path.join(os.path.dirname(__file__), 'external_project_example.zip')



# upload data
print('uploading  a file with "external_id/AcAg/vasp.xml" inside ...')
with open(upload_file, 'rb') as f:
    upload = client.uploads.upload(file=f).response().result

print('processing ...')
while upload.tasks_running:
    upload = client.uploads.get_upload(upload_id=upload.upload_id).response().result
    time.sleep(5)
    print('processed: %d, failures: %d' % (upload.processed_calcs, upload.failed_calcs))

# check if processing was a success
if upload.tasks_status != 'SUCCESS':
    print('something went wrong')
    print('errors: %s' % str(upload.errors))
    # delete the unsuccessful upload
    client.uploads.delete_upload(upload_id=upload.upload_id).response().result
    sys.exit(1)

# publish data
print('publishing ...')
client.uploads.exec_upload_operation(upload_id=upload.upload_id, payload={
    'operation': 'publish',
    'metadata': {
        # these metadata are applied to all calcs in the upload
        'comment': 'Data from a cool external project',
        'references': ['http://external.project.eu'],
        'calculations': [
            {
                # these metadata are only applied to the calc identified by its 'mainfile'
                'mainfile': 'external_id/AcAg/vasp.xml',

                # 'coauthors': ['sheldon.cooper@ucla.edu'],  this does not YET work with emails,
                # Currently you have to use user_ids: leonard (the uploader, who is automatically an author) is 2 and sheldon is 1.
                # Ask NOMAD developers about how to find out about user_ids.
                'coauthors': [1],

                # If users demand, we can implement a specific metadata keys (e.g. 'external_id', 'external_url') for external projects.
                # This could allow to directly search for, or even have API endpoints that work with external_ids
                # 'external_id': 'external_id',
                # 'external_url': 'http://external.project.eu/data/calc/external_id/'
            }
        ]


    }
}).response().result

while upload.process_running:
    upload = client.uploads.get_upload(upload_id=upload.upload_id).response().result
    time.sleep(1)
if upload.tasks_status != 'SUCCESS' or len(upload.errors) > 0:
    print('something went wrong')
    print('errors: %s' % str(upload.errors))
    # delete the unsuccessful upload
    client.uploads.delete_upload(upload_id=upload.upload_id).response().result
    sys.exit(1)


# search for data
result = client.repo.search(paths=['external_id']).response().result
if result.pagination.total == 0:
    print('not found')
    sys.exit(1)
elif result.pagination.total > 1:
    print('my ids are not specific enough, bummer ... or did I uploaded stuff multiple times?')
# The results key holds an array with the current page data
print('Found the following calcs for my "external_id".')
print(', '.join(calc['calc_id'] for calc in result.results))

# download data
calc = result.results[0]
client.raw.get(upload_id=calc['upload_id'], path=calc['mainfile']).response()
print('Download of first calc works.')

# download urls, e.g. for curl
print('Possible download URLs are:')
print('%s/raw/%s/%s' % (nomad_url, calc['upload_id'], calc['mainfile']))
print('%s/raw/%s/%s/*' % (nomad_url, calc['upload_id'], os.path.dirname(calc['mainfile'])))

# direct download urls without having to search before
print('%s/raw/query?paths=external_id' % nomad_url)
