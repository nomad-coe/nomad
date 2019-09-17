"""
This example shows how to read files from many sources (here .tar.gz files),
chunk the data into even sized uploads and upload/process them in parallel. The assumption
is that each source file is much smaller than the targeted upload size.
"""

from typing import Iterator, Iterable, Union, Tuple
from bravado.requests_client import RequestsClient
from bravado.client import SwaggerClient
from urllib.parse import urlparse, urlencode
import requests
import time
import os.path
import tarfile
import io
import zipfile
import zipstream

# config
nomad_url = 'http://labdev-nomad.esc.rzg.mpg.de/fairdi/nomad/testing/api'
user = 'leonard.hofstadter@nomad-fairdi.tests.de'
password = 'password'
approx_upload_size = 32 * 1024 * 1024 * 1024  # you can make it really small for testing
max_parallel_uploads = 9

# create the bravado client
host = urlparse(nomad_url).netloc.split(':')[0]
http_client = RequestsClient()
http_client.set_basic_auth(host, user, password)
client = SwaggerClient.from_url('%s/swagger.json' % nomad_url, http_client=http_client)


def source_generator() -> Iterable[Tuple[str, Union[str, None]]]:
    """
    Yields all data sources. Yields tuples (path to .tgz, prefix). Prefix denotes
    a subdirectory to put the contents in. Use None for no prefix.
    """
    yield os.path.join(os.path.dirname(__file__), 'example-1.tar.gz'), 'example_1'
    yield os.path.join(os.path.dirname(__file__), 'example-2.tar.gz'), 'example_2'
    yield os.path.join(os.path.dirname(__file__), 'example-3.tar.gz'), 'example_3'


def upload_next_data(sources: Iterator[Tuple[str, str]], upload_name='next upload'):
    """
    Reads data from the given sources iterator. Creates and uploads a .zip-stream of
    approx. size. Returns the upload, or raises StopIteration if the sources iterator
    was empty. Should be used repeatedly on the same iterator until it is empty.
    """

    # potentially raises StopIteration before being streamed
    first_source = next(sources)

    def iterator():
        """
        Yields dicts with keys arcname, iterable, as required for the zipstream
        library. Will read from generator until the zip-stream has the desired size.
        """
        size = 0
        first = True
        while(True):
            if first:
                source_file, prefix = first_source
                first = False
            else:
                try:
                    source_file, prefix = next(sources)
                except StopIteration:
                    break

            source_tar = tarfile.open(source_file)
            source = source_tar.fileobj
            bufsize = source_tar.copybufsize
            for source_member in source_tar.getmembers():
                if not source_member.isfile():
                    continue

                target = io.BytesIO()
                source.seek(source_member.offset_data)
                tarfile.copyfileobj(  # type: ignore
                    source, target, source_member.size, tarfile.ReadError, bufsize)

                size += source_member.size
                target.seek(0)

                def iter_content():
                    while True:
                        data = target.read(io.DEFAULT_BUFFER_SIZE)
                        if not data:
                            break
                        yield data

                name = source_member.name
                if prefix is not None:
                    name = os.path.join(prefix, name)

                yield dict(arcname=source_member.name, iterable=iter_content())

            if size > approx_upload_size:
                break

    # create the zip-stream from the iterator above
    zip_stream = zipstream.ZipFile(mode='w', compression=zipfile.ZIP_STORED, allowZip64=True)
    zip_stream.paths_to_write = iterator()

    zip_stream

    user = client.auth.get_user().response().result
    token = user.token
    url = nomad_url + '/uploads/?%s' % urlencode(dict(name=upload_name))

    def content():
        for chunk in zip_stream:
            if len(chunk) != 0:
                yield chunk

    # stream .zip to nomad
    response = requests.put(url=url, headers={'X-Token': token}, data=content())

    if response.status_code != 200:
        raise Exception('nomad return status %d' % response.status_code)

    upload_id = response.json()['upload_id']

    return client.uploads.get_upload(upload_id=upload_id).response().result


def publish_upload(upload):
    client.uploads.exec_upload_operation(upload_id=upload.upload_id, payload={
        'operation': 'publish',
        'metadata': {
            # these metadata are applied to all calcs in the upload
            'comment': 'Data from a cool external project',
            'references': ['http://external.project.eu']
        }
    }).response()


if __name__ == '__main__':
    source_iter = iter(source_generator())
    all_uploaded = False
    processing_completed = False

    # run until there are no more uploads and everything is processed (and published)
    while not (all_uploaded and processing_completed):
        # process existing uploads
        while True:
            uploads = client.uploads.get_uploads().response().result
            for upload in uploads.results:
                if not upload.process_running:
                    if upload.tasks_status == 'SUCCESS':
                        print('publish %s(%s)' % (upload.name, upload.upload_id))
                        publish_upload(upload)
                    elif upload.tasks_status == 'FAILURE':
                        print('could not process %s(%s)' % (upload.name, upload.upload_id))
                        client.uploads.delete_upload(upload_id=upload.upload_id).response().result

            if uploads.pagination.total < max_parallel_uploads:
                # still processing some, but there is room for more uploads
                break
            else:
                # wait for processing
                time.sleep(10)

        # add a new upload
        if all_uploaded:
            processing_completed = uploads.pagination.total == 0

        try:
            upload = upload_next_data(source_iter)
            processing_completed = False
            print('uploaded %s(%s)' % (upload.name, upload.upload_id))
        except StopIteration:
            all_uploaded = True
        except Exception as e:
            print('could not upload next upload: %s' % str(e))
