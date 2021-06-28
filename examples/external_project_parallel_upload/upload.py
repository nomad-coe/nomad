"""
This example shows how to read files from many sources (here .tar.gz files),
chunk the data into even sized uploads and upload/process them in parallel. The assumption
is that each source file is much smaller than the targeted upload size.
"""

from typing import Iterator, Iterable, Union, Tuple, Dict, Any
from bravado.requests_client import RequestsClient
from bravado.client import SwaggerClient
from urllib.parse import urlparse, urlencode
import requests
import re
import time
import os
import os.path
import tarfile
import io
import zipfile
import zipstream
import uuid

# config
nomad_url = 'http://labdev-nomad.esc.rzg.mpg.de/fairdi/nomad/mp/api'
user = 'leonard.hofstadter@nomad-fairdi.tests.de'
password = 'password'
approx_upload_size = 32 * 1024 * 1024 * 1024  # you can make it really small for testing
max_parallel_uploads = 9
direct_stream = False

# create the bravado client
host = urlparse(nomad_url).netloc.split(':')[0]
http_client = RequestsClient()
http_client.set_basic_auth(host, user, password)
client = SwaggerClient.from_url('%s/swagger.json' % nomad_url, http_client=http_client)


def source_generator() -> Iterable[Tuple[str, Union[str, None], Union[str, None]]]:
    """
    Yields all data sources. Yields tuples (path to .tgz, prefix, external_id). Prefix denotes
    a subdirectory to put the contents in. Use None for no prefix. The external_id is
    used when "publishing" the data to populate the external_id field.
    """
    yield os.path.join(os.path.dirname(__file__), 'example-1.tar.gz'), 'example_1', 'external_1'
    yield os.path.join(os.path.dirname(__file__), 'example-2.tar.gz'), 'example_2', 'external_2'
    yield os.path.join(os.path.dirname(__file__), 'example-3.tar.gz'), 'example_3', 'external_3'


def upload_next_data(sources: Iterator[Tuple[str, str, str]], upload_name='next upload'):
    """
    Reads data from the given sources iterator. Creates and uploads a .zip-stream of
    approx. size. Returns the upload and corresponding metadata, or raises StopIteration
    if the sources iterator was empty. Should be used repeatedly on the same iterator
    until it is empty.
    """

    # potentially raises StopIteration before being streamed
    first_source = next(sources)
    calc_metadata = []

    def iterator():
        """
        Yields dicts with keys arcname, iterable, as required for the zipstream
        library. Will read from generator until the zip-stream has the desired size.
        """
        size = 0
        first = True
        while(True):
            if first:
                source_file, prefix, external_id = first_source
                first = False
            else:
                try:
                    source_file, prefix, external_id = next(sources)
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

                if re.search(r'vasp(run)?\.xml(.gz)?$', name):
                    calc_metadata.append(dict(
                        mainfile=name,
                        external_id=external_id
                        # references=['http://external.project.eu/.../<your_id>']
                    ))

                # there was bug using the wrong name (source_member.name) here
                yield dict(arcname=name, iterable=iter_content(), buffer_size=source_member.size)

            if size > approx_upload_size:
                break

    # create the zip-stream from the iterator above
    zip_stream = zipstream.ZipFile(mode='w', compression=zipfile.ZIP_STORED, allowZip64=True)
    zip_stream.paths_to_write = iterator()

    user = client.auth.get_user().response().result
    token = user.token
    url = nomad_url + '/uploads/?%s' % urlencode(dict(name=upload_name))

    def content():
        for chunk in zip_stream:
            if len(chunk) != 0:
                yield chunk

    if direct_stream:
        # stream .zip to nomad
        response = requests.put(url=url, headers={'X-Token': token, 'Content-type': 'application/octet-stream'}, data=content())

    else:
        # save .zip and upload file to nomad
        zipfile_name = '/tmp/%s.zip' % str(uuid.uuid4())
        with open(zipfile_name, 'wb') as f:
            for c in content():
                f.write(c)
        try:
            with open(zipfile_name, 'rb') as f:
                response = requests.put(url=url, headers={'X-Token': token, 'Content-type': 'application/octet-stream'}, data=f)
        finally:
            os.remove(zipfile_name)

    if response.status_code != 200:
        raise Exception('nomad return status %d' % response.status_code)

    upload_id = response.json()['upload_id']

    return client.uploads.get_upload(upload_id=upload_id).response().result, calc_metadata


def publish_upload(upload, calc_metadata):
    metadata = {
        # these metadata are applied to all calcs in the upload
        'comment': 'Data from a cool external project',
        'references': ['http://external.project.eu'],
        # 'coauthors': [<nomad_user_id>, <nomad_user_id>, <nomad_user_id>]
        # these are calc specific metadata that supercede any upload metadata
        'calculations': calc_metadata}

    client.uploads.exec_upload_operation(upload_id=upload.upload_id, payload={
        'operation': 'publish',
        'metadata': metadata
    }).response()


if __name__ == '__main__':
    source_iter = iter(source_generator())
    all_uploaded = False
    processing_completed = False
    all_calc_metadata: Dict[str, Any] = {}

    # run until there are no more uploads and everything is processed (and published)
    while not (all_uploaded and processing_completed):
        # process existing uploads
        while True:
            uploads = client.uploads.get_uploads().response().result

            for upload in uploads.results:
                calc_metadata = all_calc_metadata.get(upload.upload_id, None)
                if calc_metadata is None:
                    continue

                if not upload.process_running:
                    if upload.process_status == 'SUCCESS':
                        print('publish %s(%s)' % (upload.name, upload.upload_id))
                        publish_upload(upload, calc_metadata)
                    elif upload.process_status == 'FAILURE':
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
            upload, calc_metadata = upload_next_data(source_iter)
            all_calc_metadata[upload.upload_id] = calc_metadata
            processing_completed = False
            print('uploaded %s(%s)' % (upload.name, upload.upload_id))
        except StopIteration:
            all_uploaded = True
        except Exception as e:
            print('could not upload next upload: %s' % str(e))
