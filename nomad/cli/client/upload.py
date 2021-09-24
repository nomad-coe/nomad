#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import os
import time
import click
import urllib.parse
import requests

from nomad import config
from nomad.processing import ProcessStatus

from .client import client


def stream_upload_with_client(client, stream, upload_name=None):
    user = client.auth.get_auth().response().result
    token = user.access_token
    url = config.client.url + '/uploads/'
    if upload_name is not None:
        url += '?upload_name=%s' % urllib.parse.quote(upload_name)

    response = requests.put(url, headers={'Authorization': 'Bearer %s' % token}, data=stream)
    if response.status_code != 200:
        raise Exception('nomad return status %d' % response.status_code)
    upload_id = response.json()['upload_id']

    return client.uploads.get_upload(upload_id=upload_id).response().result


def upload_file(file_path: str, upload_name: str = None, offline: bool = False, publish: bool = False, client=None):
    '''
    Upload a file to nomad.

    Arguments:
        file_path: path to the file, absolute or relative to call directory
        upload_name: optional name, default is the file_path's basename
        offline: allows to process data without upload, requires client to be run on the server
        publish: automatically publish after successful processing

    Returns: The upload_id
    '''
    if client is None:
        from nomad.cli.client import create_client
        client = create_client()
    if offline:
        upload = client.uploads.upload(
            local_path=os.path.abspath(file_path), upload_name=upload_name).response().result
        click.echo('process offline: %s' % file_path)
    else:
        # bravado does not seem to support streaming?
        try:
            with open(file_path, 'rb') as f:
                upload = stream_upload_with_client(client, f, upload_name=upload_name)
        except Exception as e:
            click.echo('could not upload the file: %s' % str(e))
            return

    while upload is not None and upload.process_status not in [ProcessStatus.SUCCESS, ProcessStatus.FAILURE]:
        upload = client.uploads.get_upload(upload_id=upload.upload_id).response().result
        calcs = upload.calcs.pagination
        if calcs is None:
            total, successes, failures = 0, 0, 0
        else:
            total, successes, failures = (calcs.total, calcs.successes, calcs.failures)

        ret = '\n' if upload.process_status in (ProcessStatus.SUCCESS, ProcessStatus.FAILURE) else '\r'

        print(
            'status: %s; process: %s; parsing: %d/%d/%d                %s' %
            (upload.process_status, upload.current_process, successes, failures, total, ret), end='')

        time.sleep(1)

    if upload.process_status == ProcessStatus.FAILURE:
        click.echo('There have been errors:')
        for error in upload.errors:
            click.echo('    %s' % error)
    elif publish:
        client.uploads.exec_upload_operation(upload_id=upload.upload_id, payload=dict(operation='publish')).response()

    return upload.upload_id


@client.command(
    help='Upload files to nomad. The given path can be a single file or a directory. '
    'All .zip files in a directory will be uploaded.')
@click.argument('PATH', nargs=-1, required=True, type=click.Path(exists=True))
@click.option(
    '--upload_name',
    help='Optional name for the upload of a single file. Will be ignored on directories.')
@click.option(
    '--offline', is_flag=True, default=False,
    help='Upload files "offline": files will not be uploaded, but processed were they are. '
    'Only works when run on the nomad host.')
@click.option(
    '--publish', is_flag=True, default=False,
    help='Automatically move upload out of the staging area after successful processing')
def upload(path, upload_name: str, offline: bool, publish: bool):
    paths = path
    click.echo('uploading files from %s paths' % len(paths))
    for path in paths:
        click.echo('uploading %s' % path)
        if os.path.isfile(path):
            upload_name = upload_name if upload_name is not None else os.path.basename(path)
            upload_file(path, upload_name, offline, publish)

        elif os.path.isdir(path):
            for (dirpath, _, filenames) in os.walk(path):
                for filename in filenames:
                    if filename.endswith('.zip'):
                        file_path = os.path.abspath(os.path.join(dirpath, filename))
                        upload_name = os.path.basename(file_path)
                        upload_file(file_path, upload_name, offline, publish)

        else:
            click.echo('Unknown path type %s.' % path)
