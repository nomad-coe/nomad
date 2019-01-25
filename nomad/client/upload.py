# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os.path
import os
import time
import click

from nomad import utils
from nomad.processing import FAILURE, SUCCESS

from .main import cli, create_client


def upload_file(file_path: str, name: str = None, offline: bool = False, commit: bool = False, client=None):
    """
    Upload a file to nomad.

    Arguments:
        file_path: path to the file, absolute or relative to call directory
        name: optional name, default is the file_path's basename
        offline: allows to process data without upload, requires client to be run on the server
        commit: automatically commit after successful processing

    Returns: The upload_id
    """
    if client is None:
        client = create_client()
    if offline:
        upload = client.uploads.upload(
            local_path=os.path.abspath(file_path), name=name).reponse().result
        click.echo('process offline: %s' % file_path)
    else:
        with open(file_path, 'rb') as f:
            upload = client.uploads.upload(file=f, name=name).response().result
        click.echo('process online: %s' % file_path)

    while upload.tasks_status not in [SUCCESS, FAILURE]:
        upload = client.uploads.get_upload(upload_id=upload.upload_id).response().result
        calcs = upload.calcs.pagination
        if calcs is None:
            total, successes, failures = 0, 0, 0
        else:
            total, successes, failures = (calcs.total, calcs.successes, calcs.failures)

        ret = '\n' if upload.tasks_status in (SUCCESS, FAILURE) else '\r'

        print(
            'status: %s; task: %s; parsing: %d/%d/%d                %s' %
            (upload.tasks_status, upload.current_task, successes, failures, total, ret), end='')

        time.sleep(3)

    if upload.tasks_status == FAILURE:
        click.echo('There have been errors:')
        for error in upload.errors:
            click.echo('    %s' % error)
    elif commit:
        client.uploads.exec_upload_command(upload_id=upload.upload_id, command='commit').reponse()

    return upload.upload_id


@cli.command(
    help='Upload files to nomad. The given path can be a single file or a directory. '
    'All .zip files in a directory will be uploaded.')
@click.argument('PATH', nargs=-1, required=True, type=click.Path(exists=True))
@click.option(
    '--name',
    help='Optional name for the upload of a single file. Will be ignored on directories.')
@click.option(
    '--offline', is_flag=True, default=False,
    help='Upload files "offline": files will not be uploaded, but processed were they are. '
    'Only works when run on the nomad host.')
@click.option(
    '--commit', is_flag=True, default=False,
    help='Automatically move upload out of the staging area after successful processing')
def upload(path, name: str, offline: bool, commit: bool):
    utils.configure_logging()
    paths = path
    click.echo('uploading files from %s paths' % len(paths))
    for path in paths:
        click.echo('uploading %s' % path)
        if os.path.isfile(path):
            name = name if name is not None else os.path.basename(path)
            upload_file(path, name, offline, commit)

        elif os.path.isdir(path):
            for (dirpath, _, filenames) in os.walk(path):
                for filename in filenames:
                    if filename.endswith('.zip'):
                        file_path = os.path.abspath(os.path.join(dirpath, filename))
                        name = os.path.basename(file_path)
                        upload_file(file_path, name, offline, commit)

        else:
            click.echo('Unknown path type %s.' % path)
