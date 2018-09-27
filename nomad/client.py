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

"""
Simple client library for the nomad api that allows to bulk upload files via shell command.
"""

import os.path
import os
import sys
import subprocess
import shlex
import time
import requests
from requests.auth import HTTPBasicAuth
import click

api_base = 'http://localhost/nomad/api'


def handle_common_errors(func):
    def wrapper(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except requests.exceptions.ConnectionError:
            click.echo(
                '\nCould not connect to nomad at %s. '
                'Check connection and host/port options.' % api_base)
            sys.exit(0)
    return wrapper


@handle_common_errors
def upload_file(file_path, name=None, offline=False, user='other@gmail.com', pw='nomad'):
    """
    Upload a file to nomad.

    Arguments:
        file_path: Path to the file, absolute or relative to call directory.
        name: Optional name, default is the file_path's basename
    """
    auth = HTTPBasicAuth(user, pw)

    if name is None:
        name = os.path.basename(file_path)

    post_data = dict(name=name)
    if offline:
        post_data.update(dict(local_path=os.path.abspath(file_path)))
        click.echo('process offline: %s' % file_path)

    upload = requests.post('%s/uploads' % api_base, json=post_data, auth=auth).json()

    if not offline:
        upload_cmd = upload['upload_command']
        upload_cmd = upload_cmd.replace('your_file', file_path)

        subprocess.call(shlex.split(upload_cmd))

        click.echo('uploaded: %s' % file_path)

    while True:
        upload = requests.get('%s/uploads/%s' % (api_base, upload['upload_id']), auth=auth).json()
        status = upload['status']
        calcs_pagination = upload['calcs'].get('pagination')
        if calcs_pagination is None:
            total, successes, failures = 0, 0, 0
        else:
            total, successes, failures = (
                calcs_pagination[key] for key in ('total', 'successes', 'failures'))

        ret = '\n' if status in ('SUCCESS', 'FAILURE') else '\r'

        print(
            'status: %s; task: %s; parsing: %d/%d/%d                %s' %
            (status, upload['current_task'], successes, failures, total, ret), end='')

        if status in ('SUCCESS', 'FAILURE'):
            break

        time.sleep(5)

    if status == 'FAILURE':
        click.echo('There have been errors:')
        for error in upload['errors']:
            click.echo('    %s' % error)


def walk_through_files(path, extension='.zip'):
    """
    Returns all abs path of all files in a sub tree of the given path that match
    the given extension.

    Arguments:
        path (str): the directory
        extension (str): the extension, incl. '.', e.g. '.zip' (default)
    """

    for (dirpath, _, filenames) in os.walk(path):
        for filename in filenames:
            if filename.endswith(extension):
                yield os.path.abspath(os.path.join(dirpath, filename))


@click.group()
@click.option('--host', default='localhost', help='The host nomad runs on, default is "localhost".')
@click.option('--port', default=80, help='the port nomad runs with, default is 80.')
def cli(host: str, port: int):
    global api_base
    api_base = 'http://%s:%d/nomad/api' % (host, port)


@cli.command(
    help='Upload files to nomad. The given path can be a single file or a directory. '
    'All .zip files in a directory will be uploaded.')
@click.argument('PATH', type=click.Path(exists=True))
@click.option(
    '--name',
    help='Optional name for the upload of a single file. Will be ignored on directories.')
@click.option(
    '--offline', is_flag=True, default=False,
    help='Upload files "offline": files will not be uploaded, but processed were they are. '
    'Only works when run on the nomad host.')
def upload(path, name: str, offline: bool):
    if os.path.isfile(path):
        name = name if name is not None else os.path.basename(path)
        upload_file(path, name, offline)

    elif os.path.isdir(path):
        for file in walk_through_files(path):
            name = os.path.basename(file)
            upload_file(path, name, offline)

    else:
        sys.exit(0)


if __name__ == '__main__':
    cli()  # pylint: disable=E1120
