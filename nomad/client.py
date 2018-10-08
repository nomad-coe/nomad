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

from nomad import config
from nomad.files import UploadFile


api_base = 'http://localhost/nomad/api'
user = 'other@gmail.com'
pw = 'nomad'


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
def upload_file(file_path: str, name: str = None, offline: bool = False):
    """
    Upload a file to nomad.

    Arguments:
        file_path: path to the file, absolute or relative to call directory
        name: optional name, default is the file_path's basename
        offline: allows to process data without upload, requires client to be run on the server

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
        upload_cmd = upload_cmd.replace('local_file', file_path)

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

        time.sleep(3)

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


class CalcProcReproduction(UploadFile):
    """
    Instances represent a local reproduction of the processing for a single calculation.
    It allows to download raw data from a nomad server and reproduce its processing
    (parsing, normalizing) with the locally installed parsers and normalizers.

    The use-case is error/warning reproduction. Use ELK to identify errors, use
    the upload, archive ids/hashes to given by ELK, and reproduce and fix the error
    in your development environment.

    This is a class of :class:`UploadFile` the downloaded raw data will be treated as
    an fake 'upload' that only contains the respective calculation data. This allows us
    to locally run processing code that is very similar to the one used on the server.
    """
    def __init__(self, archive_id: str) -> CalcProcReproduction:
        local_path = os.path.join(config.fs.tmp, '%s.zip' % archive_id)
        if not os.path.exists(local_path):
            # download raw if not already downloaded
            req = requests.get('%s/raw/%s?all=1' % (api_base, archive_id), stream=True)
            with open(local_path, 'wb') as f:
                for chunk in req.iter_content():
                    f.write(chunk)

        super().__init__(upload_id='tmp_%s' % archive_id, local_path=local_path)

    def parse(parser_name: str = None):
        """
        Run the given parser on the downloaded calculation. If no parser is given,
        do parser matching and use the respective parser.
        """
        pass

    def normalize(normalizer_name: str):
        """
        Parse the downloaded calculation and run the given normalizer.
        """
        pass

    def normalize_all():
        """
        Parse the downloaded calculation and run the whole normalizer chain.
        """
        pass


@click.group()
@click.option('--host', default='localhost', help='The host nomad runs on, default is "localhost".')
@click.option('--port', default=80, help='the port nomad runs with, default is 80.')
def cli(host: str, port: int):
    global api_base
    api_base = 'http://%s:%d/nomad/api' % (host, port)


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
def upload(path, name: str, offline: bool):
    paths = path
    click.echo('uploading files from %s paths' % len(paths))
    for path in paths:
        click.echo('uploading %s' % path)
        if os.path.isfile(path):
            name = name if name is not None else os.path.basename(path)
            upload_file(path, name, offline)

        elif os.path.isdir(path):
            for file_path in walk_through_files(path):
                name = os.path.basename(file_path)
                upload_file(file_path, name, offline)

        else:
            click.echo('Unknown path type %s.' % path)


@cli.command(help='Attempts to reset the nomad.')
def reset():
    response = requests.post('%s/admin/reset' % api_base, auth=HTTPBasicAuth(user, pw))
    if response.status_code != 200:
        click.echo('API return %s' % str(response.status_code))
        click.echo(response.text)
        sys.exit(1)


@cli.group(help='Run a nomad service locally (outside docker).')
def run():
    pass


@run.command(help='Run the nomad development worker.')
def worker():
    from nomad import processing
    processing.app.worker_main(['worker', '--loglevel=INFO'])


@run.command(help='Run the nomad development api.')
def api():
    from nomad import infrastructure, api
    infrastructure.setup()
    api.app.run(debug=True, port=8000)


@cli.command(help='Runs tests and linting. Useful before commit code.')
def qa():
    os.chdir(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    ret_code = 0
    click.echo('Run tests ...')
    ret_code += os.system('python -m pytest tests')
    click.echo('Run code style checks ...')
    ret_code += os.system('python -m pycodestyle --ignore=E501,E701 nomad tests')
    click.echo('Run linter ...')
    ret_code += os.system('python -m pylint --load-plugins=pylint_mongoengine nomad tests')
    click.echo('Run static type checks ...')
    ret_code += os.system('python -m mypy --ignore-missing-imports --follow-imports=silent --no-strict-optional nomad tests')

    sys.exit(ret_code)


if __name__ == '__main__':
    cli()  # pylint: disable=E1120
