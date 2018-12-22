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
Swagger/bravado based python client library for the API and various usefull shell commands.
"""

import os.path
import os
import sys
import time
import requests
import click
from typing import Union, Callable, cast
import logging
from bravado.requests_client import RequestsClient
from bravado.client import SwaggerClient

from nomad import config, utils
from nomad.files import UploadFile
from nomad.parsing import parsers, parser_dict, LocalBackend
from nomad.normalizing import normalizers


api_base = 'http://%s:%d/%s' % (config.services.api_host, config.services.api_port, config.services.api_base_path)
user = 'leonard.hofstadter@nomad-fairdi.tests.de'
pw = 'password'


def _cli_client():
    return create_client()


def create_client(
        host: str = config.services.api_host,
        port: int = config.services.api_port,
        base_path: str = config.services.api_base_path,
        user: str = user, password: str = None):
    """ A factory method to create the client. """

    if user is not None:
        http_client = RequestsClient()
        http_client.set_basic_auth(host, user, pw)
    else:
        http_client = None

    client = SwaggerClient.from_url(
        'http://%s:%d%s/swagger.json' % (host, port, base_path),
        http_client=http_client)

    return client


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


def upload_file(file_path: str, name: str = None, offline: bool = False, unstage: bool = False):
    """
    Upload a file to nomad.

    Arguments:
        file_path: path to the file, absolute or relative to call directory
        name: optional name, default is the file_path's basename
        offline: allows to process data without upload, requires client to be run on the server
        unstage: automatically unstage after successful processing
    """
    client = _cli_client()
    if offline:
        upload = client.uploads.upload(
            local_path=os.path.abspath(file_path), name=name).reponse().result
        click.echo('process offline: %s' % file_path)
    else:
        with open(file_path, 'rb') as f:
            upload = client.uploads.upload(file=f, name=name).response().result
        click.echo('process online: %s' % file_path)

    while upload.status not in ['SUCCESS', 'FAILURE']:
        upload = client.uploads.get_upload(upload_id=upload.upload_id).response().result
        calcs = upload.calcs.pagination
        if calcs is None:
            total, successes, failures = 0, 0, 0
        else:
            total, successes, failures = (calcs.total, calcs.successes, calcs.failures)

        ret = '\n' if upload.status in ('SUCCESS', 'FAILURE') else '\r'

        print(
            'status: %s; task: %s; parsing: %d/%d/%d                %s' %
            (upload.status, upload.current_task, successes, failures, total, ret), end='')

        time.sleep(3)

    if upload.status == 'FAILURE':
        click.echo('There have been errors:')
        for error in upload.errors:
            click.echo('    %s' % error)
    elif unstage:
        client.uploads.exec(upload_id=upload.upload_id, operation='unstage').reponse()


class CalcProcReproduction:
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

    Arguments:
        archive_id: The archive_id of the calculation to locally process.
        override: Set to true to override any existing local calculation data.
    """
    def __init__(self, archive_id: str, override: bool = False) -> None:
        self.calc_hash = utils.archive.calc_hash(archive_id)
        self.upload_hash = utils.archive.upload_hash(archive_id)
        self.mainfile = None
        self.parser = None
        self.logger = utils.get_logger(__name__, archive_id=archive_id)

        local_path = os.path.join(config.fs.tmp, 'repro_%s.zip' % archive_id)
        if not os.path.exists(os.path.dirname(local_path)):
            os.makedirs(os.path.dirname(local_path))
        if not os.path.exists(local_path) or override:
            # download raw if not already downloaded or if override is set
            # download with request, since bravado does not support streaming
            # TODO currently only downloads mainfile
            self.logger.info('Downloading calc.')
            req = requests.get('%s/raw/%s/%s' % (api_base, self.upload_hash, os.path.dirname(self.mainfile)), stream=True)
            with open(local_path, 'wb') as f:
                for chunk in req.iter_content(chunk_size=1024):
                    f.write(chunk)
        else:
            self.logger.info('Calc already downloaded.')

        self.upload_file = UploadFile(upload_id='tmp_%s' % archive_id, local_path=local_path)

    def __enter__(self):
        # open/extract upload file
        self.logger.info('Extracting calc data.')
        self.upload_file.__enter__()

        # find mainfile matching calc_hash
        self.mainfile = next(
            filename for filename in self.upload_file.filelist
            if utils.hash(filename) == self.calc_hash)

        assert self.mainfile is not None, 'The mainfile could not be found.'
        self.logger = self.logger.bind(mainfile=self.mainfile)
        self.logger.info('Identified mainfile.')

        return self

    def __exit__(self, *args):
        self.upload_file.__exit__(*args)

    def parse(self, parser_name: str = None) -> LocalBackend:
        """
        Run the given parser on the downloaded calculation. If no parser is given,
        do parser matching and use the respective parser.
        """
        mainfile = self.upload_file.get_file(self.mainfile)
        if parser_name is not None:
            parser = parser_dict.get(parser_name)
        else:
            for potential_parser in parsers:
                with mainfile.open() as mainfile_f:
                    if potential_parser.is_mainfile(self.mainfile, lambda fn: mainfile_f):
                        parser = potential_parser
                        break

        assert parser is not None, 'there is not parser matching %s' % self.mainfile
        self.logger = self.logger.bind(parser=parser.name)  # type: ignore
        self.logger.info('identified parser')

        parser_backend = parser.run(mainfile.os_path, logger=self.logger)
        self.logger.info('ran parser')
        return parser_backend

    def normalize(self, normalizer: Union[str, Callable], parser_backend: LocalBackend = None):
        """
        Parse the downloaded calculation and run the given normalizer.
        """
        if parser_backend is None:
            parser_backend = self.parse()

        if isinstance(normalizer, str):
            normalizer = next(
                normalizer_instance for normalizer_instance in normalizers
                if normalizer_instance.__class__.__name__ == normalizer)

        assert normalizer is not None, 'there is no normalizer %s' % str(normalizer)
        normalizer_instance = cast(Callable, normalizer)(parser_backend)
        logger = self.logger.bind(normalizer=normalizer_instance.__class__.__name__)
        self.logger.info('identified normalizer')

        normalizer_instance.normalize(logger=logger)
        self.logger.info('ran normalizer')
        return parser_backend

    def normalize_all(self, parser_backend: LocalBackend = None):
        """
        Parse the downloaded calculation and run the whole normalizer chain.
        """
        for normalizer in normalizers:
            parser_backend = self.normalize(normalizer, parser_backend=parser_backend)

        return parser_backend


@click.group()
@click.option('-h', '--host', default=config.services.api_host, help='The host nomad runs on, default is "%s".' % config.services.api_host)
@click.option('-p', '--port', default=config.services.api_port, help='the port nomad runs with, default is %d.' % config.services.api_port)
@click.option('-u', '--user', default=None, help='the user name to login, default no login.')
@click.option('-w', '--password', default=None, help='the password use to login.')
@click.option('-v', '--verbose', help='sets log level to debug', is_flag=True)
def cli(host: str, port: int, verbose: bool, user: str, password: str):
    if verbose:
        config.console_log_level = logging.DEBUG
    else:
        config.console_log_level = logging.WARNING

    global api_base
    api_base = 'http://%s:%d/nomad/api' % (host, port)

    global _cli_client

    def _cli_client():  # pylint: disable=W0612
        if user is not None:
            return create_client(host=host, port=port, user=user, password=password)
        else:
            return create_client(host=host, port=port)


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
    '--unstage', is_flag=True, default=False,
    help='Automatically move upload out of the staging area after successful processing')
def upload(path, name: str, offline: bool, unstage: bool):
    utils.configure_logging()
    paths = path
    click.echo('uploading files from %s paths' % len(paths))
    for path in paths:
        click.echo('uploading %s' % path)
        if os.path.isfile(path):
            name = name if name is not None else os.path.basename(path)
            upload_file(path, name, offline, unstage)

        elif os.path.isdir(path):
            for (dirpath, _, filenames) in os.walk(path):
                for filename in filenames:
                    if filename.endswith('.zip'):
                        file_path = os.path.abspath(os.path.join(dirpath, filename))
                        name = os.path.basename(file_path)
                        upload_file(file_path, name, offline, unstage)

        else:
            click.echo('Unknown path type %s.' % path)


@cli.command(help='Attempts to reset the nomad.')
def reset():
    _cli_client().admin.exec(operation='reset').reponse()


@cli.command(help='Run processing locally.')
@click.argument('ARCHIVE_ID', nargs=1, required=True, type=str)
@click.option(
    '--override', is_flag=True, default=False,
    help='Override existing local calculation data.')
def local(archive_id, **kwargs):
    utils.configure_logging()
    utils.get_logger(__name__).info('Using %s' % api_base)
    with CalcProcReproduction(archive_id, **kwargs) as local:
        backend = local.parse()
        local.normalize_all(parser_backend=backend)
        # backend.write_json(sys.stdout, pretty=True)


@cli.group(help='Run a nomad service locally (outside docker).')
def run():
    pass


@run.command(help='Run the nomad development worker.')
def worker():
    config.service = 'nomad_worker'
    from nomad import processing
    processing.app.worker_main(['worker', '--loglevel=INFO'])


@run.command(help='Run the nomad development api.')
def api():
    config.service = 'nomad_api'
    from nomad import infrastructure
    from nomad.api.__main__ import run_dev_server
    infrastructure.setup()
    run_dev_server(debug=True, port=8000)


@cli.command(help='Runs tests and linting. Useful before commit code.')
@click.option('--skip-tests', help='Do not test, just do code checks.', is_flag=True)
def qa(skip_tests: bool):
    os.chdir(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    ret_code = 0
    if not skip_tests:
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
