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

import os
import io
import requests
import click
import typing
import sys
import bravado.exception
import json

from nomad import config, utils
from nomad import files
from nomad import datamodel
from nomad.cli import parse as cli_parse

from .client import client


class CalcProcReproduction:
    '''
    Instances represent a local reproduction of the processing for a single calculation.
    It allows to download raw data from a nomad server and reproduce its processing
    (parsing, normalizing) with the locally installed parsers and normalizers.

    The use-case is error/warning reproduction. Use ELK to identify errors, use
    the upload, archive ids to given by ELK, and reproduce and fix the error
    in your development environment.

    Arguments:
        calc_id: The calc_id of the calculation to locally process.
        override: Set to true to override any existing local calculation data.
    '''
    def __init__(self, archive_id: str, override: bool = False, mainfile: str = None) -> None:
        if '/' in archive_id:
            self.calc_id = utils.archive.calc_id(archive_id)
            self.upload_id = utils.archive.upload_id(archive_id)
        else:
            self.calc_id = archive_id
            self.upload_id = 'unknown'

        self.logger = utils.get_logger(__name__, upload_id=self.upload_id, calc_id=self.calc_id)

        self.mainfile = mainfile
        self.parser = None

        from nomad.cli.client import create_client
        client = create_client()
        if self.mainfile is None:
            try:
                calc = client.repo.get_repo_calc(
                    upload_id=self.upload_id, calc_id=self.calc_id).response().result
            except bravado.exception.HTTPNotFound:
                self.logger.error(
                    'could not find calculation, maybe it was deleted or the ids are wrong')
                sys.exit(1)
            self.mainfile = calc['mainfile']
            self.upload_id = calc['upload_id']
            self.logger = self.logger.bind(upload_id=self.upload_id, mainfile=self.mainfile)
        else:
            self.logger = self.logger.bind(mainfile=self.mainfile)
            self.logger.info('Using provided mainfile', mainfile=self.mainfile)

        local_path = os.path.join(config.fs.tmp, 'repro_%s.zip' % archive_id)
        if not os.path.exists(os.path.dirname(local_path)):
            os.makedirs(os.path.dirname(local_path))

        if not os.path.exists(local_path) or override:
            # download raw if not already downloaded or if override is set
            # download with request, since bravado does not support streaming
            self.logger.info('Downloading calc.', mainfile=self.mainfile)
            try:
                token = client.auth.get_auth().response().result.signature_token
                dir_name = os.path.dirname(self.mainfile)
                req = requests.get(
                    '%s/raw/%s/%s/*?signature_token=%s' % (config.client.url, self.upload_id, dir_name, token),
                    stream=True, headers={})
                with open(local_path, 'wb') as f:
                    for chunk in req.iter_content(chunk_size=io.DEFAULT_BUFFER_SIZE):
                        f.write(chunk)
            except bravado.exception.HTTPNotFound:
                self.logger.error(
                    'could not find calculation, maybe it was deleted or the ids are wrong')
                sys.exit(1)
            except Exception as e:
                self.logger.error('could not download calculation', exc_info=e)
                sys.exit(1)
        else:
            self.logger.info('Calc already downloaded.')

        self.upload_files = files.ArchiveBasedStagingUploadFiles(
            upload_id='tmp_%s' % archive_id, upload_path=local_path, create=True,
            is_authorized=lambda: True)

    def __enter__(self):
        # open/extract upload file
        self.logger.info('Extracting calc data.')
        try:
            self.upload_files.extract()
        except AssertionError:
            # already extracted
            pass

        assert self.mainfile is not None
        self.logger = self.logger.bind(mainfile=self.mainfile)
        self.logger.info('Identified mainfile.')

        return self

    def __exit__(self, *args):
        self.upload_files.delete()

    def parse(self, parser_name: str = None, **kwargs):
        '''
        Run the given parser on the downloaded calculation. If no parser is given,
        do parser matching and use the respective parser.
        '''
        return cli_parse.parse(
            self.upload_files.raw_file_object(self.mainfile).os_path,
            parser_name=parser_name, logger=self.logger, **kwargs)

    def normalize(self, normalizer: typing.Union[str, typing.Callable], parser_backend=None):
        '''
        Parse the downloaded calculation and run the given normalizer.
        '''
        if parser_backend is None:
            parser_backend = self.parse()

        return cli_parse.normalize(parser_backend=parser_backend, normalizer=normalizer, logger=self.logger)

    def normalize_all(self, parser_backend=None):
        '''
        Parse the downloaded calculation and run the whole normalizer chain.
        '''
        return cli_parse.normalize_all(parser_backend=parser_backend, logger=self.logger)


@client.command(help='Run processing locally.')
@click.argument('CALC_ID', nargs=1, required=True, type=str)
@click.option('--override', is_flag=True, help='Override existing local calculation data.')
@click.option('--show-backend', is_flag=True, help='Print the backend data.')
@click.option('--show-metadata', is_flag=True, help='Print the extracted repo metadata.')
@click.option('--mainfile', default=None, type=str, help='Use this mainfile (in case mainfile cannot be retrived via API.')
@click.option('--skip-normalizers', is_flag=True, help='Do not normalize.')
@click.option('--not-strict', is_flag=True, help='Also match artificial parsers.')
def local(calc_id, show_backend, show_metadata, skip_normalizers, not_strict, **kwargs):
    utils.get_logger(__name__).info('Using %s' % config.client.url)

    with CalcProcReproduction(calc_id, **kwargs) as local:
        if local.upload_id != 'unknown':
            print(
                'Data being saved to .volumes/fs/tmp/repro_'
                '%s if not already there' % local.upload_id)
        backend = local.parse(strict=not not_strict)

        if not skip_normalizers:
            local.normalize_all(parser_backend=backend)

        if show_backend:
            json.dump(backend.resource.m_to_dict(), sys.stdout, indent=2)

        if show_metadata:
            metadata = datamodel.EntryMetadata(domain='dft')  # TODO take domain from matched parser
            metadata.apply_domain_metadata(backend)
            json.dump(metadata.m_to_dict(), sys.stdout, indent=4)
