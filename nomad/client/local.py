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
import io
import requests
import click
from typing import Union, Callable, cast
import sys
import ujson

from nomad import config, utils
from nomad.files import ArchiveBasedStagingUploadFiles
from nomad.parsing import parser_dict, LocalBackend, match_parser
from nomad.normalizing import normalizers

from .main import cli, get_nomad_url


class CalcProcReproduction:
    """
    Instances represent a local reproduction of the processing for a single calculation.
    It allows to download raw data from a nomad server and reproduce its processing
    (parsing, normalizing) with the locally installed parsers and normalizers.

    The use-case is error/warning reproduction. Use ELK to identify errors, use
    the upload, archive ids to given by ELK, and reproduce and fix the error
    in your development environment.

    Arguments:
        archive_id: The archive_id of the calculation to locally process.
        override: Set to true to override any existing local calculation data.
    """
    def __init__(self, archive_id: str, override: bool = False, mainfile: str = None) -> None:
        self.calc_id = utils.archive.calc_id(archive_id)
        self.upload_id = utils.archive.upload_id(archive_id)
        self.mainfile = mainfile
        self.parser = None
        self.logger = utils.get_logger(__name__, archive_id=archive_id)

        from .main import create_client
        client = create_client()
        if self.mainfile is None:
            calc = client.repo.get_repo_calc(upload_id=self.upload_id, calc_id=self.calc_id).response().result
            self.mainfile = calc['mainfile']
        else:
            self.logger.info('Using provided mainfile', mainfile=self.mainfile)

        local_path = os.path.join(config.fs.tmp, 'repro_%s.zip' % archive_id)
        if not os.path.exists(os.path.dirname(local_path)):
            os.makedirs(os.path.dirname(local_path))
        if not os.path.exists(local_path) or override:
            # download raw if not already downloaded or if override is set
            # download with request, since bravado does not support streaming
            # TODO currently only downloads mainfile
            self.logger.info('Downloading calc.', mainfile=self.mainfile)
            token = client.auth.get_user().response().result.token
            req = requests.get('%s/raw/%s/%s' % (get_nomad_url(), self.upload_id, os.path.dirname(self.mainfile)) + '/*', stream=True, headers={'X-Token': token})
            with open(local_path, 'wb') as f:
                for chunk in req.iter_content(chunk_size=io.DEFAULT_BUFFER_SIZE):
                    f.write(chunk)
        else:
            self.logger.info('Calc already downloaded.')

        self.upload_files = ArchiveBasedStagingUploadFiles(upload_id='tmp_%s' % archive_id, local_path=local_path, create=True, is_authorized=lambda: True)

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

    def parse(self, parser_name: str = None) -> LocalBackend:
        """
        Run the given parser on the downloaded calculation. If no parser is given,
        do parser matching and use the respective parser.
        """
        if parser_name is not None:
            parser = parser_dict.get(parser_name)
        else:
            parser = match_parser(self.mainfile, self.upload_files)

        assert parser is not None, 'there is not parser matching %s' % self.mainfile
        self.logger = self.logger.bind(parser=parser.name)  # type: ignore
        self.logger.info('identified parser')

        parser_backend = parser.run(self.upload_files.raw_file_object(self.mainfile).os_path, logger=self.logger)

        if not parser_backend.status[0] == 'ParseSuccess':
            self.logger.error('parsing was not successful', status=parser_backend.status)

        parser_backend.openNonOverlappingSection('section_calculation_info')
        parser_backend.addValue('upload_id', self.upload_id)
        parser_backend.addValue('calc_id', self.calc_id)
        parser_backend.addValue('calc_hash', "no hash")
        parser_backend.addValue('main_file', self.mainfile)
        parser_backend.addValue('parser_name', parser.__class__.__name__)
        parser_backend.closeNonOverlappingSection('section_calculation_info')

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


@cli.command(help='Run processing locally.')
@click.argument('ARCHIVE_ID', nargs=1, required=True, type=str)
@click.option('--override', is_flag=True, default=False, help='Override existing local calculation data.')
@click.option('--show-backend', is_flag=True, default=False, help='Print the backend data.')
@click.option('--show-metadata', is_flag=True, default=False, help='Print the extracted repo metadata.')
@click.option('--mainfile', default=None, type=str, help='Use this mainfile (in case mainfile cannot be retrived via API.')
def local(archive_id, show_backend=False, show_metadata=False, **kwargs):
    print(kwargs)
    utils.configure_logging()
    utils.get_logger(__name__).info('Using %s' % get_nomad_url())
    with CalcProcReproduction(archive_id, **kwargs) as local:
        backend = local.parse()
        local.normalize_all(parser_backend=backend)
        if show_backend:
            backend.write_json(sys.stdout, pretty=True)
        if show_metadata:
            metadata = backend.to_calc_with_metadata()
            ujson.dump(metadata.to_dict(), sys.stdout, indent=4)
