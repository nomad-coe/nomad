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
import io
import typing
import sys

from nomad import config, utils, datamodel

from .api import Auth


def parse(
        mainfile_path: str,
        parser_name: str = None,
        backend_factory: typing.Callable = None,
        strict: bool = True,
        logger=None,
        server_context: bool = False,
        username: str = None,
        password: str = None) -> typing.List[datamodel.EntryArchive]:
    '''
    Run the given parser on the provided mainfile. If parser_name is given, we only try
    to match this parser, otherwise we try to match all parsers.
    '''
    from nomad import parsing
    from nomad.parsing import parsers
    mainfile = os.path.basename(mainfile_path)

    if logger is None:
        logger = utils.get_logger(__name__)

    mainfile_path = os.path.abspath(mainfile_path)
    parser, mainfile_keys = parsers.match_parser(mainfile_path, strict=strict, parser_name=parser_name)
    if isinstance(parser, parsing.MatchingParser):
        parser_name = parser.name
    else:
        parser_name = parser.__class__.__name__

    assert parser is not None, 'there is no parser matching %s' % mainfile
    logger = logger.bind(parser=parser.name)  # type: ignore
    logger.info('identified parser')
    if hasattr(parser, 'backend_factory'):
        setattr(parser, 'backend_factory', backend_factory)

    entry_archives = parsers.run_parser(
        mainfile_path, parser, mainfile_keys, logger, server_context, username, password)

    logger.info('ran parser')
    return entry_archives


def normalize(
        normalizer: typing.Union[str, typing.Callable], entry_archive, logger=None):
    from nomad import normalizing

    if logger is None:
        logger = utils.get_logger(__name__)

    if isinstance(normalizer, str):
        normalizer = next(
            normalizer_instance for normalizer_instance in normalizing.normalizers
            if normalizer_instance.__class__.__name__ == normalizer)

    assert normalizer is not None, 'there is no normalizer %s' % str(normalizer)
    normalizer_instance = typing.cast(typing.Callable, normalizer)(entry_archive)
    logger = logger.bind(normalizer=normalizer_instance.__class__.__name__)
    logger.info('identified normalizer')

    normalizer_instance.normalize(logger=logger)
    logger.info('ran normalizer')


def normalize_all(entry_archive, logger=None):
    '''
    Parse the downloaded entry and run the whole normalizer chain.
    '''
    from nomad import normalizing

    for normalizer in normalizing.normalizers:
        if normalizer.domain is None or normalizer.domain == entry_archive.metadata.domain:
            normalize(normalizer, entry_archive, logger=logger)


class LocalEntryProcessing:
    '''
    Instances represent a local reproduction of the processing for a single entry.
    It allows to download raw data from a nomad server and reproduce its processing
    (parsing, normalizing) with the locally installed parsers and normalizers.

    The use-case is error/warning reproduction. Use ELK to identify errors, use
    the upload, archive ids to given by ELK, and reproduce and fix the error
    in your development environment.

    Arguments:
        entry_id: The entry_id of the entry to locally process.
        override: Set to true to override any existing local entry data.
        auth: Optional Auth object to download private data.
    '''
    def __init__(self, entry_id: str, override: bool = False, auth: Auth = None) -> None:
        from nomad import files
        from nomad.client import api

        self.entry_id = entry_id
        response = self.__handle_response(
            api.get(f'entries/{self.entry_id}', auth=auth))
        self.mainfile = response.json()['data']['mainfile']
        self.upload_id = response.json()['data']['upload_id']

        self.local_path = os.path.join(config.fs.tmp, f'repro_{self.entry_id}.zip')
        if not os.path.exists(os.path.dirname(self.local_path)):
            os.makedirs(os.path.dirname(self.local_path))

        if not os.path.exists(self.local_path) or override:
            # download raw if not already downloaded or if override is set
            print('Downloading', self.entry_id)
            response = self.__handle_response(
                api.get(f'entries/{self.entry_id}/raw', auth=auth))

            iter_content = getattr(response, 'iter_content', None)
            with open(self.local_path, 'wb') as f:
                if iter_content:
                    for chunk in iter_content(chunk_size=io.DEFAULT_BUFFER_SIZE):
                        f.write(chunk)
                else:
                    # Fallback for clients that don't support iterating the content
                    f.write(response.content)
        else:
            print('Entry already downloaded.')

        self.upload_files = files.StagingUploadFiles(upload_id=f'tmp_{self.entry_id}', create=True)

    def __handle_response(self, response):
        if response.status_code == 404:
            print('Could not find entry, maybe it was deleted or the ids are wrong:', file=sys.stderr)
            print(response.text, file=sys.stderr)
            sys.exit(1)
        elif response.status_code != 200:
            print(f'API error: {response.text}', file=sys.stderr)
            sys.exit(1)

        return response

    def __enter__(self):
        # open/extract upload file
        print('Extracting entry data.')
        if self.upload_files.is_empty():  # Only add the files once
            self.upload_files.add_rawfiles(self.local_path)

        for raw_file in self.upload_files.raw_directory_list(recursive=True):
            if raw_file.path == f'{self.upload_id}/{self.mainfile}':
                print('Identified mainfile.')
                return self

        print('Count not find file in download.', file=sys.stderr)
        sys.exit(1)

    def __exit__(self, exception_type, exception, trace):
        if exception:
            import traceback
            traceback.print_exc()
        try:
            self.upload_files.delete()
        except Exception:
            pass
        if exception:
            sys.exit(1)

    def parse(self, parser_name: str = None, **kwargs) -> typing.List[datamodel.EntryArchive]:
        '''
        Run the given parser on the downloaded entry. If no parser is given,
        do parser matching and use the respective parser.
        '''
        return parse(
            self.upload_files.raw_file_object(f'{self.upload_id}/{self.mainfile}').os_path,
            parser_name=parser_name, logger=utils.get_logger(__name__), **kwargs)

    def normalize(self, normalizer: typing.Union[str, typing.Callable], entry_archive=None):
        '''
        Parse the downloaded entry and run the given normalizer.
        '''
        if entry_archive is None:
            entry_archive = self.parse()

        return normalize(
            entry_archive=entry_archive, normalizer=normalizer, logger=utils.get_logger(__name__))

    def normalize_all(self, entry_archive=None):
        '''
        Parse the downloaded entry and run the whole normalizer chain.
        '''
        return normalize_all(
            entry_archive=entry_archive, logger=utils.get_logger(__name__))
