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

import typing
import os
import json
import click
import sys

from nomad import utils, parsing, normalizing, datamodel

from .cli import cli


def parse(
        mainfile_path: str,
        parser_name: str = None,
        backend_factory: typing.Callable = None,
        strict: bool = True, logger=None):
    '''
    Run the given parser on the downloaded calculation. If no parser is given,
    do parser matching and use the respective parser.
    '''
    from nomad.parsing import parsers
    mainfile = os.path.basename(mainfile_path)

    if logger is None:
        logger = utils.get_logger(__name__)
    if parser_name is not None:
        parser = parsers.parser_dict.get(parser_name)
        assert parser is not None, 'the given parser must exist'
    else:
        parser = parsers.match_parser(mainfile_path, strict=strict)
        if isinstance(parser, parsing.MatchingParser):
            parser_name = parser.name
        else:
            parser_name = parser.__class__.__name__

    assert parser is not None, 'there is no parser matching %s' % mainfile
    logger = logger.bind(parser=parser.name)  # type: ignore
    logger.info('identified parser')
    if hasattr(parser, 'backend_factory'):
        setattr(parser, 'backend_factory', backend_factory)

    entry_archive = datamodel.EntryArchive()
    metadata = entry_archive.m_create(datamodel.EntryMetadata)
    try:
        cwd = os.getcwd()
        mainfile_path = os.path.abspath(mainfile_path)
        os.chdir(os.path.abspath(os.path.dirname(mainfile_path)))
        parser.parse(mainfile_path, entry_archive, logger=logger)
        os.chdir(cwd)
    except Exception as e:
        logger.error('parsing was not successful', exc_info=e)

    if metadata.domain is None:
        metadata.domain = parser.domain

    logger.info('ran parser')
    return entry_archive


def normalize(
        normalizer: typing.Union[str, typing.Callable], entry_archive, logger=None):

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
    Parse the downloaded calculation and run the whole normalizer chain.
    '''
    for normalizer in normalizing.normalizers:
        if normalizer.domain == entry_archive.section_metadata.domain:
            normalize(normalizer, entry_archive, logger=logger)


@cli.command(help='Run parsing and normalizing locally.', name='parse')
@click.argument('MAINFILE', nargs=1, required=True, type=str)
@click.option('--show-archive', is_flag=True, default=False, help='Print the archive data.')
@click.option('--show-metadata', is_flag=True, default=False, help='Print the extracted repo metadata.')
@click.option('--skip-normalizers', is_flag=True, default=False, help='Do not run the normalizer.')
@click.option('--not-strict', is_flag=True, help='Do also match artificial parsers.')
@click.option('--parser', help='Skip matching and use the provided parser')
def _parse(mainfile, show_archive, show_metadata, skip_normalizers, not_strict, parser):
    kwargs = dict(strict=not not_strict, parser_name=parser)

    entry_archive = parse(mainfile, **kwargs)

    if not skip_normalizers:
        normalize_all(entry_archive)
        entry_archive.section_metadata.apply_domain_metadata(entry_archive)

    if show_archive:
        json.dump(entry_archive.m_to_dict(), sys.stdout, indent=2)

    if show_metadata:
        metadata = entry_archive.section_metadata
        metadata.apply_domain_metadata(entry_archive)
        json.dump(metadata.m_to_dict(), sys.stdout, indent=4)
