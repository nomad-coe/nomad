from typing import Union, Callable, cast
import os.path
import json
import click
import sys

from nomad import config, utils, files
from nomad.parsing import LocalBackend, parser_dict, match_parser, MatchingParser, MetainfoBackend
from nomad.metainfo.legacy import convert
from nomad.normalizing import normalizers
from nomad.datamodel import EntryMetadata

from nomadcore import simple_parser

from .cli import cli


def parse(
        mainfile: str, upload_files: Union[str, files.StagingUploadFiles],
        parser_name: str = None,
        backend_factory: Callable = None,
        strict: bool = True, logger=None) -> LocalBackend:
    '''
    Run the given parser on the downloaded calculation. If no parser is given,
    do parser matching and use the respective parser.
    '''
    if logger is None:
        logger = utils.get_logger(__name__)
    if parser_name is not None:
        parser = parser_dict.get(parser_name)
        assert parser is not None, 'the given parser must exist'
    else:
        parser = match_parser(mainfile, upload_files, strict=strict)
        if isinstance(parser, MatchingParser):
            parser_name = parser.name
        else:
            parser_name = parser.__class__.__name__

    assert parser is not None, 'there is no parser matching %s' % mainfile
    logger = logger.bind(parser=parser.name)  # type: ignore
    logger.info('identified parser')
    if hasattr(parser, 'backend_factory'):
        setattr(parser, 'backend_factory', backend_factory)

    if isinstance(upload_files, str):
        mainfile_path = os.path.join(upload_files, mainfile)
    else:
        mainfile_path = upload_files.raw_file_object(mainfile).os_path

    parser_backend = parser.run(mainfile_path, logger=logger)

    if not parser_backend.status[0] == 'ParseSuccess':
        logger.error('parsing was not successful', status=parser_backend.status)

    parser_backend.openNonOverlappingSection('section_entry_info')
    parser_backend.addValue('upload_id', config.services.unavailable_value)
    parser_backend.addValue('calc_id', config.services.unavailable_value)
    parser_backend.addValue('calc_hash', "no hash")
    parser_backend.addValue('mainfile', mainfile)
    parser_backend.addValue('parser_name', parser_name)
    parser_backend.closeNonOverlappingSection('section_entry_info')

    logger.info('ran parser')
    return parser_backend


def normalize(
        normalizer: Union[str, Callable], parser_backend: LocalBackend = None,
        logger=None) -> LocalBackend:

    if logger is None:
        logger = utils.get_logger(__name__)

    if isinstance(normalizer, str):
        normalizer = next(
            normalizer_instance for normalizer_instance in normalizers
            if normalizer_instance.__class__.__name__ == normalizer)

    assert normalizer is not None, 'there is no normalizer %s' % str(normalizer)
    normalizer_instance = cast(Callable, normalizer)(parser_backend)
    logger = logger.bind(normalizer=normalizer_instance.__class__.__name__)
    logger.info('identified normalizer')

    normalizer_instance.normalize(logger=logger)
    logger.info('ran normalizer')
    return parser_backend


def normalize_all(parser_backend: LocalBackend = None, logger=None) -> LocalBackend:
    '''
    Parse the downloaded calculation and run the whole normalizer chain.
    '''
    for normalizer in normalizers:
        if normalizer.domain == parser_backend.domain:
            parser_backend = normalize(
                normalizer, parser_backend=parser_backend, logger=logger)

    return parser_backend


@cli.command(help='Run parsing and normalizing locally.', name='parse')
@click.argument('MAINFILE', nargs=1, required=True, type=str)
@click.option('--show-backend', is_flag=True, default=False, help='Print the backend data.')
@click.option('--show-metadata', is_flag=True, default=False, help='Print the extracted repo metadata.')
@click.option('--skip-normalizers', is_flag=True, default=False, help='Do not run the normalizer.')
@click.option('--not-strict', is_flag=True, help='Do also match artificial parsers.')
@click.option('--parser', help='Skip matching and use the provided parser')
@click.option('--metainfo', is_flag=True, help='Use the new metainfo instead of the legacy metainfo.')
@click.option('--annotate', is_flag=True, help='Sub-matcher based parsers will create a .annotate file.')
def _parse(
        mainfile, show_backend, show_metadata, skip_normalizers, not_strict, parser,
        metainfo, annotate):

    simple_parser.annotate = annotate

    utils.configure_logging()
    kwargs = dict(strict=not not_strict, parser_name=parser)

    if metainfo:

        def backend_factory(env, logger):
            # from vaspparser.metainfo import m_env
            # from nomad.metainfo import Section
            # m_env.resolve_definition('section_basis_set_atom_centered', Section)
            # return MetainfoBackend(m_env, logger=logger)

            return MetainfoBackend(convert(env), logger=logger)

        kwargs.update(backend_factory=backend_factory)

    backend = parse(mainfile, '.', **kwargs)

    if not skip_normalizers:
        normalize_all(backend)

    if show_backend:
        backend.write_json(sys.stdout, pretty=True)
    if show_metadata:
        metadata = EntryMetadata(domain='dft')  # TODO take domain from matched parser
        metadata.apply_domain_metadata(backend)
        json.dump(metadata.m_to_dict(), sys.stdout, indent=4)
