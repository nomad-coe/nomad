from typing import Union, Callable, cast
import os.path
import ujson
import click
import sys

from nomad import config, utils, files
from nomad.parsing import LocalBackend, parser_dict, match_parser, MatchingParser
from nomad.normalizing import normalizers
from nomad.datamodel import CalcWithMetadata

from .cli import cli


def parse(
        mainfile: str, upload_files: Union[str, files.StagingUploadFiles],
        parser_name: str = None, strict: bool = True, logger=None) -> LocalBackend:
    """
    Run the given parser on the downloaded calculation. If no parser is given,
    do parser matching and use the respective parser.
    """
    if logger is None:
        logger = utils.get_logger(__name__)
    if parser_name is not None:
        parser = parser_dict.get(parser_name)
    else:
        parser = match_parser(mainfile, upload_files, strict=strict)
        if isinstance(parser, MatchingParser):
            parser_name = parser.name
        else:
            parser_name = parser.__class__.__name__

    assert parser is not None, 'there is not parser matching %s' % mainfile
    logger = logger.bind(parser=parser.name)  # type: ignore
    logger.info('identified parser')

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
    """
    Parse the downloaded calculation and run the whole normalizer chain.
    """
    for normalizer in normalizers:
        parser_backend = normalize(normalizer, parser_backend=parser_backend, logger=logger)

    return parser_backend


@cli.command(help='Run parsing and normalizing locally.', name='parse')
@click.argument('MAINFILE', nargs=1, required=True, type=str)
@click.option('--show-backend', is_flag=True, default=False, help='Print the backend data.')
@click.option('--show-metadata', is_flag=True, default=False, help='Print the extracted repo metadata.')
@click.option('--skip-normalizers', is_flag=True, default=False, help='Do not run the normalizer.')
@click.option('--not-strict', is_flag=True, help='Do also match artificial parsers.')
def _parse(mainfile, show_backend, show_metadata, skip_normalizers, not_strict):
    utils.configure_logging()

    backend = parse(mainfile, '.', strict=not not_strict)

    if not skip_normalizers:
        normalize_all(backend)

    if show_backend:
        backend.write_json(sys.stdout, pretty=True)
    if show_metadata:
        metadata = CalcWithMetadata()
        metadata.apply_domain_metadata(backend)
        ujson.dump(metadata.to_dict(), sys.stdout, indent=4)
