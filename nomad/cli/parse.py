import typing
import os
import json
import click
import sys

from nomad import utils, parsing, normalizing, datamodel

import nomadcore

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

    parser_backend = parser.run(mainfile_path, logger=logger)

    if not parser_backend.status[0] == 'ParseSuccess':
        logger.error('parsing was not successful', status=parser_backend.status)

    logger.info('ran parser')
    return parser_backend


def normalize(
        normalizer: typing.Union[str, typing.Callable], parser_backend=None,
        logger=None):

    if logger is None:
        logger = utils.get_logger(__name__)

    if isinstance(normalizer, str):
        normalizer = next(
            normalizer_instance for normalizer_instance in normalizing.normalizers
            if normalizer_instance.__class__.__name__ == normalizer)

    assert normalizer is not None, 'there is no normalizer %s' % str(normalizer)
    normalizer_instance = typing.cast(typing.Callable, normalizer)(parser_backend.entry_archive)
    logger = logger.bind(normalizer=normalizer_instance.__class__.__name__)
    logger.info('identified normalizer')

    normalizer_instance.normalize(logger=logger)
    logger.info('ran normalizer')
    return parser_backend


def normalize_all(parser_backend=None, logger=None):
    '''
    Parse the downloaded calculation and run the whole normalizer chain.
    '''
    for normalizer in normalizing.normalizers:
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
@click.option('--annotate', is_flag=True, help='Sub-matcher based parsers will create a .annotate file.')
def _parse(
        mainfile, show_backend, show_metadata, skip_normalizers, not_strict, parser,
        annotate):
    nomadcore.simple_parser.annotate = annotate
    kwargs = dict(strict=not not_strict, parser_name=parser)

    backend = parse(mainfile, **kwargs)

    if not skip_normalizers:
        normalize_all(backend)

    if show_backend:
        json.dump(backend.resource.m_to_dict(), sys.stdout, indent=2)

    if show_metadata:
        metadata = datamodel.EntryMetadata(domain='dft')  # TODO take domain from matched parser
        metadata.apply_domain_metadata(backend)
        json.dump(metadata.m_to_dict(), sys.stdout, indent=4)
