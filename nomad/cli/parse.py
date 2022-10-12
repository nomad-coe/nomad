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

import click

from .cli import cli


@cli.command(help='Run parsing and normalizing locally.', name='parse')
@click.argument('MAINFILE', nargs=1, required=True, type=str)
@click.option('--show-archive', is_flag=True, default=False, help='Print the archive data.')
@click.option('--archive-with-meta', is_flag=True, default=False, help='If the archvie is printed, it would include metadata like m_def and m_annotations.')
@click.option('--show-metadata', is_flag=True, default=False, help='Print the extracted repo metadata.')
@click.option('--skip-normalizers', is_flag=True, default=False, help='Do not run the normalizer.')
@click.option('--not-strict', is_flag=True, help='Do also match artificial parsers.')
@click.option('--parser', help='Skip matching and use the provided parser')
@click.option('--server-context', is_flag=True, default=False, help='Whether to use server context.')
@click.option('--username', default=None, help='Username for authentication.')
@click.option('--password', default=None, help='Password for authentication.')
def _parse(
        mainfile, show_archive, archive_with_meta, show_metadata, skip_normalizers, not_strict, parser, server_context,
        username, password):
    import sys
    import json

    from nomad.client import parse, normalize_all

    kwargs = dict(
        strict=not not_strict,
        parser_name=parser,
        server_context=server_context,
        username=username,
        password=password)

    entry_archives = parse(mainfile, **kwargs)

    for entry_archive in entry_archives:
        if not skip_normalizers:
            normalize_all(entry_archive)
            entry_archive.metadata.apply_archive_metadata(entry_archive)

        if show_archive:
            json.dump(entry_archive.m_to_dict(with_meta=archive_with_meta), sys.stdout, indent=2)

        if show_metadata:
            metadata = entry_archive.metadata
            metadata.apply_archive_metadata(entry_archive)
            json.dump(metadata.m_to_dict(), sys.stdout, indent=4)
