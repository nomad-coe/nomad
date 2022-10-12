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


from nomad.cli.cli import cli
from nomad import config


@cli.group(help='Commands that use the nomad API to do useful things')
@click.option('-n', '--url', default=config.client.url, help='The URL where nomad is running, default is "%s".' % config.client.url)
@click.option('-u', '--user', default=None, help='the user name to login, default is no login.')
@click.option('-w', '--password', default=config.client.password, help='the password used to login.')
@click.option('--token-via-api', is_flag=True, help='retrieve the access token from the api not keycloak.')
@click.option('--no-ssl-verify', help='disables SSL verificaton when talking to nomad.', is_flag=True)
@click.pass_context
def client(ctx, url: str, user: str, password: str, no_ssl_verify: bool, token_via_api: bool):
    config.client.url = url

    ctx.obj.user = user
    ctx.obj.password = password
    ctx.obj.token_via_api = token_via_api


def _create_auth(ctx):
    print(f'Used nomad is {config.client.url}')
    print(f'Used user is {ctx.obj.user}')

    from nomad.client import Auth
    if ctx.obj.user is None:
        return None

    return Auth(user=ctx.obj.user, password=ctx.obj.password, from_api=ctx.obj.token_via_api)


@client.command(help='Runs a few example operations as a test.')
@click.option(
    '--skip-parsers', is_flag=True,
    help='Skip extensive upload and parser tests.')
@click.option(
    '--skip-publish', is_flag=True,
    help='Skip publish the upload. Should not be done on an production environment.')
@click.option(
    '--skip-doi', is_flag=True,
    help='Skip assigning a doi to a dataset.')
@click.pass_context
def integrationtests(ctx, skip_parsers, skip_publish, skip_doi):
    import sys

    if ctx.obj.user is None:
        print('Authorization is required. Please provide credentials via --user, --password.')
        sys.exit(1)

    from .integrationtests import integrationtests
    integrationtests(_create_auth(ctx), skip_parsers, skip_publish, skip_doi)


@client.command(help='Metainfo compatibility tests against data in a NOMAD installation.')
@click.pass_context
def datatests(ctx):
    from .datatests import datatests
    datatests(_create_auth(ctx))


@client.command(
    help='Upload files to nomad. The given path can be a single file or a directory. '
    'All .zip files in a directory will be uploaded.')
@click.argument('PATH', nargs=-1, required=True, type=click.Path(exists=True))
@click.option(
    '--upload-name',
    help='Optional name for the upload of a single file. Will be ignored on directories.')
@click.option(
    '--local-path', is_flag=True, default=False,
    help='Upload files "offline": files will not be uploaded, but processed were they are. '
    'Only works when run on the nomad host.')
@click.option(
    '--publish', is_flag=True, default=False,
    help='Automatically move upload out of the staging area after successful processing')
@click.pass_context
def upload(ctx, path, upload_name: str, local_path: bool, publish: bool):
    import os
    import sys

    from nomad.client import upload_file as client_upload_file

    auth = _create_auth(ctx)
    paths = path

    def upload_file(path, upload_name):
        result = client_upload_file(
            path, auth, upload_name=upload_name, local_path=local_path, publish=publish)
        if result is None:
            sys.exit(1)

    click.echo('uploading files from %s paths' % len(paths))
    for path in paths:
        click.echo('uploading %s' % path)
        if os.path.isfile(path):
            upload_name = upload_name if upload_name is not None else os.path.basename(path)
            upload_file(path, upload_name)

        elif os.path.isdir(path):
            for (dirpath, _, filenames) in os.walk(path):
                for filename in filenames:
                    if filename.endswith('.zip'):
                        file_path = os.path.abspath(os.path.join(dirpath, filename))
                        upload_name = os.path.basename(file_path)
                        upload_file(file_path, upload_name)

        else:
            click.echo('Unknown path type %s.' % path)


@client.command(help='Run processing locally.')
@click.argument('ENTRY_ID', nargs=1, required=True, type=str)
@click.option('--override', is_flag=True, help='Override existing local entry data.')
@click.option('--show-archive', is_flag=True, help='Print the archive data.')
@click.option('--show-metadata', is_flag=True, help='Print the extracted repo metadata.')
@click.option('--skip-normalizers', is_flag=True, help='Do not normalize.')
@click.option('--not-strict', is_flag=True, help='Also match artificial parsers.')
@click.pass_context
def local(ctx, entry_id, show_archive, show_metadata, skip_normalizers, not_strict, **kwargs):
    import sys
    import json

    from nomad.client import LocalEntryProcessing

    print('Using %s' % config.client.url)
    auth = _create_auth(ctx)

    with LocalEntryProcessing(entry_id, auth=auth, **kwargs) as local:
        print(f'Data being saved to .volumes/fs/tmp/repro_{entry_id} if not already there')
        entry_archives = local.parse(strict=not not_strict)

        for entry_archive in entry_archives:
            if not skip_normalizers:
                local.normalize_all(entry_archive=entry_archive)

            if show_archive:
                json.dump(entry_archive.m_to_dict(), sys.stdout, indent=2)

            if show_metadata:
                metadata = entry_archive.metadata
                metadata.apply_archive_metadata(entry_archive)
                json.dump(metadata.m_to_dict(), sys.stdout, indent=4)


@client.command(
    help='Synchronizes the NOMAD database with the given external database.')
@click.argument('db_name', nargs=1, required=True)
@click.argument('root_url', nargs=1, required=True)
@click.option(
    '--outfile', default=None,
    help='File to read/write files missing in NOMAD database')
@click.option(
    '--nomadfile', default=None,
    help='File to read/write files in NOMAD database')
@click.option(
    '--dbfile', default=None,
    help='File to read/write files in given database')
@click.option(
    '--local_path', default='/nomad/fairdi/external',
    help='Directory to which the files will be downloaded')
@click.option(
    '--parallel', default=2,
    help='Number of processes to spawn to download/upload files')
@click.option(
    '--do-download', is_flag=True, default=False,
    help='Flag to automatically download downloaded files')
@click.option(
    '--do-upload', is_flag=True, default=False,
    help='Flag to automatically upload downloaded files')
@click.option(
    '--do-publish', is_flag=True, default=False,
    help='Flag to automatically publish upload')
@click.option(
    '--cleanup', is_flag=True, default=False,
    help='Flag to clean up downloaded files')
@click.pass_context
def synchdb(ctx, **kwargs):
    from nomad.cli.aflow import DbUpdater

    db = DbUpdater(auth=_create_auth(ctx), **kwargs)
    db.update()
