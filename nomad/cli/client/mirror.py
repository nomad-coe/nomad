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

import click
import sys
import shutil
import json
import os
import os.path

from nomad import utils, processing as proc, search, config, files, infrastructure

from .client import client


__in_test = False
""" Will be monkeypatched by tests to alter behavior for testing. """


class Mapping:

    def __init__(self, mapping: str):
        if mapping is None:
            self.source_len = 0
            self.target = ''
        else:
            mapping_split = mapping.split(':')
            assert len(mapping_split) == 2, 'Mapping format is dir:mapping'

            source, self.target = mapping_split
            self.source_len = len(source)

    def apply(self, path):
        return self.target + path[self.source_len:]


@client.command(
    help='''
        Mirror data from another nomad deployment.

        It uses its 'client' nature to export data from the other nomad (source) via
        REST API. I.e., use the -n (--url) client parameter to specify the source deployment.

        The data will be added directly to the databases, filesystem, etc. of
        'this' nomad deployment (target), i.e. it be haves like an 'admin' command.
        This means you either run it in the environement of the target deployment
        or use the --config nomad parameter.''')
@click.argument('QUERY', nargs=1, required=False)
@click.option(
    '--move', is_flag=True, default=False,
    help='Instead of copying the underlying upload files, we move it and replace it with a symlink.')
@click.option(
    '--source-mapping', type=str, default=None,
    help=(
        'A mapping in the form "dir:mapped" that replaces path prefix "dir" with '
        '"mapped" in all paths provided by the API of the source. Allows to handle local mounts '
        'paths in source deployment. E.g. use ".volumes/fs:/nomad/fairdi/<source>/fs".'))
@click.option(
    '--target-mapping', type=str, default=None,
    help=(
        'A mapping in the form "dir:mapped" that replaces path prefix "mirror" with '
        '"mapped" in all paths used for the target. Allows to handle local mounts '
        'paths in target deployment. E.g. use ".volumes/fs:/nomad/fairdi/<target>/fs".'))
@click.option(
    '--dry', is_flag=True, default=False,
    help='Do not actually mirror data, just fetch data and report.')
def mirror(query, move: bool, dry: bool, source_mapping, target_mapping):
    infrastructure.setup_mongo()
    infrastructure.setup_elastic()

    source_mapping = Mapping(source_mapping)
    target_mapping = Mapping(target_mapping)

    if query is not None:
        try:
            query = json.loads(query)
        except Exception as e:
            print('Cannot parse the given query %s: %s' % (query, str(e)))
            sys.exit(1)
    else:
        query = {}

    query.update(owner='admin')

    utils.configure_logging()

    from nomad.cli.client import create_client
    client = create_client()

    uploads = client.mirror.get_uploads_mirror(payload=dict(query={})).response().result

    for upload_with_id in uploads:
        upload_id = upload_with_id['upload_id']
        upload_data = client.mirror.get_upload_mirror(upload_id=upload_id).response().result

        try:
            upload = proc.Upload.get(upload_id)
            if __in_test:
                proc.Calc.objects(upload_id=upload_id).delete()
                proc.Upload.objects(upload_id=upload_id).delete()
                search.delete_upload(upload_id)

                raise KeyError()

            print(
                'Upload %s already exists, updating existing uploads is not implemented yet. '
                'Skip upload.' % upload_id)
            continue
        except KeyError:
            pass

        if dry:
            print(
                'Need to mirror %s with %d calcs at %s' %
                (upload_data.upload_id, len(upload_data.calcs), upload_data.upload_files_path))
            continue

        # copy/mv file
        upload_files_path = upload_data.upload_files_path
        if __in_test:
            tmp = os.path.join(config.fs.tmp, 'to_mirror')
            os.rename(upload_files_path, tmp)
            upload_files_path = tmp

        upload_files_path = source_mapping.apply(upload_files_path)

        target_upload_files_path = files.PathObject(
            config.fs.public, upload_id, create_prefix=True, prefix=True).os_path
        target_upload_files_path = target_mapping.apply(target_upload_files_path)
        if not os.path.exists(target_upload_files_path):
            os.makedirs(target_upload_files_path)

        if move:
            os.rename(upload_files_path, target_upload_files_path)
            os.symlink(os.path.abspath(target_upload_files_path), upload_files_path)
        else:
            for to_copy in os.listdir(upload_files_path):
                shutil.copyfile(
                    os.path.join(upload_files_path, to_copy),
                    os.path.join(target_upload_files_path, to_copy))

        # create mongo
        upload = proc.Upload.from_json(upload_data.upload, created=True).save()
        for calc_data in upload_data.calcs:
            proc.Calc.from_json(calc_data, created=True).save()

        # index es
        search.index_all(upload.to_upload_with_metadata().calcs)

        print(
            'Mirrored %s with %d calcs at %s' %
            (upload_data.upload_id, len(upload_data.calcs), upload_data.upload_files_path))
