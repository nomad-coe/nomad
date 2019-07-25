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

from nomad import utils, processing as proc, search, config, files

from .client import client


__in_test = False
""" Will be monkeypatched by tests to alter behavior for testing. """


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
    '--dry', is_flag=True, default=False,
    help='Do not actually mirror data, just fetch data and report.')
def mirror(query, move: bool, dry: bool):
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

    while True:
        query_results = client.repo.quantity_search(quantity='upload_id', **query).response().result
        upload_ids = query_results.quantities['upload_id']

        for upload_id in upload_ids['values'].keys():
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
                    (upload_data.upload_id, upload_ids['values'][upload_id], upload_data.upload_files_path))
                continue

            # create mongo
            upload = proc.Upload.from_json(upload_data.upload, created=True).save()
            for calc_data in upload_data.calcs:
                proc.Calc.from_json(calc_data, created=True).save()

            # index es
            search.index_all(upload.to_upload_with_metadata().calcs)

            # copy/mv file
            upload_files_path = upload_data.upload_files_path
            if __in_test:
                tmp = os.path.join(config.fs.tmp, 'to_mirror')
                os.rename(upload_files_path, tmp)
                upload_files_path = tmp

            target_upload_files_path = files.PathObject(config.fs.public, upload.upload_id, create_prefix=True, prefix=True).os_path
            if move:
                os.rename(upload_files_path, target_upload_files_path)
                os.symlink(os.path.abspath(target_upload_files_path), upload_files_path)
            else:
                shutil.copytree(upload_files_path, target_upload_files_path)

            print(
                'Mirrored %s with %d calcs at %s' %
                (upload_data.upload_id, upload_ids['values'][upload_id], upload_data.upload_files_path))

        if 'after' not in upload_ids:
            break

        query.update(after=upload_ids['after'])
