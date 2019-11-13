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
from bravado.exception import HTTPBadRequest

from nomad import utils, processing as proc, search, config, files, infrastructure
from nomad.datamodel import Dataset, User

from .client import client


__in_test = False
""" Will be monkeypatched by tests to alter behavior for testing. """

__Dataset = Dataset.m_def.m_x('me')
__logger = utils.get_logger(__name__)


def v0Dot6(upload_data):
    """ Inplace transforms v0.6.x upload data into v0.7.x upload data. """

    def tarnsform_user_id(source_user_id):
        target_user = User.repo_users.get(source_user_id)
        if target_user is None:
            __logger.error('user does not exist in target', source_user_id=source_user_id)
            raise KeyError

        return target_user.user_id

    def transform_dataset(source_dataset):
        legacy_id = source_dataset['id']
        target_dataset = __Dataset.objects(legacy_id=legacy_id).first()
        if target_dataset is not None:
            return target_dataset.dataset_id

        target_dataset = __Dataset(
            dataset_id=utils.create_uuid(),
            legacy_id=source_dataset['id'],
            name=source_dataset['name'])
        if 'doi' in source_dataset and len(source_dataset['doi']) > 0:
            source_doi = source_dataset['doi'][0]
            target_dataset.doi = source_doi.replace('http://dx.doi.org/', '')
        target_dataset.save()

        return target_dataset.dataset_id

    upload = json.loads(upload_data.upload)
    upload['user_id'] = tarnsform_user_id(upload['user_id'])
    upload_data.upload = json.dumps(upload)

    for calc_data_json, i in enumerate(upload_data.calcs):
        calc_data = json.loads(calc_data_json)
        metadata = calc_data['metadata']

        # transform users
        metadata['uploader'] = tarnsform_user_id(metadata['uploader']['id'])
        metadata['coauthors'] = [tarnsform_user_id(user['id']) for user in metadata['coauthors']]
        metadata['shared_with'] = [tarnsform_user_id(user['id']) for user in metadata['shared_with']]

        # transform datasets
        metadata['datasets'] = [transform_dataset(dataset) for dataset in metadata['datasets']]

        upload_data.calcs[i] = json.dumps(calc_data)
    return upload_data


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
    '--link', is_flag=True, default=False,
    help='Instead of copying the underlying upload files, we create symlinks in the target.'
)
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
@click.option(
    '--files-only', is_flag=True, default=False,
    help=(
        'Will only copy/move files and not even look at the calculations. Useful, '
        'when moving metadata via mongo dump/restore.'))
@click.option(
    '--migration', type=str, default=None,
    help='The name of a migration script used to transform the metadata.'
)
def mirror(
        query, move: bool, link: bool, dry: bool, files_only: bool, source_mapping: str,
        target_mapping: str, migration: str):

    migration_func = None
    if migration is not None:
        if migration == 'v0.6.x':
            migration_func = v0Dot6
        else:
            print('Migration %s does not exist.' % migration)
            sys.exit(1)

    infrastructure.setup_mongo()
    infrastructure.setup_elastic()

    source_mapping_obj = Mapping(source_mapping)
    target_mapping_obj = Mapping(target_mapping)

    if query is not None:
        try:
            query = json.loads(query)
        except Exception as e:
            print('Cannot parse the given query %s: %s' % (query, str(e)))
            sys.exit(1)
    else:
        query = dict(published=True)

    utils.configure_logging()

    from nomad.cli.client import create_client
    client = create_client()

    uploads = client.mirror.get_uploads_mirror(payload=dict(query=query)).response().result

    for upload_data in uploads:
        upload_id = upload_data.upload_id

        if not files_only:
            try:
                upload = proc.Upload.get(upload_id)
                if __in_test:
                    # In tests, we mirror from our selves, fake that the upload does not exist
                    raise KeyError()

                if len(query) > 0:
                    print(
                        'Upload %s already exists, updating existing uploads is not implemented yet. '
                        'Skip upload.' % upload_id)
                continue
            except KeyError:
                pass

            try:
                upload_data = client.mirror.get_upload_mirror(upload_id=upload_id).response().result
                n_calcs = len(upload_data.calcs)
            except HTTPBadRequest:
                print('Could not mirror %s, it is probably not published.' % upload_id)
                n_calcs = 0
                continue

            if __in_test:
                # In tests, we mirror from our selves, remove it so it is not there for import
                proc.Calc.objects(upload_id=upload_id).delete()
                proc.Upload.objects(upload_id=upload_id).delete()
                search.delete_upload(upload_id)
        else:
            n_calcs = 0

        if dry:
            print(
                'Need to mirror %s with %d calcs at %s' %
                (upload_id, n_calcs, upload_data.upload_files_path))
            continue

        if not files_only:
            # migrate
            if migration_func is not None:
                try:
                    upload_data = migration_func(upload_data)
                except Exception as e:
                    __logger.error('could not migrate upload_data', exc_info=e)
                    continue

        # copy/link/mv file
        upload_files_path = upload_data.upload_files_path
        if __in_test:
            tmp = os.path.join(config.fs.tmp, 'to_mirror')
            os.rename(upload_files_path, tmp)
            upload_files_path = tmp

        upload_files_path = source_mapping_obj.apply(upload_files_path)

        target_upload_files_path = files.PathObject(
            config.fs.public, upload_id, create_prefix=False, prefix=True).os_path
        target_upload_files_path = target_mapping_obj.apply(target_upload_files_path)

        if not os.path.exists(target_upload_files_path):
            if move:
                os.rename(upload_files_path, target_upload_files_path)
                os.symlink(os.path.abspath(target_upload_files_path), upload_files_path)

            elif link:
                os.symlink(os.path.abspath(upload_files_path), target_upload_files_path)

            else:
                os.makedirs(target_upload_files_path)
                for to_copy in os.listdir(upload_files_path):
                    shutil.copyfile(
                        os.path.join(upload_files_path, to_copy),
                        os.path.join(target_upload_files_path, to_copy))

        if not files_only:
            # create mongo
            upload = proc.Upload.from_json(upload_data.upload, created=True).save()
            for calc_data in upload_data.calcs:
                proc.Calc.from_json(calc_data, created=True).save()

            # index es
            search.index_all(upload.to_upload_with_metadata().calcs)

        print(
            'Mirrored %s with %d calcs at %s' %
            (upload_id, n_calcs, upload_data.upload_files_path))
