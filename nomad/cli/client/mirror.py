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
import bravado.exception
import datetime
import traceback

from nomad import utils, processing as proc, search, config, files, infrastructure
from nomad import datamodel
from nomad import doi as nomad_doi
from nomad.cli.admin import uploads as admin_uploads

from .client import client


__in_test = False
''' Will be monkeypatched by tests to alter behavior for testing. '''


# will be initialized in mirror command func
_Dataset = None
__logger = None


def fix_time(data, keys):
    for key in keys:
        time = data.get(key)
        if isinstance(time, int):
            data[key] = datetime.datetime.utcfromtimestamp(time)


def tarnsform_user_id(source_user_id):
    target_user = datamodel.User.repo_users().get(str(source_user_id))
    if target_user is None:
        __logger.error('user does not exist in target', source_user_id=source_user_id)
        raise KeyError

    return target_user.user_id


def transform_dataset(source_dataset):
    pid = str(source_dataset['id'])
    target_dataset = _Dataset.objects(pid=pid).first()
    if target_dataset is not None:
        return target_dataset.dataset_id

    target_dataset = _Dataset(
        dataset_id=utils.create_uuid(),
        pid=pid,
        name=source_dataset['name'])

    if 'doi' in source_dataset and source_dataset['doi'] is not None:
        source_doi = source_dataset['doi']

        if isinstance(source_doi, dict):
            source_doi = source_doi['value']

        if source_doi is not None:
            target_dataset.doi = source_doi.replace('http://dx.doi.org/', '')

    target_dataset.save()

    return target_dataset.dataset_id


def transform_reference(reference):
    return reference['value']


def v0Dot7(upload_data):
    ''' Inplace transforms v0.7.x upload data into v0.8.x upload data. '''
    __mongo_properties = set(d.name for d in datamodel.MongoMetadata.m_def.definitions)
    for calc in upload_data['calcs']:
        calc_metadata = calc['metadata']
        if 'pid' in calc_metadata:
            calc_metadata['pid'] = str(calc_metadata['pid'])
        metadata = {
            key: value
            for key, value in calc_metadata.items()
            if key in __mongo_properties
        }
        entry_metadata = datamodel.EntryMetadata(**metadata)
        calc['metadata'] = entry_metadata.m_to_dict(
            include_defaults=True,
            categories=[datamodel.MongoMetadata])

    return upload_data


def v0Dot6(upload_data):
    ''' Inplace transforms v0.6.x upload data into v0.7.x upload data. '''
    upload = json.loads(upload_data.upload)
    upload['user_id'] = tarnsform_user_id(upload['user_id'])
    upload_data.upload = json.dumps(upload)

    for i, calc_data_json in enumerate(upload_data.calcs):
        calc_data = json.loads(calc_data_json)
        metadata = calc_data['metadata']

        # transform users
        metadata['uploader'] = tarnsform_user_id(metadata['uploader']['id'])
        metadata['coauthors'] = [tarnsform_user_id(user['id']) for user in metadata['coauthors']]
        metadata['shared_with'] = [tarnsform_user_id(user['id']) for user in metadata['shared_with']]

        # transform datasets
        metadata['datasets'] = [transform_dataset(dataset) for dataset in metadata['datasets']]

        # transform references
        metadata['references'] = [transform_reference(reference) for reference in metadata['references']]

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
    '--move', is_flag=True,
    help='Instead of copying the underlying upload files, we move it and replace it with a symlink.')
@click.option(
    '--link', is_flag=True,
    help='Instead of copying the underlying upload files, we create symlinks in the target.')
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
    '--dry', is_flag=True, help='Do not actually mirror data, just fetch data and report.')
@click.option(
    '--files-only', is_flag=True,
    help=(
        'Will only copy/move files and not even look at the calculations. Useful, '
        'when moving metadata via mongo dump/restore.'))
@click.option(
    '--skip-files', is_flag=True, help='Will not copy/move/link any files.')
@click.option(
    '--migration', type=str, default=None,
    help='The name of a migration script used to transform the metadata.')
@click.option(
    '--skip-es', is_flag=True, help='Do not add mirrored data to elastic search')
@click.option(
    '--staging', is_flag=True, help='Mirror non published uploads. Only works with --move or --link.')
@click.option(
    '--replace', is_flag=True, help='Replace existing uploads.')
def mirror(
        query, move: bool, link: bool, dry: bool, files_only: bool, skip_files: bool,
        source_mapping: str, target_mapping: str, migration: str, staging: bool,
        skip_es: bool, replace: bool):

    # init global vars
    global _Dataset
    global __logger
    _Dataset = datamodel.Dataset.m_def.a_mongo.mongo_cls
    __logger = utils.get_logger(__name__)

    if staging and not (move or link):
        print('--with-staging requires either --move or --link')
        sys.exit(1)

    migration_func = None
    if migration is not None:
        if migration == 'v0.6.x':
            migration_func = v0Dot6
        if migration == 'v0.7.x':
            migration_func = v0Dot7
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
        if staging:
            query = dict(published=False)
        else:
            query = dict(published=True)

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

                if replace and not dry:
                    admin_uploads.delete_upload(upload=upload, skip_files=True)

                else:
                    if len(query) > 0:
                        print(
                            'Upload %s already exists, updating existing uploads is not '
                            'implemented yet. Skip upload.' % upload_id)
                    continue
            except KeyError:
                pass

            try:
                upload_data = client.mirror.get_upload_mirror(upload_id=upload_id).response().result
                n_calcs = len(upload_data.calcs)
            except bravado.exception.HTTPBadRequest:
                print('Could not mirror %s, it is probably not published.' % upload_id)
                n_calcs = 0
                continue

            if __in_test:
                # In tests, we mirror from our selves, remove it so it is not there for import
                proc.Calc.objects(upload_id=upload_id).delete()
                proc.Upload.objects(upload_id=upload_id).delete()
                _Dataset.objects().delete()
                nomad_doi.DOI.objects().delete()
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
        if not skip_files:
            upload_files_path = upload_data.upload_files_path
            if __in_test:
                tmp = os.path.join(config.fs.tmp, 'to_mirror')
                os.rename(upload_files_path, tmp)
                upload_files_path = tmp

            upload_files_path = source_mapping_obj.apply(upload_files_path)

            target_upload_files_path = files.PathObject(
                config.fs.public if not staging else config.fs.staging,
                upload_id, create_prefix=False, prefix=True).os_path
            target_upload_files_path = target_mapping_obj.apply(target_upload_files_path)

            if not os.path.exists(target_upload_files_path):
                if move:
                    os.rename(upload_files_path, target_upload_files_path)
                    os.symlink(os.path.abspath(target_upload_files_path), upload_files_path)

                elif link:
                    os.makedirs(os.path.dirname(target_upload_files_path.rstrip('/')), exist_ok=True)
                    os.symlink(os.path.abspath(upload_files_path), target_upload_files_path)

                else:
                    os.makedirs(target_upload_files_path)
                    for to_copy in os.listdir(upload_files_path):
                        shutil.copyfile(
                            os.path.join(upload_files_path, to_copy),
                            os.path.join(target_upload_files_path, to_copy))

        if not files_only:
            try:
                # create mongo
                upload = proc.Upload.from_json(upload_data.upload, created=True)
                if upload_data.datasets is not None:
                    for dataset in upload_data.datasets.values():
                        fix_time(dataset, ['created'])
                        _Dataset._get_collection().update(dict(_id=dataset['_id']), dataset, upsert=True)
                if upload_data.dois is not None:
                    for doi in upload_data.dois.values():
                        if doi is not None and nomad_doi.DOI.objects(doi=doi).first() is None:
                            fix_time(doi, ['create_time'])
                            nomad_doi.DOI._get_collection().update(dict(_id=doi['_id']), doi, upsert=True)
                if len(upload_data.calcs) > 0:
                    for calc in upload_data.calcs:
                        fix_time(calc, ['create_time', 'complete_time'])
                        fix_time(calc['metadata'], ['upload_time', 'last_processing'])
                    proc.Calc._get_collection().insert(upload_data.calcs)
                upload.save()
            except Exception as e:
                traceback.print_exc()

                print(
                    'Could not mirror %s with %d calcs at %s' %
                    (upload_id, n_calcs, upload_data.upload_files_path))

                print(
                    'Rolling back %s, files might need to be removed manually' % upload_id)

                # rollback
                try:
                    if upload:
                        upload.delete
                        proc.Calc.objects(upload_id=upload.upload_id).delete()
                except Exception:
                    pass

                continue

            # index es
            if not skip_es:
                with upload.entries_metadata() as entries:
                    search.index_all(entries)

        print(
            'Mirrored %s with %d calcs at %s' %
            (upload_id, n_calcs, upload_data.upload_files_path))
