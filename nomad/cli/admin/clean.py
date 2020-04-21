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
import os
import shutil
import tabulate
import elasticsearch_dsl

from nomad import config as nomad_config, infrastructure, processing
from nomad import search as nomad_search

from .admin import admin


@admin.command(help='Checks consistency of files and es vs mongo and deletes orphan entries.')
@click.option('--dry', is_flag=True, help='Do not delete anything, just check.')
@click.option('--skip-calcs', is_flag=True, help='Skip cleaning calcs with missing uploads.')
@click.option('--skip-fs', is_flag=True, help='Skip cleaning the filesystem.')
@click.option('--skip-es', is_flag=True, help='Skip cleaning the es index.')
@click.option('--staging-too', is_flag=True, help='Also clean published entries in staging, make sure these files are not due to reprocessing')
@click.option('--force', is_flag=True, help='Do not ask for confirmation.')
def clean(dry, skip_calcs, skip_fs, skip_es, staging_too, force):
    mongo_client = infrastructure.setup_mongo()
    infrastructure.setup_elastic()

    if not skip_calcs:
        uploads_for_calcs = mongo_client[nomad_config.mongo.db_name]['calc'].distinct('upload_id')
        uploads = {}
        for upload in mongo_client[nomad_config.mongo.db_name]['upload'].distinct('_id'):
            uploads[upload] = True

        missing_uploads = []
        for upload_for_calc in uploads_for_calcs:
            if upload_for_calc not in uploads:
                missing_uploads.append(upload_for_calc)

        if not dry and len(missing_uploads) > 0:
            if not force:
                input('Will delete calcs (mongo + es) for %d missing uploads. Press any key to continue ...' % len(missing_uploads))

            for upload in missing_uploads:
                mongo_client[nomad_config.mongo.db_name]['calc'].remove(dict(upload_id=upload))
                elasticsearch_dsl.Search(index=nomad_config.elastic.index_name).query('term', upload_id=upload).delete()
        else:
            print('Found %s uploads that have calcs in mongo, but there is no upload entry.' % len(missing_uploads))
            print('List first 10:')
            for upload in missing_uploads[:10]:
                print(upload)

    if not skip_fs:
        upload_dirs = []
        public_dirs, staging_dirs = {}, {}
        for bucket, dir_map in [(nomad_config.fs.public, public_dirs), (nomad_config.fs.staging, staging_dirs)]:
            for prefix in os.listdir(bucket):
                for upload in os.listdir(os.path.join(bucket, prefix)):
                    upload_dirs.append((upload, os.path.join(bucket, prefix, upload)))
                    dir_map[upload] = os.path.join(bucket, prefix, upload)

        to_delete = list(
            path for upload, path in upload_dirs
            if processing.Upload.objects(upload_id=upload).first() is None)

        if not dry and len(to_delete) > 0:
            if not force:
                input('Will delete %d upload directories. Press any key to continue ...' % len(to_delete))

            for path in to_delete:
                shutil.rmtree(path)
        else:
            print('Found %d upload directories with no upload in mongo.' % len(to_delete))
            print('List first 10:')
            for path in to_delete[:10]:
                print(path)

    if staging_too and not skip_fs:
        to_delete = list(
            path for upload, path in staging_dirs.items()
            if upload in public_dirs)

        if not dry and len(to_delete) > 0:
            if not force:
                input('Will delete %d staging upload directories. Press any key to continue ...' % len(to_delete))

            for path in to_delete:
                shutil.rmtree(path)
        else:
            print('Found %d staging upload directories with upload directory in public.' % len(to_delete))
            print('List first 10:')
            for path in to_delete[:10]:
                print(path)

    if not skip_es:
        search = nomad_search.Search(index=nomad_config.elastic.index_name)
        search.aggs.bucket('uploads', elasticsearch_dsl.A('terms', field='upload_id', size=12000))
        response = search.execute()

        to_delete = list(
            (bucket.key, bucket.doc_count)
            for bucket in response.aggregations.uploads.buckets
            if processing.Upload.objects(upload_id=bucket.key).first() is None)

        calcs = 0
        for _, upload_calcs in to_delete:
            calcs += upload_calcs

        if not dry and len(to_delete) > 0:
            if not force:
                input(
                    'Will delete %d calcs in %d uploads from ES. Press any key to continue ...' %
                    (calcs, len(to_delete)))
            for upload, _ in to_delete:
                nomad_search.Search(index=nomad_config.elastic.index_name).query('term', upload_id=upload).delete()
        else:
            print('Found %d calcs in %d uploads from ES with no upload in mongo.' % (calcs, len(to_delete)))
            print('List first 10:')
            tabulate.tabulate(to_delete, headers=['id', '#calcs'])
