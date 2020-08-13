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

from .admin import admin


@admin.group(help='Migrate data from older NOMAD versions')
@click.option('--mongo-db', help='The database name of the existing data', type=str)
def migrate(mongo_db: str, elastic_index: str):
    import pymongo
    import sys
    import json
    from nomad import config, processing as proc, doi as nomad_doi, datamodel
    from nomad.app.api.mirror import _upload_data
    from nomad.cli.client.mirror import v0Dot7, fix_time, _Dataset

    _Dataset = datamodel.Dataset.m_def.a_mongo.mongo_cls

    client = pymongo.MongoClient(config.mongo.host, config.mongo.port)
    db = getattr(client, mongo_db)
    if db is None:
        print('The given mongo database %s does not exist' % mongo_db)
        sys.exit(1)

    for upload in db.uploads.find():
        print('migrating upload with id %s' % upload['_id'])
        upload_data = _upload_data(upload['_id'], json.dumps(upload), calcs_col=db.calcs, datasets_col=db.datasets, dois_col=db.d_o_i)
        upload_data = v0Dot7(upload_data)

        proc.Upload._get_collection().insert(upload)
        for calc in db.calcs.find(dict(upload_id=upload['_id'])):
            proc.Upload.from_dict(upload).save()

        # create mongo
        try:
            upload = proc.Upload.from_json(upload_data['upload'], created=True)
            if upload_data['datasets'] is not None:
                for dataset in upload_data['datasets'].values():
                    fix_time(dataset, ['created'])
                    _Dataset._get_collection().update(dict(_id=dataset['_id']), dataset, upsert=True)
            if upload_data['dois'] is not None:
                for doi in upload_data['dois'].values():
                    if doi is not None and nomad_doi.DOI.objects(doi=doi).first() is None:
                        fix_time(doi, ['create_time'])
                        nomad_doi.DOI._get_collection().update(dict(_id=doi['_id']), doi, upsert=True)
            if len(upload_data['calcs']) > 0:
                for calc in upload_data['calcs']:
                    fix_time(calc, ['create_time', 'complete_time'])
                    fix_time(calc['metadata'], ['upload_time', 'last_processing'])
                proc.Calc._get_collection().insert(upload_data['calcs'])
            upload.save()
        except Exception as e:
            print('Could not migrate the uploads: %s' % str(e))
            print('Please reset and try again.')
            sys.exit(1)

        # reprocess
        upload.reset()
        upload.re_process_upload()
        upload.block_until_complete(interval=.5)

        if upload.tasks_status == proc.FAILURE:
            print('upload processed with failure')
