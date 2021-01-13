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

''' Temporary code needed for ops, version transitions, etc. '''

from nomad.datamodel import EntryMetadata
from nomad.normalizing.optimade import transform_to_v1 as optimade


external_dbs = {
    'AFLOW': '81b96683-7170-49d7-8c4e-e9f34906b3ea',
    'Materials Project': 'ab0305be-b3bb-41e6-a000-b004b49e475a',
    'OQMD': '9ac9db26-b268-4c22-b3f7-fc82ce35967e'
}
external_dbs_from_user_id = {value: key for key, value in external_dbs.items()}


def transform_to_v0_10(entry: EntryMetadata) -> EntryMetadata:
    entry = optimade(entry)

    # external db
    uploader_id = entry.uploader.user_id
    if uploader_id in external_dbs_from_user_id and entry.external_db is None:
        entry.external_db = external_dbs_from_user_id[uploader_id]

    return entry


def set_external_db_in_mongo(upload_id: str = None):
    from nomad import infrastructure
    from nomad.processing import Calc
    infrastructure.setup_mongo()

    calcs = Calc._get_collection()
    print(calcs)

    for external_db, user_id in external_dbs.items():
        query = {
            'metadata.uploader': user_id,
            'metadata.external_db': {
                '$exists': False
            }}
        if upload_id is not None:
            query['upload_id'] = upload_id

        update = {
            '$set': {
                'metadata.external_db': external_db
            }
        }

        result = calcs.update_many(query, update)
        print(
            f'updated for {external_db}, update acknowledged ({result.acknowledged}), '
            f'matches = {result.matched_count}')
