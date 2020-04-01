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


import typing
import mongoengine
import pymongo
import time
import datetime
import json

from nomad import utils, processing as proc, search
from nomad import datamodel
from nomad.cli.client import mirror


__logger = utils.get_logger(__name__)


class SourceCalc(mongoengine.Document):
    '''
    Mongo document used as a calculation, upload, and metadata db and index
    build from a given source db. Each :class:`SourceCacl` entry relates
    a pid, mainfile, upload "id" with each other for a corressponding calculation.
    It might alos contain the user metadata. The uploads are "id"ed via the
    specific path segment that identifies an upload on the CoE repo FS(s) without
    any prefixes (e.g. $EXTRACTED, /data/upload, etc.)
    '''
    pid = mongoengine.IntField(primary_key=True)
    mainfile = mongoengine.StringField()
    upload = mongoengine.StringField()
    metadata = mongoengine.DictField()

    migration_version = mongoengine.IntField(default=-1)

    extracted_prefix = '$EXTRACTED/'
    sites = ['/data/nomad/extracted/', '/nomad/repository/extracted/']
    prefixes = [extracted_prefix] + sites

    meta = dict(indexes=['upload', 'mainfile', 'migration_version'])

    _dataset_cache: dict = {}


def update_user_metadata(bulk_size: int = 1000, update_index: bool = False, **kwargs):
    ''' Goes through the whole source index to sync differences between repo user metadata
    and metadata in fairdi.

    It goes through the source index calc by calc, working in bulks. Getting the samedata
    for fairdi and updating the different calcs in mongo. Will only update user metadata.

    Uses kwargs as filters for the used source index query.
    '''
    logger = utils.get_logger(__name__)
    start_time = time.time()

    # iterate the source index in bulk
    size = SourceCalc.objects(**kwargs).count()
    count = 0
    important_changes: typing.Dict[str, typing.Any] = dict(missing_calcs=dict(), replaced=dict(), lifted_embargo=list())

    try:
        for start in range(0, size, bulk_size):
            source_bulk = SourceCalc.objects(**kwargs)[start: start + bulk_size]
            count += bulk_size

            # retrieve fairdi data for bulk (by pid)
            pids = [int(calc.pid) for calc in source_bulk]
            target_bulk = proc.Calc.objects(metadata__pid__in=pids)
            target_bulk_dict = {
                str(target.metadata['pid']): target
                for target in target_bulk}

            # comparing entries and preparing mongo update
            updates = []
            updated_calcs = []
            for source in source_bulk:
                target = target_bulk_dict.get(str(source.pid))
                if target is None:
                    # missing calc (maybe we find it another way)
                    potential_replacements = proc.Calc.objects(mainfile=source.mainfile, metadata__pid=None)
                    if potential_replacements.count() == 1:
                        target = potential_replacements.first()
                        important_changes['replaced'].setdefault(source.upload, []).append(source.pid)
                    else:
                        important_changes['missing_calcs'].setdefault(source.upload, []).append(source.pid)
                        continue

                target_metadata = datamodel.EntryMetadata(**target.metadata)
                source_metadata_normalized: typing.Dict[str, typing.Any] = dict(
                    comment=source.metadata.get('comment'),
                    references={mirror.transform_reference(ref) for ref in source.metadata['references']},
                    coauthors={mirror.tarnsform_user_id(user['id']) for user in source.metadata['coauthors']},
                    shared_with={mirror.tarnsform_user_id(user['id']) for user in source.metadata['shared_with']},
                    datasets={mirror.transform_dataset(ds) for ds in source.metadata['datasets']},
                    with_embargo=source.metadata['with_embargo'])

                target_metadata_normalized: typing.Dict[str, typing.Any] = dict(
                    comment=target_metadata.comment,
                    references=set(target_metadata.references),
                    coauthors=set(target_metadata.coauthors),
                    shared_with=set(target_metadata.shared_with),
                    datasets=set(target_metadata.datasets),
                    with_embargo=target_metadata.with_embargo)

                if source_metadata_normalized != target_metadata_normalized:
                    # do a full update of all metadata!
                    update = pymongo.UpdateOne(
                        dict(_id=target.calc_id),
                        {
                            "$set": {
                                "metadata.comment": source.metadata.get('comment'),
                                "metadata.references": list(source_metadata_normalized['references']),
                                "metadata.coauthors": list(source_metadata_normalized['coauthors']),
                                "metadata.shared_with": list(source_metadata_normalized['shared_with']),
                                "metadata.datasets": list(source_metadata_normalized['datasets']),
                                "metadata.with_embargo": source.metadata['with_embargo']
                            }
                        }
                    )
                    updates.append(update)
                    updated_calcs.append(target_metadata)

                    if target_metadata_normalized['with_embargo'] != source_metadata_normalized['with_embargo']:
                        important_changes['lifted_embargo'].append(source.pid)

            # execute mongo update
            if len(updates) > 0:
                result = proc.Calc._get_collection().bulk_write(updates)
                if result.bulk_api_result['nModified'] != len(updates):
                    logger.error('incomplete update in syncing user metadata')

                if update_index:
                    search.index_all(updated_calcs, do_refresh=False)

            # log
            eta = ((time.time() - start_time) / float(count)) * (size - count)
            print('Synced calcs %d through %d of %d with %d diffs, %s' % (
                start, start + bulk_size, size, len(updates), datetime.timedelta(seconds=eta)), flush=True)

    finally:
        with open('sync_important_changes.json', 'wt') as f:
            json.dump(important_changes, f, indent=2)

    print('done')
