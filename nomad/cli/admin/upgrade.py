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

from .admin import admin

import json
import time
from typing import List, Dict, Set, Any
from pydantic import BaseModel

from pymongo import ReplaceOne
from pymongo.database import Database
from pymongo.cursor import Cursor
from nomad import config, utils, infrastructure
from nomad.processing import ProcessStatus, Upload, Calc
from nomad.processing.data import generate_entry_id
from nomad.datamodel import Dataset


_metadata_keys_to_flatten_v0 = (
    'calc_hash', 'pid', 'external_id', 'nomad_version', 'nomad_commit', 'comment',
    'references', 'datasets')
_metadata_keys_to_remove_v0 = (
    'upload_name', 'upload_time', 'uploader', 'published', 'license', 'with_embargo',
    'external_db', 'shared_with', 'coauthors', 'processed', 'domain', 'raw_id')


@admin.group(help='Commands for upgrading to a newer NOMAD version')
def upgrade():
    pass


@upgrade.command(
    help='''Converts (upgrades) records from one mongodb and migrates to another.
            Note, it is strongly recommended to run this command with loglevel verbosed, i.e.

                nomad -v upgrade migrate-mongo ...

         ''')
@click.option(
    '--host', type=str, default=config.mongo.host,
    help='The mongodb host. By default same as the configureed db.')
@click.option(
    '--port', type=int, default=config.mongo.port,
    help='The mongodb port. By default same as the configured db.')
@click.option(
    '--src-db-name', type=str, required=True,
    help='The name of the source database.')
@click.option(
    '--dst-db-name', type=str, default=config.mongo.db_name,
    help='The name of the destination database. By default same as the configured db.')
@click.option(
    '--query', type=str,
    help='An mongo query, for selecting only some of the uploads in the source database for import.')
@click.option(
    '--ids-from-file', type=str,
    help='''Reads upload IDs from the specified file. Cannot be used together with the --query
            option.
            This can for example be used to retry just the uploads that has previously failed
            (as these ids can be exported to file using --failed-ids-to-file). You can specify both
            --ids-from-file and --failed-ids-to-file at the same time with the same file name.''')
@click.option(
    '--failed-ids-to-file', type=str,
    help='''Write the IDs of failed and skipped uploads to the specified file.
            This can for example be used to subsequently retry just the uploads that failed
            (as these ids can be loaded from file using --ids-from-file). You can specify both
            --ids-from-file and --failed-ids-to-file at the same time with the same file name.''')
@click.option(
    '--fix-problems', is_flag=True,
    help='''If a minor, fixable problem is encountered, fixes it automaticall; otherwise fail.''')
@click.option(
    '--dry', is_flag=True,
    help='Dry run (not writing anything to the destination database).')
def migrate_mongo(
        host, port, src_db_name, dst_db_name, query, ids_from_file, failed_ids_to_file,
        fix_problems, dry):
    logger = utils.get_logger('migrate-mongo')
    config.mongo.host = host
    config.mongo.port = port
    config.mongo.db_name = dst_db_name
    infrastructure.setup_mongo()

    db_src: Database = infrastructure.mongo_client.get_database(src_db_name)
    db_dst: Database = infrastructure.mongo_client.get_database(dst_db_name)

    if not dry:
        _create_collections_if_needed(db_dst)

    if ids_from_file:
        if query is not None:
            print('Cannot specify a query when using --ids-from-file.')
            return -1
        try:
            with open(ids_from_file, 'r') as f:
                upload_ids = [line.strip() for line in f.readlines() if line.strip()]
        except FileNotFoundError:
            logger.error(f'Could not open file {ids_from_file}')
            return -1
        query = {'_id': {'$in': upload_ids}}
    elif query:
        query = json.loads(query)

    logger.info('Quering uploads...')
    uploads = db_src.upload.find(query)

    _migrate_mongo_uploads(
        db_src, db_dst, uploads, failed_ids_to_file, fix_problems, dry, logger)


class _CollectionStatistics(BaseModel):
    collection_name: str
    total: int = 0
    converted: int = 0
    migrated: int = 0


class _UpgradeStatistics(BaseModel):
    uploads = _CollectionStatistics(collection_name='Uploads')
    entries = _CollectionStatistics(collection_name='Entries')
    datasets = _CollectionStatistics(collection_name='Datasets')


def _create_collections_if_needed(db_dst: Database):
    '''
    If the collections haven't yet been created, create them by calling .objects() on the
    MongoDocument class.
    '''
    if 'upload' not in db_dst.collection_names():
        Upload.objects()
    if 'calc' not in db_dst.collection_names():
        Calc.objects()
    if 'dataset' not in db_dst.collection_names():
        Dataset.m_def.a_mongo.objects()


def _migrate_mongo_uploads(
        db_src: Database, db_dst: Database, uploads: Cursor, failed_ids_to_file: bool,
        fix_problems: bool, dry: bool, logger):
    ''' Converts and/or migrates an upload and all related entries and datasets. '''
    number_of_uploads = uploads.count()
    logger.info(f'Found {number_of_uploads} uploads to import.')
    migrated_dataset_ids: Set[str] = set()
    stats = _UpgradeStatistics()
    stats.uploads.total = number_of_uploads
    count_treated = count_failures = count_processing = 0
    failed_and_skipped = []
    start_time = time.time()
    next_report_time = start_time + 60
    for upload_dict in uploads:
        upload_id = upload_dict.get('_id')
        try:
            if _is_processing(upload_dict):
                logger.warn(f'Skipping upload because it is currently processing', upload_id=upload_id)

                count_processing += 1
                failed_and_skipped.append(upload_id)
            else:
                entry_dicts, dataset_dicts = _convert_mongo_upload(
                    db_src, upload_dict, fix_problems, migrated_dataset_ids, stats, logger)
                if not dry:
                    _commit_upload(upload_dict, entry_dicts, dataset_dicts, db_dst, stats)
                del entry_dicts, dataset_dicts  # To free up memory immediately
        except Exception as e:
            logger.error(f'Failed to migrate upload: {str(e)}', upload_id=upload_id)
            count_failures += 1
            failed_and_skipped.append(upload_id)
        count_treated += 1
        # log progress report and estimated remaining time at regular intervals
        t = time.time()
        if t > next_report_time:
            elapsed = t - start_time
            progress = count_treated / number_of_uploads
            est_remain = elapsed * number_of_uploads / count_treated - elapsed
            logger.info(
                f'Elapsed: {_seconds_to_h_m(elapsed)}, progress: {progress * 100.0:6.2f} %, '
                f'est. remain: {_seconds_to_h_m(est_remain)}')
            next_report_time += 60
    # Write failed ids to file, if requested
    if failed_ids_to_file:
        with open(failed_ids_to_file, 'w') as f:
            for upload_id in failed_and_skipped:
                f.write(upload_id + '\n')
    # Log a summary
    summary = f'Summary:\n\nElapsed time: {_seconds_to_h_m(time.time() - start_time)}\n\n'
    summary += f'{"Doc type":<10}{"found":>20}{"converted":>20}{"migrated":>20}\n'
    summary += '-' * 70 + '\n'
    for sub_stats in (stats.uploads, stats.entries, stats.datasets):
        summary += (
            f'{sub_stats.collection_name:<10}{sub_stats.total:>20}'  # pylint: disable=E1101
            f'{sub_stats.converted:>20}{sub_stats.migrated:>20}\n')
    if count_failures:
        summary += f'\n!!! THERE WERE ERRORS: {count_failures} upload(s) could not be converted/migrated\n\n'
    else:
        summary += '\nNo errors occurred :-)\n\n'
    if count_processing:
        summary += f'{count_processing} uploads were skipped since they are currently processing\n\n'
    if dry:
        summary += 'Dry run - nothing written to the destination db\n\n'
    logger.info(summary)
    if count_failures:
        return -1


def _convert_mongo_upload(
        db_src: Database, upload_dict: Dict[str, Any], fix_problems: bool,
        migrated_dataset_ids: Set[str], stats: _UpgradeStatistics, logger):
    '''
    Converts (upgrades) an upload_dict and all related records. If successful,
    returns two lists: one with converted entry dicts, and one with converted dataset dicts.
    Datasets whose IDs are in `migrated_dataset_ids` are skipped (as they are assumed to
    have been migrated, or at least attempted to migrated already).
    '''
    upload_id = upload_dict['_id']
    published = upload_dict.get('publish_time') is not None
    # Determine which version we are migrating from
    is_v0 = 'user_id' in upload_dict

    _convert_mongo_proc(upload_dict)

    if is_v0:
        _rename_key(upload_dict, 'name', 'upload_name')
        _rename_key(upload_dict, 'create_time', 'upload_create_time')
        _rename_key(upload_dict, 'user_id', 'main_author')
        # Verify and then remove redundant legacy field 'published'.
        if 'published' in upload_dict:
            assert upload_dict['published'] == published, 'Inconsistency: published flag vs publish_time'
            upload_dict.pop('published')
        # Set license if missing
        if 'license' not in upload_dict:
            upload_dict['license'] = 'CC BY 4.0'

    # Fetch all entries as a dictionary
    entry_dicts = list(db_src.calc.find({'upload_id': upload_id}))
    stats.entries.total += len(entry_dicts)
    if is_v0 and entry_dicts:
        # Need to pass through all entries to check consistency and calculate common_coauthos
        first_metadata = entry_dicts[0].get('metadata', {})
        first_with_embargo = first_metadata.get('with_embargo')
        first_shared_with = first_metadata.get('shared_with')
        first_shared_with_set = set(first_shared_with or ())
        first_entry_uploader = first_metadata.get('uploader')
        first_external_db = first_metadata.get('external_db')
        first_entry_coauthors = first_metadata.get('coauthors', ())
        common_coauthors = set(_wrap_author(ca) for ca in first_entry_coauthors)

        fixed_external_db = False
        for entry_dict in entry_dicts:
            assert 'metadata' in entry_dict, 'Entry dict has no metadata key'
            entry_metadata_dict: Dict[str, Any] = entry_dict['metadata']
            with_embargo = entry_metadata_dict.get('with_embargo')
            assert with_embargo == first_with_embargo, 'Inconsistent embargo settings for entries'
            shared_with_set = set(entry_metadata_dict.get('shared_with') or ())
            assert shared_with_set == first_shared_with_set, 'Inconsistent shared_with settings for entries'
            uploader = entry_metadata_dict.get('uploader')
            assert uploader == first_entry_uploader, 'Inconsistent uploader for entries'
            external_db = entry_metadata_dict.get('external_db')
            if external_db != first_external_db:
                if external_db and first_external_db:
                    # Problem is unfixable (two non-empty values encountered)
                    assert False, 'Inconsistent external_db for entries - unfixable'
                elif not fix_problems:
                    assert False, 'Inconsistent external_db for entries - use --fix-problems to fix'
                first_external_db = first_external_db or external_db  # Fix it
                fixed_external_db = True
            common_coauthors.intersection_update(
                _wrap_author(ca) for ca in entry_metadata_dict.get('coauthors', ()))

        if fixed_external_db:
            logger.warn('Fixed inconsistent external_db for upload', upload_id=upload_id)

        # Check embargo setting on the upload
        if published:
            if first_with_embargo:
                assert upload_dict.get('embargo_length'), (
                    f'Embargo flag set on entries, but no embargo_length specified on upload')
            else:
                upload_dict['embargo_length'] = 0
        if 'embargo_length' not in upload_dict:
            upload_dict['embargo_length'] = 0
        # Fields moving to the upload level
        upload_dict['reviewers'] = first_shared_with
        upload_dict['external_db'] = first_external_db
        # main_author
        if first_entry_uploader != upload_dict['main_author']:
            assert first_external_db, 'Different uploader on entry and upload, but external_db not set'
            # It's ok, but we should fetch main_author from entries.uploader instead of uplad.user_id
            upload_dict['main_author'] = first_entry_uploader
        # coauthors - put as many as possible in Upload.coauthors, the rest in entry.entry_coauthors
        if common_coauthors:
            upload_dict['coauthors'] = [ca for ca in first_entry_coauthors if _wrap_author(ca) in common_coauthors]
    else:
        common_coauthors = set(_wrap_author(ca) for ca in upload_dict.get('coauthors', ()))

    # Check that all required fields are there
    for field in ('_id', 'upload_create_time', 'main_author', 'embargo_length', 'license'):
        assert upload_dict.get(field) is not None, f'Missing required upload field: {field}'

    # migrate entries
    new_dataset_ids: Set[str] = set()
    for entry_dict in entry_dicts:
        assert not _is_processing(entry_dict), (
            f'the entry {entry_dict["_id"]} has status processing, but the upload is not processing.')
        _convert_mongo_entry(entry_dict, common_coauthors, fix_problems, logger)
        stats.entries.converted += 1
        # Add encountered datasets
        datasets = entry_dict.get('datasets')
        if datasets:
            for dataset_id in datasets:
                if dataset_id not in migrated_dataset_ids:
                    new_dataset_ids.add(dataset_id)

    # Migrate newly discovered datasets
    stats.datasets.total += len(new_dataset_ids)
    migrated_dataset_ids.update(new_dataset_ids)
    dataset_dicts: List[Dict[str, Any]] = []
    for dataset_id in new_dataset_ids:
        dataset_dict = db_src.dataset.find_one({'_id': dataset_id})
        assert dataset_dict is not None, f'Missing dataset reference: {dataset_id}'
        _convert_mongo_dataset(dataset_dict)
        dataset_dicts.append(dataset_dict)
        stats.datasets.converted += 1
    # All conversion successful!
    stats.uploads.converted += 1
    return entry_dicts, dataset_dicts


def _convert_mongo_entry(entry_dict: Dict[str, Any], common_coauthors: Set, fix_problems: bool, logger):
    _convert_mongo_proc(entry_dict)
    # Validate the id and possibly fix problems
    generated_entry_id = generate_entry_id(entry_dict['upload_id'], entry_dict['mainfile'])
    if entry_dict['_id'] != generated_entry_id:
        if not fix_problems:
            assert False, f'Entry id {entry_dict["_id"]} does not match generated value - use --fix-problems to fix'
        logger.warn('Fixing bad id for entry', old_id=entry_dict['_id'], new_id=generated_entry_id)
        entry_dict['_id'] = generated_entry_id
    # Convert old v0 metadata
    if 'metadata' in entry_dict:
        _rename_key(entry_dict, 'parser', 'parser_name')
        _rename_key(entry_dict, 'create_time', 'entry_create_time')
        _rename_key(entry_dict, 'metadata.last_processing', 'last_processing_time')
        _rename_key(entry_dict, 'metadata.last_edit', 'last_edit_time')

        entry_metadata = entry_dict['metadata']
        # Entry coauthors
        entry_coauthors = entry_metadata.get('coauthors')
        if entry_coauthors:
            entry_coauthors = [ca for ca in entry_coauthors if _wrap_author(ca) not in common_coauthors]
            entry_dict['entry_coauthors'] = entry_coauthors or None
        # Flatten keys
        for key in _metadata_keys_to_flatten_v0:
            if key in entry_metadata:
                entry_dict[key] = entry_metadata.pop(key)
        # Remove redunant keys from the metadata dict; then it should be empty
        for key in _metadata_keys_to_remove_v0:
            if key in entry_metadata:
                entry_metadata.pop(key)
        assert not entry_metadata, f'Unexpected fields in Calc.metadata: {repr(entry_metadata)}'

    # Check that all required fields are populated
    for field in ('_id', 'upload_id', 'entry_create_time', 'parser_name'):
        assert entry_dict.get(field) is not None, f'Missing required entry field: {field}'


def _convert_mongo_proc(proc_dict: Dict[str, Any]):
    if 'tasks_status' in proc_dict:
        # Old v0 version
        process_status = proc_dict['tasks_status']
        if process_status == 'CREATED':
            process_status = ProcessStatus.READY
        assert process_status in ProcessStatus.STATUSES_VALID_IN_DB, f'Invalid tasks_status: {process_status}'
        proc_dict['process_status'] = process_status
        last_status_message = proc_dict.get('last_status_message')
        if not last_status_message:
            # Generate a nicer last_status_message
            current_process: str = proc_dict.get('current_process')
            errors: List[str] = proc_dict.get('errors')
            if errors:
                last_status_message = f'Process {current_process} failed: {errors[-1]}'
            elif current_process and process_status == ProcessStatus.SUCCESS:
                last_status_message = f'Process {current_process} completed successfully'
            else:
                assert False, f'Unexpected proc state: {current_process} / {process_status}'
            proc_dict['last_status_message'] = last_status_message

        for field in ('tasks', 'current_task', 'tasks_status'):
            proc_dict.pop(field, None)


def _convert_mongo_dataset(dataset_dict: Dict[str, Any]):
    _rename_key(dataset_dict, 'name', 'dataset_name')
    _rename_key(dataset_dict, 'created', 'dataset_create_time')
    _rename_key(dataset_dict, 'modified', 'dataset_modified_time')
    # Check that all required fields are there
    for field in ('dataset_name',):
        assert dataset_dict.get(field) is not None, (
            f'Dataset {dataset_dict["_id"]} missing required field {field}')


def _is_processing(proc_dict: Dict[str, Any]) -> bool:
    process_status = proc_dict.get('tasks_status')  # Used in v0
    if not process_status:
        process_status = proc_dict['process_status']
    return process_status in ProcessStatus.STATUSES_PROCESSING


def _commit_upload(
        upload_dict: Dict[str, Any], entry_dicts: List[Dict[str, Any]], dataset_dicts: List[Dict[str, Any]],
        db_dst: Database, stats: _UpgradeStatistics):
    # Commit datasets
    if dataset_dicts:
        dataset_writes = []
        for dataset_dict in dataset_dicts:
            dataset_writes.append(ReplaceOne({'_id': dataset_dict['_id']}, dataset_dict, upsert=True))
        db_dst.dataset.bulk_write(dataset_writes)
        stats.datasets.migrated += len(dataset_dicts)
    # Commit upload
    db_dst.upload.replace_one({'_id': upload_dict['_id']}, upload_dict, upsert=True)
    stats.uploads.migrated += 1
    # Commit entries
    if entry_dicts:
        entry_writes = []
        for entry_dict in entry_dicts:
            entry_writes.append(ReplaceOne({'_id': entry_dict['_id']}, entry_dict, upsert=True))
        db_dst.calc.bulk_write(entry_writes)
        stats.entries.migrated += len(entry_dicts)


def _rename_key(d: Dict[str, Any], old_name: str, new_name: str):
    '''
    Renames a key in the provided dictionary `d`, from `old_name` to `new_name`. We may use
    "point notation" in `old_name`, i.e. "metadata.external_id" will look for a
    key 'external_id' in d['metadata'].
    '''
    parts = old_name.split('.')
    d_current = d
    for i, part in enumerate(parts):
        if i == len(parts) - 1:
            if d_current and part in d_current:
                d[new_name] = d_current.pop(part)
        else:
            if d_current and part in d_current:
                d_current = d_current[part]
            else:
                break  # Not found


def _seconds_to_h_m(seconds: float):
    ''' Converst a time in seconds to a nicer string with hour and minute parts. '''
    return f'{seconds // 3600:>3.0f} h {(seconds % 3600) / 60:2.0f} min'


def _wrap_author(author):
    '''
    Takes an author reference (str or dict) and if it is a dict, converts it into a tuple.
    If the author is a str, it is instead returned as it is. This is used to get an object
    which is hashable and can be used in sets.
    '''
    if type(author) == str:
        return author
    return tuple((k, author[k]) for k in sorted(author.keys()))
