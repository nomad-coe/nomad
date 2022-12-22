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

import typing
import click
import json
import os
import traceback

from nomad import config

from .admin import admin


def _run_parallel(uploads, parallel: int, callable, label: str, print_progress: int = 0):
    import threading
    import time

    from nomad import utils, processing as proc

    if isinstance(uploads, (tuple, list)):
        uploads_count = len(uploads)

    else:
        uploads_count = uploads.count()
        uploads = list(uploads)  # copy the whole mongo query set to avoid cursor timeouts

    cv = threading.Condition()
    threads: typing.List[threading.Thread] = []

    state = dict(
        completed_count=0,
        skipped_count=0,
        available_threads_count=parallel)

    logger = utils.get_logger(__name__)

    print('%d uploads selected, %s ...' % (uploads_count, label))

    def process_upload(upload: proc.Upload):
        logger.info('%s started' % label, upload_id=upload.upload_id)

        completed = False
        try:
            if callable(upload, logger):
                completed = True
        except Exception as e:
            completed = True
            logger.error('%s failed' % label, upload_id=upload.upload_id, exc_info=e)

        with cv:
            state['completed_count'] += 1 if completed else 0
            state['skipped_count'] += 1 if not completed else 0
            state['available_threads_count'] += 1

            print(
                '   %s %s and skipped %s of %s uploads' %
                (label, state['completed_count'], state['skipped_count'], uploads_count))

            cv.notify()

    for upload in uploads:
        logger.info(
            'cli schedules parallel %s processing for upload' % label,
            current_process=upload.current_process,
            last_status_message=upload.last_status_message, upload_id=upload.upload_id)
        with cv:
            cv.wait_for(lambda: state['available_threads_count'] > 0)
            state['available_threads_count'] -= 1
            thread = threading.Thread(target=lambda: process_upload(upload))
            threads.append(thread)
            thread.start()

    def print_progress_lines():
        while True:
            time.sleep(print_progress)
            print('.', flush=True)

    if print_progress > 0:
        progress_thread = threading.Thread(target=print_progress_lines)
        progress_thread.daemon = True
        progress_thread.start()

    for thread in threads:
        thread.join()


def _run_processing(
        uploads, parallel: int, process, label: str, process_running: bool = False,
        wait_until_complete: bool = True, reset_first: bool = False, **kwargs):

    from nomad import processing as proc

    def run_process(upload, logger):
        logger.info(
            'cli calls %s processing' % label,
            current_process=upload.current_process,
            last_status_message=upload.last_status_message, upload_id=upload.upload_id)
        if upload.process_running and not process_running:
            logger.warn(
                'cannot trigger %s, since the upload is already/still processing' % label,
                current_process=upload.current_process,
                last_status_message=upload.last_status_message, upload_id=upload.upload_id)
            return False

        if reset_first:
            upload.reset(force=True)
        elif upload.process_running:
            upload.reset(force=True, process_status=proc.ProcessStatus.FAILURE)

        process(upload)
        if wait_until_complete:
            upload.block_until_complete(interval=.5)
        else:
            upload.block_until_complete_or_waiting_for_result(interval=.5)

        if upload.process_status == proc.ProcessStatus.FAILURE:
            logger.info('%s with failure' % label, upload_id=upload.upload_id)

        logger.info('%s complete' % label, upload_id=upload.upload_id)
        return True

    _run_parallel(uploads, parallel=parallel, callable=run_process, label=label, **kwargs)


@admin.group(help='Upload related commands')
@click.option('--uploads-mongo-query', type=str, help='A query')
@click.option('--entries-mongo-query', type=str, help='A query')
@click.option('--entries-es-query', type=str, help='A query')
@click.option('--unpublished', help='Select only uploads in staging', is_flag=True)
@click.option('--published', help='Select only uploads that are publised', is_flag=True)
@click.option('--outdated', help='Select published uploads with older nomad version', is_flag=True)
@click.option('--processing', help='Select only processing uploads', is_flag=True)
@click.option('--processing-failure-uploads', is_flag=True, help='Select uploads with failed processing')
@click.option('--processing-failure-entries', is_flag=True, help='Select uploads with entries with failed processing')
@click.option('--processing-failure', is_flag=True, help='Select uploads where the upload or any entry has failed processing')
@click.option('--processing-incomplete-uploads', is_flag=True, help='Select uploads that have not yet been processed')
@click.option('--processing-incomplete-entries', is_flag=True, help='Select uploads where any entry has net yot been processed')
@click.option('--processing-incomplete', is_flag=True, help='Select uploads where the upload or any entry has not yet been processed')
@click.option('--processing-necessary', is_flag=True, help='Select uploads where the upload or any entry has either not been processed or processing has failed in the past')
@click.option('--unindexed', is_flag=True, help='Select uploads that have no entries in the elastic search index.')
@click.pass_context
def uploads(ctx, **kwargs):
    ctx.obj.uploads_kwargs = kwargs


def _query_uploads(
        uploads,
        unpublished: bool, published: bool, processing: bool, outdated: bool,
        uploads_mongo_query: str, entries_mongo_query: str, entries_es_query: str,
        processing_failure_uploads: bool, processing_failure_entries: bool,
        processing_failure: bool, processing_incomplete_uploads: bool,
        processing_incomplete_entries: bool, processing_incomplete: bool,
        processing_necessary: bool, unindexed: bool):

    '''
    Produces a list of uploads (mongoengine proc.Upload objects) based on a given
    list of upoad ids and further filter parameters.
    '''

    from typing import Set, cast
    import json
    from mongoengine import Q

    from nomad import infrastructure, processing as proc, search
    from nomad.app.v1 import models

    infrastructure.setup_mongo()
    infrastructure.setup_elastic()

    if uploads is not None and len(uploads) == 0:
        uploads = None  # None meaning all uploads
    else:
        uploads = set(uploads)

    entries_mongo_query_q = Q()
    if entries_mongo_query:
        entries_mongo_query_q = Q(**json.loads(entries_mongo_query))

    entries_query_uploads: Set[str] = None

    if entries_es_query is not None:
        entries_es_query_dict = json.loads(entries_es_query)
        results = search.search(
            owner='admin',
            query=entries_es_query_dict,
            pagination=models.MetadataPagination(page_size=0),
            user_id=config.services.admin_user_id,
            aggregations={
                'uploads': models.Aggregation(
                    terms=models.TermsAggregation(
                        quantity='upload_id',
                        pagination=models.AggregationPagination(
                            page_size=10000
                        )
                    )
                )
            })

        entries_query_uploads = set([
            cast(str, bucket.value)
            for bucket in results.aggregations['uploads'].terms.data])  # pylint: disable=no-member

    if outdated:
        entries_mongo_query_q &= Q(nomad_version={'$ne': config.meta.version})

    if processing_failure_entries or processing_failure or processing_necessary:
        entries_mongo_query_q &= Q(process_status=proc.ProcessStatus.FAILURE)

    if processing_incomplete_entries or processing_incomplete or processing_necessary:
        entries_mongo_query_q &= Q(process_status__in=proc.ProcessStatus.STATUSES_PROCESSING)

    if entries_mongo_query_q == Q():
        # If there is no entry based query, we get the list of all uploads from the upload
        # and not the entry collection. This ensures that we will also catch uploads that
        # do not have an entry.
        mongo_entry_based_uploads = set(proc.Upload.objects().distinct(field="upload_id"))
    else:
        mongo_entry_based_uploads = set(proc.Entry.objects(entries_mongo_query_q).distinct(field="upload_id"))

    if entries_query_uploads is not None:
        entries_query_uploads = entries_query_uploads.intersection(mongo_entry_based_uploads)
    else:
        entries_query_uploads = mongo_entry_based_uploads

    if entries_query_uploads:
        uploads_mongo_query_q = Q(upload_id__in=list(entries_query_uploads))
    else:
        uploads_mongo_query_q = Q()

    if uploads_mongo_query:
        uploads_mongo_query_q &= Q(**json.loads(uploads_mongo_query))

    if published:
        uploads_mongo_query_q &= Q(publish_time__exists=True)

    if unpublished:
        uploads_mongo_query_q &= Q(publish_time__exists=False)

    if processing:
        uploads_mongo_query_q &= Q(process_status__in=proc.ProcessStatus.STATUSES_PROCESSING)

    if processing_failure_uploads or processing_failure or processing_necessary:
        uploads_mongo_query_q &= Q(process_status=proc.ProcessStatus.FAILURE)

    if processing_incomplete_uploads or processing_incomplete or processing_necessary:
        uploads_mongo_query_q &= Q(process_status__in=proc.ProcessStatus.STATUSES_PROCESSING)

    final_query = uploads_mongo_query_q
    if uploads is not None:
        final_query &= Q(upload_id__in=list(uploads))

    return final_query, proc.Upload.objects(final_query)


@uploads.command(help='List selected uploads')
@click.argument('UPLOADS', nargs=-1)
@click.option('-e', '--entries', is_flag=True, help='Show details about entries.')
@click.option('--ids', is_flag=True, help='Only show a list of ids.')
@click.option('--json', is_flag=True, help='Output a JSON array of ids.')
@click.pass_context
def ls(ctx, uploads, entries, ids, json):
    import tabulate

    _, uploads = _query_uploads(uploads, **ctx.obj.uploads_kwargs)

    def row(upload):
        row = [
            upload.upload_id,
            upload.upload_name,
            upload.main_author,
            upload.process_status,
            upload.published]

        if entries:
            row += [
                upload.total_entries_count,
                upload.failed_entries_count,
                upload.total_entries_count - upload.processed_entries_count]

        return row

    headers = ['id', 'upload_name', 'user', 'process', 'published']
    if entries:
        headers += ['entries', 'failed', 'processing']

    if ids:
        for upload in uploads:
            print(upload.upload_id)
        return

    if json:
        print('[%s]' % ','.join(['"%s"' % upload.upload_id for upload in uploads]))
        return

    print('%d uploads selected, showing no more than first 10' % uploads.count())
    print(tabulate.tabulate(
        [row(upload) for upload in uploads[:10]],
        headers=headers))


@uploads.command(help='Change the owner of the upload and all its entries.')
@click.argument('USERNAME', nargs=1)
@click.argument('UPLOADS', nargs=-1)
@click.pass_context
def chown(ctx, username, uploads):
    from nomad import datamodel

    _, uploads = _query_uploads(uploads, **ctx.obj.uploads_kwargs)

    print('%d uploads selected, changing owner ...' % uploads.count())

    user = datamodel.User.get(username=username)
    for upload in uploads:
        upload.edit_upload_metadata(
            edit_request_json=dict(metadata={'main_author': user.user_id}),
            user_id=config.services.admin_user_id)


@uploads.command(help='Reset the processing state.')
@click.argument('UPLOADS', nargs=-1)
@click.option('--with-entries', is_flag=True, help='Also reset all entries.')
@click.option('--success', is_flag=True, help='Set the process status to success instead of pending')
@click.option('--failure', is_flag=True, help='Set the process status to failure instead of pending.')
@click.pass_context
def reset(ctx, uploads, with_entries, success, failure):
    from nomad import processing as proc

    _, uploads = _query_uploads(uploads, **ctx.obj.uploads_kwargs)
    uploads_count = uploads.count()

    print('%d uploads selected, resetting their processing ...' % uploads_count)

    i = 0
    for upload in uploads:
        if with_entries:
            entry_update = proc.Entry.reset_pymongo_update()
            if success:
                entry_update['process_status'] = proc.ProcessStatus.SUCCESS
            if failure:
                entry_update['process_status'] = proc.ProcessStatus.FAILURE

            proc.Entry._get_collection().update_many(
                dict(upload_id=upload.upload_id), {'$set': entry_update})

        upload.reset(force=True)
        if success:
            upload.process_status = proc.ProcessStatus.SUCCESS
        if failure:
            upload.process_status = proc.ProcessStatus.FAILURE
        upload.save()
        i += 1
        print('resetted %d of %d uploads' % (i, uploads_count))


@uploads.command(help='(Re-)index all entries of the given uploads.')
@click.argument('UPLOADS', nargs=-1)
@click.option('--parallel', default=1, type=int, help='Use the given amount of parallel processes. Default is 1.')
@click.option('--transformer', help='Qualified name to a Python function that should be applied to each EntryMetadata.')
@click.option('--skip-materials', is_flag=True, help='Only update the entries index.')
@click.option('--print-progress', default=0, type=int, help='Prints a dot every given seconds. Can be used to keep terminal open that have an i/o-based timeout.')
@click.pass_context
def index(ctx, uploads, parallel, transformer, skip_materials, print_progress):
    from nomad import search

    transformer_func = None
    if transformer is not None:
        import importlib
        module_name, func_name = transformer.rsplit('.', 1)
        module = importlib.import_module(module_name)
        transformer_func = getattr(module, func_name)

    _, uploads = _query_uploads(uploads, **ctx.obj.uploads_kwargs)

    def transform(entries):
        for entry in entries:
            try:
                entry = transformer_func(entry)
            except Exception as e:
                import traceback
                traceback.print_exc()
                print(f'   ERROR failed to transform entry (stop transforming for upload): {str(e)}')
                break

    def index_upload(upload, logger):
        with upload.entries_metadata() as entries:
            if transformer is not None:
                transform(entries)
            archives = [entry.m_parent for entry in entries]
            search.index(archives, update_materials=not skip_materials, refresh=True)

        return True

    _run_parallel(uploads, parallel, index_upload, 'index', print_progress=print_progress)


def delete_upload(upload, skip_es: bool = False, skip_files: bool = False, skip_mongo: bool = False):
    from nomad import search, files, utils, processing as proc

    # delete elastic
    if not skip_es:
        search.delete_upload(upload_id=upload.upload_id, update_materials=True, refresh=True)

    # delete files
    if not skip_files:
        # do it twice to get the two potential versions 'public' and 'staging'
        for _ in range(0, 2):
            upload_files = files.UploadFiles.get(upload_id=upload.upload_id)

            try:
                if upload_files is not None:
                    upload_files.delete()
            except Exception as e:
                logger = utils.get_logger(__name__)
                logger.error('could not delete files', exc_info=e)
                break

    # delete mongo
    if not skip_mongo:
        proc.Entry.objects(upload_id=upload.upload_id).delete()
        upload.delete()


@uploads.command(help='Delete selected upload')
@click.argument('UPLOADS', nargs=-1)
@click.option('--skip-es', help='Keep the elastic index version of the data.', is_flag=True)
@click.option('--skip-mongo', help='Keep uploads and entries in mongo.', is_flag=True)
@click.option('--skip-files', help='Keep all related files.', is_flag=True)
@click.pass_context
def rm(ctx, uploads, skip_es, skip_mongo, skip_files):
    _, uploads = _query_uploads(uploads, **ctx.obj.uploads_kwargs)

    print('%d uploads selected, deleting ...' % uploads.count())

    for upload in uploads:
        delete_upload(upload, skip_es=skip_es, skip_mongo=skip_mongo, skip_files=skip_files)


@uploads.command(help='Reprocess selected uploads.')
@click.argument('UPLOADS', nargs=-1)
@click.option('--parallel', default=1, type=int, help='Use the given amount of parallel processes. Default is 1.')
@click.option('--process-running', is_flag=True, help='Also reprocess already running processes.')
@click.option('--setting', type=str, multiple=True, help='key=value to overwrite a default reprocess config setting.')
@click.option('--print-progress', default=0, type=int, help='Prints a dot every given seconds. Can be used to keep terminal open that have an i/o-based timeout.')
@click.pass_context
def process(ctx, uploads, parallel: int, process_running: bool, setting: typing.List[str], print_progress: int):
    _, uploads = _query_uploads(uploads, **ctx.obj.uploads_kwargs)
    settings: typing.Dict[str, bool] = {}
    for settings_str in setting:
        key, value = settings_str.split('=')
        settings[key] = bool(value)
    _run_processing(
        uploads, parallel, lambda upload: upload.process_upload(reprocess_settings=settings),
        'processing', process_running=process_running, reset_first=True, print_progress=print_progress)


@uploads.command(help='Repack selected uploads.')
@click.argument('UPLOADS', nargs=-1)
@click.pass_context
def re_pack(ctx, uploads):
    _, uploads = _query_uploads(uploads, **ctx.obj.uploads_kwargs)

    for upload in uploads:
        if not upload.published:
            print(f'Cannot repack unpublished upload {upload.upload_id}')
            continue

        upload.upload_files.re_pack(upload.with_embargo)
        print(f'successfully re-packed {upload.upload_id}')


@uploads.command(help='Attempt to abort the processing of uploads.')
@click.argument('UPLOADS', nargs=-1)
@click.option('--entries', is_flag=True, help='Only stop entries processing.')
@click.option('--kill', is_flag=True, help='Use the kill signal and force task failure.')
@click.option('--no-celery', is_flag=True, help='Do not attempt to stop the actual celery tasks')
@click.pass_context
def stop(ctx, uploads, entries: bool, kill: bool, no_celery: bool):
    import mongoengine

    from nomad import utils, processing as proc
    query, _ = _query_uploads(uploads, **ctx.obj.uploads_kwargs)

    logger = utils.get_logger(__name__)

    def stop_all(query):
        for process in query:
            logger_kwargs = dict(upload_id=process.upload_id)
            if isinstance(process, proc.Entry):
                logger_kwargs.update(entry_id=process.entry_id)

            if not no_celery:
                logger.info(
                    'send terminate celery task', celery_task_id=process.celery_task_id,
                    kill=kill, **logger_kwargs)

            kwargs = {}
            if kill:
                kwargs.update(signal='SIGKILL')
            try:
                if not no_celery:
                    proc.app.control.revoke(process.celery_task_id, terminate=True, **kwargs)
            except Exception as e:
                logger.warning(
                    'could not revoke celery task', exc_info=e,
                    celery_task_id=process.celery_task_id, **logger_kwargs)

            if kill:
                logger.info(
                    'fail proc', celery_task_id=process.celery_task_id, kill=kill,
                    **logger_kwargs)

                process.fail('process terminate via nomad cli')

    running_query = query & mongoengine.Q(process_status__in=proc.ProcessStatus.STATUSES_PROCESSING)
    stop_all(proc.Entry.objects(running_query))
    if not entries:
        stop_all(proc.Upload.objects(running_query))


@uploads.group(help='Check certain integrity criteria')
@click.pass_context
def integrity(ctx):
    pass


@integrity.command(help='Uploads that have more entries in mongo than in ES.')
@click.argument('UPLOADS', nargs=-1)
@click.pass_context
def entry_index(ctx, uploads):
    from nomad.search import search
    from nomad.processing import Upload
    from nomad.app.v1.models import Pagination

    _, uploads = _query_uploads(uploads, **ctx.obj.uploads_kwargs)

    upload: Upload = None
    for upload in uploads:
        search_results = search(
            owner='admin',
            query=dict(upload_id=upload.upload_id),
            pagination=Pagination(page_size=0),
            user_id=config.services.admin_user_id)

        if search_results.pagination.total != upload.total_entries_count:
            print(upload.upload_id)


@uploads.command(help='''Export one or more uploads as bundles.''')
@click.argument('UPLOADS', nargs=-1)
@click.option(
    '--out-dir', type=str,
    help=f'Output folder. Default value is "{config.bundle_export.default_cli_bundle_export_path}" (defined in config)')
@click.option(
    '--uncompressed', is_flag=True,
    help='Specify to export each bundle as an uncompressed folder, instead of a zip-file.')
@click.option(
    '--overwrite', is_flag=True,
    help='Specify to, for each bundle, overwrite the destination file/folder if it already exists.')
@click.option(
    '--settings', '-s', type=str,
    help='''The export settings, specified as json. Settings not specified in the dictionary
            will be set to the default values.''')
@click.option(
    '--ignore-errors', '-i', is_flag=True,
    help='''Specify to ignore errors on individual uploads, and continue exporting (the default
            behaviour is to abort on first failing upload).''')
@click.pass_context
def export_bundle(ctx, uploads, out_dir, uncompressed, overwrite, settings, ignore_errors):
    from nomad.bundles import BundleExporter

    _, uploads = _query_uploads(uploads, **ctx.obj.uploads_kwargs)

    default_export_settings = config.bundle_export.default_settings.customize(
        config.bundle_export.default_settings_cli)
    if settings:
        settings = json.loads(settings)
        try:
            export_settings = default_export_settings.customize(settings)
            BundleExporter.check_export_settings(export_settings)
        except Exception as e:
            # Invalid setting provided
            print(e)
            print('\nAvailable settings and their configured default values:')
            for k, v in default_export_settings.dict().items():
                print(f'    {k:<40}: {v}')
            return -1
    else:
        export_settings = default_export_settings

    out_dir = out_dir or config.bundle_export.default_cli_bundle_export_path

    count = 0
    count_failed = 0
    try:
        for upload in uploads:
            try:
                count += 1
                print(f'Exporting upload {count} of {len(uploads)}: {upload.upload_id}')
                bundle_name = f'bundle_{upload.upload_id}' + ('.zip' if not uncompressed else '')
                export_path = os.path.abspath(os.path.join(out_dir, bundle_name))

                BundleExporter(
                    upload,
                    export_as_stream=False,
                    export_path=export_path,
                    zipped=not uncompressed,
                    overwrite=overwrite,
                    export_settings=export_settings).export_bundle()

            except Exception:
                count_failed += 1
                print(f'ERROR: Failed to export bundle: {upload.upload_id}')
                traceback.print_exc()
                if not ignore_errors:
                    print('Aborting export ...')
                    return -1
    finally:
        print('-' * 80 + '\nSummary:\n' + '-' * 80)
        print(f'Successfully exported: {count - count_failed} out of {len(uploads)}')
        if count_failed:
            print(f'FAILED to export: {count_failed}')
            if not ignore_errors:
                print(f'Aborted export for {len(uploads) - count} subsequent uploads.')


@uploads.command(help='''Import one or more uploads from bundles. Unless specified by the user,
                         the configured default import settings are used.''')
@click.option(
    '--in', 'input_path', type=str,
    help='The input path, specifying a bundle or a folder containing multiple bundles.')
@click.option(
    '--multi', '-m', is_flag=True,
    help=f'''Specify this flag if the input_path is a folder containing multiple bundles, and
             all these should be imported. If this option is specified without specifying --in, we
             will default the input path to {config.bundle_import.default_cli_bundle_import_path}''')
@click.option(
    '--settings', '-s', type=str,
    help='''The import settings, specified as json. Settings not specified in the dictionary
            will be set to the default values.''')
@click.option(
    '--embargo_length', '-e', type=int,
    help='''The embargo length (0-36 months). 0 means no embargo. If unspecified, the embargo
            period defined in the bundle will be used.''')
@click.option(
    '--use-celery', '-c', is_flag=True,
    help='''If specified, uses celery and the worker pool to do the main part of the import.
            NOTE: this requires that the workers can access the bundle via the exact same path.''')
@click.option(
    '--ignore-errors', '-i', is_flag=True,
    help='''Specify this flag to ignore errors on individual bundles, and continue importing
            (the default behaviour is to abort on first failing bundle).''')
@click.pass_context
def import_bundle(ctx, input_path, multi, settings, embargo_length, use_celery, ignore_errors):
    from nomad.bundles import BundleImporter
    from nomad import infrastructure

    for key, value in ctx.obj.uploads_kwargs.items():
        if value:
            print(f'Bad argument: "{key}" (query args are not applicable for bundle-import)')
            return -1
    if not input_path and not multi:
        print('Need to specify a bundle source, using --in, --multi, or both.')
        return -1
    if multi and not input_path:
        input_path = config.bundle_import.default_cli_bundle_import_path

    if multi:
        if not os.path.isdir(input_path):
            print(f'No such folder: "{input_path}"')
            return -1
        bundle_paths = [
            os.path.abspath(os.path.join(input_path, element)) for element in sorted(os.listdir(input_path))]
    else:
        if not os.path.exists(input_path):
            print(f'Path not found: {input_path}')
        bundle_paths = [os.path.abspath(input_path)]

    default_import_settings = config.bundle_import.default_settings.customize(
        config.bundle_import.default_settings_cli)

    if settings:
        settings = json.loads(settings)
        try:
            import_settings = default_import_settings.customize(settings)
        except Exception as e:
            # Invalid setting provided
            print(e)
            print('\nAvailable settings and their configured default values:')
            for k, v in default_import_settings.dict().items():
                print(f'    {k:<40}: {v}')
            return -1
    else:
        import_settings = default_import_settings

    infrastructure.setup()

    count = count_failed = 0
    try:
        for bundle_path in bundle_paths:
            if BundleImporter.looks_like_a_bundle(bundle_path):
                count += 1
                print(f'Importing bundle: {bundle_path}')
                bundle_importer: BundleImporter = None
                try:
                    bundle_importer = BundleImporter(None, import_settings, embargo_length)
                    bundle_importer.open(bundle_path)
                    upload = bundle_importer.create_upload_skeleton()
                    if use_celery:
                        # Run using celery (as a @process)
                        bundle_importer.close()
                        upload.import_bundle(bundle_path, import_settings, embargo_length)
                    else:
                        # Run in same thread (as a @process_local)
                        upload.import_bundle_local(bundle_importer)
                    if upload.errors:
                        raise RuntimeError(f'Import failed: {upload.errors[0]}')
                except Exception:
                    count_failed += 1
                    print(f'ERROR: Failed to import bundle: {bundle_path}')
                    traceback.print_exc()
                    if not ignore_errors:
                        print('Aborting import ...')
                        return -1
                finally:
                    if bundle_importer:
                        bundle_importer.close()
            else:
                print(f'Skipping, does not look like a bundle: {bundle_path}')
    finally:
        print('-' * 80 + '\nSummary:\n' + '-' * 80)
        print(f'Number of bundles successfully {"sent to worker" if use_celery else "imported"}: {count - count_failed}')
        if count_failed:
            print(f'FAILED to import: {count_failed}')
            if not ignore_errors:
                print('Aborted import after first failure.')
