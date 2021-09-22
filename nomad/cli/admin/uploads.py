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
import typing

from nomad import config

from .admin import admin


def _run_parallel(uploads, parallel: int, callable, label: str):
    import threading

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
        if callable(upload, logger):
            completed = True

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
            current_process_step=upload.current_process_step, upload_id=upload.upload_id)
        with cv:
            cv.wait_for(lambda: state['available_threads_count'] > 0)
            state['available_threads_count'] -= 1
            thread = threading.Thread(target=lambda: process_upload(upload))
            threads.append(thread)
            thread.start()

    for thread in threads:
        thread.join()


def _run_processing(
        uploads, parallel: int, process, label: str, reprocess_running: bool = False,
        wait_until_complete: bool = True, reset_first: bool = False):

    from nomad import processing as proc

    def run_process(upload, logger):
        logger.info(
            'cli calls %s processing' % label,
            current_process=upload.current_process,
            current_process_step=upload.current_process_step, upload_id=upload.upload_id)
        if upload.process_running and not reprocess_running:
            logger.warn(
                'cannot trigger %s, since the upload is already/still processing' % label,
                current_process=upload.current_process,
                current_process_step=upload.current_process_step, upload_id=upload.upload_id)
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

    _run_parallel(uploads, parallel=parallel, callable=run_process, label=label)


@admin.group(help='Upload related commands')
@click.option('--user', help='Select uploads of user with given id', type=str)
@click.option('--unpublished', help='Select only uploads in staging', is_flag=True)
@click.option('--published', help='Select only uploads that are publised', is_flag=True)
@click.option('--outdated', help='Select published uploads with older nomad version', is_flag=True)
@click.option('--code', multiple=True, type=str, help='Select only uploads with calcs of given codes')
@click.option('--query-mongo', is_flag=True, help='Select query mongo instead of elastic search.')
@click.option('--processing', help='Select only processing uploads', is_flag=True)
@click.option('--processing-failure-uploads', is_flag=True, help='Select uploads with failed processing')
@click.option('--processing-failure-calcs', is_flag=True, help='Select uploads with calcs with failed processing')
@click.option('--processing-failure', is_flag=True, help='Select uploads where the upload or any calc has failed processing')
@click.option('--processing-incomplete-uploads', is_flag=True, help='Select uploads that have not yet been processed')
@click.option('--processing-incomplete-calcs', is_flag=True, help='Select uploads where any calc has net yot been processed')
@click.option('--processing-incomplete', is_flag=True, help='Select uploads where the upload or any calc has not yet been processed')
@click.option('--processing-necessary', is_flag=True, help='Select uploads where the upload or any calc has either not been processed or processing has failed in the past')
@click.option('--unindexed', is_flag=True, help='Select uploads that have no calcs in the elastic search index.')
@click.pass_context
def uploads(ctx, **kwargs):
    ctx.obj.uploads_kwargs = kwargs


def _query_uploads(
        uploads,
        user: str, unpublished: bool, published: bool, processing: bool, outdated: bool,
        code: typing.List[str], query_mongo: bool,
        processing_failure_uploads: bool, processing_failure_calcs: bool,
        processing_failure: bool, processing_incomplete_uploads: bool,
        processing_incomplete_calcs: bool, processing_incomplete: bool,
        processing_necessary: bool, unindexed: bool):

    import mongoengine
    import json
    import elasticsearch_dsl as es

    from nomad import infrastructure, processing as proc

    mongo_client = infrastructure.setup_mongo()
    infrastructure.setup_elastic()

    query = mongoengine.Q()
    calc_query = None
    if user is not None:
        query |= mongoengine.Q(user_id=user)
    if unpublished:
        query |= mongoengine.Q(published=False)
    if published:
        query |= mongoengine.Q(published=True)
    if processing:
        query |= mongoengine.Q(process_status__in=proc.ProcessStatus.STATUSES_PROCESSING)

    if outdated:
        uploads = proc.Calc._get_collection().distinct(
            'upload_id',
            {'metadata.nomad_version': {'$ne': config.meta.version}})
        query |= mongoengine.Q(upload_id__in=uploads)

    if code is not None and len(code) > 0:
        code_queries = [es.Q('match', **{'dft.code_name': code_name}) for code_name in code]
        code_query = es.Q('bool', should=code_queries, minimum_should_match=1)

        code_search = es.Search(index=config.elastic.index_name)
        code_search = code_search.query(code_query)
        code_search.aggs.bucket('uploads', es.A(
            'terms', field='upload_id', size=10000, min_doc_count=1))
        uploads = [
            upload['key']
            for upload in code_search.execute().aggs['uploads']['buckets']]

        query |= mongoengine.Q(upload_id__in=uploads)

    if processing_failure_calcs or processing_failure or processing_necessary:
        if calc_query is None:
            calc_query = mongoengine.Q()
        calc_query |= mongoengine.Q(process_status=proc.ProcessStatus.FAILURE)
    if processing_failure_uploads or processing_failure or processing_necessary:
        query |= mongoengine.Q(process_status=proc.ProcessStatus.FAILURE)
    if processing_incomplete_calcs or processing_incomplete or processing_necessary:
        if calc_query is None:
            calc_query = mongoengine.Q()
        calc_query |= mongoengine.Q(process_status__in=proc.ProcessStatus.STATUSES_PROCESSING)
    if processing_incomplete_uploads or processing_incomplete or processing_necessary:
        query |= mongoengine.Q(process_status__in=proc.ProcessStatus.STATUSES_PROCESSING)

    if unindexed:
        from nomad.search.v1 import quantity_values
        uploads_in_es = set(quantity_values('upload_id', page_size=1000, owner='all'))

        uploads_in_mongo = mongo_client[config.mongo.db_name]['calc'].distinct('upload_id')

        uploads_not_in_es = []
        for upload_id in uploads_in_mongo:
            if upload_id not in uploads_in_es:
                uploads_not_in_es.append(upload_id)

        query |= mongoengine.Q(
            upload_id__in=uploads_not_in_es)

    try:
        json_query = json.loads(' '.join(uploads))
        if query_mongo:
            uploads = proc.Calc.objects(**json_query).distinct(field="upload_id")
        else:
            from nomad.search.v1 import quantity_values
            uploads = list(quantity_values(
                'upload_id', query=es.Q(json_query), page_size=1000, owner='all'))
    except Exception:
        pass

    if calc_query is not None:
        query |= mongoengine.Q(
            upload_id__in=proc.Calc.objects(calc_query).distinct(field="upload_id"))
    if len(uploads) > 0:
        query |= mongoengine.Q(upload_id__in=uploads)

    return query, proc.Upload.objects(query)


@uploads.command(help='List selected uploads')
@click.argument('UPLOADS', nargs=-1)
@click.option('-c', '--calculations', is_flag=True, help='Show details about calculations.')
@click.option('--ids', is_flag=True, help='Only show a list of ids.')
@click.option('--json', is_flag=True, help='Output a JSON array of ids.')
@click.pass_context
def ls(ctx, uploads, calculations, ids, json):
    import tabulate

    _, uploads = _query_uploads(uploads, **ctx.obj.uploads_kwargs)

    def row(upload):
        row = [
            upload.upload_id,
            upload.name,
            upload.user_id,
            upload.process_status,
            upload.published]

        if calculations:
            row += [
                upload.total_calcs,
                upload.failed_calcs,
                upload.total_calcs - upload.processed_calcs]

        return row

    headers = ['id', 'name', 'user', 'process', 'published']
    if calculations:
        headers += ['calcs', 'failed', 'processing']

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


@uploads.command(help='Change the owner of the upload and all its calcs.')
@click.argument('USERNAME', nargs=1)
@click.argument('UPLOADS', nargs=-1)
@click.pass_context
def chown(ctx, username, uploads):
    from nomad import datamodel

    _, uploads = _query_uploads(uploads, **ctx.obj.uploads_kwargs)

    print('%d uploads selected, changing owner ...' % uploads.count())

    user = datamodel.User.get(username=username)
    upload_metadata = datamodel.UploadMetadata(uploader=user)
    for upload in uploads:
        upload.set_upload_metadata(upload_metadata.m_to_dict())


@uploads.command(help='Reset the processing state.')
@click.argument('UPLOADS', nargs=-1)
@click.option('--with-calcs', is_flag=True, help='Also reset all calculations.')
@click.option('--success', is_flag=True, help='Set the process status to success instead of pending')
@click.option('--failure', is_flag=True, help='Set the process status to failure instead of pending.')
@click.pass_context
def reset(ctx, uploads, with_calcs, success, failure):
    from nomad import processing as proc

    _, uploads = _query_uploads(uploads, **ctx.obj.uploads_kwargs)
    uploads_count = uploads.count()

    print('%d uploads selected, resetting their processing ...' % uploads_count)

    i = 0
    for upload in uploads:
        if with_calcs:
            calc_update = proc.Calc.reset_pymongo_update()
            if success:
                calc_update['process_status'] = proc.ProcessStatus.SUCCESS
            if failure:
                calc_update['process_status'] = proc.ProcessStatus.FAILURE

            proc.Calc._get_collection().update_many(
                dict(upload_id=upload.upload_id), {'$set': calc_update})

        upload.process_status = None
        upload.reset()
        if success:
            upload.process_status = proc.ProcessStatus.SUCCESS
        if failure:
            upload.process_status = proc.ProcessStatus.FAILURE
        upload.save()
        i += 1
        print('resetted %d of %d uploads' % (i, uploads_count))


@uploads.command(help='(Re-)index all calcs of the given uploads.')
@click.argument('UPLOADS', nargs=-1)
@click.option('--parallel', default=1, type=int, help='Use the given amount of parallel processes. Default is 1.')
@click.option('--transformer', help='Qualified name to a Python function that should be applied to each EntryMetadata.')
@click.pass_context
def index(ctx, uploads, parallel, transformer):
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
                print(f'   ERROR failed to transform calc (stop transforming for upload): {str(e)}')
                break

    def index_upload(upload, logger):
        with upload.entries_metadata() as entries:
            if transformer is not None:
                transform(entries)
            search.index([entry.m_parent for entry in entries], update_materials=True, refresh=True)

        return True

    _run_parallel(uploads, parallel, index_upload, 'index')


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
        proc.Calc.objects(upload_id=upload.upload_id).delete()
        upload.delete()


@uploads.command(help='Delete selected upload')
@click.argument('UPLOADS', nargs=-1)
@click.option('--skip-es', help='Keep the elastic index version of the data.', is_flag=True)
@click.option('--skip-mongo', help='Keep uploads and calcs in mongo.', is_flag=True)
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
@click.option('--reprocess-running', is_flag=True, help='Also reprocess already running processes.')
@click.pass_context
def re_process(ctx, uploads, parallel: int, reprocess_running: bool):
    _, uploads = _query_uploads(uploads, **ctx.obj.uploads_kwargs)
    _run_processing(
        uploads, parallel, lambda upload: upload.process_upload(), 'processing',
        reprocess_running=reprocess_running, reset_first=True)


@uploads.command(help='Repack selected uploads.')
@click.argument('UPLOADS', nargs=-1)
@click.option('--parallel', default=1, type=int, help='Use the given amount of parallel processes. Default is 1.')
@click.pass_context
def re_pack(ctx, uploads, parallel: int):
    _, uploads = _query_uploads(uploads, **ctx.obj.uploads_kwargs)
    _run_processing(
        uploads, parallel, lambda upload: upload.re_pack(), 're-packing',
        wait_until_complete=False)


@uploads.command(help='Attempt to abort the processing of uploads.')
@click.argument('UPLOADS', nargs=-1)
@click.option('--calcs', is_flag=True, help='Only stop calculation processing.')
@click.option('--kill', is_flag=True, help='Use the kill signal and force task failure.')
@click.option('--no-celery', is_flag=True, help='Do not attempt to stop the actual celery tasks')
@click.pass_context
def stop(ctx, uploads, calcs: bool, kill: bool, no_celery: bool):
    import mongoengine

    from nomad import utils, processing as proc
    query, _ = _query_uploads(uploads, **ctx.obj.uploads_kwargs)

    logger = utils.get_logger(__name__)

    def stop_all(query):
        for process in query:
            logger_kwargs = dict(upload_id=process.upload_id)
            if isinstance(process, proc.Calc):
                logger_kwargs.update(calc_id=process.calc_id)

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
    stop_all(proc.Calc.objects(running_query))
    if not calcs:
        stop_all(proc.Upload.objects(running_query))
