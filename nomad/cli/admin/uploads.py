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
import click
import tabulate
import mongoengine
import pymongo
import elasticsearch_dsl as es
import json

from nomad import processing as proc, config, infrastructure, utils, search, files, datamodel, archive

from .admin import admin, __run_processing, __run_parallel


@admin.group(help='Upload related commands')
@click.option('--user', help='Select uploads of user with given id', type=str)
@click.option('--staging', help='Select only uploads in staging', is_flag=True)
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
@click.pass_context
def uploads(
        ctx, user: str, staging: bool, processing: bool, outdated: bool,
        code: typing.List[str], query_mongo: bool,
        processing_failure_uploads: bool, processing_failure_calcs: bool, processing_failure: bool,
        processing_incomplete_uploads: bool, processing_incomplete_calcs: bool, processing_incomplete: bool,
        processing_necessary: bool):
    infrastructure.setup_mongo()
    infrastructure.setup_elastic()

    query = mongoengine.Q()
    calc_query = None
    if user is not None:
        query |= mongoengine.Q(user_id=user)
    if staging:
        query |= mongoengine.Q(published=False)
    if processing:
        query |= mongoengine.Q(process_status=proc.PROCESS_RUNNING) | mongoengine.Q(tasks_status=proc.RUNNING)

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
        calc_query |= mongoengine.Q(tasks_status=proc.FAILURE)
    if processing_failure_uploads or processing_failure or processing_necessary:
        query |= mongoengine.Q(tasks_status=proc.FAILURE)
    if processing_incomplete_calcs or processing_incomplete or processing_necessary:
        if calc_query is None:
            calc_query = mongoengine.Q()
        calc_query |= mongoengine.Q(process_status__ne=proc.PROCESS_COMPLETED)
    if processing_incomplete_uploads or processing_incomplete or processing_necessary:
        query |= mongoengine.Q(process_status__ne=proc.PROCESS_COMPLETED)

    ctx.obj.query = query
    ctx.obj.calc_query = calc_query
    ctx.obj.uploads = proc.Upload.objects(query)
    ctx.obj.query_mongo = query_mongo


def query_uploads(ctx, uploads):
    try:
        json_query = json.loads(' '.join(uploads))
        if ctx.obj.query_mongo:
            uploads = proc.Calc.objects(**json_query).distinct(field="upload_id")
        else:
            request = search.SearchRequest()
            request.q = es.Q(json_query)
            request.quantity('upload_id', size=10000)
            uploads = list(request.execute()['quantities']['upload_id']['values'])
    except Exception:
        pass

    query = ctx.obj.query
    if ctx.obj.calc_query is not None:
        query |= mongoengine.Q(
            upload_id__in=proc.Calc.objects(ctx.obj.calc_query).distinct(field="upload_id"))
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
    _, uploads = query_uploads(ctx, uploads)

    def row(upload):
        row = [
            upload.upload_id,
            upload.name,
            upload.user_id,
            upload.process_status,
            upload.tasks_status,
            upload.published]

        if calculations:
            row += [
                upload.total_calcs,
                upload.failed_calcs,
                upload.total_calcs - upload.processed_calcs]

        return row

    headers = ['id', 'name', 'user', 'process', 'tasks', 'published']
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
    _, uploads = query_uploads(ctx, uploads)

    print('%d uploads selected, changing its owner ...' % uploads.count())

    user = datamodel.User.get(username=username)

    for upload in uploads:
        upload.user_id = user.user_id
        calcs = upload.entries_metadata()

        def create_update(calc_id):
            return pymongo.UpdateOne(
                {'_id': calc_id},
                {'$set': {'metadata.uploader': user.user_id}})

        proc.Calc._get_collection().bulk_write(
            [create_update(calc_id) for calc_id in upload.entry_ids()])
        upload.save()

        with upload.entries_metadata() as calcs:
            search.index_all(calcs, do_refresh=False)
        search.refresh()


@uploads.command(help='Reset the processing state.')
@click.argument('UPLOADS', nargs=-1)
@click.option('--with-calcs', is_flag=True, help='Also reset all calculations.')
@click.pass_context
def reset(ctx, uploads, with_calcs):
    _, uploads = query_uploads(ctx, uploads)
    uploads_count = uploads.count()

    print('%d uploads selected, resetting their processing ...' % uploads_count)

    i = 0
    for upload in uploads:
        proc.Calc._get_collection().update_many(
            dict(upload_id=upload.upload_id),
            {'$set': proc.Calc.reset_pymongo_update()})

        upload.process_status = None
        upload.reset()
        upload.save()
        i += 1
        print('resetted %d of %d uploads' % (i, uploads_count))


@uploads.command(help='(Re-)index all calcs of the given uploads.')
@click.argument('UPLOADS', nargs=-1)
@click.option('--parallel', default=1, type=int, help='Use the given amount of parallel processes. Default is 1.')
@click.pass_context
def index(ctx, uploads, parallel):
    _, uploads = query_uploads(ctx, uploads)

    def index_upload(upload, logger):
        with upload.entries_metadata() as calcs:
            # This is just a temporary fix to update the group hash without re-processing
            try:
                for calc in calcs:
                    if calc.dft is not None:
                        calc.dft.update_group_hash()
            except Exception:
                pass
            failed = search.index_all(calcs)
            if failed > 0:
                print('    WARNING failed to index %d entries' % failed)

        return True

    __run_parallel(uploads, parallel, index_upload, 'index')


def delete_upload(upload, skip_es: bool = False, skip_files: bool = False, skip_mongo: bool = False):
    # delete elastic
    if not skip_es:
        search.delete_upload(upload_id=upload.upload_id)

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
    _, uploads = query_uploads(ctx, uploads)

    print('%d uploads selected, deleting ...' % uploads.count())

    for upload in uploads:
        delete_upload(upload, skip_es=skip_es, skip_mongo=skip_mongo, skip_files=skip_files)


@uploads.command(help='Create msgpack file for upload')
@click.argument('UPLOADS', nargs=-1)
@click.pass_context
def msgpack(ctx, uploads):
    _, uploads = query_uploads(ctx, uploads)

    for upload in uploads:
        upload_files = files.UploadFiles.get(upload_id=upload.upload_id)

        if isinstance(upload_files, files.PublicUploadFiles):
            def iterator(zf, names):
                for name in names:
                    calc_id = name.strip('.json')
                    with zf.open(name) as f:
                        yield (calc_id, json.load(f))

            for access in ['public', 'restricted']:
                with upload_files._open_zip_file('archive', access, 'json') as zf:
                    archive_path = upload_files._file_object('archive', access, 'msg', 'msg').os_path
                    names = [name for name in zf.namelist() if name.endswith('json')]
                    archive.write_archive(archive_path, len(names), iterator(zf, names))
                print('wrote msgpack archive %s' % archive_path)


@uploads.command(help='Reprocess selected uploads.')
@click.argument('UPLOADS', nargs=-1)
@click.option('--parallel', default=1, type=int, help='Use the given amount of parallel processes. Default is 1.')
@click.option('--reprocess-running', is_flag=True, help='Also reprocess already running processes.')
@click.pass_context
def re_process(ctx, uploads, parallel: int, reprocess_running: bool):
    _, uploads = query_uploads(ctx, uploads)
    __run_processing(
        uploads, parallel, lambda upload: upload.re_process_upload(), 're-processing',
        reprocess_running=reprocess_running)


@uploads.command(help='Repack selected uploads.')
@click.argument('UPLOADS', nargs=-1)
@click.option('--parallel', default=1, type=int, help='Use the given amount of parallel processes. Default is 1.')
@click.pass_context
def re_pack(ctx, uploads, parallel: int):
    _, uploads = query_uploads(ctx, uploads)
    __run_processing(uploads, parallel, lambda upload: upload.re_pack(), 're-packing')


@uploads.command(help='Attempt to abort the processing of uploads.')
@click.argument('UPLOADS', nargs=-1)
@click.option('--calcs', is_flag=True, help='Only stop calculation processing.')
@click.option('--kill', is_flag=True, help='Use the kill signal and force task failure.')
@click.option('--no-celery', is_flag=True, help='Do not attempt to stop the actual celery tasks')
@click.pass_context
def stop(ctx, uploads, calcs: bool, kill: bool, no_celery: bool):
    query, _ = query_uploads(ctx, uploads)

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
                process.process_status = proc.PROCESS_COMPLETED
                process.on_process_complete(None)
                process.save()

    running_query = query & (mongoengine.Q(process_status=proc.PROCESS_RUNNING) | mongoengine.Q(tasks_status=proc.RUNNING))
    stop_all(proc.Calc.objects(running_query))
    if not calcs:
        stop_all(proc.Upload.objects(running_query))
