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

from typing import List
import click
from tabulate import tabulate
from mongoengine import Q
from pymongo import UpdateOne
import threading
import elasticsearch_dsl as es

from nomad import processing as proc, config, infrastructure, utils, search, files, coe_repo
from .admin import admin


@admin.group(help='Upload related commands')
@click.option('--user', help='Select uploads of user with given id', type=str)
@click.option('--staging', help='Select only uploads in staging', is_flag=True)
@click.option('--processing', help='Select only processing uploads', is_flag=True)
@click.option('--outdated', help='Select published uploads with older nomad version', is_flag=True)
@click.option('--code', multiple=True, type=str, help='Select only uploads with calcs of given codes')
@click.pass_context
def uploads(ctx, user: str, staging: bool, processing: bool, outdated: bool, code: List[str]):
    infrastructure.setup_mongo()
    infrastructure.setup_elastic()

    query = Q()
    if user is not None:
        query &= Q(user_id=user)
    if staging:
        query &= Q(published=False)
    if processing:
        query &= Q(process_status=proc.PROCESS_RUNNING) | Q(tasks_status=proc.RUNNING)

    if outdated:
        uploads = proc.Calc._get_collection().distinct(
            'upload_id',
            {'metadata.nomad_version': {'$ne': config.version}})
        query &= Q(upload_id__in=uploads)

    if code is not None and len(code) > 0:
        code_queries = [es.Q('match', code_name=code_name) for code_name in code]
        code_query = es.Q('bool', should=code_queries, minimum_should_match=1)

        code_search = es.Search(index=config.elastic.index_name)
        code_search = code_search.query(code_query)
        code_search.aggs.bucket('uploads', es.A(
            'terms', field='upload_id', size=10000, min_doc_count=1))
        uploads = [
            upload['key']
            for upload in code_search.execute().aggs['uploads']['buckets']]

        query &= Q(upload_id__in=uploads)

    ctx.obj.query = query
    ctx.obj.uploads = proc.Upload.objects(query)


def query_uploads(ctx, uploads):
    query = ctx.obj.query
    if len(uploads) > 0:
        query &= Q(upload_id__in=uploads)

    return query, proc.Upload.objects(query)


@uploads.command(help='List selected uploads')
@click.argument('UPLOADS', nargs=-1)
@click.pass_context
def ls(ctx, uploads):
    _, uploads = query_uploads(ctx, uploads)

    print('%d uploads selected, showing no more than first 10' % uploads.count())
    print(tabulate(
        [
            [upload.upload_id, upload.name, upload.user_id, upload.process_status, upload.published]
            for upload in uploads[:10]],
        headers=['id', 'name', 'user', 'status', 'published']))


@uploads.command(help='Change the owner of the upload and all its calcs.')
@click.argument('USER', nargs=1)
@click.argument('UPLOADS', nargs=-1)
@click.pass_context
def chown(ctx, user, uploads):
    infrastructure.setup_repository_db()
    _, uploads = query_uploads(ctx, uploads)

    print('%d uploads selected, changing its owner ...' % uploads.count())

    user_id = user
    user = coe_repo.User.from_user_id(int(user_id))

    for upload in uploads:
        upload.user_id = user_id
        upload_with_metadata = upload.to_upload_with_metadata()
        calcs = upload_with_metadata.calcs

        def create_update(calc):
            return UpdateOne(
                {'_id': calc.calc_id},
                {'$set': {'metadata.uploader': user.to_popo()}})

        proc.Calc._get_collection().bulk_write([create_update(calc) for calc in calcs])
        upload.save()

        upload_with_metadata = upload.to_upload_with_metadata()
        calcs = upload_with_metadata.calcs
        search.publish(calcs)


@uploads.command(help='(Re-)index all calcs of the given uploads.')
@click.argument('UPLOADS', nargs=-1)
@click.pass_context
def index(ctx, uploads):
    infrastructure.setup_repository_db()
    _, uploads = query_uploads(ctx, uploads)
    uploads_count = uploads.count()

    print('%d uploads selected, indexing ...' % uploads_count)

    i, failed = 0, 0
    for upload in uploads:
        upload_with_metadata = upload.to_upload_with_metadata()
        calcs = upload_with_metadata.calcs
        failed += search.index_all(calcs)
        i += 1

        print('   indexed %d of %d uploads, failed to index %d entries' % (i, uploads_count, failed))


@uploads.command(help='Delete selected upload')
@click.argument('UPLOADS', nargs=-1)
@click.option('--with-coe-repo', help='Also attempt to delete from repository db', is_flag=True)
@click.option('--skip-es', help='Keep the elastic index version of the data.', is_flag=True)
@click.option('--skip-mongo', help='Keep uploads and calcs in mongo.', is_flag=True)
@click.option('--skip-files', help='Keep all related files.', is_flag=True)
@click.pass_context
def rm(ctx, uploads, with_coe_repo, skip_es, skip_mongo, skip_files):
    _, uploads = query_uploads(ctx, uploads)

    logger = utils.get_logger(__name__)
    print('%d uploads selected, deleting ...' % uploads.count())

    if with_coe_repo:
        from nomad import coe_repo
        infrastructure.setup_repository_db()

    for upload in uploads:
        # delete repository db entry
        if with_coe_repo:
            coe_repo.Upload.delete(upload.upload_id)

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
                    logger.error('could not delete files', exc_info=e)
                    break

        # delete mongo
        if not skip_mongo:
            proc.Calc.objects(upload_id=upload.upload_id).delete()
            upload.delete()


@uploads.command(help='Reprocess selected uploads.')
@click.argument('UPLOADS', nargs=-1)
@click.option('--parallel', default=1, type=int, help='Use the given amount of parallel processes. Default is 1.')
@click.pass_context
def re_process(ctx, uploads, parallel: int):
    _, uploads = query_uploads(ctx, uploads)
    uploads_count = uploads.count()
    uploads = list(uploads)  # copy the whole mongo query set to avoid cursor timeouts

    cv = threading.Condition()
    threads: List[threading.Thread] = []

    state = dict(
        completed_count=0,
        skipped_count=0,
        available_threads_count=parallel)

    logger = utils.get_logger(__name__)

    print('%d uploads selected, re-processing ...' % uploads_count)

    def re_process_upload(upload: proc.Upload):
        logger.info('re-processing started', upload_id=upload.upload_id)

        completed = False
        if upload.process_running:
            logger.warn(
                'cannot trigger re-process, since the upload is already/still processing',
                current_process=upload.current_process,
                current_task=upload.current_task, upload_id=upload.upload_id)

        else:
            upload.reset()
            upload.re_process_upload()
            upload.block_until_complete(interval=.5)

            if upload.tasks_status == proc.FAILURE:
                logger.info('re-processing with failure', upload_id=upload.upload_id)

            completed = True

            logger.info('re-processing complete', upload_id=upload.upload_id)

        with cv:
            state['completed_count'] += 1 if completed else 0
            state['skipped_count'] += 1 if not completed else 0
            state['available_threads_count'] += 1

            print(
                '   re-processed %s and skipped %s of %s uploads' %
                (state['completed_count'], state['skipped_count'], uploads_count))

            cv.notify()

    for upload in uploads:
        with cv:
            cv.wait_for(lambda: state['available_threads_count'] > 0)
            state['available_threads_count'] -= 1
            thread = threading.Thread(target=lambda: re_process_upload(upload))
            threads.append(thread)
            thread.start()

    for thread in threads:
        thread.join()


@uploads.command(help='Attempt to abort the processing of uploads.')
@click.argument('UPLOADS', nargs=-1)
@click.option('--calcs', is_flag=True, help='Only stop calculation processing.')
@click.option('--kill', is_flag=True, help='Use the kill signal and force task failure.')
@click.option('--no-celery', is_flag=True, help='Do not attempt to stop the actual celery tasks')
@click.pass_context
def stop(ctx, uploads, calcs: bool, kill: bool, no_celery: bool):
    infrastructure.setup_repository_db()
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

    running_query = query & (Q(process_status=proc.PROCESS_RUNNING) | Q(tasks_status=proc.RUNNING))
    stop_all(proc.Calc.objects(running_query))
    if not calcs:
        stop_all(proc.Upload.objects(running_query))
