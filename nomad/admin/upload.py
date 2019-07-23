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
from tabulate import tabulate
from mongoengine import Q
from pymongo import UpdateOne

from nomad import processing as proc, config, infrastructure, utils, search, files, coe_repo
from .__main__ import cli


@cli.group(help='Upload related commands')
@click.option('--user', help='Select uploads of user with given id', type=str)
@click.option('--staging', help='Select only uploads in staging', is_flag=True)
@click.option('--processing', help='Select only processing uploads', is_flag=True)
@click.option('--outdated', help='Select published uploads with older nomad version', is_flag=True)
@click.pass_context
def upload(ctx, user: str, staging: bool, processing: bool, outdated: bool):
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
        uploads = proc.Calc.objects(
            {'metadata.nomad_version__ne': config.version}).distinct(field='upload_id')
        query &= Q(uploads__in=uploads)

    ctx.obj.query = query
    ctx.obj.uploads = proc.Upload.objects(query)


def query_uploads(ctx, uploads):
    query = ctx.obj.query
    if len(uploads) > 0:
        query &= Q(upload_id__in=uploads)

    return query, proc.Upload.objects(query)


@upload.command(help='List selected uploads')
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


@upload.command(help='Change the owner of the upload and all its calcs.')
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


@upload.command(help='Delete selected upload')
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


@upload.command(help='Reprocess selected uploads.')
@click.argument('UPLOADS', nargs=-1)
@click.pass_context
def re_process(ctx, uploads):
    _, uploads = query_uploads(ctx, uploads)

    logger = utils.get_logger(__name__)
    print('%d uploads selected, re-processing ...' % uploads.count())

    def re_process_upload(upload_id: str):
        logger.info('re-processing started', upload_id=upload_id)
        upload = proc.Upload.objects(upload_id=upload_id).first()
        if upload is None:
            logger.error('upload for re-processing does not exist', upload_id=upload_id)

        upload.re_process_upload()
        upload.block_until_complete(interval=.1)

        logger.info('re-processing complete', upload_id=upload_id)

    count = 0
    for upload_id in uploads:
        re_process_upload(upload_id)
        count += 1
        print('   re-processed %s of %s uploads' % (count, len(uploads)))


@upload.command(help='Attempt to abort the processing of uploads.')
@click.argument('UPLOADS', nargs=-1)
@click.option('--calcs', is_flag=True, help='Only stop calculation processing.')
@click.option('--kill', is_flag=True, help='Use the kill signal and force task failure.')
@click.pass_context
def stop(ctx, uploads, calcs: bool, kill: bool):
    query, _ = query_uploads(ctx, uploads)

    logger = utils.get_logger(__name__)

    def stop_all(query):
        for proc in query:
            logger_kwargs = dict(upload_id=proc.upload_id)
            if isinstance(proc, proc.Calc):
                logger_kwargs.update(calc_id=proc.calc_id)

            logger.info(
                'send terminate celery task', celery_task_id=proc.celery_task_id,
                kill=kill, **logger_kwargs)

            kwargs = {}
            if kill:
                kwargs.update(signal='SIGKILL')
            try:
                proc.app.control.revoke(proc.celery_task_id, terminate=True, **kwargs)
            except Exception as e:
                logger.warning(
                    'could not revoke celery task', exc_info=e,
                    celery_task_id=proc.celery_task_id, **logger_kwargs)
            if kill:
                logger.info(
                    'fail proc', celery_task_id=proc.celery_task_id, kill=kill,
                    **logger_kwargs)

                proc.fail('process terminate via nomad cli')
                proc.process_status = proc.PROCESS_COMPLETED
                proc.on_process_complete(None)
                proc.save()

    stop_all(proc.Calc.objects(query))
    if not calcs:
        stop_all(proc.Upload.objects(query))
