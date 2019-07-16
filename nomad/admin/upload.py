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

from nomad import processing as proc, infrastructure, utils, search, files
from .__main__ import cli


uploads = None
query = None


@cli.group(help='Upload related commands')
@click.option('--upload', help='Select upload of with given id', type=str)
@click.option('--user', help='Select uploads of user with given id', type=str)
@click.option('--staging', help='Select only uploads in staging', is_flag=True)
@click.option('--processing', help='Select only processing uploads', is_flag=True)
def upload(upload: str, user: str, staging: bool, processing: bool):
    infrastructure.setup_mongo()
    infrastructure.setup_elastic()

    global query
    query = Q()
    if upload is not None:
        query &= Q(upload_id=upload)
    if user is not None:
        query &= Q(user_id=user)
    if staging:
        query &= Q(published=False)
    if processing:
        query &= Q(process_status=proc.PROCESS_RUNNING) | Q(tasks_status=proc.RUNNING)

    global uploads
    uploads = proc.Upload.objects(query)


@upload.command(help='List selected uploads')
def ls():
    print('%d uploads selected, showing no more than first 10' % uploads.count())
    print(tabulate(
        [
            [upload.upload_id, upload.name, upload.user_id, upload.process_status, upload.published]
            for upload in uploads[:10]],
        headers=['id', 'name', 'user', 'status', 'published']))


@upload.command(help='Delete selected upload')
@click.option('--with-coe-repo', help='Also attempt to delete from repository db', is_flag=True)
@click.option('--skip-es', help='Keep the elastic index version of the data.', is_flag=True)
@click.option('--skip-mongo', help='Keep uploads and calcs in mongo.', is_flag=True)
@click.option('--skip-files', help='Keep all related files.', is_flag=True)
def rm(with_coe_repo, skip_es, skip_mongo, skip_files):
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


@upload.command(help='Attempt to abort the processing of uploads.')
@click.option('--calcs', is_flag=True, help='Only stop calculation processing.')
@click.option('--kill', is_flag=True, help='Use the kill signal and force task failure.')
def stop(calcs: bool, kill: bool):
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
