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

import os.path
import os
import shutil
import sys
import click
import asyncio
from concurrent.futures import ProcessPoolExecutor
from mongoengine import Q
from elasticsearch_dsl import Search, A

from nomad import config, infrastructure, processing, utils, files, search

from .main import cli


@cli.group(help='Processing related functions')
def proc():
    pass


@proc.command(help='List processing tasks')
def ls():
    infrastructure.setup_logging()
    infrastructure.setup_mongo()

    def ls(query):
        for proc in query:
            print(proc)

    query = Q(process_status=processing.PROCESS_RUNNING) | Q(tasks_status=processing.RUNNING)

    ls(processing.Calc.objects(query))
    ls(processing.Upload.objects(query))


@proc.command(help='Remove uploads')
@click.argument('upload-ids', nargs=-1)
@click.option('--still-processing', is_flag=True, help='Target all uploads that are still processing')
@click.option('--stop-processing', is_flag=True, help='Attempt to stop processing on to delete uploads')
def rm(upload_ids, still_processing, stop_processing):
    infrastructure.setup_logging()
    infrastructure.setup_mongo()
    infrastructure.setup_elastic()

    upload_ids = list(upload_ids)

    if still_processing:
        query = Q(process_status=processing.PROCESS_RUNNING) | Q(tasks_status=processing.RUNNING)
        for upload in processing.Upload.objects(query):
            upload_ids.append(upload.upload_id)

    logger = utils.get_logger(__name__)
    logger.info('start deleting uploads', count=len(upload_ids))

    print('Will delete: ')
    for upload_id in upload_ids:
        print('   %s' % upload_id)
    input("Press Enter to continue...")

    for upload_id in upload_ids:
        logger = logger.bind(upload_id=upload_id)

        # stop processing
        if stop_processing:
            upload = processing.Upload.objects(upload_id=upload_id).first()
            try:
                processing.app.control.revoke(upload.celery_task_id, terminate=True, signal='SIGKILL')
            except Exception as e:
                logger.debug(
                    'could not revoke celery task', exc_info=e, celery_task_id=upload.celery_task_id)
            query = Q(upload_id=upload_id)
            query &= Q(process_status=processing.PROCESS_RUNNING) | Q(tasks_status=processing.RUNNING)
            for calc in processing.Calc.objects(query):
                try:
                    processing.app.control.revoke(calc.celery_task_id, terminate=True, signal='SIGKILL')
                except Exception as e:
                    logger.debug(
                        'could not revoke celery task', exc_info=e,
                        celery_task_id=calc.celery_task_id, calc_id=calc.calc_id)

                calc.fail('process terminate via nomad cli')
                calc.process_status = processing.PROCESS_COMPLETED
                calc.on_process_complete(None)
                calc.save()

            upload.fail('process terminate via nomad cli')
            upload.process_status = processing.PROCESS_COMPLETED
            upload.on_process_complete(None)
            upload.save()

            logger.info('stopped upload processing')

        # delete elastic
        search.delete_upload(upload_id=upload_id)

        # delete files
        for _ in range(0, 2):
            upload_files = files.UploadFiles.get(upload_id=upload_id)

            try:
                if upload_files is not None:
                    upload_files.delete()
            except Exception as e:
                logger.error('could not delete files', exc_info=e)
                break

        # delete mongo
        processing.Calc.objects(upload_id=upload_id).delete()
        processing.Upload.objects(upload_id=upload_id).delete()


@proc.command(help='Stop all running processing')
@click.option('--calcs', is_flag=True, help='Only stop calculation processing')
@click.option('--kill', is_flag=True, help='Use the kill signal and force task failure')
def stop_all(calcs: bool, kill: bool):
    infrastructure.setup_logging()
    infrastructure.setup_mongo()

    logger = utils.get_logger(__name__)

    def stop_all(query):
        for proc in query:
            logger_kwargs = dict(upload_id=proc.upload_id)
            if isinstance(proc, processing.Calc):
                logger_kwargs.update(calc_id=proc.calc_id)

            logger.info(
                'send terminate celery task', celery_task_id=proc.celery_task_id,
                kill=kill, **logger_kwargs)

            kwargs = {}
            if kill:
                kwargs.update(signal='SIGKILL')
            try:
                processing.app.control.revoke(proc.celery_task_id, terminate=True, **kwargs)
            except Exception as e:
                logger.warning(
                    'could not revoke celery task', exc_info=e,
                    celery_task_id=proc.celery_task_id, **logger_kwargs)
            if kill:
                logger.info(
                    'fail proc', celery_task_id=proc.celery_task_id, kill=kill,
                    **logger_kwargs)

                proc.fail('process terminate via nomad cli')
                proc.process_status = processing.PROCESS_COMPLETED
                proc.on_process_complete(None)
                proc.save()

    query = Q(process_status=processing.PROCESS_RUNNING) | Q(tasks_status=processing.RUNNING)

    stop_all(processing.Calc.objects(query))
    if not calcs:
        stop_all(processing.Upload.objects(query))


@cli.command(help='Removes all public (non staging) files that do not belong to an upload.')
@click.option('--dry', is_flag=True, help='Dont delete anything, just check.')
def clean_public(dry):
    infrastructure.setup_logging()
    infrastructure.setup_mongo()

    def uploads():
        for prefix in os.listdir(config.fs.public):
            for upload in os.listdir(os.path.join(config.fs.public, prefix)):
                yield upload, os.path.join(config.fs.public, prefix, upload)

    to_delete = list(
        path for upload, path in uploads()
        if processing.Upload.objects(upload_id=upload).first() is None)

    input('Will delete %d uploads. Press any key to continue ...' % len(to_delete))

    if not dry:
        for path in to_delete:
            shutil.rmtree(path)


@cli.command(help='Removes data that does not belong to an upload.')
@click.option('--dry', is_flag=True, help='Dont delete anything, just check.')
def clean_es(dry):
    infrastructure.setup_logging()
    infrastructure.setup_mongo()
    infrastructure.setup_elastic()

    def uploads():
        search = Search(index=config.elastic.index_name)
        search.aggs.bucket('uploads', A('terms', field='upload_id', size=12000))

        response = search.execute()

        print('Found %d uploads in ES.' % len(response.aggregations.uploads.buckets))
        for bucket in response.aggregations.uploads.buckets:
            yield bucket.key, bucket.doc_count

    to_delete = list(
        (upload, calcs) for upload, calcs in uploads()
        if processing.Upload.objects(upload_id=upload).first() is None)

    calcs = 0
    for _, upload_calcs in to_delete:
        calcs += upload_calcs

    input(
        'Will delete %d calcs in %d uploads from ES. Press any key to continue ...' %
        (calcs, len(to_delete)))

    if not dry:
        for upload, _ in to_delete:
            Search(index=config.elastic.index_name).query('term', upload_id=upload).delete()


@cli.command(help='Attempts to reset the nomad.')
def reset():
    from .main import create_client
    create_client().admin.exec_reset_command().response()


@cli.group(help='Run a nomad service locally (outside docker).')
def run():
    pass


@run.command(help='Run the nomad development worker.')
def worker():
    run_worker()


@run.command(help='Run the nomad development api.')
@click.option('--debug', help='Does run flask in debug.', is_flag=True)
def api(debug: bool):
    run_api(debug=debug)


def run_api(**kwargs):
    config.service = 'api'
    from nomad import infrastructure
    from nomad.api.__main__ import run_dev_server
    infrastructure.setup()
    run_dev_server(port=8000, **kwargs)


def run_worker():
    config.service = 'worker'
    from nomad import processing
    processing.app.worker_main(['worker', '--loglevel=INFO', '-Q', 'celery,uploads,calcs'])


@run.command(help='Run both api and worker.')
def apiworker():
    executor = ProcessPoolExecutor(2)
    loop = asyncio.get_event_loop()
    loop.run_in_executor(executor, run_api)
    loop.run_in_executor(executor, run_worker)


@cli.command(help='Runs tests and linting. Useful before commit code.')
@click.option('--skip-tests', help='Do not test, just do code checks.', is_flag=True)
def qa(skip_tests: bool):
    os.chdir(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
    ret_code = 0
    if not skip_tests:
        click.echo('Run tests ...')
        ret_code += os.system('python -m pytest -svx tests')
    click.echo('Run code style checks ...')
    ret_code += os.system('python -m pycodestyle --ignore=E501,E701 nomad tests')
    click.echo('Run linter ...')
    ret_code += os.system('python -m pylint --load-plugins=pylint_mongoengine nomad tests')
    click.echo('Run static type checks ...')
    ret_code += os.system('python -m mypy --ignore-missing-imports --follow-imports=silent --no-strict-optional nomad tests')

    sys.exit(ret_code)
