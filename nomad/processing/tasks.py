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
from celery import Task, chord, group
from celery.canvas import Signature
from datetime import datetime

from nomad import files, utils, search
from nomad.parsing import parsers, parser_dict
from nomad.normalizing import normalizers
import nomad.patch  # pylint: disable=unused-import

from nomad.processing.app import app
from nomad.processing.state import UploadProc, CalcProc

import json


@app.task(bind=True, name='extracting')
def extracting_task(task: Task, proc: UploadProc) -> UploadProc:
    logger = utils.get_logger(__name__, task=task.name, upload_id=proc.upload_id)
    if not proc.continue_with(task.name):
        return proc

    try:
        upload = files.Upload(proc.upload_id)
        upload.open()
    except KeyError as e:
        logger.debug('Process request for non existing upload')
        proc.fail(e)
        return proc
    except files.UploadError as e:
        logger.debug('Could not open upload, %s' % e)
        proc.fail(e)
        return proc
    except Exception as e:
        logger.error('Unknown exception %s', exc_info=e)
        proc.fail(e)
        return proc

    logger.debug('Upload opened')

    try:
        proc.upload_hash = upload.hash()
    except files.UploadError as e:
        logger.error('Could not create upload hash', exc_info=e)
        proc.fail(e)
        return proc

    try:
        # TODO: deal with multiple possible parser specs
        for filename in upload.filelist:
            for parser in parsers:
                if parser.is_mainfile(upload, filename):
                    tmp_mainfile = upload.get_path(filename)
                    calc_processing = CalcProc(
                        proc.upload_hash, filename, parser.name, tmp_mainfile)

                    proc.calc_procs.append(calc_processing)
    except files.UploadError as e:
        logger.warn('Could find parse specs in open upload', exc_info=e)
        proc.fail(e)
        return proc

    return proc


@app.task(bind=True, name='cleanup')
def cleanup_task(task, calc_procs: List[CalcProc], upload_proc: UploadProc) -> UploadProc:
    logger = utils.get_logger(__name__, task=task.name, upload_id=upload_proc.upload_id)
    if upload_proc.continue_with(task.name):
        try:
            upload = files.Upload(upload_proc.upload_id)
        except KeyError as e:
            logger.warn('Upload does not exist')
            upload_proc.fail(e)
            return upload_proc

        try:
            upload.close()
        except Exception as e:
            logger.error('Could not close upload', exc_info=e)
            upload_proc.fail(e)
            return upload_proc

        logger.debug('Closed upload')
        upload_proc.success()

    return upload_proc


def _report_progress(task, dct):
    if not task.request.called_directly:
        task.update_state(state='PROGRESS', meta=dct)


@app.task(bind=True, name='parse_all')
def parse_all_task(task: Task, upload_proc: UploadProc, cleanup: Signature) -> UploadProc:
    if not upload_proc.continue_with(task.name):
        chord([])(cleanup.clone(args=(upload_proc,)))
        return upload_proc

    # prepare the group of parallel calc processings
    parses = group(parse_task.s(calc_proc) for calc_proc in upload_proc.calc_procs)

    # save the calc processing task ids to the overall processing
    i = 0
    for child in parses.freeze().children:
        upload_proc.calc_procs[i].celery_task_id = child.task_id
        i = i + 1

    # initiate the chord that runs calc processings first, and close_upload afterwards
    chord(parses)(cleanup.clone(args=(upload_proc,)))

    return upload_proc


@app.task(bind=True, name='parse')
def parse_task(self, proc: CalcProc) -> CalcProc:
    assert proc.upload_hash is not None

    upload_hash, parser, mainfile = proc.upload_hash, proc.parser_name, proc.mainfile
    logger = utils.get_logger(__name__, upload_hash=upload_hash, mainfile=mainfile)

    # parsing
    proc.continue_with(parser)
    logger.debug('Start parsing with %s' % parser)
    try:
        parser_backend = parser_dict[parser].run(proc.tmp_mainfile)
        if parser_backend.status[0] != 'ParseSuccess':
            proc.fail(parser_backend.status[1])
            return proc
    except Exception as e:
        logger.warn('Exception wile parsing', exc_info=e)
        proc.fail(e)
        return proc
    _report_progress(self, proc)

    # normalization
    for normalizer in normalizers:
        normalizer_name = normalizer.__name__
        proc.continue_with(normalizer_name)
        logger.debug('Start %s.' % normalizer)
        try:
            normalizer(parser_backend).normalize()
            if parser_backend.status[0] != 'ParseSuccess':
                proc.fail(parser_backend.status[1])
                return proc
        except Exception as e:
            logger.warn('Exception wile normalizing with %s' % normalizer, exc_info=e)
            proc.fail(e)
            return proc
    _report_progress(self, proc)

    # update search
    proc.continue_with('index')
    try:
        search.Calc.add_from_backend(
            parser_backend,
            upload_hash=upload_hash,
            calc_hash=proc.calc_hash,
            mainfile=mainfile,
            upload_time=datetime.now())
    except Exception as e:
        logger.error('Could not index', exc_info=e)
        proc.fail(e)
        return proc
    _report_progress(self, proc)

    # calc data persistence
    proc.continue_with('archiving')
    archive_id = proc.archive_id
    try:
        with files.write_archive_json(archive_id) as out:
            parser_backend.write_json(out, pretty=True)
    except Exception as e:
        logger.error('Could not write archive', exc_info=e)
        proc.fail(e)
        return proc

    logger.debug('Completed')

    proc.success()
    return proc
