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

from celery import Celery, group, subtask
from celery.result import result_from_tuple
import logging
import time

import nomad.config as config
import nomad.files as files
from nomad.parsers import parsers, parser_dict

broker_url = 'pyamqp://%s:%s@localhost//' % (config.rabbitmq.user, config.rabbitmq.password)
backend_url = 'redis://localhost/0'
app = Celery('nomad.processing', backend=backend_url, broker=broker_url)
app.conf.update(
    accept_content=['pickle'],
    task_serializer='pickle',
    result_serializer='pickle',
)

LOGGER = logging.getLogger(__name__)


@app.task()
def open_upload(upload_id):
    try:
        upload = files.upload(upload_id)
        upload.open()
        return upload
    except Exception as e:
        return e


@app.task()
def find_mainfiles(upload):
    if isinstance(upload, Exception):
        return list()

    mainfile_specs = list()
    for filename in upload.filelist:
        for parser in parsers:
            if parser.is_mainfile(upload, filename):
                mainfile_specs.append((upload, filename, parser.name))

    return mainfile_specs


@app.task()
def close_upload(parse_results, upload_id):
    try:
        upload = files.upload(upload_id)
    except KeyError as e:
        return e

    upload.close()
    return parse_results


@app.task()
def parse(mainfile_spec):
    upload, mainfile, parser = mainfile_spec
    LOGGER.debug('Start parsing mainfile %s/%s with %s.' % (upload, mainfile, parser))
    parser_dict[parser].run(upload.get_path(mainfile))

    return True


@app.task()
def dmap(it, callback):
    callback = subtask(callback)
    return group(callback.clone([arg, ]) for arg in it)()


def start_process_upload(upload_id):
    """
    Starts the processing of uploaded data and returns the celery
    :class:`~celery.results.AsyncResults` instance in serialized *tuple* form.

    The serialized form is understood by all respective methods in this
    module. We use the serialized form to allow storage. Keep in mind
    that the sheer `task_id` is not enough, because it does not contain
    the parent tasks. See [third comment](https://github.com/celery/celery/issues/1328)
    for details.
    """
    parsing_workflow = (
        open_upload.s(upload_id) |
        find_mainfiles.s() |
        dmap.s(parse.s()) |
        close_upload.s(upload_id)
    )

    async_result = parsing_workflow.delay()
    return async_result.as_tuple()


def get_process_upload_state(async_result):
    """
    Extract the current state from the various tasks involved in upload processing.
    """
    async_result = result_from_tuple(async_result)

    close = async_result
    parse = close.parent
    find_mainfiles = parse.parent
    open_task = find_mainfiles.parent

    return {
        'open': open_task.state,
        'find_mainfiles': find_mainfiles.state,
        'parse': parse.state,
        'close': close.state
    }


if __name__ == '__main__':
    upload_id = 'examples_vasp.zip'

    task = start_process_upload(upload_id)

    result = None
    while(True):
        time.sleep(0.0001)
        new_result = get_process_upload_state(task)
        if result != new_result:
            result = new_result
            print(result)
            if result['close'] == 'SUCCESS' or result['close'] == 'FAILURE':
                break
