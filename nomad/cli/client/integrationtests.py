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

"""
A command that runs some example operations on a working nomad@FAIRDI installation
as a final integration test.
"""

import time
import os
import click

from .client import client


multi_code_example_file = 'tests/data/integration/multi_code_data.zip'
simple_example_file = 'tests/data/integration/examples_vasp.zip'


@client.command(help='Runs a few example operations as a test.')
@click.option(
    '--with-publish', is_flag=True,
    help='Also publish the upload. Should not be done on an production environment.')
def integrationtests(with_publish):
    from .client import create_client
    client = create_client()

    print('get the upload command')
    command = client.uploads.get_upload_command().response().result.upload_command_with_name

    print('upload multi code test data with curl')
    command = command.replace('<local_file>', multi_code_example_file)
    command = command.replace('<name>', 'integration_test_upload')
    command += ' -k'
    code = os.system(command)
    assert code == 0, 'curl command must be successful'
    uploads = client.uploads.get_uploads(name='integration_test_upload').response().result.results
    assert len(uploads) == 1, 'exactly one test upload must be on the server'
    upload = uploads[0]

    def get_upload(upload):
        upload = client.uploads.get_upload(
            upload_id=upload.upload_id, per_page=100).response().result

        while upload.tasks_running:
            time.sleep(0.3)
            upload = client.uploads.get_upload(
                upload_id=upload.upload_id, per_page=100).response().result

        return upload

    print('observe the upload process to be finished')
    upload = get_upload(upload)

    assert upload.tasks_status == 'SUCCESS'
    total = upload.calcs.pagination.total
    assert 100 > total > 0
    assert len(upload.calcs.results) == total

    try:
        print('get repo data')
        for calc in upload.calcs.results:
            repo = client.repo.get_repo_calc(
                upload_id=upload.upload_id, calc_id=calc.calc_id).response().result
            repo['calc_id'] == calc.calc_id

        print('get archive data')
        for calc in upload.calcs.results:
            client.archive.get_archive_calc(
                upload_id=upload.upload_id, calc_id=calc.calc_id).response()

        print('get archive logs')
        for calc in upload.calcs.results:
            client.archive.get_archive_logs(
                upload_id=upload.upload_id, calc_id=calc.calc_id).response()

        print('perform search on data')
        search = client.repo.search(owner='staging', per_page=100).response().result
        assert search.pagination.total >= total
        assert len(search.results) <= search.pagination.total

    finally:
        print('delete the upload again')
        client.uploads.delete_upload(upload_id=upload.upload_id).response()
        while upload.process_running:
            upload = client.uploads.get_upload(
                upload_id=upload.upload_id).response().result

    print('upload simple data with API')
    with open(simple_example_file, 'rb') as f:
        upload = client.uploads.upload(
            name='integration test upload', file=f).response().result

    print('observe the upload process to be finished')
    upload = get_upload(upload)

    if with_publish:
        try:
            print('publish upload')
            client.uploads.exec_upload_operation(
                upload_id=upload.upload_id,
                payload=dict(operation='publish')).response()

            while upload.process_running:
                upload = client.uploads.get_upload(
                    upload_id=upload.upload_id).response().result

            assert upload.tasks_status == 'SUCCESS', 'publish must be successful'

        except Exception as e:
            print('delete the upload after exception')
            client.uploads.delete_upload(upload_id=upload.upload_id).response()
            while upload.process_running:
                upload = client.uploads.get_upload(
                    upload_id=upload.upload_id).response().result
            raise e

    else:
        print('delete the upload again')
        client.uploads.delete_upload(upload_id=upload.upload_id).response()
        while upload.process_running:
            upload = client.uploads.get_upload(
                upload_id=upload.upload_id).response().result
