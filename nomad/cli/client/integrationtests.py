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
@click.option(
    '--with-edit', is_flag=True,
    help='Edits upload i.e. by adding comment, references, coauthors, shared_with.')
@click.option(
    '--with-doi', is_flag=True,
    help='Assign a doi to a dataset.')
@click.option(
    '--with-mirror', is_flag=True,
    help='Get mirror uploads.')
@click.option(
    '--with-query', is_flag=True,
    help='Perform query operations.')
def integrationtests(with_publish, with_edit, with_doi, with_mirror, with_query):
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

    if with_edit:
        # add comment, references, coauthors, shared_with, datasets
        dataset = 'test_dataset'
        actions = {
            'comment': {'value': 'Test comment'},
            'references': [{'value': 'http;//test_reference.com'}],
            'coauthors': [{'value': 'author1-id'}, {'value': 'author2-id'}],
            'shared_with': [{'value': 'author3-id'}],
            'datasets': [{'value': dataset}]}

        payload = dict(actions=actions, query=dict(upload_id=upload.upload_id))
        print('editing upload')
        result = client.repo.edit_repo(payload=payload).response().result
        assert result.success

        # check if dataset was created
        page = 1
        found = False
        while not found:
            result = client.datasets.list_datasets(page=page, per_page=10).response().result
            results = result.results
            if len(results) == 0:
                break
            found = dataset in [res.name for res in results]
            page += 1
        assert found
        print('successfully created dataset')

        # assign a doi
        if with_doi:
            print('assigning a DOI')
            result = client.datasets.assign_doi(name=dataset).response().result
            doi = result.doi
            assert doi

        # delete dataset
        print('deleting dataset')
        result = client.datasets.delete_dataset(name=dataset).response().result

    if with_mirror:
        print('getting upload mirror')
        # get_upload_mirror gives 404
        payload = dict(query=dict())
        result = client.mirror.get_uploads_mirror(payload=payload).response().result
        assert len(result) > 0

    if with_query:
        # paginated search
        print('performing paginated search')
        result = client.archive.archive_query(upload_id=[upload.upload_id], page=1, per_page=10).response().result
        # why no result?

        # scrolled search
        print('performing scrolled search')
        result = client.archive.archive_query(upload_id=[upload.upload_id], scroll=True).response().result
