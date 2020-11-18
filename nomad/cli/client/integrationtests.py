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

'''
A command that runs some example operations on a working nomad@FAIRDI installation
as a final integration test.
'''

import time
import os
import click

from .client import client


multi_code_example_file = 'tests/data/integration/multi_code_data.zip'
simple_example_file = 'tests/data/integration/examples_vasp.zip'


@client.command(help='Runs a few example operations as a test.')
@click.option(
    '--skip-parsers', is_flag=True,
    help='Skip extensive upload and parser tests.')
@click.option(
    '--skip-publish', is_flag=True,
    help='Skip publish the upload. Should not be done on an production environment.')
@click.option(
    '--skip-doi', is_flag=True,
    help='Skip assigning a doi to a dataset.')
@click.option(
    '--skip-mirror', is_flag=True,
    help='Skip get mirror tests.')
@click.pass_context
def integrationtests(ctx, skip_parsers, skip_publish, skip_doi, skip_mirror):
    from .client import create_client
    client = create_client()
    has_doi = False
    published = False

    print('get the upload command')
    command = client.uploads.get_upload_command().response().result.upload_command_with_name

    def get_upload(upload):
        upload = client.uploads.get_upload(
            upload_id=upload.upload_id, per_page=100).response().result

        while upload.tasks_running:
            time.sleep(0.3)
            upload = client.uploads.get_upload(
                upload_id=upload.upload_id, per_page=100).response().result

        return upload

    if not skip_parsers:
        print('upload multi code test data with curl')
        command = command.replace('<local_file>', multi_code_example_file)
        command = command.replace('<name>', 'integration_test_upload')
        command += ' -k'
        code = os.system(command)
        assert code == 0, 'curl command must be successful'
        uploads = client.uploads.get_uploads(name='integration_test_upload').response().result.results
        assert len(uploads) == 1, 'exactly one test upload must be on the server'
        upload = uploads[0]

        print('observe the upload process to be finished')
        upload = get_upload(upload)

        assert upload.tasks_status == 'SUCCESS'
        total = upload.calcs.pagination.total
        assert 100 > total > 0
        assert len(upload.calcs.results) == total

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
    total = upload.calcs.pagination.total
    assert total > 0
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

        query = dict(owner='staging', upload_id=[upload.upload_id])

        print('perform repo search on data')
        search = client.repo.search(per_page=100, **query).response().result
        assert search.pagination.total >= total
        assert len(search.results) <= search.pagination.total

        print('performing archive paginated search')
        result = client.archive.post_archive_query(payload={
            'pagination': {
                'page': 1,
                'per_page': 10
            },
            'query': query
        }).response().result
        assert len(result.results) > 0

        print('performing archive scrolled search')
        result = client.archive.post_archive_query(payload={
            'scroll': {
                'scroll': True
            },
            'query': query
        }).response().result
        assert len(result.results) > 0

        print('performing download')
        client.raw.raw_files_from_query(**query)

        if not skip_publish:
            print('publish upload')
            client.uploads.exec_upload_operation(
                upload_id=upload.upload_id,
                payload=dict(operation='publish')).response()

            while upload.process_running:
                upload = client.uploads.get_upload(
                    upload_id=upload.upload_id).response().result

            assert upload.tasks_status == 'SUCCESS', 'publish must be successful'
            published = True

        print('editing upload')
        dataset = 'test_dataset'
        actions = {
            'comment': {'value': 'Test comment'},
            'references': [{'value': 'http;//test_reference.com'}],
            'coauthors': [{'value': 'author1-id'}, {'value': 'author2-id'}],
            'shared_with': [{'value': 'author3-id'}],
            'datasets': [{'value': dataset}]}

        payload = dict(actions=actions, query=dict(upload_id=[upload.upload_id]))
        result = client.repo.edit_repo(payload=payload).response().result
        assert result.success
        assert client.datasets.get_dataset(name=dataset).response().result['name'] == dataset

        print('list datasets')
        result = client.datasets.list_datasets(page=1, per_page=10).response().result
        results = result.results
        assert len(results) > 0

        if not skip_doi and published:
            print('assigning a DOI')
            result = client.datasets.assign_doi(name=dataset).response().result
            doi = result.doi
            assert doi
            has_doi = True

        if not has_doi or ctx.obj.user == 'admin':
            print('deleting dataset')
            result = client.datasets.delete_dataset(name=dataset).response().result

        if not skip_mirror and ctx.obj.user == 'admin':
            print('getting upload mirror')
            # get_upload_mirror gives 404
            payload = dict(query=dict(upload_id=upload.upload_id))
            result = client.mirror.get_uploads_mirror(payload=payload).response().result
            assert len(result) == 1
            assert len(client.mirror.get_upload_mirror(upload_id=upload.upload_id).response().result.calcs) > 0

    finally:
        if not published or ctx.obj.user == 'admin':
            print('delete the upload again')
            client.uploads.delete_upload(upload_id=upload.upload_id).response()
            while upload.process_running:
                upload = client.uploads.get_upload(
                    upload_id=upload.upload_id).response().result
