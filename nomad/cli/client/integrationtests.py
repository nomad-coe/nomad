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
import json

from nomad.client import api


def integrationtests(auth: api.Auth, skip_parsers: bool, skip_publish: bool, skip_doi: bool):
    multi_code_example_file = 'tests/data/integration/multi_code_data.zip'
    simple_example_file = 'tests/data/integration/examples_vasp.zip'
    has_doi = False
    published = False

    print('get the upload command')
    response = api.get('uploads/command-examples', auth=auth)
    assert response.status_code == 200, response.text
    command = response.json()['upload_command_with_name']

    def get_upload(upload):
        first = True
        while first or upload['process_running']:
            first = False
            response = api.get(f'uploads/{upload["upload_id"]}', auth=auth)
            if response.status_code == 404:
                return None
            assert response.status_code == 200, response.text
            upload = response.json()['data']
            time.sleep(0.3)

        return upload

    response = api.get('uploads', params=dict(upload_name='integration_test_upload'), auth=auth)
    assert response.status_code == 200, response.text
    uploads = response.json()['data']
    assert len(uploads) == 0, 'the test upload must not exist before'

    if not skip_parsers:
        print('upload multi code test data with curl')
        command = command.replace('<local_file>', multi_code_example_file)
        command = command.replace('<name>', 'integration_test_upload')
        command += ' -k'
        code = os.system(command)
        assert code == 0, 'curl command must be successful'
        response = api.get('uploads', params=dict(upload_name='integration_test_upload'), auth=auth)
        assert response.status_code == 200, response.text
        response_json = response.json()
        assert len(response_json['data']) == 1, 'exactly one test upload must be on the server'
        upload = response_json['data'][0]

        print('observe the upload process to be finished')
        upload = get_upload(upload)

        assert upload['process_status'] == 'SUCCESS'

        print('delete the upload again')
        upload = api.delete(f'uploads/{upload["upload_id"]}', auth=auth).json()['data']
        upload = get_upload(upload)

    print('upload simple data with API')
    with open(simple_example_file, 'rb') as f:
        response = api.post(
            'uploads', files=dict(file=f), params=dict(upload_name='integration_test_upload'),
            auth=auth, headers={'Accept': 'application/json'})
        assert response.status_code == 200, response.text
        upload = response.json()['data']

    print('observe the upload process to be finished')
    upload = get_upload(upload)
    response = api.get(f'uploads/{upload["upload_id"]}/entries', auth=auth)
    assert response.status_code == 200, response.text
    entries = response.json()['data']
    assert upload['entries'] == len(entries)

    try:
        print('get repo data')
        for entry in entries:
            response = api.get(f'entries/{entry["entry_id"]}', auth=auth)
            assert response.status_code == 200, response.text
            entry_metadata = response.json()['data']
            entry_metadata['entry_id'] == entry['entry_id']

        print('get archive data')
        for entry in entries:
            api.get(f'entries/{entry["entry_id"]}/archive/download', auth=auth)
            assert response.status_code == 200, response.text

        print('get archive logs')
        for entry in entries:
            response = api.post(
                f'entries/{entry["entry_id"]}/archive/query',
                data=json.dumps({
                    'required': {
                        'processing_logs': '*'
                    }
                }), auth=auth)
            assert response.status_code == 200, response.text
            response_archive_keys = [
                key for key in response.json()['data']['archive'].keys()
                if not key.startswith('m_')]
            assert response_archive_keys == ['processing_logs'], \
                f'Archive keys {response_archive_keys} should only contain processing_logs'

        query_request_params = dict(
            owner='staging',
            query={
                'upload_id': upload['upload_id']
            })

        print('perform repo search on data')
        response = api.post('entries/query', data=json.dumps(query_request_params), auth=auth)
        assert response.status_code == 200, response.text
        response_json = response.json()
        assert response_json['pagination']['total'] == 2
        assert response_json['pagination']['total'] == len(response_json['data'])

        print('performing archive paginated search')
        response = api.post('entries/archive/query', data=json.dumps(dict(
            pagination=dict(page_size=1, page_offset=1),
            **query_request_params)), auth=auth)
        assert response.status_code == 200, response.text
        response_json = response.json()
        assert response_json['pagination']['total'] == 2
        assert len(response_json['data']) == 1

        print('performing archive scrolled search')
        response = api.post('entries/archive/query', data=json.dumps(dict(
            pagination=dict(page_size=1),
            **query_request_params)), auth=auth)
        response_json = response.json()
        response = api.post('entries/archive/query', data=json.dumps(dict(
            pagination=dict(page_size=1, page_after_value=response_json['pagination']['next_page_after_value']),
            **query_request_params)), auth=auth)
        assert response.status_code == 200, response.text
        response_json = response.json()
        assert response_json['pagination']['total'] == 2
        assert len(response_json['data']) == 1

        print('performing download')
        response = api.get(
            'entries/raw',
            params=dict(upload_id=upload['upload_id'], owner='visible'), auth=auth)
        assert response.status_code == 200, response.text

        if not skip_publish:
            print('publish upload')
            api.post(f'uploads/{upload["upload_id"]}/action/publish')

            upload = get_upload(upload)
            assert upload['process_status'] == 'SUCCESS', 'publish must be successful'
            published = True

        print('editing upload')
        response = api.get('users', params=dict(prefix='Markus Scheidgen'))
        assert response.status_code == 200, response.text
        user = response.json()['data'][0]
        dataset = 'test_dataset'
        actions = {
            'comment': {'value': 'Test comment'},
            'references': [{'value': 'http;//test_reference.com'}],
            'entry_coauthors': [{'value': user['user_id']}],
            'datasets': [{'value': dataset}]}

        response = api.post(
            'entries/edit_v0',
            data=json.dumps(dict(actions=actions, **query_request_params)),
            auth=auth)
        assert response.status_code == 200, response.text

        print('list datasets')
        response = api.get('datasets', auth=auth, params=dict(dataset_name=dataset))
        assert response.status_code == 200, response.text
        response_json = response.json()
        assert len(response_json['data']) == 1, response.text
        dataset_id = response_json['data'][0]['dataset_id']

        if not skip_doi and published:
            print('assigning a DOI')
            response = api.post(f'datasets/{dataset_id}/action/doi', auth=auth)
            assert response.status_code == 200, response.text
            has_doi = True

        if not has_doi or auth.user == 'admin':
            print('deleting dataset')
            response = api.delete(f'datasets/{dataset_id}', auth=auth)
            assert response.status_code == 200, response.text

    finally:
        if not published or auth.user == 'admin':
            print('delete the upload again')
            upload = api.delete(f'uploads/{upload["upload_id"]}', auth=auth).json()['data']
            assert get_upload(upload) is None
