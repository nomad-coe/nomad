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

import os.path
import time

from .api import Auth


def upload_file(
        file_path: str, auth: Auth,
        upload_name: str = None, local_path: bool = False, publish: bool = False):
    '''
    Upload a file to nomad. Will print status information to stdout.

    Arguments:
        file_path: path to the file, absolute or relative to current working directory
        auth: mandatory Auth instance
        upload_name: optional name, default is the file_path's basename
        local_path: allows to process data without upload, requires client to be run on the server
        publish: automatically publish after successful processing

    Returns: The upload_id if successful or None if not.
    '''
    from nomad.processing import ProcessStatus
    from nomad.client import api
    if local_path:
        response = api.post(
            'uploads',
            params=dict(local_path=os.path.abspath(file_path), upload_name=upload_name),
            headers={'Accept': 'application/json'},
            auth=auth)
        print('process offline: %s' % file_path)
    else:
        with open(file_path, 'rb') as f:
            response = api.post(
                'uploads',
                params=dict(upload_name=upload_name), data=f,
                headers={'Accept': 'application/json'},
                auth=auth)
    if response.status_code != 200:
        print('Could not create upload: %s' % response.text)
        return None
    upload = response.json()['data']

    while upload['process_status'] not in [ProcessStatus.SUCCESS, ProcessStatus.FAILURE]:
        response = api.get(f'uploads/{upload["upload_id"]}/entries', auth=auth)
        response_json = response.json()

        upload = response_json['upload']
        total = response_json['pagination']['total']
        failures = response_json['processing_failed']
        successes = response_json['processing_successful']

        ret = '\n' if upload['process_status'] in (ProcessStatus.SUCCESS, ProcessStatus.FAILURE) else '\r'

        print(
            'status: %s; process: %s; parsing: %d/%d/%d                %s' %
            (upload['process_status'], upload['current_process'], successes, failures, total, ret), end='')

        time.sleep(1)

    if upload['process_status'] == ProcessStatus.FAILURE:
        print('There have been errors:')
        for error in upload['errors']:
            print('    %s' % error)
    elif publish:
        response = api.post(f'uploads/{upload["upload_id"]}/action/publish', auth=auth)
        if response.status_code != 200:
            print('Could not publish upload: %s' % response.text)
            return None

    return upload['upload_id']
