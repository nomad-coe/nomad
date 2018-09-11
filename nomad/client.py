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
Simple client library for the nomad api that allows to bulk upload files via shell command.
"""

import os.path
import subprocess
import shlex
import time
import sys
import requests
from requests.auth import HTTPBasicAuth

api_base = 'http://localhost/nomadxt/api'


def upload_file(file_path, name=None, user='other@gmail.com', pw='nomad'):
    """
    Upload a file to nomad.

    Arguments:
        file_path: Path to the file, absolute or relative to call directory.
        name: Optional name, default is the file_path's basename
    """
    auth = HTTPBasicAuth(user, pw)

    if name is None:
        name = os.path.basename(file_path)

    upload = requests.post('%s/uploads' % api_base, data={name: name}, auth=auth).json()

    upload_cmd = upload['upload_command']
    upload_cmd = upload_cmd.replace('your_file', file_path)

    subprocess.call(shlex.split(upload_cmd))

    print('File uploaded')

    while True:
        upload = requests.get('%s/uploads/%s' % (api_base, upload['upload_id']), auth=auth).json()
        status = upload['status']
        calcs_pagination = upload['calcs'].get('pagination')
        if calcs_pagination is None:
            total, successes, failures = 0, 0, 0
        else:
            total, successes, failures = (
                calcs_pagination[key] for key in ('total', 'successes', 'failures'))

        print(
            'status: %s; task: %s; parsing: %d/%d/%d' %
            (status, upload['current_task'], successes, failures, total))

        if status in ('SUCCESS', 'FAILURE'):
            break

        time.sleep(5)


if __name__ == '__main__':
    if len(sys.argv) > 3 or len(sys.argv) == 1:
        print('usage is: <client> filte_to_upload [upload_name]')
    else:
        upload_file(sys.argv[1], sys.argv[2] if len(sys.argv) == 3 else None)
