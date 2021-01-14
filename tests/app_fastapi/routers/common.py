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

from devtools import debug


def assert_response(response, status_code=None):
    ''' General assertions for status_code and error messages '''
    if status_code and response.status_code != status_code:
        try:
            debug(response.json())
        except Exception:
            pass

    if status_code is not None:
        if response.status_code != status_code and response.status_code == 422:
            print(response.json()['detail'])
        assert response.status_code == status_code

    if status_code == 422:
        response_json = response.json()
        details = response_json['detail']
        assert len(details) > 0
        for detail in details:
            assert 'loc' in detail
            assert 'msg' in detail
        return

    if 400 <= status_code < 500:
        response_json = response.json()
        assert 'detail' in response_json
