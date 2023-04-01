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
import pytest


# try:
#     from rich.pretty import pprint
# except ImportError:
#     def pprint(x):
#         print(x)


def assert_path_exists(path, response):
    for segment in path:
        response = response[segment]
    # empty container implies wrong id
    if isinstance(response, dict):
        response.pop('m_errors', None)
    if not response:
        raise KeyError


@pytest.mark.skip
@pytest.mark.parametrize('upload_id,entry_id,user,status_code', [
    pytest.param('id_embargo', 'id_embargo_1', 'test_user', 200, id='ok'),
    pytest.param('id_child_entries', 'id_child_entries_child1', 'test_user', 200, id='child-entry'),
    pytest.param('id_embargo', 'id_embargo_1', 'admin_user', 200, id='admin-access'),
    pytest.param('id_embargo', 'id_embargo_1', None, 401, id='no-credentials'),
    pytest.param('id_embargo', 'id_embargo_1', 'invalid', 401, id='invalid-credentials'),
    pytest.param('id_embargo', 'id_embargo_1', 'other_test_user', 404, id='no-access'),
    pytest.param('silly_value', 'id_embargo_1', 'test_user', 404, id='invalid-upload_id'),
    pytest.param('id_embargo', 'silly_value', 'test_user', 404, id='invalid-entry_id')
])
def test_graph_query(
        client, mongo_module, test_auth_dict, example_data,
        upload_id, entry_id, user, status_code
):
    user_auth, _ = test_auth_dict[user]
    response = client.post(
        'graph/query',
        json={'m_uploads': {upload_id: {'m_entries': {entry_id: '*'}}}},
        headers={'Accept': 'application/json'} | (user_auth if user_auth else {}))
    target_path = ('m_response', 'm_uploads', upload_id, 'm_entries', entry_id)
    if 200 == status_code:
        assert_path_exists(target_path, response.json())
    elif 401 == status_code:
        assert response.status_code == 401
    else:
        with pytest.raises(KeyError):
            assert_path_exists(target_path, response.json())
