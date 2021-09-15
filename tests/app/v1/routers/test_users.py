
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

def test_me(client, test_user_auth):
    response = client.get('users/me', headers=test_user_auth)
    assert response.status_code == 200


def test_me_auth_required(client):
    response = client.get('users/me')
    assert response.status_code == 401


def test_me_auth_bad_token(client):
    response = client.get('users/me', headers={'Authentication': 'Bearer NOTATOKEN'})
    assert response.status_code == 401


def test_invite(client, test_user_auth, no_warn):
    rv = client.put(
        'users/invite', headers=test_user_auth, json={
            'first_name': 'John',
            'last_name': 'Doe',
            'affiliation': 'Affiliation',
            'email': 'john.doe@affiliation.edu'
        })
    assert rv.status_code == 200
    data = rv.json()
    keys = data.keys()
    required_keys = ['name', 'email', 'user_id']
    assert all(key in keys for key in required_keys)


def test_users(client):
    rv = client.get('users?prefix=Sheldon')
    assert rv.status_code == 200

    data = rv.json()
    assert len(data['data']) == 1
    user = data['data'][0]

    for key in ['name', 'user_id']:
        assert key in user

    for value in user.values():
        assert value is not None

    assert 'email' not in user
