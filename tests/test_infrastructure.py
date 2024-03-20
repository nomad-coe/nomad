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

from nomad.infrastructure import UserManagement
from tests.fixtures.users import test_user_uuid as create_test_user_uuid


@pytest.fixture(scope='function')
def user_management(api_v1):
    from nomad.infrastructure import OasisUserManagement

    return OasisUserManagement('users')


@pytest.mark.parametrize(
    'query,count',
    [
        pytest.param('Sheldon', 1, id='exists'),
        pytest.param('Does not exist $%&#', 0, id='does-not-exist'),
    ],
)
def test_search_user(user_management: UserManagement, query, count):
    users = user_management.search_user(query)
    assert len(users) == count


@pytest.mark.parametrize(
    'key,value',
    [
        pytest.param('username', 'scooper', id='username'),
        pytest.param('email', 'sheldon.cooper@nomad-coe.eu', id='email'),
        pytest.param('user_id', create_test_user_uuid(1), id='user_id'),
    ],
)
def test_get_user(user_management: UserManagement, key, value):
    user = user_management.get_user(**{key: value})
    assert user is not None
    assert getattr(user, key) == (value if key != 'email' else None)


def test_get_admin_user(monkeypatch, user_management: UserManagement):
    user = user_management.get_user(username='scooper')
    assert user is not None
    monkeypatch.setattr('nomad.config.services.admin_user_id', user.user_id)
    assert user.is_admin
