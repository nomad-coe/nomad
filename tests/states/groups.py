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

from collections.abc import Iterable

from nomad.app.v1.models.groups import UserGroupMember
from nomad.auth import user_management
from nomad.mongo.groups import create_mongo_user_group_for_test, get_mongo_user_group
from tests.utils import list_without


def _create(group_id: str, group_name: str, owner: str, members: Iterable[str]):
    no_owner_ids = list_without(members, owner)
    members = [owner] + no_owner_ids
    members_info = [UserGroupMember(user_id=owner, role='owner')]
    members_info.extend(
        UserGroupMember(user_id=uid, role='member') for uid in no_owner_ids
    )

    return create_mongo_user_group_for_test(
        group_id=group_id,
        group_name=group_name,
        owner=owner,
        members=members,
        members_info=members_info,
    )


def delete_group(group_id):
    get_mongo_user_group(group_id).delete()


def init_gui_test_groups():
    user0 = user_management.user_management.get_user(username='admin').user_id
    user1 = user_management.user_management.get_user(username='test').user_id
    user2 = user_management.user_management.get_user(username='scooper').user_id
    user3 = user_management.user_management.get_user(username='ttester').user_id

    groups = {
        'group0': (
            'group0',
            'Group Admin',
            user0,
        ),
        'group1': ('group1', 'Group Test', user1),
        'group2': ('group2', 'Group Cooper', user2),
        'group3': ('group3', 'Group Tester', user3),
        'group23': ('group23', 'Group 23', user2, [user3]),
    }
    groups = {k: _create(*args) for k, args in groups.items()}

    return groups
