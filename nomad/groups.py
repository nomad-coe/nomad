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
from typing import Iterable, Optional
from mongoengine import Document, StringField, ListField
from mongoengine.queryset.visitor import Q

from nomad.utils import create_uuid


class UserGroup(Document):
    """
    A group of users. One user is the owner, all others are members.
    """

    id_field = 'group_id'

    group_id = StringField(primary_key=True)
    group_name = StringField()
    owner = StringField(required=True)
    members = ListField(StringField())

    meta = {'indexes': ['group_name', 'owner', 'members']}

    @classmethod
    def get_by_user_id(cls, user_id):
        group_query = Q(owner=user_id) | Q(members=user_id)
        user_groups = cls.objects(group_query)
        return user_groups


def create_user_group(
    *,
    group_id: Optional[str] = None,
    group_name: Optional[str] = None,
    owner: Optional[str] = None,
    members: Optional[Iterable[str]] = None,
) -> UserGroup:
    user_group = UserGroup(
        group_id=group_id, group_name=group_name, owner=owner, members=members
    )
    if user_group.group_id is None:
        user_group.group_id = create_uuid()
    if user_group.group_name is None:
        user_group.group_name = user_group.group_id
    user_group.save()

    return user_group


def get_user_ids_by_group_ids(group_ids: list[str]) -> set[str]:
    user_ids = set()

    groups = UserGroup.objects(group_id__in=group_ids)
    for group in groups:
        user_ids.add(group.owner)
        user_ids.update(group.members)

    return user_ids
