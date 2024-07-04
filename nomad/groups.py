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
import operator
from functools import reduce
from typing import Iterable, List, Optional

from mongoengine import Document, ListField, StringField
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
    def get_by_ids(cls, group_ids: Iterable[str]):
        """
        Returns UserGroup objects with group_ids.
        """
        user_groups = cls.objects(group_id__in=group_ids)
        return user_groups

    @classmethod
    def get_by_user_id(cls, user_id: Optional[str]):
        """
        Returns UserGroup objects where user_id is owner or member, or None.
        Does not include special group 'all' because it has no UserGroup object.
        """
        group_query = Q(owner=user_id) | Q(members=user_id)
        user_groups = cls.objects(group_query)
        return user_groups

    @classmethod
    def get_ids_by_user_id(cls, user_id: Optional[str], include_all=True) -> List[str]:
        """
        Returns ids of all user groups where user_id is owner or member.
        Does include special group 'all', even if user_id is missing or not a user.
        """
        group_ids = ['all'] if include_all else []
        if user_id is not None:
            group_ids.extend(group.group_id for group in cls.get_by_user_id(user_id))
        return group_ids

    @classmethod
    def get_by_search_terms(cls, search_terms: str):
        """
        Returns UserGroup objects where group_name includes search_terms (no case).
        """
        search_terms = str(search_terms).split()
        if not search_terms:
            return []

        query = (Q(group_name__icontains=term) for term in search_terms)
        query = reduce(operator.and_, query)
        user_groups = cls.objects(query)
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


def get_user_group(group_id: str) -> Optional[UserGroup]:
    return UserGroup.objects(group_id=group_id).first()


def user_group_exists(group_id: str, *, include_all=True) -> bool:
    if include_all and group_id == 'all':
        return True
    return get_user_group(group_id) is not None


def get_group_ids(user_id, include_all=True):
    return UserGroup.get_ids_by_user_id(user_id, include_all=include_all)
