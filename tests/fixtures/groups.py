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

"""
Group fixtures:
- groupO: group owned by userO without members
- groupOMN…: group owned by userO with members userM, userN, …
- groupOMmN…: group owned by userO with maintainer userM and members userN, …
- special (unfinished) groups, often only partly defined
"""

import pytest

from nomad.app.v1.models.groups import UserGroupMember
from nomad.mongo.groups import create_mongo_user_group_for_test
from tests.utils import (
    fake_group_uuid,
    fake_user_uuid,
    generate_convert_label,
    list_without,
)


def make_member(user_id, role):
    return dict(user_id=fake_user_uuid(user_id), role=role)


def make_members_info(owner=None, maintainers=None, members=None):
    members_info = []
    if owner is not None:
        members_info.append(make_member(owner, 'owner'))
    if maintainers is not None:
        members_info.extend(make_member(id, 'maintainer') for id in maintainers)
    if members is not None:
        rest = list_without(members, owner, *(maintainers or []))
        members_info.extend(make_member(id, 'member') for id in rest)

    return members_info


def custom_group(
    name=None,
    owner=None,
    maintainers=None,
    members=None,
    *,
    fill=None,
    group_id=None,
):
    if fill is None:
        fill = {'owner', 'members', 'members_info'}

    mold = {}
    if group_id is not None:
        mold['group_id'] = fake_group_uuid(group_id)
    if name is not None:
        mold['group_name'] = name
    if 'owner' in fill and owner is not None:
        mold['owner'] = fake_user_uuid(owner)
    if 'members_info' in fill:
        mold['members_info'] = make_members_info(owner, maintainers, members)
    if 'members' in fill:
        mold['members'] = [fake_user_uuid(m) for m in members or []]

    return mold


def default_group(owner, members=None, *, maintainers=None, label=None):
    """Return group mold with conforming id, name and members_info."""
    if maintainers is None:
        maintainers = []
    if members is None:
        members = []

    if label is None:
        label = str(owner)
        if len(maintainers) > 0:
            label += ''.join(str(m) for m in maintainers) + 'm'
        label += ''.join(str(m) for m in list_without(members, owner, *maintainers))

    return custom_group(
        name=f'Group {label}',
        owner=owner,
        maintainers=maintainers,
        members=members,
        group_id=label,
    )


def full_group(*args, **kwargs):
    return custom_group(*args, **kwargs, fill=None)


def info_group(*args, **kwargs):
    return custom_group(*args, **kwargs, fill={'members_info'})


def members_group(*args, **kwargs):
    return custom_group(*args, **kwargs, fill={'owner', 'members'})


@pytest.fixture(scope='session')
def group_molds():
    """Return a dict: group label -> group data (dict)."""

    # lots of groups are also used for upload tests
    default_groups = {
        'group0': default_group(0, [0]),
        'group1': default_group(1, [1]),
        'group2': default_group(2, [2]),
        'group3': default_group(3, [3]),
        'group6': default_group(6, [6]),
        'group8': default_group(8, [8]),
        'group9': default_group(9, [9]),
        'group14': default_group(1, [1, 4]),
        'group15': default_group(1, [1, 5]),
        'group18': default_group(1, [1, 8]),
        'group19': default_group(1, [1, 9]),
        'group123': default_group(1, [1, 2, 3]),
        'group12m3': default_group(1, [1, 2, 3], maintainers=[2]),
        'uniq': default_group(0, [0], label='Uniq'),
        'twin1': default_group(0, [0], label='Twin One'),
        'twin2': default_group(0, [0], label='Twin Two'),
        'numerals': default_group(0, [0], label='One Two Three'),
        'dirty234': custom_group(
            'Group Dirty 234',
            owner=2,
            members=[3, 4],
            group_id=fake_group_uuid('dirty234'),
            fill={'owner', 'members'},
        ),
    }

    custom_groups = {
        'dirty234_ref': full_group('Group Dirty 234', owner=2, members=[2, 3, 4]),
        'edit_x23': info_group('New Group X23', members=[2, 3]),
        'edit_123': info_group('New Group 123', owner=1, members=[2, 3]),
        'edit_123_ref': full_group('New Group 123', owner=1, members=[1, 2, 3]),
        'edit_14m5': info_group(
            'New Group 14m5', owner=1, maintainers=[4], members=[5]
        ),
        'edit_14m5_ref': full_group(
            'New Group 14m5', owner=1, maintainers=[4], members=[1, 4, 5]
        ),
        'edit_1232': info_group('Group 1232', members=[2, 3, 2]),
        'edit_1232_ref': info_group('Group 1232', members=[1, 2, 3]),
        'empty_name': info_group('', owner=1),
        'short_name': info_group('G', owner=1),
        'long_name': info_group('G' * 65, owner=1),
        'whitespace_name': info_group(' \t ', owner=1),
        'special_chars_name': info_group(
            '!@#$%^&*()_+-=[]}\r\t\n{:; "\'|\\<,>.?/…😀', owner=1
        ),
        'no_name': info_group(owner=1),
        'no_name_ref': full_group(owner=1, members=[1]),
        'only_name': custom_group('Only Name Group', fill=set()),
        'only_name_ref': full_group('Only Name Group', owner=1, members=[1]),
        'both_fields': full_group('Both Fields Group', owner=1, members=[1, 2, 3]),
    }

    old_custom_groups = {
        'old_edit_x23': members_group('Group X23', members=[2, 3]),
        'old_edit_x23_ref1': full_group('Group X23', owner=1, members=[1, 2, 3]),
        'old_edit_x23_ref2': full_group('Group X23', owner=2, members=[2, 3]),
        'old_edit_x123': members_group('Group X123', members=[1, 2, 3]),
        'old_edit_x123_ref': full_group('Group X123', owner=1, members=[1, 2, 3]),
        'old_edit_1232': members_group('Group 1232', members=[2, 3, 2]),
        'old_edit_1232_ref': full_group('Group 1232', owner=1, members=[1, 2, 3]),
        'old_edit_x24': members_group('Group x24', members=[2, 4]),
        'old_edit_x24_ref': full_group(
            'Group x24', owner=1, maintainers=[2], members=[1, 2, 4]
        ),
    }

    is_disjoint = set(default_groups).isdisjoint(old_custom_groups | custom_groups)
    assert is_disjoint, 'Duplicate group labels in default and custom groups.'

    return {**default_groups, **custom_groups, **old_custom_groups}


@pytest.fixture(scope='session')
def group_label_id_mapping(group_molds):
    """Return a dict: group label -> group id."""
    return {label: value.get('group_id') for label, value in group_molds.items()}


@pytest.fixture(scope='session')
def convert_group_labels_to_ids(group_label_id_mapping):
    """Returned function converts group labels to ids, also in lists and dicts."""
    return generate_convert_label(group_label_id_mapping)


@pytest.fixture(scope='session')
def convert_agent_labels_to_ids(user_label_id_mapping, group_label_id_mapping):
    """Returned function converts agent labels to ids, also in lists and dicts."""
    is_disjoint = set(user_label_id_mapping).isdisjoint(group_label_id_mapping)
    assert is_disjoint, 'Duplicate labels in users and groups.'
    mapping = {**user_label_id_mapping, **group_label_id_mapping}
    return generate_convert_label(mapping)


@pytest.fixture(scope='session')
def create_user_groups(group_molds):
    """Returned function creates and returns predefined user groups for testing."""

    def create():
        groups_with_id = {k: g for k, g in group_molds.items() if 'group_id' in g}
        user_groups = {
            k: create_mongo_user_group_for_test(
                group_id=g.get('group_id'),
                group_name=g.get('group_name'),
                owner=g.get('owner'),
                members=g.get('members'),
                members_info=[UserGroupMember(**m) for m in g['members_info']]
                if 'members_info' in g
                else None,
            )
            for k, g in groups_with_id.items()
        }
        return user_groups

    return create


@pytest.fixture(scope='module')
def groups_module(mongo_module, create_user_groups):
    """Create and return predefined user groups for testing (module scope)."""
    return create_user_groups()


@pytest.fixture
def groups_function(mongo_function, create_user_groups):
    """Create and return predefined user groups for testing (function scope)."""
    return create_user_groups()


@pytest.fixture
def group_owner_not_member(mongo_function, group_molds):
    """Create and return a group where owner is not a member (old behavior)."""
    mold = group_molds['dirty234']
    group = create_mongo_user_group_for_test(**mold)
    return group
