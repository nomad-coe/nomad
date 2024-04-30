"""
Group fixtures:
- groupO: group owned by userO without members
- groupOMN…: group owned by userO with members userM, userN, …
- special (unfinished) groups, often only partly defined
"""

import pytest

from nomad.groups import create_user_group
from tests.utils import fake_group_uuid, fake_user_uuid, generate_convert_label


@pytest.fixture(scope='session')
def group_molds():
    """Return a dict: group label -> group data (dict)."""

    def old_group(owner, members):
        group_str = str(owner) + ''.join(str(m) for m in members)
        return dict(
            group_id=fake_group_uuid(group_str),
            group_name=f'Group {group_str}',
            owner=fake_user_uuid(owner),
            members=[fake_user_uuid(member) for member in members],
        )

    def new_group(group_name, members):
        return dict(
            group_name=group_name,
            members=[fake_user_uuid(member) for member in members],
        )

    old_groups = {
        'group0': old_group(0, []),
        'group1': old_group(1, []),
        'group2': old_group(2, []),
        'group3': old_group(3, []),
        'group6': old_group(6, []),
        'group8': old_group(8, []),
        'group9': old_group(9, []),
        'group14': old_group(1, [4]),
        'group15': old_group(1, [5]),
        'group18': old_group(1, [8]),
        'group19': old_group(1, [9]),
        'group123': old_group(1, [2, 3]),
    }

    new_groups = {
        'new_group': new_group('New Group X23', [2, 3]),
        'short_name': new_group('GG', []),
        'long_name': new_group('G' * 33, []),
        'double_member': new_group('Double Member', [2, 3, 2]),
        'double_member_ref': new_group('Double Member', [2, 3]),
        'special_char': new_group('G!G', []),
    }

    return {**old_groups, **new_groups}


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
        groups_with_id = {k: v for k, v in group_molds.items() if 'group_id' in v}
        user_groups = {k: create_user_group(**v) for k, v in groups_with_id.items()}

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
