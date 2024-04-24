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
    """Returns mapping from group label to field value dictionary."""

    def old_group(group_id, group_name, owner, members):
        return dict(
            group_id=fake_group_uuid(group_id),
            group_name=group_name,
            owner=fake_user_uuid(owner),
            members=[fake_user_uuid(member) for member in members],
        )

    def new_group(group_name, members):
        return dict(
            group_name=group_name,
            members=[fake_user_uuid(member) for member in members],
        )

    return {
        'group0': old_group(0, 'User0 Group', 0, []),
        'group1': old_group(1, 'User1 Group', 1, []),
        'group2': old_group(2, 'User2 Group', 2, []),
        'group3': old_group(3, 'User3 Group', 3, []),
        'group012': old_group(9012, 'Mixed Group', 0, [1, 2]),
        'new_group': new_group('New Group', [2, 3]),
        'short_name': new_group('GG', []),
        'long_name': new_group('G' * 33, []),
        'double_member': new_group('Double Member', [2, 3, 2]),
        'double_member_ref': new_group('Double Member', [2, 3]),
        'special_char': new_group('G!G', []),
    }


@pytest.fixture(scope='session')
def group_label_id_mapping(group_molds):
    return {label: value.get('group_id') for label, value in group_molds.items()}


@pytest.fixture(scope='session')
def convert_group_labels_to_ids(group_label_id_mapping):
    """Returned function converts group labels to ids, also in lists and dicts."""
    return generate_convert_label(group_label_id_mapping)


@pytest.fixture(scope='session')
def create_user_groups(group_molds):
    def create():
        user_groups = {}
        for label in [
            'group0',
            'group1',
            'group2',
            'group3',
            'group012',
        ]:
            group = group_molds[label]
            user_group = create_user_group(**group)
            user_groups[label] = user_group

        return user_groups

    return create


@pytest.fixture(scope='module')
def groups_module(mongo_module, create_user_groups):
    return create_user_groups()


@pytest.fixture(scope='function')
def groups_function(mongo_function, create_user_groups):
    return create_user_groups()
