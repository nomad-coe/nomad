"""
Group fixtures:
- groupO: group owned by userO without members
- groupOMN…: group owned by userO with members userM, userN, …
- special (unfinished) groups, often only partly defined
"""

import pytest

from nomad.groups import UserGroup, create_user_group
from nomad.utils.exampledata import ExampleData
from tests.utils import fake_group_uuid, fake_user_uuid


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
def convert_group_labels_to_ids(group_molds):
    mapping = {label: group.get('group_id') for label, group in group_molds.items()}

    def convert(raw):
        if isinstance(raw, str):
            return mapping.get(raw, raw)

        if isinstance(raw, list):
            return [convert(v) for v in raw]

        if isinstance(raw, dict):
            return {k: convert(v) for k, v in raw.items()}

        return raw

    return convert


@pytest.fixture(scope='session')
def group1(group_molds):
    return UserGroup(**group_molds['group1'])


@pytest.fixture(scope='session')
def group2(group_molds):
    return UserGroup(**group_molds['group2'])


@pytest.fixture(scope='session')
def group012(group_molds):
    return UserGroup(**group_molds['group012'])


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
def user_groups_module(mongo_module, create_user_groups):
    return create_user_groups()


@pytest.fixture(scope='function')
def user_groups_function(mongo_function, create_user_groups):
    return create_user_groups()


@pytest.fixture('module')
def fill_group_data(convert_group_labels_to_ids):
    def fill(data: ExampleData, id_label, c_groups, r_groups, **kwargs):
        upload_id = f'id_{id_label}'
        entry_id = f'{upload_id}_1'

        d = kwargs.copy()
        d['upload_id'] = upload_id
        d['coauthor_groups'] = convert_group_labels_to_ids(c_groups)
        d['reviewer_groups'] = convert_group_labels_to_ids(r_groups)
        d.setdefault('published', d.get('embargo_length') is not None)
        if d.get('embargo_length') is None:
            d['embargo_length'] = 0

        data.create_upload(**d)
        data.create_entry(upload_id=upload_id, entry_id=entry_id)

    return fill


@pytest.fixture(scope='module')
def example_data_groups(
    elastic_module, mongo_module, user_groups_module, user1, fill_group_data
):
    data = ExampleData(main_author=user1)
    fill_group_data(data, 'no_group', [], [])
    fill_group_data(data, 'coauthor_group2', ['group2'], [])
    fill_group_data(data, 'reviewer_group2', [], ['group2'])
    fill_group_data(data, 'coauthor_group012', ['group012'], [])
    fill_group_data(data, 'reviewer_group012', [], ['group012'])
    fill_group_data(data, 'reviewer_all', [], ['all'])

    data.save(with_files=False)

    yield data

    data.delete()


@pytest.fixture(scope='function')
def upload_no_group(mongo_function, user1):
    data = ExampleData(main_author=user1)
    data.create_upload(upload_id='id_no_group')
    data.save()

    yield data

    data.delete()


@pytest.fixture(scope='function')
def upload_coauthor_group2_and_group012(
    mongo_function,
    user1,
    group2,
    group012,
):
    data = ExampleData(main_author=user1)
    data.create_upload(
        upload_id='id_coauthor_ogroup_mgroup',
        coauthor_groups=[group2.group_id, group012.group_id],
    )
    data.save()

    yield data

    data.delete()


@pytest.fixture(scope='function')
def upload_reviewer_all_group(mongo_function, user1):
    data = ExampleData(main_author=user1)
    data.create_upload(upload_id='id_reviewer_all', reviewer_groups=['all'])
    data.save()

    yield data

    data.delete()
