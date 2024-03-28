import pytest

from nomad.groups import UserGroup, create_user_group
from nomad.utils.exampledata import ExampleData
from tests.utils import test_user_group_uuid, test_user_uuid


@pytest.fixture(scope='session')
def test_user_group_molds():
    """Returns mapping from group label to field value dictionary."""

    def old_group(group_id, group_name, owner, members):
        return dict(
            group_id=test_user_group_uuid(group_id),
            group_name=group_name,
            owner=test_user_uuid(owner),
            members=[test_user_uuid(member) for member in members],
        )

    def new_group(group_name, members):
        return dict(
            group_name=group_name,
            members=[test_user_uuid(member) for member in members],
        )

    return {
        'admin_owner_group': old_group(0, 'Admin Owner Group', 0, []),
        'user_owner_group': old_group(1, 'User Owner Group', 1, []),
        'other_owner_group': old_group(2, 'Other Owner Group', 2, []),
        'yet_owner_group': old_group(3, 'Yet Owner Group', 3, []),
        'mixed_group': old_group(4, 'Mixed Group', 0, [1, 2]),
        'new_group': new_group('New Group', [2, 3]),
        'short_name': new_group('GG', []),
        'long_name': new_group('G' * 33, []),
        'double_member': new_group('Double Member', [2, 3, 2]),
        'double_member_ref': new_group('Double Member', [2, 3]),
        'special_char': new_group('G!G', []),
    }


@pytest.fixture(scope='session')
def convert_group_labels_to_ids(test_user_group_molds):
    mapping = {
        label: group.get('group_id') for label, group in test_user_group_molds.items()
    }

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
def user_owner_group(test_user_group_molds):
    return UserGroup(**test_user_group_molds['user_owner_group'])


@pytest.fixture(scope='session')
def other_owner_group(test_user_group_molds):
    return UserGroup(**test_user_group_molds['other_owner_group'])


@pytest.fixture(scope='session')
def mixed_group(test_user_group_molds):
    return UserGroup(**test_user_group_molds['mixed_group'])


@pytest.fixture(scope='session')
def create_user_groups(test_user_group_molds):
    def create():
        user_groups = {}
        for label in [
            'admin_owner_group',
            'user_owner_group',
            'other_owner_group',
            'yet_owner_group',
            'mixed_group',
        ]:
            group = test_user_group_molds[label]
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
    elastic_module, mongo_module, user_groups_module, test_user, fill_group_data
):
    data = ExampleData(main_author=test_user)
    fill_group_data(data, 'no_group', [], [])
    fill_group_data(data, 'coauthor_other', ['other_owner_group'], [])
    fill_group_data(data, 'reviewer_other', [], ['other_owner_group'])
    fill_group_data(data, 'coauthor_mixed', ['mixed_group'], [])
    fill_group_data(data, 'reviewer_mixed', [], ['mixed_group'])
    fill_group_data(data, 'reviewer_all', [], ['all'])

    data.save(with_files=False)

    yield data

    data.delete()


@pytest.fixture(scope='function')
def upload_no_group(mongo_function, test_user):
    data = ExampleData(main_author=test_user)
    data.create_upload(upload_id='id_no_group')
    data.save()

    yield data

    data.delete()


@pytest.fixture(scope='function')
def upload_coauthor_other_and_mixed_group(
    mongo_function,
    test_user,
    other_owner_group,
    mixed_group,
):
    data = ExampleData(main_author=test_user)
    data.create_upload(
        upload_id='id_coauthor_ogroup_mgroup',
        coauthor_groups=[other_owner_group.group_id, mixed_group.group_id],
    )
    data.save()

    yield data

    data.delete()


@pytest.fixture(scope='function')
def upload_reviewer_all_group(mongo_function, test_user):
    data = ExampleData(main_author=test_user)
    data.create_upload(upload_id='id_reviewer_all', reviewer_groups=['all'])
    data.save()

    yield data

    data.delete()
