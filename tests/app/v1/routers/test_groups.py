import pytest
from .common import assert_response, perform_get, perform_post
from nomad.app.v1.routers.groups import UserGroup, UserGroups
from nomad.groups import get_user_group, user_group_exists


base_url = 'groups'


def get_val(obj, key):
    if isinstance(obj, dict):
        return obj[key]

    return getattr(obj, key)


def assert_unordered_lists(list1, list2):
    assert sorted(list1) == sorted(list2)


def assert_group(group, ref_group, keys=None):
    if keys is None:
        keys = UserGroup.__fields__

    excluded_fields = {'members'}
    fields = set(keys) - excluded_fields
    for field in fields:
        assert get_val(group, field) == get_val(ref_group, field)

    if 'members' in keys:
        val = get_val(group, 'members')
        ref_val = get_val(ref_group, 'members')
        assert_unordered_lists(val, ref_val)


@pytest.mark.parametrize(
    'user_label, expected_status_code',
    [
        pytest.param('user1', 200, id='user1'),
        pytest.param('user2', 200, id='user2'),
        pytest.param('invalid', 200, id='invalid-user'),
        pytest.param(None, 200, id='guest-user'),
    ],
)
def test_get_groups(
    client,
    mongo_module,
    auth_dict,
    user_groups_module,
    user_label,
    expected_status_code,
):
    user_auth, _ = auth_dict[user_label]

    response = perform_get(client, base_url, user_auth)
    assert_response(response, expected_status_code)

    response_groups = UserGroups.parse_raw(response.content)
    for response_group in response_groups.data:
        group = get_user_group(response_group.group_id)
        assert_group(group, response_group)


@pytest.mark.parametrize(
    'user_label, expected_status_code',
    [
        pytest.param('user1', 200, id='user1'),
        pytest.param('user2', 200, id='user2'),
        pytest.param('invalid', 200, id='invalid-user'),
        pytest.param(None, 200, id='guest-user'),
    ],
)
def test_get_group(
    client,
    auth_dict,
    user_groups_module,
    user_label,
    expected_status_code,
):
    user_auth, _ = auth_dict[user_label]
    ref_group = user_groups_module['group2']

    response = perform_get(client, f'{base_url}/{ref_group.group_id}', user_auth)
    assert_response(response, expected_status_code)

    response_group = UserGroup.parse_raw(response.content)
    group = get_user_group(response_group.group_id)
    assert_group(group, response_group)
    assert_group(group, ref_group)


@pytest.mark.parametrize(
    'user_label, expected_status_code',
    [
        pytest.param('user1', 404, id='user1'),
        pytest.param('user2', 404, id='user2'),
        pytest.param('invalid', 404, id='invalid-user'),
        pytest.param(None, 404, id='guest-user'),
    ],
)
def test_get_group_invalid(
    client,
    mongo_module,
    auth_dict,
    user_groups_module,
    user_label,
    expected_status_code,
):
    user_auth, _ = auth_dict[user_label]

    response = perform_get(client, f'{base_url}/invalid-group-id', user_auth)
    assert_response(response, expected_status_code)


@pytest.mark.parametrize(
    'user_label, new_group_label, ref_group_label, expected_status_code',
    [
        pytest.param('user1', 'new_group', 'new_group', 201, id='user1'),
        pytest.param('user2', 'new_group', 'new_group', 201, id='user2'),
        pytest.param('invalid', 'new_group', None, 401, id='invalid-user'),
        pytest.param(None, 'new_group', None, 401, id='guest-user'),
        pytest.param('user1', 'short_name', None, 422, id='short-name-fails'),
        pytest.param('user1', 'long_name', None, 422, id='long-name-fails'),
        pytest.param(
            'user1',
            'double_member',
            'double_member_ref',
            201,
            id='double-member-skipped',
        ),
    ],
)
def test_create_group(
    client,
    mongo_function,
    request,
    auth_dict,
    group_molds,
    user_label,
    new_group_label,
    ref_group_label,
    expected_status_code,
):
    user_auth, _ = auth_dict[user_label]
    new_group = group_molds[new_group_label]

    response = perform_post(client, base_url, user_auth, json=new_group)
    assert_response(response, expected_status_code)

    if response.status_code != 201:
        return

    response_group = UserGroup.parse_raw(response.content)
    group = get_user_group(response_group.group_id)
    assert_group(group, response_group)
    ref_group = group_molds[ref_group_label]
    assert_group(group, ref_group, ref_group.keys())


@pytest.mark.parametrize(
    'user_label, group_edit_label, ref_group_label, expected_status_code',
    [
        pytest.param(None, 'new_group', None, 401, id='guest-fails'),
        pytest.param('invalid', 'new_group', None, 401, id='faker-fails'),
        pytest.param('user2', 'new_group', None, 401, id='other-fails'),
        pytest.param('user1', 'new_group', 'new_group', 200, id='edit-ok'),
        pytest.param('user1', 'short_name', None, 422, id='short-name-fails'),
        pytest.param('user1', 'long_name', None, 422, id='long-name-fails'),
        pytest.param('user1', 'special_char', None, 422, id='special-chars-fails'),
        pytest.param(
            'user1',
            'double_member',
            'double_member_ref',
            200,
            id='double-member-skipped',
        ),
    ],
)
def test_update_user_group(
    client,
    mongo_function,
    auth_dict,
    group_molds,
    user_groups_function,
    group1,
    user_label,
    group_edit_label,
    ref_group_label,
    expected_status_code,
):
    user_auth, _ = auth_dict[user_label]
    group_before = get_user_group(group1.group_id)
    group_edit = group_molds[group_edit_label]

    url = f'{base_url}/{group_before.group_id}/edit'
    response = perform_post(client, url, user_auth, json=group_edit)
    assert_response(response, expected_status_code)
    group_after = get_user_group(group_before.group_id)

    if response.status_code != 200:
        assert_group(group_after, group_before)
        return

    response_group = UserGroup.parse_raw(response.content)
    assert_group(group_after, response_group)
    ref_group = group_molds[ref_group_label]
    assert_group(group_after, ref_group, ref_group.keys())


@pytest.mark.parametrize(
    'user_label, expected_status_code',
    [
        pytest.param('user1', 204, id='user1'),
        pytest.param('user2', 401, id='user2'),
        pytest.param('invalid', 401, id='invalid-user'),
        pytest.param(None, 401, id='guest-user'),
    ],
)
def test_delete_group(
    auth_dict,
    client,
    group1,
    group2,
    user_groups_function,
    user_label,
    expected_status_code,
):
    user_auth, _ = auth_dict[user_label]

    response = client.delete(f'{base_url}/{group1.group_id}', headers=user_auth)
    assert_response(response, expected_status_code)

    if response.status_code != 204:
        assert user_group_exists(group1.group_id)
        return

    assert not user_group_exists(group1.group_id)
    assert user_group_exists(group2.group_id)


@pytest.mark.parametrize(
    'user_label, expected_status_code',
    [
        pytest.param('user1', 404, id='user1'),
        pytest.param('user2', 404, id='user2'),
        pytest.param('invalid', 401, id='invalid-user'),
        pytest.param(None, 401, id='guest-user'),
    ],
)
def test_delete_group_invalid(
    client,
    auth_dict,
    user_groups_function,
    group1,
    user_label,
    expected_status_code,
):
    user_auth, _ = auth_dict[user_label]

    response = client.delete(f'{base_url}/invalid-group-id', headers=user_auth)
    assert_response(response, expected_status_code)
    assert user_group_exists(group1.group_id)
