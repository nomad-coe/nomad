import pytest
from .common import assert_response, perform_get, perform_post
from nomad.app.v1.routers.groups import UserGroup, UserGroups
from nomad.groups import get_user_group, user_group_exists


base_url = 'groups'


@pytest.fixture
def new_group(test_user, other_test_user):
    return {
        'group_name': 'New Group',
        'members': [test_user.user_id, other_test_user.user_id],
    }


@pytest.mark.parametrize(
    'user_label, expected_status_code',
    [
        pytest.param('test_user', 200, id='test-user'),
        pytest.param('other_test_user', 200, id='other-test-user'),
        pytest.param('invalid', 200, id='invalid-user'),
        pytest.param(None, 200, id='guest-user'),
    ],
)
def test_get_groups(
    client,
    mongo_module,
    user_label,
    test_auth_dict,
    user_groups_module,
    expected_status_code,
):
    user_auth, __token = test_auth_dict[user_label]

    response = perform_get(client, base_url, user_auth)
    assert_response(response, expected_status_code)

    groups = UserGroups.parse_raw(response.content)
    for group, ex_group in zip(groups.data, user_groups_module):
        ex_group = UserGroup.from_orm(ex_group)
        assert group == ex_group


@pytest.mark.parametrize(
    'user_label, expected_status_code',
    [
        pytest.param('test_user', 200, id='test-user'),
        pytest.param('other_test_user', 200, id='other-test-user'),
        pytest.param('invalid', 200, id='invalid-user'),
        pytest.param(None, 200, id='guest-user'),
    ],
)
def test_get_group(
    client,
    user_label,
    test_auth_dict,
    user_groups_module,
    user_owner_group,
    expected_status_code,
):
    user_auth, __token = test_auth_dict[user_label]

    response = perform_get(client, f'{base_url}/{user_owner_group.group_id}', user_auth)
    assert_response(response, expected_status_code)

    group = UserGroup.parse_raw(response.content)
    ex_group = UserGroup.from_orm(user_owner_group)
    assert group == ex_group


@pytest.mark.parametrize(
    'user_label, expected_status_code',
    [
        pytest.param('test_user', 404, id='test-user'),
        pytest.param('other_test_user', 404, id='other-test-user'),
        pytest.param('invalid', 404, id='invalid-user'),
        pytest.param(None, 404, id='guest-user'),
    ],
)
def test_get_group_invalid(
    client,
    mongo_module,
    user_label,
    test_auth_dict,
    expected_status_code,
):
    user_auth, __token = test_auth_dict[user_label]
    response = perform_get(client, f'{base_url}/invalid-group-id', user_auth)
    assert_response(response, expected_status_code)


@pytest.mark.parametrize(
    'user_label, expected_status_code',
    [
        pytest.param('test_user', 201, id='test-user'),
        pytest.param('other_test_user', 201, id='other-test-user'),
        pytest.param('invalid', 401, id='invalid-user'),
        pytest.param(None, 401, id='guest-user'),
    ],
)
def test_create_group(
    client, mongo_function, user_label, test_auth_dict, new_group, expected_status_code
):
    user_auth, __token = test_auth_dict[user_label]
    response = perform_post(client, base_url, user_auth, json=new_group)
    assert_response(response, expected_status_code)

    if response.status_code != 201:
        return

    groups = UserGroup.parse_raw(response.content)
    assert groups.group_name == new_group['group_name']
    assert set(groups.members) == set(new_group['members'])


@pytest.mark.parametrize(
    'user_label, expected_status_code',
    [
        pytest.param('test_user', 200, id='test-user'),
        pytest.param('other_test_user', 401, id='other-test-user'),
        pytest.param('invalid', 401, id='invalid-user'),
        pytest.param(None, 401, id='guest-user'),
    ],
)
def test_update_user_group(
    client,
    mongo_function,
    user_label,
    test_auth_dict,
    user_groups_function,
    user_owner_group,
    new_group,
    expected_status_code,
):
    user_auth, __token = test_auth_dict[user_label]
    group_before = get_user_group(user_owner_group.group_id)

    response = perform_post(
        client,
        f'{base_url}/{user_owner_group.group_id}/edit',
        user_auth,
        json=new_group,
    )
    assert_response(response, expected_status_code)
    group_after = get_user_group(user_owner_group.group_id)

    if response.status_code != 200:
        assert group_before == group_after
        return

    assert group_after.group_name == new_group['group_name']
    assert set(group_after.members) == set(new_group['members'])


@pytest.mark.parametrize(
    'user_label, expected_status_code',
    [
        pytest.param('test_user', 204, id='test-user'),
        pytest.param('other_test_user', 401, id='other-test-user'),
        pytest.param('invalid', 401, id='invalid-user'),
        pytest.param(None, 401, id='guest-user'),
    ],
)
def test_delete_group(
    client,
    user_label,
    test_auth_dict,
    user_groups_function,
    user_owner_group,
    other_owner_group,
    expected_status_code,
):
    user_auth, __token = test_auth_dict[user_label]

    response = client.delete(
        f'{base_url}/{user_owner_group.group_id}', headers=user_auth
    )
    assert_response(response, expected_status_code)

    if response.status_code != 204:
        assert user_group_exists(user_owner_group.group_id)
        return

    assert not user_group_exists(user_owner_group.group_id)
    assert user_group_exists(other_owner_group.group_id)


@pytest.mark.parametrize(
    'user_label, expected_status_code',
    [
        pytest.param('test_user', 404, id='test-user'),
        pytest.param('other_test_user', 404, id='other-test-user'),
        pytest.param('invalid', 401, id='invalid-user'),
        pytest.param(None, 401, id='guest-user'),
    ],
)
def test_delete_group_invalid(
    client,
    user_label,
    test_auth_dict,
    user_groups_function,
    user_owner_group,
    expected_status_code,
):
    user_auth, __token = test_auth_dict[user_label]

    response = client.delete(f'{base_url}/invalid-group-id', headers=user_auth)
    assert_response(response, expected_status_code)
    assert user_group_exists(user_owner_group.group_id)
