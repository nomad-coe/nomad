import pytest

from nomad.processing.data import Upload
from ..common import assert_response, perform_get, perform_post
from .common import assert_upload


@pytest.mark.parametrize(
    'user, query_params, expected_upload_ids, expected_status_code',
    [
        pytest.param(None, {}, [], 401, id='guest-user'),
        pytest.param('invalid', {}, [], 401, id='invalid-user'),
        pytest.param(
            'other_test_user',
            {},
            [
                'id_coauthor_other',
                'id_reviewer_other',
                'id_coauthor_mixed',
                'id_reviewer_mixed',
                'id_reviewer_all',
            ],
            200,
            id='no-args',
        ),
        pytest.param(
            'other_test_user',
            {'roles': 'coauthor'},
            ['id_coauthor_other', 'id_coauthor_mixed'],
            200,
            id='coauthor',
        ),
        pytest.param(
            'other_test_user',
            {'roles': 'reviewer'},
            [
                'id_reviewer_other',
                'id_reviewer_mixed',
                'id_reviewer_all',
            ],
            200,
            id='reviewer',
        ),
    ],
)
def test_get_group_uploads(
    client,
    example_data_groups,
    test_auth_dict,
    user,
    query_params,
    expected_upload_ids,
    expected_status_code,
):
    user_auth, __token = test_auth_dict[user]

    response = perform_get(client, 'uploads', user_auth=user_auth, **query_params)
    assert_response(response, expected_status_code)
    if expected_status_code != 200:
        return

    response_json = response.json()
    response_data = response_json['data']

    assert len(response_data) == len(expected_upload_ids)
    found_upload_ids = [upload['upload_id'] for upload in response_data]
    assert found_upload_ids == expected_upload_ids


@pytest.mark.parametrize(
    'user, upload_id, expected_status_code',
    [
        pytest.param('other_test_user', 'id_coauthor_other', 200, id='coauthor-other'),
        pytest.param('invalid', 'id_coauthor_other', 401, id='coauthor-other-invalid'),
        pytest.param(None, 'id_coauthor_other', 401, id='coauthor-other-guest'),
        pytest.param('other_test_user', 'id_reviewer_other', 200, id='reviewer-other'),
        pytest.param('invalid', 'id_reviewer_other', 401, id='reviewer-other-invalid'),
        pytest.param(None, 'id_reviewer_other', 401, id='reviewer-other-guest'),
        pytest.param('other_test_user', 'id_coauthor_mixed', 200, id='coauthor-mixed'),
        pytest.param('other_test_user', 'id_reviewer_mixed', 200, id='reviewer-mixed'),
        pytest.param('other_test_user', 'id_reviewer_all', 200, id='reviewer-all'),
    ],
)
def test_get_group_upload(
    client,
    example_data_groups,
    test_auth_dict,
    user,
    upload_id,
    expected_status_code,
):
    user_auth, __token = test_auth_dict[user]
    response = perform_get(client, f'uploads/{upload_id}', user_auth)
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        assert_upload(response.json())


@pytest.mark.parametrize(
    'user, expected_status_code, group_quantity, new_groups',
    [
        pytest.param(
            'test_user',
            200,
            'coauthor_groups',
            ['other_owner_group'],
            id='coauthor-other-group',
        ),
        pytest.param(
            'test_user',
            200,
            'coauthor_groups',
            ['user_owner_group', 'other_owner_group', 'mixed_group'],
            id='coauthor-multiple-groups',
        ),
        pytest.param(
            'test_user',
            200,
            'reviewer_groups',
            ['other_owner_group'],
            id='reviewer-other-group',
        ),
        pytest.param('test_user', 200, 'reviewer_groups', ['all'], id='reviewer-all'),
        pytest.param(
            'test_user', 422, 'coauthor_groups', ['all'], id='coauthor-all-fails'
        ),
        pytest.param(
            'other_test_user',
            422,
            'reviewer_groups',
            ['other_owner_group'],
            id='other-user-reviewer-other-group-fails',
        ),
    ],
)
def test_add_groups_to_upload(
    client,
    user_groups_module,
    proc_infra,
    convert_group_labels_to_ids,
    upload_no_group,
    test_auth_dict,
    user,
    expected_status_code,
    group_quantity,
    new_groups,
):
    user_auth, __token = test_auth_dict[user]
    upload_id = list(upload_no_group.uploads)[0]
    new_group_ids = convert_group_labels_to_ids(new_groups)

    url = f'uploads/{upload_id}/edit'
    metadata = {group_quantity: new_group_ids}
    edit_request = dict(metadata=metadata)
    response = perform_post(client, url, user_auth, json=edit_request)

    assert_response(response, expected_status_code)
    upload = Upload.get(upload_id)
    upload.block_until_complete()

    if expected_status_code == 200:
        assert getattr(upload, group_quantity) == new_group_ids
    else:
        assert getattr(upload, group_quantity) == []


@pytest.mark.parametrize(
    'user, action, expected_groups, expected_status_code',
    [
        pytest.param('test_user', [], [], 200, id='test-user-empty'),
        pytest.param('test_user', {'set': []}, [], 200, id='test-user-set-empty'),
        pytest.param(
            'test_user',
            {'remove': 'mixed_group'},
            ['other_owner_group'],
            200,
            id='test-user-remove-other',
        ),
        pytest.param('other_test_user', [], None, 422, id='other-user-fail'),
    ],
)
def test_remove_groups_from_upload(
    client,
    user_groups_module,
    proc_infra,
    convert_group_labels_to_ids,
    upload_coauthor_other_and_mixed_group,
    test_auth_dict,
    user,
    action,
    expected_groups,
    expected_status_code,
):
    user_auth, __token = test_auth_dict[user]
    upload_id = list(upload_coauthor_other_and_mixed_group.uploads)[0]

    url = f'uploads/{upload_id}/edit'
    metadata = {'coauthor_groups': convert_group_labels_to_ids(action)}
    edit_request = dict(metadata=metadata)
    response = perform_post(client, url, user_auth, json=edit_request)

    assert_response(response, expected_status_code)
    if expected_status_code != 200:
        return

    upload = Upload.get(upload_id)
    upload.block_until_complete()
    expected_groups = convert_group_labels_to_ids(expected_groups)
    assert upload.coauthor_groups == expected_groups
