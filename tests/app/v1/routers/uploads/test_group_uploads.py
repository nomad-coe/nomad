import pytest

from nomad.processing.data import Upload
from tests.conftest import convert_group_labels_to_ids
from ..common import assert_response, perform_get, perform_post
from .common import assert_upload


@pytest.mark.parametrize(
    'kwargs',
    [
        pytest.param(
            dict(
                user=None,
                expected_upload_ids=[],
                expected_status_code=401,
            ),
            id='groups-guest',
        ),
        pytest.param(
            dict(
                user='invalid',
                expected_upload_ids=[],
                expected_status_code=401,
            ),
            id='groups-invalid-user',
        ),
        pytest.param(
            dict(
                user='other_test_user',
                expected_upload_ids=[
                    'id_coauthor_other_group',
                    'id_reviewer_other_group',
                    'id_coauthor_mixed_group',
                    'id_reviewer_mixed_group',
                ],
            ),
            id='groups-no-args',
        ),
        pytest.param(
            dict(
                user='other_test_user',
                query_params={'roles': 'coauthor'},
                expected_upload_ids=[
                    'id_coauthor_other_group',
                    'id_coauthor_mixed_group',
                ],
            ),
            id='groups-coauthor',
        ),
        pytest.param(
            dict(
                user='other_test_user',
                query_params={'roles': 'reviewer'},
                expected_upload_ids=[
                    'id_reviewer_other_group',
                    'id_reviewer_mixed_group',
                ],
            ),
            id='groups-reviewer',
        ),
    ],
)
def test_get_group_uploads(client, example_data_groups, test_auth_dict, kwargs):
    user = kwargs.get('user', 'test_user')
    query_params = kwargs.get('query_params', {})
    expected_status_code = kwargs.get('expected_status_code', 200)
    expected_upload_ids = kwargs.get('expected_upload_ids', None)
    user_auth, __token = test_auth_dict[user]

    response = perform_get(client, 'uploads', user_auth=user_auth, **query_params)
    assert_response(response, expected_status_code)
    if expected_status_code != 200:
        return

    response_json = response.json()
    response_data = response_json['data']

    if expected_upload_ids is not None:
        assert (
            len(response_data) == len(expected_upload_ids)
        ), f'Wrong number of records returned, expected {len(expected_upload_ids)}, got {len(response_data)}'
        found_upload_ids = [upload['upload_id'] for upload in response_data]
        assert (
            expected_upload_ids == found_upload_ids
        ), f'Wrong upload is list returned. Expected {repr(expected_upload_ids)}, got {repr(found_upload_ids)}.'


@pytest.mark.parametrize(
    'user, upload_id, expected_status_code',
    [
        pytest.param(
            'other_test_user', 'id_coauthor_other_group', 200, id='coauthor-other-group'
        ),
        pytest.param(
            'invalid', 'id_coauthor_other_group', 401, id='coauthor-other-group-invalid'
        ),
        pytest.param(
            None, 'id_coauthor_other_group', 401, id='coauthor-other-groups-guest'
        ),
        pytest.param(
            'other_test_user', 'id_reviewer_other_group', 200, id='reviewer-other-group'
        ),
        pytest.param(
            'invalid', 'id_reviewer_other_group', 401, id='reviewer-other-group-invalid'
        ),
        pytest.param(
            None, 'id_reviewer_other_group', 401, id='reviewer-other-group-guest'
        ),
        pytest.param(
            'other_test_user', 'id_coauthor_mixed_group', 200, id='coauthor-mixed-group'
        ),
        pytest.param(
            'other_test_user', 'id_reviewer_mixed_group', 200, id='reviewer-mixed-group'
        ),
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
        pytest.param(
            'other_test_user',
            422,
            'reviewer_groups',
            ['other_owner_group'],
            id='other-user-reviewer-other-group',
        ),
    ],
)
def test_add_groups_to_upload(
    client,
    user_groups_module,
    proc_infra,
    upload_no_group,
    test_auth_dict,
    user,
    expected_status_code,
    group_quantity,
    new_groups,
):
    user_auth, __token = test_auth_dict[user]
    upload_id = list(upload_no_group.uploads)[0]
    new_group_ids = [user_groups_module[label].group_id for label in new_groups]

    url = f'uploads/{upload_id}/edit'
    metadata = {group_quantity: convert_group_labels_to_ids(new_groups)}
    edit_request = dict(metadata=metadata)
    response = perform_post(client, url, user_auth, json=edit_request)

    assert_response(response, expected_status_code)
    if expected_status_code != 200:
        return

    upload = Upload.get(upload_id)
    upload.block_until_complete()
    assert getattr(upload, group_quantity) == new_group_ids


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
    assert upload.coauthor_groups == convert_group_labels_to_ids(expected_groups)
