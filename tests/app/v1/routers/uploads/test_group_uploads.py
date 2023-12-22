import pytest
from ..common import assert_response, perform_get
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
