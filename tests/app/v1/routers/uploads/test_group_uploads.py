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
            'user2',
            {},
            [
                'id_coauthor_group2',
                'id_reviewer_group2',
                'id_coauthor_group012',
                'id_reviewer_group012',
                'id_reviewer_all',
            ],
            200,
            id='no-args',
        ),
        pytest.param(
            'user2',
            {'roles': 'coauthor'},
            ['id_coauthor_group2', 'id_coauthor_group012'],
            200,
            id='coauthor',
        ),
        pytest.param(
            'user2',
            {'roles': 'reviewer'},
            [
                'id_reviewer_group2',
                'id_reviewer_group012',
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
    auth_dict,
    user,
    query_params,
    expected_upload_ids,
    expected_status_code,
):
    user_auth, __token = auth_dict[user]

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
        pytest.param('user2', 'id_coauthor_group2', 200, id='coauthor-group2'),
        pytest.param(
            'invalid', 'id_coauthor_group2', 401, id='coauthor-group2-invalid'
        ),
        pytest.param(None, 'id_coauthor_group2', 401, id='coauthor-group2-guest'),
        pytest.param('user2', 'id_reviewer_group2', 200, id='reviewer-group2'),
        pytest.param(
            'invalid', 'id_reviewer_group2', 401, id='reviewer-group2-invalid'
        ),
        pytest.param(None, 'id_reviewer_group2', 401, id='reviewer-group2-guest'),
        pytest.param('user2', 'id_coauthor_group012', 200, id='coauthor-group012'),
        pytest.param('user2', 'id_reviewer_group012', 200, id='reviewer-group012'),
        pytest.param('user2', 'id_reviewer_all', 200, id='reviewer-all-user2'),
        pytest.param(None, 'id_reviewer_all', 200, id='reviewer-all-guest'),
    ],
)
def test_get_group_upload(
    client,
    example_data_groups,
    auth_dict,
    user,
    upload_id,
    expected_status_code,
):
    user_auth, __token = auth_dict[user]
    response = perform_get(client, f'uploads/{upload_id}', user_auth)
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        assert_upload(response.json())


@pytest.mark.parametrize(
    'upload_fixture, user, metadata, expected_status_code, expected_groups',
    [
        pytest.param(
            'upload_no_group',
            'user1',
            {'coauthor_groups': 'group2'},
            200,
            ['group2'],
            id='coauthor-group2',
        ),
        pytest.param(
            'upload_no_group',
            'user1',
            {
                'coauthor_groups': [
                    'group1',
                    'group2',
                    'group012',
                ]
            },
            200,
            ['group1', 'group2', 'group012'],
            id='coauthor-multiple-groups',
        ),
        pytest.param(
            'upload_no_group',
            'user1',
            {'reviewer_groups': ['group2', 'group2']},
            200,
            ['group2'],
            id='reviewer-group2-double',
        ),
        pytest.param(
            'upload_coauthor_group2_and_group012',
            'user1',
            {'coauthor_groups': {'add': 'group2'}},
            200,
            ['group2', 'group012'],
            id='coauthor-add-existing',
        ),
        pytest.param(
            'upload_coauthor_group2_and_group012',
            'user1',
            {'coauthor_groups': 'unknown_group'},
            422,
            ['group2', 'group012'],
            id='coauthor-set0-unknown-fails',
        ),
        pytest.param(
            'upload_coauthor_group2_and_group012',
            'user1',
            {'coauthor_groups': {'set': 'unknown_group'}},
            422,
            ['group2', 'group012'],
            id='coauthor-set1-unknown-fails',
        ),
        pytest.param(
            'upload_coauthor_group2_and_group012',
            'user1',
            {'coauthor_groups': {'add': 'unknown_group'}},
            422,
            ['group2', 'group012'],
            id='coauthor-add-unknown-fails',
        ),
        pytest.param(
            'upload_reviewer_all_group',
            'user1',
            {'reviewer_groups': 'unknown_group'},
            422,
            ['all'],
            id='reviewer-set0-unknown-fails',
        ),
        pytest.param(
            'upload_reviewer_all_group',
            'user1',
            {'reviewer_groups': {'set': 'unknown_group'}},
            422,
            ['all'],
            id='reviewer-set1-unknown-fails',
        ),
        pytest.param(
            'upload_reviewer_all_group',
            'user1',
            {'reviewer_groups': {'add': 'unknown_group'}},
            422,
            ['all'],
            id='reviewer-add-unknown-fails',
        ),
        pytest.param(
            'upload_no_group',
            'user1',
            {'reviewer_groups': ['group2']},
            200,
            ['group2'],
            id='reviewer-group2',
        ),
        pytest.param(
            'upload_no_group',
            'user1',
            {'reviewer_groups': ['all']},
            200,
            ['all'],
            id='reviewer-all',
        ),
        pytest.param(
            'upload_no_group',
            'user1',
            {'coauthor_groups': ['all']},
            422,
            [],
            id='coauthor-all-fails',
        ),
        pytest.param(
            'upload_no_group',
            'user2',
            {'reviewer_groups': ['group2']},
            422,
            [],
            id='other-user-fails',
        ),
    ],
)
def test_add_groups_to_upload(
    request,
    client,
    user_groups_function,
    proc_infra,
    convert_group_labels_to_ids,
    upload_fixture,
    auth_dict,
    user,
    expected_status_code,
    metadata,
    expected_groups,
):
    user_auth, __token = auth_dict[user]
    upload_fixture = request.getfixturevalue(upload_fixture)
    upload_id = list(upload_fixture.uploads)[0]
    group_quantity = list(metadata)[0]

    url = f'uploads/{upload_id}/edit'
    metadata = convert_group_labels_to_ids(metadata)
    edit_request = dict(metadata=metadata)
    response = perform_post(client, url, user_auth, json=edit_request)

    assert_response(response, expected_status_code)
    upload = Upload.get(upload_id)
    upload.block_until_complete()

    group_ids = getattr(upload, group_quantity)
    expected_ids = convert_group_labels_to_ids(expected_groups)
    assert group_ids == expected_ids


@pytest.mark.parametrize(
    'user, metadata, expected_status_code, expected_groups',
    [
        pytest.param('user1', {'coauthor_groups': []}, 200, [], id='user1-empty'),
        pytest.param(
            'user1',
            {'coauthor_groups': {'set': []}},
            200,
            [],
            id='user1-set-empty',
        ),
        pytest.param(
            'user1',
            {'coauthor_groups': {'remove': 'group012'}},
            200,
            ['group2'],
            id='user1-remove-single',
        ),
        pytest.param(
            'user3',
            {'reviewer_groups': ['group3']},
            422,
            ['group2', 'group012'],
            id='user3-fails',
        ),
    ],
)
def test_remove_groups_from_upload(
    client,
    user_groups_function,
    proc_infra,
    convert_group_labels_to_ids,
    upload_coauthor_group2_and_group012,
    auth_dict,
    user,
    metadata,
    expected_status_code,
    expected_groups,
):
    user_auth, __token = auth_dict[user]
    upload_id = list(upload_coauthor_group2_and_group012.uploads)[0]

    url = f'uploads/{upload_id}/edit'
    metadata = convert_group_labels_to_ids(metadata)
    edit_request = dict(metadata=metadata)
    response = perform_post(client, url, user_auth, json=edit_request)

    assert_response(response, expected_status_code)
    if expected_status_code != 200:
        return

    upload = Upload.get(upload_id)
    upload.block_until_complete()
    expected_groups = convert_group_labels_to_ids(expected_groups)
    assert upload.coauthor_groups == expected_groups
