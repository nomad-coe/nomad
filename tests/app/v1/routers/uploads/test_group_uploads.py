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
                'id_CGg2',
                'id_RGg2',
                'id_CGg123',
                'id_RGg123',
                'id_RGall',
            ],
            200,
            id='no-args',
        ),
        pytest.param(
            'user2',
            {'roles': 'coauthor'},
            ['id_CGg2', 'id_CGg123'],
            200,
            id='coauthor',
        ),
        pytest.param(
            'user2',
            {'roles': 'reviewer'},
            [
                'id_RGg2',
                'id_RGg123',
                'id_RGall',
            ],
            200,
            id='reviewer',
        ),
    ],
)
def test_get_group_uploads(
    client,
    uploads_get_groups,
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
        pytest.param('user2', 'id_CGg2', 200, id='CGg2'),
        pytest.param('invalid', 'id_CGg2', 401, id='CGg2-invalid'),
        pytest.param(None, 'id_CGg2', 401, id='CGg2-guest'),
        pytest.param('user2', 'id_RGg2', 200, id='RGg2'),
        pytest.param('invalid', 'id_RGg2', 401, id='RGg2-invalid'),
        pytest.param(None, 'id_RGg2', 401, id='RGg2-guest'),
        pytest.param('user2', 'id_CGg123', 200, id='CGg123'),
        pytest.param('user2', 'id_RGg123', 200, id='RGg123'),
        pytest.param('user2', 'id_RGall', 200, id='RGall-user2'),
        pytest.param(None, 'id_RGall', 200, id='RGall-guest'),
    ],
)
def test_get_group_upload(
    client,
    uploads_get_groups,
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
            id='CGg2',
        ),
        pytest.param(
            'upload_no_group',
            'user1',
            {
                'coauthor_groups': [
                    'group1',
                    'group2',
                    'group123',
                ]
            },
            200,
            ['group1', 'group2', 'group123'],
            id='CGg1g2g123',
        ),
        pytest.param(
            'upload_no_group',
            'user1',
            {'reviewer_groups': ['group2', 'group2']},
            200,
            ['group2'],
            id='RGg2-double',
        ),
        pytest.param(
            'upload_CGg2g123',
            'user1',
            {'coauthor_groups': {'add': 'group2'}},
            200,
            ['group2', 'group123'],
            id='CG-add-existing',
        ),
        pytest.param(
            'upload_CGg2g123',
            'user1',
            {'coauthor_groups': 'unknown_group'},
            422,
            ['group2', 'group123'],
            id='CG-set0-unknown-fails',
        ),
        pytest.param(
            'upload_CGg2g123',
            'user1',
            {'coauthor_groups': {'set': 'unknown_group'}},
            422,
            ['group2', 'group123'],
            id='CG-set1-unknown-fails',
        ),
        pytest.param(
            'upload_CGg2g123',
            'user1',
            {'coauthor_groups': {'add': 'unknown_group'}},
            422,
            ['group2', 'group123'],
            id='CG-add-unknown-fails',
        ),
        pytest.param(
            'upload_RGall',
            'user1',
            {'reviewer_groups': 'unknown_group'},
            422,
            ['all'],
            id='RG-set0-unknown-fails',
        ),
        pytest.param(
            'upload_RGall',
            'user1',
            {'reviewer_groups': {'set': 'unknown_group'}},
            422,
            ['all'],
            id='RG-set1-unknown-fails',
        ),
        pytest.param(
            'upload_RGall',
            'user1',
            {'reviewer_groups': {'add': 'unknown_group'}},
            422,
            ['all'],
            id='RG-add-unknown-fails',
        ),
        pytest.param(
            'upload_no_group',
            'user1',
            {'reviewer_groups': ['group2']},
            200,
            ['group2'],
            id='RGg2',
        ),
        pytest.param(
            'upload_no_group',
            'user1',
            {'reviewer_groups': ['all']},
            200,
            ['all'],
            id='RGall',
        ),
        pytest.param(
            'upload_no_group',
            'user1',
            {'coauthor_groups': ['all']},
            422,
            [],
            id='CGall-fails',
        ),
        pytest.param(
            'upload_no_group',
            'user2',
            {'reviewer_groups': ['group2']},
            422,
            [],
            id='user2-fails',
        ),
    ],
)
def test_add_groups_to_upload(
    request,
    client,
    groups_function,
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
            {'coauthor_groups': {'remove': 'group123'}},
            200,
            ['group2'],
            id='user1-remove-single',
        ),
        pytest.param(
            'user4',
            {'reviewer_groups': ['group4']},
            422,
            ['group2', 'group123'],
            id='user4-fails',
        ),
    ],
)
def test_remove_groups_from_upload(
    client,
    groups_function,
    proc_infra,
    convert_group_labels_to_ids,
    upload_CGg2g123,
    auth_dict,
    user,
    metadata,
    expected_status_code,
    expected_groups,
):
    user_auth, __token = auth_dict[user]
    upload_id = list(upload_CGg2g123.uploads)[0]

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
