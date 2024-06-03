import pytest

from nomad.processing.data import Upload
from tests.utils import dict_to_params

from ..common import assert_response, perform_get, perform_post
from .common import assert_entry, assert_upload


def get_agents_from_upload(upload):
    keys = ('coauthors', 'reviewers', 'coauthor_groups', 'reviewer_groups')
    return {k: upload[k] for k in keys}


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
            ],
            200,
            id='reviewer',
        ),
        pytest.param(
            'user2',
            {'roles': 'reviewer', 'include_all': 'True'},
            [
                'id_RGg2',
                'id_RGg123',
                'id_RGall',
            ],
            200,
            id='include-all',
        ),
    ],
)
def test_get_group_uploads(
    auth_headers,
    client,
    uploads_get_groups,
    user,
    query_params,
    expected_upload_ids,
    expected_status_code,
):
    response = perform_get(client, 'uploads', auth_headers[user], **query_params)
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
def test_get_group_upload_and_entries(
    auth_headers,
    client,
    uploads_get_groups,
    user,
    upload_id,
    expected_status_code,
):
    user_auth = auth_headers[user]
    response = perform_get(client, f'uploads/{upload_id}', user_auth)
    assert_response(response, expected_status_code)

    if expected_status_code != 200:
        return

    all_entries = uploads_get_groups.entries.items()
    reference_entries = {k: v for k, v in all_entries if v['upload_id'] == upload_id}
    assert_upload(response.json(), entries=len(reference_entries))

    response = perform_get(client, f'uploads/{upload_id}/entries', user_auth)
    assert_response(response, expected_status_code)

    response_data = response.json()['data']
    assert len(response_data) == len(reference_entries)
    for entry in response_data:
        assert_entry(entry, has_metadata=False, upload_id=upload_id)


@pytest.fixture(scope='function')
def perform_edit_upload_agents_test(
    auth_headers,
    client,
    convert_agent_labels_to_ids,
    proc_infra,
    groups_function,
):
    def perform(upload_fixture, user, metadata, expected_status_code, changed_agents):
        upload = list(upload_fixture.uploads.values())[0]
        expected_agents = get_agents_from_upload(upload)
        expected_agents.update(changed_agents)
        expected_agents = convert_agent_labels_to_ids(expected_agents)
        upload_id = upload['upload_id']

        url = f'uploads/{upload_id}/edit'
        metadata = convert_agent_labels_to_ids(metadata)
        edit_request = dict(metadata=metadata)
        response = perform_post(client, url, auth_headers[user], json=edit_request)
        assert_response(response, expected_status_code)

        upload = Upload.get(upload_id)
        upload.block_until_complete()
        agents = get_agents_from_upload(upload)
        assert sorted(agents) == sorted(expected_agents)

    return perform


edit_upload_agents_C_params = {
    'set-C': ({'coauthors': {'set': []}}, {'coauthors': []}),
    'set-C9': ({'coauthors': 'user9'}, {'coauthors': ['user9']}),
    'remove-C2': ({'coauthors': {'remove': 'user2'}}, {'coauthors': ['user4']}),
    'remove-C24': ({'coauthors': {'remove': ['user2', 'user4']}}, {'coauthors': []}),
    'add-C2-noop': ({'coauthors': {'add': ['user2']}}, {}),
    'add-C24-noop': ({'coauthors': {'add': ['user2', 'user4']}}, {}),
    'add-C8': (
        {'coauthors': {'add': 'user8'}},
        {'coauthors': ['user2', 'user4', 'user8']},
    ),
    'add-C88-merged': (
        {'coauthors': {'add': ['user8', 'user8']}},
        {'coauthors': ['user2', 'user4', 'user8']},
    ),
    'add-C89': (
        {'coauthors': {'add': ['user8', 'user9']}},
        {'coauthors': ['user2', 'user4', 'user8', 'user9']},
    ),
    'add-CX-fails': ({'coauthors': {'add': 'unknown_user'}}, 422),
    'remove-C2-add-C9': (
        {'coauthors': {'remove': 'user2', 'add': 'user9'}},
        {'coauthors': ['user4', 'user9']},
    ),
}

edit_upload_agents_CG_params = {
    'set-CG': ({'coauthor_groups': {'set': []}}, {'coauthor_groups': []}),
    'set-CGg9': (
        {'coauthor_groups': {'set': 'group9'}},
        {'coauthor_groups': ['group9']},
    ),
    'set-CGg9g19': (
        {'coauthor_groups': {'set': ['group9', 'group19']}},
        {'coauthor_groups': ['group9', 'group19']},
    ),
    'remove-CGg8': (
        {'coauthor_groups': {'remove': 'group8'}},
        {'coauthor_groups': ['group2', 'group14']},
    ),
    'remove-CGg8g14': (
        {'coauthor_groups': {'remove': ['group8', 'group14']}},
        {'coauthor_groups': ['group2']},
    ),
    'remove-CGg9-noop': ({'coauthor_groups': {'remove': 'group9'}}, {}),
    'add-CGg9': (
        {'coauthor_groups': {'add': 'group9'}},
        {'coauthor_groups': ['group2', 'group14', 'group9']},
    ),
    'add-CGg9g19': (
        {'coauthor_groups': {'add': ['group9', 'group19']}},
        {'coauthor_groups': ['group2', 'group14', 'group9', 'group19']},
    ),
    'add-CGg2g14-noop': ({'coauthor_groups': {'add': ['group2', 'group14']}}, {}),
    'add-CGgX-fails': ({'coauthor_groups': {'add': 'unknown_group'}}, 422),
    'set-CGall-fails': ({'coauthor_groups': {'set': 'all'}}, 422),
    'add-CGall-fails': ({'coauthor_groups': {'add': 'all'}}, 422),
    'add-CGg9-remove-CGg8': (
        {'coauthor_groups': {'add': 'group9', 'remove': 'group8'}},
        {'coauthor_groups': ['group2', 'group14', 'group9']},
    ),
    'add-CGX-fails': ({'coauthor_groups': {'add': 'unknown'}}, 422),
}

edit_upload_agents_RG_params = {
    'set-RGall': ({'reviewer_groups': {'set': []}}, {'reviewer_groups': []}),
    'add-RGall': (
        {'reviewer_groups': {'add': 'all'}},
        {'reviewer_groups': ['all', 'group3', 'group8', 'group15']},
    ),
    'remove-RGall-noop': ({'reviewer_groups': {'remove': 'all'}}, {}),
}


edit_upload_agents_params = {
    **edit_upload_agents_C_params,
    **edit_upload_agents_CG_params,
    **edit_upload_agents_RG_params,
}
edit_upload_agents_params = dict_to_params(edit_upload_agents_params)


@pytest.mark.parametrize(
    'metadata, code_or_changed_agents',
    edit_upload_agents_params,
)
def test_edit_upload_agents(
    perform_edit_upload_agents_test,
    upload_full_agents,
    metadata,
    code_or_changed_agents,
):
    if isinstance(code_or_changed_agents, int):
        expected_status_code = code_or_changed_agents
        changed_agents = {}
    else:
        expected_status_code = 200
        changed_agents = code_or_changed_agents

    perform_edit_upload_agents_test(
        upload_full_agents, 'user1', metadata, expected_status_code, changed_agents
    )


@pytest.mark.parametrize(
    'upload_id, user, expected_status_code',
    [
        pytest.param('id_full_agents', 'user2', 200, id='full-user2'),
        pytest.param('id_full_agents', 'user3', 422, id='full-user3-fails'),
        pytest.param('id_full_agents', 'user4', 200, id='full-user4'),
        pytest.param('id_full_agents', 'user5', 422, id='full-user5-fails'),
        pytest.param('id_full_agents', 'user6', 422, id='full-user6-fails'),
        pytest.param('id_full_agents', 'user7', 422, id='full-user7-fails'),
        pytest.param('id_full_agents', 'user8', 200, id='full-user8'),
        pytest.param('id_full_agents', 'user9', 422, id='full-user9-fails'),
        pytest.param('id_C2', 'user1', 200, id='C2-user1'),
        pytest.param('id_C2', 'user2', 200, id='C2-user2'),
        pytest.param('id_C2', 'user3', 422, id='C2-user3-fails'),
        pytest.param('id_R2', 'user1', 200, id='R2-user1'),
        pytest.param('id_R2', 'user2', 422, id='R2-user2-fails'),
        pytest.param('id_R2', 'user3', 422, id='R2-user3-fails'),
        pytest.param('id_CGg2', 'user1', 200, id='CGg2-user1'),
        pytest.param('id_CGg2', 'user2', 200, id='CGg2-user2'),
        pytest.param('id_CGg2', 'user3', 422, id='CGg2-user3-fails'),
        pytest.param('id_RGg2', 'user1', 200, id='RGg2-user1'),
        pytest.param('id_RGg2', 'user2', 422, id='RGg2-user2-fails'),
        pytest.param('id_RGg2', 'user3', 422, id='RGg2-user3-fails'),
        pytest.param('id_RGall', 'user1', 200, id='RGall-user1'),
        pytest.param('id_RGall', 'user2', 422, id='RGall-user2'),
        pytest.param('id_RGall', 'user3', 422, id='RGall-user3'),
    ],
)
def test_upload_agents_write_access(
    auth_headers,
    client,
    convert_agent_labels_to_ids,
    uploads_agent_write_access,
    upload_id,
    user,
    expected_status_code,
):
    url = f'uploads/{upload_id}/edit'
    metadata = convert_agent_labels_to_ids(
        {'coauthor_groups': ['group9'], 'reviewer_groups': ['group9']}
    )
    edit_request = dict(metadata=metadata, verify_only=True)
    response = perform_post(client, url, auth_headers[user], json=edit_request)
    assert_response(response, expected_status_code)
