import os
from datetime import datetime
from unittest.mock import MagicMock, PropertyMock

import pytest
from fastapi.testclient import TestClient

from nomad.config import config
from nomad.mongo.action import ActionDocument


@pytest.fixture
def client(api_v1: TestClient) -> TestClient:
    return api_v1


@pytest.fixture
def saved_action_document(mongo_function, user1):
    action = ActionDocument(
        action_id='my-action',
        action_instance_id='workflow-1',
        status='RUNNING',
        user_id=user1.user_id,
        created_at=datetime.now(),
        updated_at=datetime.now(),
        input_data={},
        results={},
    )
    action.save()
    return action


def test_action_start(client: TestClient, auth_headers, monkeypatch):
    monkeypatch.setattr(
        'nomad.app.v1.routers.actions.validate_action_arg', lambda action_id, data: data
    )
    monkeypatch.setattr(
        'nomad.app.v1.routers.actions.start_action',
        lambda action_id, data: 'workflow-123',
    )

    response = client.post(
        '/actions/my-action/start',
        json={'data': {'arg1': 'test'}},
        headers=auth_headers['user1'],
    )

    assert response.status_code == 200
    assert response.json() == {'action_instance_id': 'workflow-123'}


def test_action_status(
    client: TestClient, auth_headers, saved_action_document, monkeypatch
):
    mock_status = MagicMock()
    type(mock_status).name = PropertyMock(return_value='RUNNING')
    monkeypatch.setattr(
        'nomad.app.v1.routers.actions.get_action_status',
        lambda action_instance_id, user_id: mock_status,
    )
    response = client.get(
        f'/actions/{saved_action_document.action_instance_id}/status',
        headers=auth_headers['user1'],
    )
    assert response.status_code == 200
    assert response.json() == {'status': 'RUNNING'}


def test_action_result(
    client: TestClient, auth_headers, saved_action_document, monkeypatch
):
    monkeypatch.setattr(
        'nomad.app.v1.routers.actions.get_action_result',
        lambda action_instance_id, user_id: {'result': 'success'},
    )
    response = client.get(
        f'/actions/{saved_action_document.action_instance_id}/result',
        headers=auth_headers['user1'],
    )
    assert response.status_code == 200
    assert response.json() == {'result': 'success'}


def test_action_input_schemas(
    client: TestClient, auth_headers, monkeypatch, fastapi_cache
):
    monkeypatch.setattr(
        'nomad.app.v1.routers.actions.get_all_action_schemas',
        lambda: [{'action_id': 'my-action', 'json_schema': {}}],
    )
    response = client.get('/actions/schemas', headers=auth_headers['user1'])
    assert response.status_code == 200
    assert response.json()[0]['action_id'] == 'my-action'


def test_actions_list(
    client: TestClient, auth_headers, saved_action_document, monkeypatch
):
    monkeypatch.setattr('nomad.actions.manager._update_status', lambda action: None)
    response = client.get('/actions/', headers=auth_headers['user1'])
    assert response.status_code == 200
    response_json = response.json()
    assert response_json[0]['action_id'] == saved_action_document.action_id
    assert 'results' not in response_json[0]
    assert 'input_data' not in response_json[0]


def test_get_action(
    client: TestClient, auth_headers, saved_action_document, monkeypatch
):
    monkeypatch.setattr('nomad.actions.manager._update_status', lambda action: None)
    response = client.get(
        f'/actions/{saved_action_document.action_instance_id}',
        headers=auth_headers['user1'],
    )
    assert response.status_code == 200
    response_json = response.json()
    assert response_json['action_id'] == saved_action_document.action_id
    assert (
        response_json['action_instance_id'] == saved_action_document.action_instance_id
    )
    assert response_json['status'] == saved_action_document.status


def test_get_action_not_found(client: TestClient, auth_headers, mongo_function):
    response = client.get('/actions/workflow-2', headers=auth_headers['user1'])
    assert response.status_code == 404


def test_action_stop(
    client: TestClient, auth_headers, saved_action_document, monkeypatch
):
    monkeypatch.setattr(
        'nomad.app.v1.routers.actions.stop_action',
        lambda action_instance_id, user_id: None,
    )
    response = client.post(
        f'/actions/{saved_action_document.action_instance_id}/stop',
        headers=auth_headers['user1'],
    )
    assert response.status_code == 200
    assert response.json() == {'status': 'stopped'}


def test_action_logs(
    client: TestClient, auth_headers, saved_action_document, monkeypatch, tmp_path
):
    # Mock get_user_action so authorization passes (and bypasses Temporal fetching)
    monkeypatch.setattr('nomad.actions.manager._update_status', lambda action: None)

    # Setup dummy log file
    log_dir = os.path.join(config.fs.actions, 'logs')
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f'{saved_action_document.action_instance_id}.log')
    with open(log_file, 'w') as f:
        f.write('Test log line 1\nTest log line 2\n')

    try:
        response = client.get(
            f'/actions/{saved_action_document.action_instance_id}/logs',
            headers=auth_headers['user1'],
        )
        assert response.status_code == 200
        assert response.text == 'Test log line 1\nTest log line 2\n'
        assert response.headers['content-type'] == 'text/plain; charset=utf-8'
    finally:
        if os.path.exists(log_file):
            os.remove(log_file)


def test_action_logs_not_found(
    client: TestClient, auth_headers, saved_action_document, monkeypatch
):
    # Mock get_user_action so authorization passes
    monkeypatch.setattr('nomad.actions.manager._update_status', lambda action: None)

    response = client.get(
        f'/actions/{saved_action_document.action_instance_id}/logs',
        headers=auth_headers['user1'],
    )
    assert response.status_code == 404
    assert response.json()['detail'] == 'Log file not found for this action.'


def test_action_logs_truncate(
    client: TestClient, auth_headers, saved_action_document, monkeypatch
):
    # Mock get_user_action so authorization passes
    monkeypatch.setattr('nomad.actions.manager._update_status', lambda action: None)

    # Setup dummy log file
    log_dir = os.path.join(config.fs.actions, 'logs')
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f'{saved_action_document.action_instance_id}.log')

    # Create a file slightly larger than 2MB
    chunk = b'A' * 1024 * 1024
    with open(log_file, 'wb') as f:
        f.write(chunk)
        f.write(chunk)
        f.write(b'B' * 10)

    try:
        response = client.get(
            f'/actions/{saved_action_document.action_instance_id}/logs',
            headers=auth_headers['user1'],
        )
        assert response.status_code == 200

        notice_len = len(
            b'[... earlier content truncated due to file size limit. Contact admin for the full log file ...]\n\n'
        )
        assert len(response.content) == notice_len + 2 * 1024 * 1024
        assert response.content.endswith(b'B' * 10)
    finally:
        if os.path.exists(log_file):
            os.remove(log_file)


def test_action_logs_stream(
    client: TestClient, auth_headers, saved_action_document, monkeypatch
):
    # Mock get_user_action so authorization passes
    monkeypatch.setattr('nomad.actions.manager._update_status', lambda action: None)

    # Setup dummy log file
    log_dir = os.path.join(config.fs.actions, 'logs')
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f'{saved_action_document.action_instance_id}.log')

    with open(log_file, 'w') as f:
        f.write('Initial line 1\nInitial line 2\n')

    status_calls = 0

    def mock_status(*args, **kwargs):
        nonlocal status_calls
        status_calls += 1
        if status_calls == 1:
            with open(log_file, 'a') as f:
                f.write('New streaming line 1\nNew streaming line 2\n')

        mock_obj = MagicMock()
        type(mock_obj).name = PropertyMock(
            return_value='RUNNING' if status_calls <= 2 else 'SUCCESS'
        )
        return mock_obj

    monkeypatch.setattr('nomad.app.v1.routers.actions.get_action_status', mock_status)

    try:
        response = client.get(
            f'/actions/{saved_action_document.action_instance_id}/logs?stream=true',
            headers=auth_headers['user1'],
        )
        assert response.status_code == 200
        assert response.headers['content-type'] == 'text/event-stream; charset=utf-8'

        # Since it streams from the end, the initial lines shouldn't be there.
        assert 'Initial line' not in response.text
        assert 'New streaming line 1' in response.text
        assert 'New streaming line 2' in response.text
    finally:
        if os.path.exists(log_file):
            os.remove(log_file)


@pytest.mark.parametrize(
    'method, endpoint',
    [
        ('GET', '/actions/workflow-1/status'),
        ('POST', '/actions/workflow-1/stop'),
        ('GET', '/actions/workflow-1/result'),
        ('GET', '/actions/workflow-1'),
        ('GET', '/actions/workflow-1/logs'),
    ],
)
def test_action_endpoints_unauthorized(client: TestClient, method: str, endpoint: str):
    response = client.request(method, endpoint)
    assert response.status_code == 401


def mock_get_action_status_raise(*args, **kwargs):
    raise Exception('Action status not found')


def mock_stop_action_raise(*args, **kwargs):
    raise Exception(
        'The action was not registered in the DB or was registered under a different user.'
    )


def mock_get_action_result_raise(*args, **kwargs):
    raise Exception('Action result not found.')


@pytest.mark.parametrize(
    'method, endpoint, mock_function_name, mock_function, expected_status_code',
    [
        (
            'GET',
            '/actions/workflow-1/status',
            'nomad.app.v1.routers.actions.get_action_status',
            mock_get_action_status_raise,
            500,
        ),
        (
            'POST',
            '/actions/workflow-1/stop',
            'nomad.app.v1.routers.actions.stop_action',
            mock_stop_action_raise,
            500,
        ),
        (
            'GET',
            '/actions/workflow-1/result',
            'nomad.app.v1.routers.actions.get_action_result',
            mock_get_action_result_raise,
            500,
        ),
        ('GET', '/actions/workflow-1', None, None, 404),
        ('GET', '/actions/workflow-1/logs', None, None, 404),
    ],
)
def test_action_endpoints_wrong_user(
    client: TestClient,
    auth_headers,
    monkeypatch,
    method: str,
    endpoint: str,
    mock_function_name: str,
    mock_function,
    expected_status_code: int,
    saved_action_document,
):
    if mock_function_name:
        monkeypatch.setattr(mock_function_name, mock_function)
    response = client.request(method, endpoint, headers=auth_headers['user2'])
    assert response.status_code == expected_status_code


# Test for GET /actions not containing other users' actions
def test_actions_list_does_not_contain_other_users_actions(
    client: TestClient, auth_headers, mongo_function, user1, user2, monkeypatch
):
    monkeypatch.setattr('nomad.actions.manager._update_status', lambda action: None)
    # user1 has one action, user2 has another
    ActionDocument(
        action_id='action1',
        action_instance_id='wf1',
        status='RUNNING',
        created_at=datetime.now(),
        updated_at=datetime.now(),
        user_id=user1.user_id,
        input_data={},
        results={},
    ).save()
    ActionDocument(
        action_id='action2',
        action_instance_id='wf2',
        status='RUNNING',
        created_at=datetime.now(),
        updated_at=datetime.now(),
        user_id=user2.user_id,
        input_data={},
        results={},
    ).save()

    # request as user1
    response = client.get('/actions/', headers=auth_headers['user1'])
    assert response.status_code == 200
    response_json = response.json()
    assert len(response_json) == 1
    assert response_json[0]['action_id'] == 'action1'

    # request as user2
    response = client.get('/actions/', headers=auth_headers['user2'])
    assert response.status_code == 200
    response_json = response.json()
    assert len(response_json) == 1
    assert response_json[0]['action_id'] == 'action2'
