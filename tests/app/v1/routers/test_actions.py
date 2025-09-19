from datetime import datetime
from unittest.mock import MagicMock, PropertyMock

import pytest
from fastapi.testclient import TestClient

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
    monkeypatch.setattr('nomad.actions.utils._update_status', lambda action: None)
    response = client.get('/actions/', headers=auth_headers['user1'])
    assert response.status_code == 200
    response_json = response.json()
    assert response_json[0]['action_id'] == saved_action_document.action_id
    assert 'results' not in response_json[0]
    assert 'input_data' not in response_json[0]


def test_get_action(
    client: TestClient, auth_headers, saved_action_document, monkeypatch
):
    monkeypatch.setattr('nomad.actions.utils._update_status', lambda action: None)
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


@pytest.mark.parametrize(
    'method, endpoint',
    [
        ('GET', '/actions/workflow-1/status'),
        ('POST', '/actions/workflow-1/stop'),
        ('GET', '/actions/workflow-1/result'),
        ('GET', '/actions/workflow-1'),
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
    monkeypatch.setattr('nomad.actions.utils._update_status', lambda action: None)
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
