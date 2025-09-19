from importlib.metadata import EntryPoint
from unittest.mock import MagicMock

import pytest
from pydantic import BaseModel
from temporalio.client import WorkflowExecutionStatus

from nomad.actions.action import Action
from nomad.actions.utils import (
    _get_param_schema,
    _validate_with_pydantic,
    get_action_result,
    get_action_status,
    get_all_action_schemas,
    get_all_user_actions,
    start_action,
    validate_action_arg,
)
from nomad.mongo.action import ActionDocument


class MyActionArgs(BaseModel):
    arg1: str
    arg2: int
    user_id: str | None = None


async def my_workflow_run(args: MyActionArgs):
    pass


@pytest.fixture
def mock_action_entry_point():
    mock_action = MagicMock(spec=Action)
    mock_action.name = 'My Action'
    mock_action.description = 'My action description'
    mock_action.task_queue = 'my-task-queue'

    mock_workflow = MagicMock()
    mock_workflow.run = my_workflow_run
    mock_action.workflow = mock_workflow

    mock_entry_point = MagicMock(spec=EntryPoint)
    mock_entry_point.load.return_value = mock_action
    mock_entry_point.name = 'My Action'
    mock_entry_point.description = 'My action description'
    mock_entry_point.task_queue = 'my-task-queue'
    mock_entry_point.plugin_package = 'my-plugin'

    return mock_entry_point


def test_validate_with_pydantic():
    def my_func(args: MyActionArgs):
        pass

    validated = _validate_with_pydantic(my_func, {'arg1': 'test', 'arg2': 123})
    assert isinstance(validated, MyActionArgs)
    assert validated.arg1 == 'test'
    assert validated.arg2 == 123

    with pytest.raises(Exception):
        _validate_with_pydantic(my_func, {'arg1': 'test'})


def test_get_param_schema():
    def my_func(args: MyActionArgs):
        pass

    schema = _get_param_schema(my_func)
    assert 'arg1' in schema['properties']
    assert 'arg2' in schema['properties']


def test_validate_action_arg(monkeypatch, mock_action_entry_point):
    monkeypatch.setattr(
        'nomad.actions.utils.get_actions',
        lambda: {'my-action': mock_action_entry_point},
    )

    validated = validate_action_arg('my-action', {'arg1': 'test', 'arg2': 123})
    assert isinstance(validated, MyActionArgs)

    with pytest.raises(ValueError):
        validate_action_arg('nonexistent-action', {})


def test_get_all_action_schemas(monkeypatch, mock_action_entry_point):
    monkeypatch.setattr(
        'nomad.actions.utils.get_actions',
        lambda: {'my-action': mock_action_entry_point},
    )
    schemas = get_all_action_schemas()
    assert len(schemas) == 1
    assert schemas[0].action_id == 'my-action'
    assert 'arg1' in schemas[0].json_schema['properties']


def test_start_action(monkeypatch, mongo_function, user1, mock_action_entry_point):
    monkeypatch.setattr(
        'nomad.actions.utils.get_actions',
        lambda: {'my-action': mock_action_entry_point},
    )

    async def mock_async_start_workflow(*args, **kwargs):
        return 'workflow-id-start-test'

    monkeypatch.setattr(
        'nomad.actions.utils._async_start_workflow', mock_async_start_workflow
    )

    args = MyActionArgs(arg1='test', arg2=123, user_id=user1.user_id)
    action_instance_id = start_action('my-action', args)

    assert action_instance_id is not None


@pytest.fixture
def mock_temporal_client(monkeypatch):
    mock_client = MagicMock()
    mock_handle = MagicMock()

    async def mock_describe(*args, **kwargs):
        mock_status = MagicMock()
        mock_status.status = WorkflowExecutionStatus.RUNNING
        return mock_status

    async def mock_result(*args, **kwargs):
        return {'result': 'success'}

    mock_handle.describe = mock_describe
    mock_handle.result = mock_result

    mock_client.get_workflow_handle.return_value = mock_handle

    async def get_client():
        return mock_client

    monkeypatch.setattr('nomad.actions.utils.get_client', get_client)
    return mock_client


def test_get_action_status(mongo_function, user1, mock_temporal_client):
    action_doc = ActionDocument(
        action_id='my-action',
        action_instance_id='workflow-123',
        user_id=user1.user_id,
        status='PENDING',
        input_data={},
    )
    action_doc.save()

    status = get_action_status('workflow-123', user1.user_id)
    assert status == WorkflowExecutionStatus.RUNNING

    action_doc.reload()
    assert action_doc.status == 'RUNNING'

    with pytest.raises(Exception):
        get_action_status('nonexistent-workflow', user1.user_id)


def test_get_action_result(mongo_function, user1, mock_temporal_client):
    action_doc = ActionDocument(
        action_id='my-action',
        action_instance_id='workflow-123',
        user_id=user1.user_id,
        status='RUNNING',
        input_data={},
    )
    action_doc.save()

    result = get_action_result('workflow-123', user1.user_id)
    assert result == {'result': 'success'}

    action_doc.reload()
    assert action_doc.status == 'COMPLETED'
    assert action_doc.results == {'result': 'success'}

    with pytest.raises(Exception):
        get_action_result('nonexistent-workflow', user1.user_id)


def test_get_all_user_actions(monkeypatch, mongo_function, user1):
    ActionDocument(
        action_id='my-action-1',
        action_instance_id='workflow-1',
        user_id=user1.user_id,
        status='PENDING',
        input_data={},
    ).save()
    ActionDocument(
        action_id='my-action-2',
        action_instance_id='workflow-2',
        user_id=user1.user_id,
        status='COMPLETED',
        input_data={},
    ).save()

    def mock_update_status(action):
        pass

    monkeypatch.setattr('nomad.actions.utils._update_status', mock_update_status)

    actions = get_all_user_actions(user1.user_id)
    assert len(actions) == 2

    # Test no actions for user
    actions = get_all_user_actions('other-user')
    assert len(actions) == 0
