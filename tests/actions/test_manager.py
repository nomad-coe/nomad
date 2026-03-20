from datetime import datetime, timedelta, timezone
from importlib.metadata import EntryPoint
from unittest.mock import MagicMock

import pytest
from pydantic import BaseModel
from temporalio.client import WorkflowExecutionStatus

from nomad.actions.action import Action
from nomad.actions.manager import (
    _get_param_schema,
    _validate_with_pydantic,
    get_action_result,
    get_action_status,
    get_all_action_schemas,
    list_user_actions,
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

    class DummyWorkflow:
        async def run(self, args: MyActionArgs):
            pass

    mock_action.workflow = DummyWorkflow

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
        'nomad.actions.manager.get_actions',
        lambda: {'my-action': mock_action_entry_point},
    )

    validated = validate_action_arg('my-action', {'arg1': 'test', 'arg2': 123})
    assert isinstance(validated, MyActionArgs)

    with pytest.raises(ValueError):
        validate_action_arg('nonexistent-action', {})


def test_get_all_action_schemas(monkeypatch, mock_action_entry_point):
    monkeypatch.setattr(
        'nomad.actions.manager.get_actions',
        lambda: {'my-action': mock_action_entry_point},
    )
    schemas = get_all_action_schemas()
    assert len(schemas) == 1
    assert schemas[0].action_id == 'my-action'
    assert 'arg1' in schemas[0].json_schema['properties']


@pytest.mark.asyncio
async def test_start_action(
    monkeypatch, mongo_function, async_mongo_function, user1, mock_action_entry_point
):
    monkeypatch.setattr(
        'nomad.actions.manager.get_actions',
        lambda: {'my-action': mock_action_entry_point},
    )

    async def mock_async_start_workflow(*args, **kwargs):
        return 'workflow-id-start-test'

    monkeypatch.setattr(
        'nomad.actions.manager._async_start_workflow', mock_async_start_workflow
    )

    args = MyActionArgs(arg1='test', arg2=123, user_id=user1.user_id)
    action_instance_id = await start_action('my-action', args)

    assert action_instance_id is not None

    # Verify the document was actually created in the DB
    doc = await ActionDocument.find_one(
        ActionDocument.action_instance_id == action_instance_id
    )
    assert doc is not None
    assert doc.action_id == 'my-action'
    assert doc.status == 'PENDING'


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

    monkeypatch.setattr('nomad.actions.manager.get_client', get_client)
    return mock_client


@pytest.mark.asyncio
async def test_get_action_status(
    mongo_function, async_mongo_function, user1, mock_temporal_client
):
    action_doc = ActionDocument(
        action_id='my-action',
        action_instance_id='workflow-123',
        user_id=user1.user_id,
        status='PENDING',
        input_data={},
    )
    await action_doc.insert()

    status = await get_action_status('workflow-123', user1.user_id)
    assert status == WorkflowExecutionStatus.RUNNING

    # Verify the status was updated in the DB
    updated_doc = await ActionDocument.find_one(
        ActionDocument.action_instance_id == 'workflow-123'
    )
    assert updated_doc.status == 'RUNNING'

    with pytest.raises(Exception):
        await get_action_status('nonexistent-workflow', user1.user_id)


@pytest.mark.asyncio
async def test_get_action_result(
    mongo_function, async_mongo_function, user1, mock_temporal_client
):
    action_doc = ActionDocument(
        action_id='my-action',
        action_instance_id='workflow-123',
        user_id=user1.user_id,
        status='RUNNING',
        input_data={},
    )
    await action_doc.insert()

    result = await get_action_result('workflow-123', user1.user_id)
    assert result == {'result': 'success'}

    # Verify results were saved to DB
    updated_doc = await ActionDocument.find_one(
        ActionDocument.action_instance_id == 'workflow-123'
    )
    assert updated_doc.status == 'COMPLETED'
    assert updated_doc.results == {'result': 'success'}

    with pytest.raises(Exception):
        await get_action_result('nonexistent-workflow', user1.user_id)


@pytest.mark.asyncio
async def test_list_user_actions(
    monkeypatch, mongo_function, async_mongo_function, user1
):
    await ActionDocument(
        action_id='my-action-1',
        action_instance_id='workflow-1',
        user_id=user1.user_id,
        status='PENDING',
        input_data={},
    ).insert()
    await ActionDocument(
        action_id='my-action-2',
        action_instance_id='workflow-2',
        user_id=user1.user_id,
        status='COMPLETED',
        input_data={},
    ).insert()

    async def mock_update_status(action):
        pass

    monkeypatch.setattr('nomad.actions.manager._update_status', mock_update_status)

    page = await list_user_actions(user1.user_id)
    assert len(page.items) == 2
    assert page.total == 2

    # Test no actions for user
    empty_page = await list_user_actions('other-user')
    assert len(empty_page.items) == 0
    assert empty_page.total == 0


@pytest.mark.asyncio
async def test_list_user_actions_page_size_one_cursor_chain(
    monkeypatch, mongo_function, async_mongo_function, user1
):
    base_time = datetime(2024, 1, 1, tzinfo=timezone.utc)
    for i in range(3):
        await ActionDocument(
            action_id=f'my-action-{i}',
            action_instance_id=f'workflow-{i}',
            user_id=user1.user_id,
            status='COMPLETED',
            input_data={},
            created_at=base_time + timedelta(seconds=i),
            updated_at=base_time + timedelta(seconds=i),
        ).insert()

    async def mock_update_status(action):
        pass

    monkeypatch.setattr('nomad.actions.manager._update_status', mock_update_status)

    first_page = await list_user_actions(user1.user_id, page_size=1)
    assert first_page.total == 3
    assert len(first_page.items) == 1
    assert first_page.next_cursor is not None
    assert first_page.items[0].action_id == 'my-action-2'

    second_page = await list_user_actions(
        user1.user_id, page_size=1, cursor=first_page.next_cursor
    )
    assert second_page.total == 3
    assert len(second_page.items) == 1
    assert second_page.next_cursor is not None
    assert second_page.items[0].action_id == 'my-action-1'

    third_page = await list_user_actions(
        user1.user_id, page_size=1, cursor=second_page.next_cursor
    )
    assert third_page.total == 3
    assert len(third_page.items) == 1
    assert third_page.next_cursor is None
    assert third_page.items[0].action_id == 'my-action-0'


@pytest.mark.asyncio
async def test_list_user_actions_filters_by_upload_id(
    monkeypatch, mongo_function, async_mongo_function, user1
):
    await ActionDocument(
        action_id='my-action-1',
        action_instance_id='workflow-1',
        user_id=user1.user_id,
        upload_id='upload-a',
        status='COMPLETED',
        input_data={},
    ).insert()
    await ActionDocument(
        action_id='my-action-2',
        action_instance_id='workflow-2',
        user_id=user1.user_id,
        upload_id='upload-b',
        status='COMPLETED',
        input_data={},
    ).insert()
    await ActionDocument(
        action_id='my-action-3',
        action_instance_id='workflow-3',
        user_id=user1.user_id,
        upload_id='upload-a',
        status='COMPLETED',
        input_data={},
    ).insert()

    async def mock_update_status(action):
        pass

    monkeypatch.setattr('nomad.actions.manager._update_status', mock_update_status)

    filtered_page = await list_user_actions(user1.user_id, upload_id='upload-a')
    assert filtered_page.total == 2
    assert len(filtered_page.items) == 2
    assert {item.upload_id for item in filtered_page.items} == {'upload-a'}
