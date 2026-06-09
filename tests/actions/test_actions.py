from unittest.mock import MagicMock, patch

import pytest

from nomad.actions import TaskQueue
from nomad.actions.action import Action, get_actions
from nomad.actions.activities.utils import get_all_activities
from nomad.config.models.plugins import ActionEntryPoint, ParserEntryPoint


class MockWorkflow:
    """A mock workflow class with a run method."""

    def run(self):
        pass


def cpu_activity():
    pass


def gpu_activity():
    pass


def create_mock_action_entry_point():
    mock_entry_point = MagicMock(spec=ActionEntryPoint)
    mock_entry_point.id = 'test-action'
    mock_entry_point.load.return_value = Action(
        task_queue='test-queue',
        activities=[lambda: None],
        workflow=MockWorkflow(),
    )
    return mock_entry_point


def create_mock_parser_entry_point():
    mock_entry_point = MagicMock(spec=ParserEntryPoint)
    mock_entry_point.id = 'test-parser'
    return mock_entry_point


def create_mock_action_wrong_return_type():
    mock_entry_point = MagicMock(spec=ActionEntryPoint)
    mock_entry_point.id = 'test-action-invalid-return'
    mock_entry_point.load.return_value = 'not-an-action'
    return mock_entry_point


def create_mock_action_load_exception():
    mock_entry_point = MagicMock(spec=ActionEntryPoint)
    mock_entry_point.id = 'test-action-load-exception'
    mock_entry_point.load.side_effect = Exception('Failed to load action')
    return mock_entry_point


@pytest.fixture(autouse=True)
def clear_get_actions_cache():
    """Clear the lru_cache before and after each test."""
    get_actions.cache_clear()
    yield
    get_actions.cache_clear()


mock_action_entry_point = create_mock_action_entry_point()
mock_parser_entry_point = create_mock_parser_entry_point()
mock_action_wrong_return_type = create_mock_action_wrong_return_type()
mock_action_load_exception = create_mock_action_load_exception()


@pytest.mark.parametrize(
    'entry_points,expected_result,exception',
    [
        pytest.param(
            [],
            {},
            None,
            id='empty-entry-points',
        ),
        pytest.param(
            [mock_action_entry_point],
            {'test-action': mock_action_entry_point},
            None,
            id='valid-action-entry-point',
        ),
        pytest.param(
            [mock_parser_entry_point],
            {},
            None,
            id='only-loads-actions',
        ),
        pytest.param(
            [mock_action_load_exception],
            {},
            'Failed to load action from entry point "test-action-load-exception"',
            id='raise-error-during-load',
        ),
        pytest.param(
            [mock_action_wrong_return_type],
            {},
            'The following entry point did not return an Action: "test-action-invalid-return"',
            id='raise-error-for-wrong-return-type',
        ),
    ],
)
def test_get_actions(
    entry_points,
    expected_result,
    exception,
):
    """Test get_actions with various entry point configurations."""
    with patch('nomad.config.config') as mock_config:
        mock_config.plugins.entry_points.filtered_values.return_value = entry_points

        if exception:
            with pytest.raises(Exception) as exc_info:
                get_actions()
            assert exception in str(exc_info.value)
        else:
            result = get_actions()
            assert result == expected_result


def test_action_accepts_task_queue_activities():
    action = Action(
        task_queue=TaskQueue.CPU,
        activities=[cpu_activity],
        workflow=MockWorkflow(),
        task_queue_activities={TaskQueue.GPU: [gpu_activity]},
    )

    assert action.task_queue_activities == {TaskQueue.GPU: [gpu_activity]}


def test_get_all_activities_includes_action_extra_task_queue_activities():
    action = Action(
        task_queue=TaskQueue.CPU,
        activities=[cpu_activity],
        workflow=MockWorkflow(),
        task_queue_activities={TaskQueue.GPU: [gpu_activity]},
    )
    mock_entry_point = MagicMock(spec=ActionEntryPoint)
    mock_entry_point.task_queue = TaskQueue.CPU
    mock_entry_point.load.return_value = action

    with patch(
        'nomad.actions.activities.utils.get_actions',
        return_value={'test-action': mock_entry_point},
    ):
        activities = get_all_activities(TaskQueue.GPU)

    assert gpu_activity in activities
    assert cpu_activity not in activities
