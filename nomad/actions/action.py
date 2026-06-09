import functools
from collections.abc import Callable
from typing import TYPE_CHECKING, Protocol

from nomad.actions import TaskQueue

if TYPE_CHECKING:
    from nomad.config.models.plugins import ActionEntryPoint


class _HasRun(Protocol):
    run: ...  # type: ignore


class Action:
    """
    Base class for handling a set of activities and workflows associated with a task queue.
    """

    task_queue: TaskQueue
    activities: list[Callable]
    workflow: _HasRun
    child_workflows: list
    task_queue_activities: dict[TaskQueue | str, list[Callable]]

    def __init__(
        self,
        task_queue: str,
        activities: list[Callable],
        workflow: _HasRun,
        child_workflows: list | None = None,
        task_queue_activities: dict[TaskQueue | str, list[Callable]] | None = None,
    ):
        """
        Initializes the Action with a task queue and a list of activities and a workflow, and child workflows.

        Args:
            task_queue: An instance of the TaskQueue associated with this action.
            activities: A list of callable functions that represent the activities
                        this handler can perform.
            workflow: The main workflow for this Action.
            child_workflows: Optionally, any child workflows of of the main workflow.
            task_queue_activities: Additional activities that should be registered
                                   on other task queues for this action.
        Raises:
            TypeError: If task_queue is not a TaskQueue or activities is not a list
                       of callables.
        """
        if not isinstance(activities, list) or not all(callable(a) for a in activities):
            raise TypeError('activities must be a list of callable functions')
        if not isinstance(workflow, object):
            raise TypeError('workflow must be a list')
        child_workflows = child_workflows if child_workflows is not None else []
        if not isinstance(child_workflows, list):
            raise TypeError('child_workflow must be a list')
        task_queue_activities = (
            task_queue_activities if task_queue_activities is not None else {}
        )
        if not isinstance(task_queue_activities, dict):
            raise TypeError('task_queue_activities must be a dict')
        if not all(
            isinstance(queue_activities, list)
            and all(callable(activity) for activity in queue_activities)
            for queue_activities in task_queue_activities.values()
        ):
            raise TypeError(
                'task_queue_activities values must be lists of callable functions'
            )
        self.task_queue = task_queue  # type: ignore
        self.activities = activities
        self.workflow = workflow
        self.child_workflows = child_workflows
        self.task_queue_activities = task_queue_activities


@functools.lru_cache
def get_actions() -> dict[str, 'ActionEntryPoint']:
    """
    Loads and returns all valid and available Actions from the nomad plugin entry points.

    Raises:
        Exception: If an entry point fails to load.
        TypeError: If a loaded entry point doesn't return an Action.

    Returns:
        Dictionary containing the entry point ID as key and ActionEntryPoint as value.
    """
    from nomad.config import config
    from nomad.config.models.plugins import ActionEntryPoint

    config.load_plugins()
    nomad_entry_points = config.plugins.entry_points.filtered_values()
    actions: dict[str, ActionEntryPoint] = {}

    for entry_point in nomad_entry_points:
        if not isinstance(entry_point, ActionEntryPoint):
            continue

        try:
            action = entry_point.load()
        except Exception as e:
            raise Exception(
                f'Failed to load action from entry point "{entry_point.id}"'
            ) from e

        if not isinstance(action, Action):
            raise TypeError(
                f'The following entry point did not return an Action: "{entry_point.id}"'
            )
        actions[entry_point.id] = entry_point

    return actions
