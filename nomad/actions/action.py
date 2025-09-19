import functools
from collections.abc import Callable
from importlib.metadata import entry_points
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

    def __init__(
        self,
        task_queue: str,
        activities: list[Callable],
        workflow: _HasRun,
        child_workflows: list | None = None,
    ):
        """
        Initializes the Action with a task queue and a list of activities and a workflow, and child workflows.

        Args:
            task_queue: An instance of the TaskQueue associated with this action.
            activities: A list of callable functions that represent the activities
                        this handler can perform.
            workflow: The main workflow for this Action.
            child_workflows: Optionally, any child workflows of of the main workflow.
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
        self.task_queue = task_queue  # type: ignore
        self.activities = activities
        self.workflow = workflow
        self.child_workflows = child_workflows


@functools.lru_cache
def get_actions() -> dict[str, 'ActionEntryPoint']:
    """
    Loads and returns all valid Actions from the 'nomad.plugin' entry points.

    Raises:
        Exception: If an entry point fails to load.
        TypeError: If a loaded entry point doesn't return an Action.

    Returns:
        list[Action]: Loaded Action instances.
    """
    from nomad.config.models.plugins import ActionEntryPoint

    nomad_entry_points = entry_points(group='nomad.plugin')

    actions: dict[str, ActionEntryPoint] = {}
    invalid_entrypoints: list = []

    for plugin_entry_point in nomad_entry_points:
        entry_point = plugin_entry_point.load()
        if not isinstance(entry_point, ActionEntryPoint):
            continue

        try:
            handler = entry_point.load()
        except Exception as e:
            raise Exception(f'Failed to load entry point {entry_point}: {e}')

        if not isinstance(handler, Action):
            invalid_entrypoints.append(str(plugin_entry_point))
        else:
            actions[plugin_entry_point.value] = entry_point

    if invalid_entrypoints:
        raise TypeError(
            'The following entry points did not return an Action:\n'
            + '\n'.join(invalid_entrypoints)
        )

    return actions
