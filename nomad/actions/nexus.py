from typing import Any

from nomad.actions import TaskQueue
from nomad.actions.action import get_actions


def get_all_nexus_service_handlers(task_queue: TaskQueue) -> list[Any]:
    """Return plugin-provided Nexus handlers registered on a task queue."""
    handlers: list[Any] = []

    for action_entry_point in get_actions().values():
        if action_entry_point.task_queue == task_queue:
            handlers.extend(action_entry_point.load().nexus_service_handlers)

    return handlers
