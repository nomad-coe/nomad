from typing import Any

from nomad.actions import TaskQueue
from nomad.actions.action import get_actions
from nomad.workflows.workflows import (
    BatchProcessEntriesWorkflow,
    DeleteUploadWorkflow,
    EditUploadMetadataWorkflow,
    ImportBundleWorkflow,
    ProcessEntryWorkflow,
    ProcessExampleUploadWorkflow,
    ProcessUploadWorkflow,
    PublishExternallyWorkflow,
    PublishUploadWorkflow,
)


def get_nomad_internal_workflows() -> list:
    return [
        BatchProcessEntriesWorkflow,
        DeleteUploadWorkflow,
        ProcessUploadWorkflow,
        ProcessEntryWorkflow,
        ProcessExampleUploadWorkflow,
        EditUploadMetadataWorkflow,
        ImportBundleWorkflow,
        PublishUploadWorkflow,
        PublishExternallyWorkflow,
    ]


def get_all_workflows(task_queue: TaskQueue) -> list:
    workflows: list[Any] = []

    for action_entry_point in get_actions().values():
        if action_entry_point.task_queue == task_queue:
            action = action_entry_point.load()
            workflows.append(action.workflow)
            workflows.extend(action.child_workflows)

    if task_queue == TaskQueue.NOMAD_INTERNAL_WORKFLOWS:
        workflows.extend(get_nomad_internal_workflows())

    return workflows
