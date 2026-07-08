#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

from typing import Any

from nomad.actions import TaskQueue
from nomad.actions.action import get_actions
from nomad.workflows.workflows import (
    BatchCleanupEntriesWorkflow,
    DeleteUploadWorkflow,
    EditUploadMetadataWorkflow,
    ImportBundleWorkflow,
    ProcessEntryWorkflow,
    ProcessExampleUploadWorkflow,
    PublishExternallyWorkflow,
    PublishUploadWorkflow,
    TransferUploadOwnershipWorkflow,
    UpdateUploadWorkflow,
)


def get_nomad_internal_workflows() -> list:
    return [
        BatchCleanupEntriesWorkflow,
        DeleteUploadWorkflow,
        UpdateUploadWorkflow,
        ProcessEntryWorkflow,
        ProcessExampleUploadWorkflow,
        EditUploadMetadataWorkflow,
        ImportBundleWorkflow,
        PublishUploadWorkflow,
        PublishExternallyWorkflow,
        TransferUploadOwnershipWorkflow,
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

    return list(set(workflows))
