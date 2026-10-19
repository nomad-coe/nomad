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

import tempfile
import uuid

import pytest

from nomad.actions import TaskQueue
from nomad.processing.base import ProcessStatus
from nomad.processing.data import Entry, Upload
from nomad.workflows.shared_objects import (
    DeleteUploadWorkflowInput,
    PublishUploadWorkflowInput,
    UploadProcessingWorkflowInput,
)


@pytest.mark.asyncio
async def test_process_upload_workflow_integration(
    non_empty_uploaded,
    user1,
    temporal_worker,
):
    """
    Integration test: run the full ProcessUploadWorkflow on a real upload and check side effects.
    """
    upload_id, upload_path = non_empty_uploaded
    upload_id = f'{upload_id}-{uuid.uuid4()}'
    upload = Upload.create(upload_id=upload_id, main_author=user1)
    upload.save()
    assert Upload.get(upload_id) is not None
    input_data = UploadProcessingWorkflowInput(
        upload_id=upload_id,
        file_operations=[
            dict(
                op='ADD',
                path=upload_path,
                target_dir='',
                temporary=False,
            )
        ],
        reprocess_settings=None,
        path_filter=None,
        only_updated_files=False,
        publish_directly_after_processing=False,
        workflow_id=f'integration-test-workflow-{upload_id}',
        workflow_tmp_dir=tempfile.mkdtemp(),
    )

    async with temporal_worker() as env:
        await env.client.execute_workflow(
            'UpdateUploadWorkflow',
            input_data,
            id=input_data.workflow_id,
            task_queue=TaskQueue.NOMAD_INTERNAL_WORKFLOWS,
        )

    # Check upload is processed
    upload.reload()
    assert upload.process_status == ProcessStatus.SUCCESS
    assert not upload.process_running
    assert len(upload.errors) == 0
    # Check entries are processed
    for entry in Entry.objects(upload_id=upload_id):
        assert entry.process_status == ProcessStatus.SUCCESS
        assert len(entry.errors) == 0


@pytest.mark.asyncio
async def test_publish_upload_workflow_integration(
    non_empty_uploaded,
    user1,
    temporal_worker,
):
    """
    Integration test: run the full PublishUploadWorkflow on a real upload and check side effects.
    """
    upload_id, upload_path = non_empty_uploaded
    upload_id = f'{upload_id}-{uuid.uuid4()}'
    upload = Upload.create(upload_id=upload_id, main_author=user1)
    upload.save()
    assert Upload.get(upload_id) is not None
    process_input = UploadProcessingWorkflowInput(
        upload_id=upload_id,
        file_operations=[
            dict(
                op='ADD',
                path=upload_path,
                target_dir='',
                temporary=False,
            )
        ],
        reprocess_settings=None,
        path_filter=None,
        only_updated_files=False,
        publish_directly_after_processing=False,
        workflow_id=f'integration-test-workflow-process-{upload_id}',
        workflow_tmp_dir=tempfile.mkdtemp(),
    )

    async with temporal_worker() as env:
        await env.client.execute_workflow(
            'UpdateUploadWorkflow',
            process_input,
            id=process_input.workflow_id,
            task_queue=TaskQueue.NOMAD_INTERNAL_WORKFLOWS,
        )

        upload.reload()
        assert upload.process_status == ProcessStatus.SUCCESS

        publish_input = PublishUploadWorkflowInput(
            upload_id=upload_id, embargo_length=0
        )
        await env.client.execute_workflow(
            'PublishUploadWorkflow',
            publish_input,
            id=f'integration-test-workflow-publish-{upload_id}',
            task_queue=TaskQueue.NOMAD_INTERNAL_WORKFLOWS,
        )

    # Check upload is published
    upload.reload()
    assert upload.published is True
    assert upload.process_status == ProcessStatus.SUCCESS
    assert not upload.process_running
    assert len(upload.errors) == 0


@pytest.mark.asyncio
async def test_delete_upload_workflow_integration(
    user1,
    temporal_worker,
):
    """
    Integration test: run the full DeleteUploadWorkflow on a real upload and check side effects.
    """
    upload_id = str(uuid.uuid4())
    upload = Upload.create(upload_id=upload_id, main_author=user1)
    upload.save()
    assert Upload.get(upload_id) is not None

    delete_input = DeleteUploadWorkflowInput(upload_id=upload_id)

    async with temporal_worker() as env:
        await env.client.execute_workflow(
            'DeleteUploadWorkflow',
            delete_input,
            id=f'integration-test-workflow-delete-{upload_id}',
            task_queue=TaskQueue.NOMAD_INTERNAL_WORKFLOWS,
        )

    with pytest.raises(Exception):  # Should raise KeyError
        Upload.get(upload_id)


@pytest.mark.asyncio
async def test_publish_with_embargo_and_lift(
    non_empty_uploaded,
    user1,
    temporal_worker,
):
    """
    Test publish with embargo, then lift embargo, and check transitions.
    """
    upload_id, upload_path = non_empty_uploaded
    upload_id = f'{upload_id}-{uuid.uuid4()}'
    upload = Upload.create(upload_id=upload_id, main_author=user1)
    upload.save()
    process_input = UploadProcessingWorkflowInput(
        upload_id=upload_id,
        file_operations=[
            dict(op='ADD', path=upload_path, target_dir='', temporary=False)
        ],
        workflow_id=f'integration-test-workflow-process-{upload_id}',
        workflow_tmp_dir=tempfile.mkdtemp(),
    )

    async with temporal_worker() as env:
        await env.client.execute_workflow(
            'UpdateUploadWorkflow',
            process_input,
            id=process_input.workflow_id,
            task_queue=TaskQueue.NOMAD_INTERNAL_WORKFLOWS,
        )
        upload.reload()
        assert upload.process_status == ProcessStatus.SUCCESS

        # Publish with embargo
        publish_input = PublishUploadWorkflowInput(
            upload_id=upload_id, embargo_length=12
        )
        await env.client.execute_workflow(
            'PublishUploadWorkflow',
            publish_input,
            id=f'integration-test-workflow-publish-{upload_id}',
            task_queue=TaskQueue.NOMAD_INTERNAL_WORKFLOWS,
        )
        upload.reload()
        assert upload.with_embargo

        # Lift embargo
        publish_input = PublishUploadWorkflowInput(
            upload_id=upload_id, embargo_length=0
        )
        await env.client.execute_workflow(
            'PublishUploadWorkflow',
            publish_input,
            id=f'integration-test-workflow-lift-embargo-{upload_id}',
            task_queue=TaskQueue.NOMAD_INTERNAL_WORKFLOWS,
        )
        upload.reload()
        assert not upload.with_embargo
