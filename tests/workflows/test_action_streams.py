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

from __future__ import annotations

from datetime import timedelta
from uuid import uuid4

import pytest
from temporalio import workflow
from temporalio.contrib.workflow_streams import WorkflowStream
from temporalio.worker import Worker

from nomad.actions import streams
from nomad.actions.models import ActionStreamEvent, ActionStreamEventType
from nomad.actions.streams import ACTION_STREAM_TOPIC


@workflow.defn(name='tests.workflows.StreamingTestWorkflow', sandboxed=False)
class StreamingTestWorkflow:
    @workflow.init
    def __init__(self, _data: dict):
        self.stream = WorkflowStream()
        self.events = self.stream.topic(ACTION_STREAM_TOPIC, type=ActionStreamEvent)

    @workflow.run
    async def run(self, _data: dict) -> str:
        self.events.publish(
            ActionStreamEvent(
                type=ActionStreamEventType.STATE,
                name='started',
                message='Workflow started',
            )
        )
        self.events.publish(
            ActionStreamEvent(
                type=ActionStreamEventType.MESSAGE,
                message='halfway there',
            )
        )
        await workflow.sleep(timedelta(milliseconds=1))
        return 'ok'


@pytest.mark.asyncio
async def test_stream_workflow_events_async_with_real_worker(
    monkeypatch, temporal_worker
):
    workflow_id = f'streaming-test-{uuid4()}'
    task_queue = f'test-stream-queue-{uuid4()}'

    async with temporal_worker() as env:

        async def _mock_get_client():
            return env.client

        monkeypatch.setattr('nomad.actions.streams.get_client', _mock_get_client)

        async with Worker(
            env.client,
            task_queue=task_queue,
            workflows=[StreamingTestWorkflow],
            activities=[],
        ):
            handle = await env.client.start_workflow(
                StreamingTestWorkflow.run,
                dict(user_id='u-1'),
                id=workflow_id,
                task_queue=task_queue,
            )

            stream = await streams._stream_workflow_events_async(
                workflow_id=workflow_id,
                topic=ACTION_STREAM_TOPIC,
                from_offset=0,
            )
            first = await anext(stream)
            second = await anext(stream)

            assert first.topic == ACTION_STREAM_TOPIC
            assert first.event.type == ActionStreamEventType.STATE
            assert first.event.name == 'started'
            assert second.event.type == ActionStreamEventType.MESSAGE

            assert await handle.result() == 'ok'
