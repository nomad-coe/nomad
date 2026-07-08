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

import asyncio
import multiprocessing
import signal
import sys
from datetime import timedelta
from typing import Any

from temporalio.worker import (
    ResourceBasedSlotConfig,
    SharedStateManager,
    Worker,
    WorkerTuner,
)

from nomad.actions import TaskQueue
from nomad.actions.activities.utils import get_all_activities
from nomad.actions.client import close_client, get_client
from nomad.actions.nexus import get_all_nexus_service_handlers
from nomad.actions.workers.health_check import (
    should_start_health_server,
    start_health_server,
)
from nomad.actions.workflows.utils import get_all_workflows
from nomad.config import config
from nomad.config.models.config import WorkerConfig
from nomad.tracing import TracingProcessPoolExecutor
from nomad.utils.structlogging import get_logger

from .utils import worker_process_initializer


def _internal_worker_process_warmup():
    # Internal activities run in ProcessPoolExecutor children; pre-warm parser imports
    # in each child so the first workflow task does not absorb that cold-start latency.
    return


async def run_worker(worker_config: WorkerConfig):
    logger = get_logger(__name__)
    loop = asyncio.get_running_loop()
    stop_event = asyncio.Event()

    def _signal_handler():
        # Handle graceful shutdown
        logger.info('Received SIGTERM. Preparing for graceful shutdown')
        stop_event.set()

    if sys.platform == 'win32':
        signal.signal(signal.SIGTERM, lambda s, f: _signal_handler())
        signal.signal(signal.SIGINT, lambda s, f: _signal_handler())
    else:
        loop.add_signal_handler(signal.SIGTERM, _signal_handler)
        loop.add_signal_handler(signal.SIGINT, _signal_handler)

    client = await get_client()
    executor_kwargs = {
        'max_workers': worker_config.pool_size,
        'initializer': worker_process_initializer,
    }
    if sys.version_info >= (3, 11) and worker_config.max_tasks_per_child is not None:
        executor_kwargs['max_tasks_per_child'] = worker_config.max_tasks_per_child

    # NOTE: internal processing is not thread safe, avoid using ThreadPoolExecutor with more than 1 worker.
    # mypy: has issues with **kwargs in this context
    with TracingProcessPoolExecutor(**executor_kwargs) as executor:  # type: ignore
        # ProcessPoolExecutor starts children lazily. Pre-start all children here so
        # the first queued workflow does not pay startup/initializer latency.
        warmup_futures = [
            executor.submit(_internal_worker_process_warmup)
            for _ in range(worker_config.pool_size)
        ]
        for future in warmup_futures:
            future.result()

        worker_kwargs: dict[str, Any] = {
            'client': client,
            'task_queue': TaskQueue.NOMAD_INTERNAL_WORKFLOWS,
            'workflows': get_all_workflows(TaskQueue.NOMAD_INTERNAL_WORKFLOWS),
            'activities': get_all_activities(TaskQueue.NOMAD_INTERNAL_WORKFLOWS),
            'nexus_service_handlers': get_all_nexus_service_handlers(
                TaskQueue.NOMAD_INTERNAL_WORKFLOWS
            ),
            'activity_executor': executor,
            'shared_state_manager': SharedStateManager.create_from_multiprocessing(
                multiprocessing.Manager()
            ),
            'graceful_shutdown_timeout': timedelta(
                seconds=config.temporal.graceful_shutdown_timeout
            ),
        }

        if worker_config.max_concurrent_activities:
            worker_kwargs['max_concurrent_activities'] = (
                worker_config.max_concurrent_activities
            )
        else:
            minimum_activity_slots = (
                worker_config.min_activity_slots
                if worker_config.min_activity_slots is not None
                else worker_config.pool_size
            )
            worker_kwargs['tuner'] = WorkerTuner.create_resource_based(
                target_memory_usage=worker_config.target_memory_usage,
                target_cpu_usage=worker_config.target_cpu_usage,
                activity_config=ResourceBasedSlotConfig(
                    minimum_slots=minimum_activity_slots,
                    maximum_slots=worker_config.max_activity_slots,
                    ramp_throttle=timedelta(
                        milliseconds=worker_config.activity_ramp_throttle
                    ),
                ),
            )

        worker = Worker(**worker_kwargs)
        health_runner = None

        if should_start_health_server(worker_config):
            health_runner = await start_health_server(
                host=worker_config.healthcheck_host,
                port=worker_config.healthcheck_port,
            )

        # Run the worker until SIGTERM
        logger.info('Starting internal processing worker.')
        worker_task = asyncio.create_task(worker.run())

        try:
            await stop_event.wait()
        finally:
            if health_runner is not None:
                await health_runner.cleanup()

        logger.info('Stopping worker.')
        worker_task.cancel()
        try:
            await worker_task
        except asyncio.CancelledError:
            logger.info('Worker shut down cleanly.')
        finally:
            await close_client()
