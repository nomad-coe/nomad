import asyncio
import logging
import signal
import sys
from concurrent.futures import ThreadPoolExecutor
from datetime import timedelta
from typing import Any

from temporalio.worker import ResourceBasedSlotConfig, Worker, WorkerTuner

from nomad.actions import TaskQueue
from nomad.actions.action_logging import (
    WorkflowLoggingInterceptor,
    WorkflowRoutingHandler,
)
from nomad.actions.activities.utils import get_all_activities
from nomad.actions.client import get_client
from nomad.actions.workflows.utils import get_all_workflows
from nomad.config import config
from nomad.config.models.config import WorkerConfig
from nomad.infrastructure import setup
from nomad.utils.structlogging import get_logger


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

    # Ensure the global root logger routes workflow logs
    root_logger = logging.getLogger()
    if not any(isinstance(h, WorkflowRoutingHandler) for h in root_logger.handlers):
        root_logger.addHandler(WorkflowRoutingHandler())

    client = await get_client()
    with ThreadPoolExecutor(max_workers=worker_config.pool_size) as executor:
        worker_kwargs: dict[str, Any] = {
            'client': client,
            'task_queue': TaskQueue.CPU.value,
            'workflows': get_all_workflows(TaskQueue.CPU),
            'activities': get_all_activities(TaskQueue.CPU),
            'activity_executor': executor,
            'interceptors': [WorkflowLoggingInterceptor()],
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
        setup()
        # Run the worker until SIGTERM
        logger.info('Starting CPU worker.')
        worker_task = asyncio.create_task(worker.run())
        await stop_event.wait()

        logger.info('Stopping worker.')
        worker_task.cancel()
        try:
            await worker_task
        except asyncio.CancelledError:
            logger.info('Worker shut down cleanly.')
