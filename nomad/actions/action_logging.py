import logging
from collections.abc import Awaitable, Callable
from logging.handlers import WatchedFileHandler
from typing import Any

import structlog
from temporalio import activity
from temporalio.worker import (
    ActivityInboundInterceptor,
    ExecuteActivityInput,
    Interceptor,
)

from nomad.actions.manager import action_log_file_path
from nomad.config import config
from nomad.utils.structlogging import ISO8601_UTC_FORMAT

_resolve_root_workflow_id: Callable[[str], Awaitable[str]]
_resolved_root_workflow_ids: dict[str, str] = {}


async def _default_resolve_root_workflow_id(workflow_id: str) -> str:
    """
    Resolve a workflow ID to its root workflow ID.

    Child workflows should write logs to the same action log file as their root action
    workflow. If lookup fails, we gracefully fall back to the current workflow ID.
    """

    if workflow_id in _resolved_root_workflow_ids:
        return _resolved_root_workflow_ids[workflow_id]

    from nomad.actions.client import get_client

    try:
        client = await get_client()
        handle = client.get_workflow_handle(workflow_id)
        description = await handle.describe()
        root_execution = (
            description.raw_description.workflow_execution_info.root_execution
        )
        root_workflow_id = root_execution.workflow_id or workflow_id
    except Exception:
        root_workflow_id = workflow_id

    _resolved_root_workflow_ids[workflow_id] = root_workflow_id
    _resolved_root_workflow_ids[root_workflow_id] = root_workflow_id
    return root_workflow_id


_resolve_root_workflow_id = _default_resolve_root_workflow_id


class WorkflowRoutingHandler(logging.Handler):
    """
    A custom logging handler that routes log records to a specific file based on the
    `workflow_id` present in `structlog.contextvars`.
    """

    def __init__(self):
        super().__init__()
        self._handlers: dict[str, WatchedFileHandler] = {}

    def emit(self, record):
        # We need to extract the workflow_id from structlog's contextvars
        context = structlog.contextvars.get_contextvars()
        workflow_id = context.get('action_instance_id') or context.get('workflow_id')

        if not workflow_id:
            # If no workflow_id is found, ignore this log for workflow routing
            return

        try:
            handler = self._get_handler(workflow_id)
            if handler and record.levelno >= handler.level:
                handler.handle(record)
        except Exception:
            self.handleError(record)

    def _get_handler(self, workflow_id: str) -> WatchedFileHandler:
        if workflow_id not in self._handlers:
            log_file = action_log_file_path(workflow_id)
            handler = WatchedFileHandler(log_file)
            handler.setLevel(config.services.actions_log_level)

            # Simple formatter to match standard Python log structures locally
            formatter = logging.Formatter(
                '{"timestamp": "%(asctime)s", "level": "%(levelname)s", "event": "%(message)s"}'
            )
            formatter.datefmt = ISO8601_UTC_FORMAT
            handler.setFormatter(formatter)
            self._handlers[workflow_id] = handler
        return self._handlers[workflow_id]

    def close(self):
        for handler in self._handlers.values():
            handler.close()
        super().close()


class _WorkflowLoggingInboundInterceptor(ActivityInboundInterceptor):
    async def execute_activity(self, input: ExecuteActivityInput) -> Any:
        info = activity.info()
        if info.workflow_id:
            action_instance_id = await _resolve_root_workflow_id(info.workflow_id)
            structlog.contextvars.bind_contextvars(
                workflow_id=info.workflow_id, action_instance_id=action_instance_id
            )
        return await super().execute_activity(input)


class WorkflowLoggingInterceptor(Interceptor):
    """
    A Temporal Interceptor that binds the `workflow_id` into `structlog.contextvars`
    for activities executed by the worker.
    """

    def intercept_activity(
        self, next: ActivityInboundInterceptor
    ) -> ActivityInboundInterceptor:
        return _WorkflowLoggingInboundInterceptor(next)
