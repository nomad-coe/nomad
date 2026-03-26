"""
This module provides utility functions for working with NOMAD actions.

It includes functions for:
- Validating action arguments.
- Retrieving action schemas.
- Managing action execution and results.
- Interacting with the Temporal workflow engine.
"""

import asyncio
import base64
import os
import threading
import uuid
from collections.abc import Coroutine
from dataclasses import asdict, is_dataclass
from datetime import datetime, timedelta, timezone
from typing import Any, get_args, get_origin, get_type_hints

from pydantic import BaseModel, Field, SecretBytes, SecretStr, TypeAdapter
from pymongo import ReturnDocument
from temporalio import activity, workflow
from temporalio.client import WorkflowExecutionStatus
from temporalio.common import RetryPolicy
from temporalio.service import RPCError, RPCStatusCode

from nomad import infrastructure
from nomad.actions.action import get_actions
from nomad.actions.client import get_client
from nomad.config import config
from nomad.files import StagingUploadFiles
from nomad.metainfo.metainfo import Callable
from nomad.mongo.action import ActionDocument
from nomad.processing.data import Upload
from nomad.utils.structlogging import get_logger


class ActionModel(BaseModel):
    action_id: str
    action_instance_id: str
    upload_id: str | None = None
    status: str
    input_data: dict[str, Any] = Field(default_factory=dict)
    signal_input_requests: list[dict[str, str | None]] = Field(default_factory=list)
    signal_inputs_submitted: list[dict[str, Any]] = Field(default_factory=list)
    results: Any = Field(default_factory=dict)
    created_at: datetime
    updated_at: datetime


class ActionModelSummary(BaseModel):
    action_id: str
    action_instance_id: str
    upload_id: str | None = None
    signal_input_requests: list[dict[str, str | None]] = Field(default_factory=list)
    status: str
    created_at: datetime
    updated_at: datetime


class ActionPage(BaseModel):
    """Paginated response for action list queries."""

    items: list[ActionModelSummary]
    next_cursor: str | None = None  # None means no further pages
    total: int  # total number of actions for the user (across all pages)


class ActionSchemaInfo(BaseModel):
    action_id: str
    json_schema: dict[str, Any]
    name: str | None = None
    plugin_package: str | None = None
    description: str | None = None
    task_queue: str | None = None
    signals: list[dict[str, Any]] | None = None


class RunThread(threading.Thread):
    def __init__(self, coro: Coroutine[Any, Any, Any]):
        self.coro = coro
        self.result = None
        self.error: BaseException | None = None
        super().__init__()

    def run(self):
        try:
            self.result = asyncio.run(self.coro)
        except BaseException as exc:
            self.error = exc


def run_async(coro: Coroutine[Any, Any, Any]) -> Any:
    async def _run_with_action_infra():
        await infrastructure.init_async_mongo()
        return await coro

    try:
        loop = asyncio.get_running_loop()
    except RuntimeError:
        loop = None
    if loop and loop.is_running():
        # In jupyter/fastapi there is already a loop running.
        thread = RunThread(_run_with_action_infra())
        thread.start()
        thread.join()
        if thread.error is not None:
            raise thread.error
        return thread.result
    else:
        # Create our own loop.
        return asyncio.run(_run_with_action_infra())


def _run_async_or_return(coro: Coroutine[Any, Any, Any]) -> Any:
    """
    Run a coroutine and return its concrete result.

    This powers synchronous compatibility facades and never returns a coroutine.
    Async call sites must use explicit ``*_async`` functions.
    """
    return run_async(coro)


def _to_dict(data: Any) -> dict:
    if isinstance(data, BaseModel):  # pydantic
        secret_types = (SecretStr, SecretBytes)
        secret_fields = {
            field_name
            for field_name, field_info in type(data).model_fields.items()
            if field_info.annotation in secret_types
            or (
                get_origin(field_info.annotation)
                and any(arg in secret_types for arg in get_args(field_info.annotation))
            )
        }
        return data.model_dump(exclude=secret_fields)
    elif is_dataclass(data) and not isinstance(data, type):
        return asdict(data)
    elif isinstance(data, dict):  # already a dict
        return data
    else:
        raise TypeError(f'Unsupported type: {type(data)}')


def _validate_with_pydantic(func: Callable, arg):
    """
    Validate the single argument of a function against its type hint using Pydantic.

    Args:
        func: The function with the argument to validate.
        arg: The argument to validate.

    Returns:
        The validated argument.
    """
    hints = get_type_hints(func)

    # get the single non-return annotation
    [(_, param_type)] = [(n, t) for n, t in hints.items() if n != 'return']

    adapter = TypeAdapter(param_type)
    return adapter.validate_python(arg)


def _get_param_schema(func: Callable) -> dict[str, Any]:
    """
    Generate a JSON Schema for the single argument of a function.

    This is useful for generating frontend forms for actions.

    Args:
        func: The function with the argument to generate the schema for.

    Returns:
        The JSON schema for the argument.
    """
    hints = get_type_hints(func)

    # get the single non-return annotation
    [(_, param_type)] = [(n, t) for n, t in hints.items() if n != 'return']

    if isinstance(param_type, type) and issubclass(param_type, BaseModel):
        schema = param_type.model_json_schema()
    else:
        adapter = TypeAdapter(param_type)
        schema = adapter.json_schema()

    # remove the user_id from the schema,
    # we rely on the user_id of the logged in user instead of form input.
    schema.get('properties', {}).pop('user_id', None)
    required = schema.get('required', [])
    if 'user_id' in required:
        required.remove('user_id')
    return schema


def _get_signal_schema(signal_fn: Callable) -> dict[str, Any]:
    """
    Generate a JSON Schema for the single argument of a signal function.
    Raises ValueError if there is more than one argument (excluding 'self').
    """
    import inspect

    sig = inspect.signature(signal_fn)  # type: ignore
    params = list(sig.parameters.values())

    # Exclude 'self'
    if params and params[0].name == 'self':
        params = params[1:]

    if len(params) > 1:
        name_str = getattr(signal_fn, '__name__', str(signal_fn))
        raise ValueError(
            f'Signal {name_str} has more than one argument. Only zero or one arguments are supported.'
        )

    if len(params) == 0:
        return {}

    hints = get_type_hints(signal_fn)
    param_type = hints.get(params[0].name, Any)

    if isinstance(param_type, type) and issubclass(param_type, BaseModel):
        return param_type.model_json_schema()
    else:
        adapter = TypeAdapter(param_type)
        return adapter.json_schema()


def validate_action_arg(action_id: str, arg: Any):
    """
    Validate the argument for a given action's `workflow.run` function
    against its type hint. Raises if the action does not exist or the
    argument is invalid.
    """
    action = get_actions().get(action_id)
    if not action:
        raise ValueError('Action not found')
    return _validate_with_pydantic(action.load().workflow.run, arg)


def get_all_action_schemas() -> list[ActionSchemaInfo]:
    """
    Return a list of JSON Schemas for all registered actions'
    `workflow.run` parameters, keyed by action_id.
    """
    data: list[ActionSchemaInfo] = []
    for action_id, action in get_actions().items():
        workflow_cls = action.load().workflow

        signals = []
        for attr_name in dir(workflow_cls):
            if attr_name.startswith('__'):
                continue
            attr = getattr(workflow_cls, attr_name, None)
            if hasattr(attr, '__temporal_signal_definition'):
                signal_fn: Callable | None = getattr(
                    getattr(attr, '__temporal_signal_definition'),
                    'fn',
                    attr,
                )
                if signal_fn is not None:
                    schema = _get_signal_schema(signal_fn)
                    # Expose Python method names as the canonical API key.
                    signals.append({attr_name: schema})

        data.append(
            ActionSchemaInfo(
                action_id=action_id,
                json_schema=_get_param_schema(workflow_cls.run),
                description=action.description,
                task_queue=action.task_queue,
                name=action.name,
                plugin_package=action.plugin_package,
                signals=signals,
            )
        )
    return data


async def _get_workflow_status_safe(
    action_instance_id: str,
) -> WorkflowExecutionStatus | None:
    """
    Safely retrieves workflow status, returning None if workflow not found.

    Args:
        action_instance_id: The unique ID of the action instance.

    Returns:
        The workflow status, or None if workflow not found.

    Raises:
        Exception: For errors other than workflow not found.
    """
    try:
        client = await get_client()
        handle = client.get_workflow_handle(action_instance_id)
        desc = await handle.describe()
        return desc.status
    except RPCError as e:
        if e.status == RPCStatusCode.NOT_FOUND:
            return None
        raise


async def _get_workflow_result_safe(action_instance_id: str) -> dict[str, Any] | None:
    """
    Safely retrieves workflow result, returning None if workflow not found.

    Args:
        action_instance_id: The unique ID of the action instance.

    Returns:
        The workflow result, or None if workflow not found.

    Raises:
        Exception: For errors other than workflow not found.
    """
    try:
        client = await get_client()
        handle = client.get_workflow_handle(action_instance_id)
        return await handle.result()
    except RPCError as e:
        if e.status == RPCStatusCode.NOT_FOUND:
            return None
        raise


async def _validate_action_ownership(
    action_instance_id: str, user_id: str
) -> ActionDocument:
    """
    Validates that an action exists and belongs to the specified user.

    Args:
        action_instance_id: The unique ID of the action instance.
        user_id: The user ID to validate ownership against.

    Returns:
        The action document if found and owned by user.

    Raises:
        Exception: If action not found or owned by different user.
    """
    action = await ActionDocument.find_one(
        ActionDocument.action_instance_id == action_instance_id,
        ActionDocument.user_id == user_id,
    )
    if not action:
        raise Exception(
            'The action was not registered in the DB or was registered under a different user.'
        )
    return action


async def _get_action_status_async(
    action_instance_id: str, user_id: str
) -> WorkflowExecutionStatus:
    """
    Retrieves the current execution status of an action.

    Args:
        action_instance_id: The unique ID of the action instance to check.
        user_id: The user who initiated the action.

    Returns:
        The current status of the action. If workflow is not found, returns
        TERMINATED.
    """
    action = await _validate_action_ownership(action_instance_id, user_id)
    logger = get_logger(__name__)

    status = await _get_workflow_status_safe(action_instance_id)

    if status is None:
        logger.warning(
            f'Workflow {action_instance_id} could not be found for user {user_id}. '
            f'Setting status to TERMINATED.'
        )
        action.status = str(WorkflowExecutionStatus.TERMINATED.name)
        await action.save()
        return WorkflowExecutionStatus.TERMINATED

    action.status = str(status.name)
    await action.save()
    return status


def get_action_status(action_instance_id: str, user_id: str) -> WorkflowExecutionStatus:
    """
    Retrieves the current execution status of an action.

    Synchronous callers can call this function directly.
    Asynchronous callers should use ``await get_action_status_async(...)``.
    """
    return _run_async_or_return(_get_action_status_async(action_instance_id, user_id))


async def get_action_status_async(
    action_instance_id: str, user_id: str
) -> WorkflowExecutionStatus:
    """
    Async-only variant of ``get_action_status`` for typed async call sites.
    """
    return await _get_action_status_async(action_instance_id, user_id)


def get_action_result(action_instance_id: str, user_id: str) -> dict[str, Any] | None:
    """
    Retrieves the result of a completed action.

    Synchronous callers can call this function directly.
    Asynchronous callers should use ``await get_action_result_async(...)``.
    """
    return _run_async_or_return(get_action_result_async(action_instance_id, user_id))


async def get_action_result_async(
    action_instance_id: str, user_id: str
) -> dict[str, Any] | None:
    """
    Retrieves the result of a completed action.
    """
    logger = get_logger(__name__)
    action = await _validate_action_ownership(action_instance_id, user_id)

    results = await _get_workflow_result_safe(action_instance_id)

    if results is None:
        logger.warning(
            f'Workflow {action_instance_id} could not be found for user {user_id}. '
            f'Result could not be retrieved.'
        )
        return None

    action.results = _to_dict(results)
    action.status = str(WorkflowExecutionStatus.COMPLETED.name)
    await action.save()
    return results


async def _update_status(action: ActionDocument):
    """
    Update the status of an action in the database.
    Silently handles workflow not found errors by setting status to UNKNOWN.

    Args:
        action: The action document to update.
    """
    status = await _get_workflow_status_safe(action.action_instance_id)
    logger = get_logger(__name__)

    if status is None:
        # Workflow not found - mark as unknown and return
        logger.warning(
            f'Workflow {action.action_instance_id} could not be found. '
            f'Setting status to UNKNOWN.'
        )
        action.status = 'UNKNOWN'
        await action.save()
        return

    action.status = str(status.name)

    if status.name == 'COMPLETED':
        results = await _get_workflow_result_safe(action.action_instance_id)
        if results:
            try:
                action.results = _to_dict(results)
            except TypeError:
                action.results = results

    await action.save()


_CURSOR_DT_FMT = '%Y-%m-%dT%H:%M:%S.%f+00:00'


def _encode_cursor(dt: datetime) -> str:
    """
    Encode a datetime as an opaque, base64url cursor string.

    The cursor encodes the ``created_at`` timestamp of the *last item on the
    current page*.  The next query will return documents whose ``created_at``
    is strictly less than this value, giving stable forward-only pagination
    even as new documents are inserted at the head of the collection.
    """
    # Beanie/Mongo can return naive datetimes when tz-awareness is disabled;
    # treat those values as UTC to avoid timezone-shifted cursors.
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)

    # Always work in UTC so the encoded string is unambiguous.
    utc_dt = dt.astimezone(timezone.utc)
    token = utc_dt.strftime(_CURSOR_DT_FMT)
    return base64.urlsafe_b64encode(token.encode()).decode()


def _decode_cursor(cursor: str) -> datetime:
    """
    Decode a cursor string produced by :func:`_encode_cursor`.

    Raises ``ValueError`` when the token is not a valid base64url string or
    does not decode to the expected timestamp format.
    """
    try:
        token = base64.urlsafe_b64decode(cursor.encode()).decode()
        return datetime.strptime(token, _CURSOR_DT_FMT).replace(tzinfo=timezone.utc)
    except Exception as exc:
        raise ValueError(f'Invalid pagination cursor: {cursor!r}') from exc


async def list_user_actions(
    user_id: str,
    page_size: int = 20,
    cursor: str | None = None,
    upload_id: str | None = None,
) -> 'ActionPage':
    """
    Get a page of actions for a given user, ordered by ``created_at`` descending
    (newest first).

    This function also updates the status of any pending or running actions
    within the returned page.

    Args:
        user_id: The ID of the user.
        page_size: Maximum number of items to return (default 20).
        cursor: Opaque pagination token returned by a previous call.  When
            supplied the query returns the next page of results after the
            cursor position.
        upload_id: Optional upload ID to filter actions by.

    Returns:
        An :class:`ActionPage` containing the items, an optional
        ``next_cursor`` for the following page, and the ``total`` count of
        all documents belonging to this user.
    """
    query_filters = [ActionDocument.user_id == user_id]
    if upload_id is not None:
        query_filters.append(ActionDocument.upload_id == upload_id)

    # Decode cursor and narrow the query to documents strictly older than it.
    cursor_filter = None
    if cursor is not None:
        cursor_dt = _decode_cursor(cursor)
        cursor_filter = ActionDocument.created_at < cursor_dt

    # Fetch one extra document to detect whether a next page exists.
    fetch_limit = page_size + 1
    if cursor_filter is not None:
        action_query = ActionDocument.find(*query_filters, cursor_filter)
    else:
        action_query = ActionDocument.find(*query_filters)
    action_documents = (
        await action_query.sort('-created_at').limit(fetch_limit).to_list()
    )

    has_next = len(action_documents) == fetch_limit
    page_docs = action_documents[:page_size]

    # Update status only for PENDING/RUNNING items in this page.
    active_actions = [a for a in page_docs if a.status in ('PENDING', 'RUNNING')]
    if active_actions:
        await asyncio.gather(*(_update_status(a) for a in active_actions))

    next_cursor: str | None = None
    if has_next and page_docs:
        next_cursor = _encode_cursor(page_docs[-1].created_at)

    # Cheap total count (uses the (user_id, created_at) compound index).
    total = await ActionDocument.find(*query_filters).count()

    return ActionPage(
        items=[
            ActionModelSummary.model_construct(**doc.model_dump()) for doc in page_docs
        ],
        next_cursor=next_cursor,
        total=total,
    )


async def get_user_action(action_instance_id: str, user_id: str) -> ActionModel | None:
    """
    Get a specific action for a given user.

    This function also updates the status of the action if it's pending or running.

    Args:
        action_instance_id: The ID of the action instance.
        user_id: The ID of the user.

    Returns:
        The action if found, otherwise None.
    """
    action_document = await ActionDocument.find_one(
        ActionDocument.action_instance_id == action_instance_id,
        ActionDocument.user_id == user_id,
    )

    if not action_document:
        return None

    if action_document.status in ('PENDING', 'RUNNING'):
        await _update_status(action_document)
    elif action_document.status == 'COMPLETED' and not action_document.results:
        # Backfill results for completed rows where results were not persisted yet.
        results = await _get_workflow_result_safe(action_instance_id)
        if results is not None:
            try:
                serialized_results = _to_dict(results)
            except TypeError:
                serialized_results = results
            action_document.results = serialized_results
            await ActionDocument.get_pymongo_collection().update_one(
                {'_id': action_document.id},
                {
                    '$set': {
                        'results': serialized_results,
                        'updated_at': datetime.now(),
                    }
                },
            )

    return ActionModel.model_construct(**action_document.model_dump())


def get_upload_files(upload_id: str, user_id: str) -> StagingUploadFiles | None:
    """
    Retrieves files for an upload after verifying user authorization.

    Checks if the user is the main author or a coauthor.

    Args:
        upload_id: The unique identifier for the upload.
        user_id: The unique identifier for the user.

    Returns:
        The UploadFiles object if found and authorized, otherwise None
        (if the upload doesn't exist or the associated files aren't found).

    Raises:
        PermissionError: If the upload exists but the user is not authorized.
    """
    if infrastructure.mongo_client is None:
        infrastructure.setup_mongo()

    upload = Upload.get(upload_id)

    if upload is None:
        return None

    # Determine if user is authorized to get the upload.
    is_coauthor = isinstance(upload.coauthors, list) and user_id in upload.coauthors
    is_authorized = upload.main_author == user_id or is_coauthor

    # Raise error if not authorized
    if not is_authorized:
        raise PermissionError(
            f'User {user_id} is not authorized to access upload {upload_id}.'
        )

    # User is authorized, retrieve and return files
    if StagingUploadFiles.exists_for(upload_id):
        return StagingUploadFiles(upload_id)

    return None


def action_artifacts_dir() -> str:
    """
    Returns the path to the action artifacts directory.

    Activities can use this directory to store artifacts that can be used
    by multiple actions, such as ML training models, global configuration,
    or reference datasets.
    """

    path = os.path.join(config.fs.actions, 'artifacts')
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)
    return path


def action_instance_artifacts_dir(action_instance_id: str) -> str:
    """
    Returns the path to the artifacts directory for a specific instance.

    Activities can use this directory to store artifacts that are generated
    by a given instance, for example a classification_result for a given input.
    """
    action_instance_dir = os.path.join(config.fs.actions, action_instance_id)
    if not os.path.exists(action_instance_dir):
        os.makedirs(action_instance_dir, exist_ok=True)
    return action_instance_dir


def action_log_file_path(action_instance_id: str) -> str:
    """
    Returns the file path for the logs of a given action instance.
    Logs are stored in config.fs.action/logs/<action_instance_id>.log.
    """
    log_dir = os.path.join(config.fs.actions, 'logs')
    if not os.path.exists(log_dir):
        os.makedirs(log_dir, exist_ok=True)
    return os.path.join(log_dir, f'{action_instance_id}.log')


async def _async_start_workflow(action, data, workflow_id) -> str:
    """
    Asynchronously starts a workflow.

    Args:
        action: The action to start.
        data: The input data for the workflow.
        workflow_id: The ID of the workflow to start.

    Returns:
        The ID of the started workflow.
    """
    client = await get_client()
    await client.start_workflow(
        action.workflow.run,
        data,
        id=workflow_id,
        task_queue=action.task_queue,
    )
    return workflow_id


async def _async_stop_workflow(workflow_id: str):
    """
    Asynchronously stops a workflow.

    Args:
        workflow_id: The ID of the workflow to stop.
    """
    client = await get_client()
    handle = client.get_workflow_handle(workflow_id)
    await handle.cancel()


async def _async_signal_workflow(
    workflow_cls,
    workflow_id: str,
    signal_fn_name: str,
    data: Any,
) -> None:
    """
    Asynchronously send a signal to a running Temporal workflow execution.

    This helper obtains a workflow handle using the provided workflow ID,
    resolves the specified signal method from the workflow class, and sends
    the signal with the supplied payload.

    Args:
        workflow_cls: The Temporal workflow class that defines the signal.
        workflow_id: The ID of the target workflow execution. This targets the
            latest run for the given ID unless a run ID is specified elsewhere.
        signal_fn_name: The name of the signal method on the workflow class.
            The method must be decorated with ``@workflow.signal``.
        data: The payload to send with the signal. Must be serializable by the
            Temporal payload converter configured for the client.
    """
    client = await get_client()
    handle = client.get_workflow_handle(workflow_id)

    try:
        signal_fn: Callable = getattr(workflow_cls, signal_fn_name)
    except AttributeError as e:
        raise ValueError(
            f"Signal '{signal_fn_name}' not found on workflow {workflow_cls.__name__}"
        ) from e

    if not callable(signal_fn):
        raise TypeError(
            f"Attribute '{signal_fn_name}' on {workflow_cls.__name__} is not callable"
        )

    await handle.signal(signal_fn, data)


async def _start_action_async(action_id: str, data: Any) -> str:
    """
    Starts a new Action with the given ID and input data.

    Args:
        action_id: The ID of the action to start.
        data: Input data for the action.

    Returns:
        The unique ID of the started action instance.
    """
    assert hasattr(data, 'user_id')
    user_id = data.user_id
    workflow_id = f'{action_id}-{user_id}-{uuid.uuid4()}'
    action_entry_point = get_actions().get(action_id)
    assert action_entry_point, f'No action data for the given {action_id} ID'
    action = action_entry_point.load()

    upload_id = getattr(data, 'upload_id', None)
    new_action = ActionDocument(
        action_id=action_id,
        action_instance_id=workflow_id,
        user_id=user_id,
        upload_id=upload_id,
        status='PENDING',
        input_data=_to_dict(data),
    )
    await new_action.insert()

    await _async_start_workflow(action, data, workflow_id)
    return workflow_id


def start_action(action_id: str, data: Any) -> str:
    """
    Starts a new Action with the given ID and input data.

    Synchronous callers can call this function directly.
    Asynchronous callers should use ``await start_action_async(...)``.
    """
    return _run_async_or_return(_start_action_async(action_id, data))


async def start_action_async(action_id: str, data: Any) -> str:
    """
    Async-only variant of ``start_action`` for typed async call sites.
    """
    return await _start_action_async(action_id, data)


async def _stop_action_async(action_instance_id: str, user_id: str):
    """
    Stops a running action.

    Args:
        action_instance_id: The unique ID of the action instance to stop.
        user_id: The user who initiated the action.
    """
    action = await ActionDocument.find_one(
        ActionDocument.action_instance_id == action_instance_id,
        ActionDocument.user_id == user_id,
    )
    if not action:
        raise Exception(
            'The action was not registered in the DB or was registered under a different user.'
        )

    if action.status not in ('PENDING', 'RUNNING'):
        raise Exception('Action is not running.')

    await _async_stop_workflow(action_instance_id)

    action.status = str(WorkflowExecutionStatus.CANCELED.name)
    await action.save()


def stop_action(action_instance_id: str, user_id: str):
    """
    Stops a running action.

    Synchronous callers can call this function directly.
    Asynchronous callers should use ``await stop_action_async(...)``.
    """
    return _run_async_or_return(_stop_action_async(action_instance_id, user_id))


async def stop_action_async(action_instance_id: str, user_id: str):
    """
    Async-only variant of ``stop_action`` for typed async call sites.
    """
    return await _stop_action_async(action_instance_id, user_id)


class RequestSignalInputActivityInput(BaseModel):
    """Input parameters for the activity that delegates to the manager's request_signal_input."""

    action_instance_id: str = Field(
        ..., description='The ID of the action instance/workflow.'
    )
    user_id: str = Field(..., description='The ID of the user who owns the action.')
    signal_fn_name: str = Field(
        ..., description='The Temporal signal name to add to pending requests.'
    )
    title: str | None = Field(
        default=None, description='Optional title to use for the signal input form.'
    )
    description: str | None = Field(
        default=None, description='Optional description for the signal input form.'
    )
    content: str | None = Field(
        default=None, description='Optional markdown content for the signal input form.'
    )
    initial_data: dict[str, Any] | None = Field(
        default=None, description='Optional initial data for the signal input form.'
    )


@activity.defn
async def request_signal_input_activity(data: RequestSignalInputActivityInput):
    """
    Activity that interacts with the main nomad API state to record that
    this workflow is waiting on a signal input signal.
    """
    request_info: dict[str, Any] = {
        'signal_fn_name': data.signal_fn_name,
    }
    if data.title is not None:
        request_info['title'] = data.title
    if data.description is not None:
        request_info['description'] = data.description
    if data.content is not None:
        request_info['content'] = data.content
    if data.initial_data is not None:
        request_info['initial_data'] = data.initial_data

    collection = ActionDocument.get_pymongo_collection()
    result = await collection.update_one(
        {
            'action_instance_id': data.action_instance_id,
            'user_id': data.user_id,
            'status': {'$in': ['PENDING', 'RUNNING']},
            'signal_input_requests.signal_fn_name': {'$ne': data.signal_fn_name},
        },
        {
            '$push': {'signal_input_requests': request_info},
            '$set': {'updated_at': datetime.now(timezone.utc)},
        },
    )
    if result.modified_count == 1:
        return {'status': 'signal_input_requested'}

    action = await ActionDocument.find_one(
        ActionDocument.action_instance_id == data.action_instance_id,
        ActionDocument.user_id == data.user_id,
    )
    if not action:
        raise Exception(
            'The action was not registered in the DB or was registered under a different user.'
        )
    if action.status not in ('PENDING', 'RUNNING'):
        raise Exception('Action is not running.')
    raise Exception(
        f"Action already has a pending signal input request for signal '{data.signal_fn_name}'."
    )


async def request_signal_input(
    action_instance_id: str,
    user_id: str,
    signal_fn_name: str,
    title: str | None = None,
    description: str | None = None,
    content: str | None = None,
    timeout: timedelta = timedelta(hours=1),
):
    """
    Record that a running action is requesting signal input for a specific signal.

    This is a helper function for workflow authors. It calls an activity that
    looks up the action instance for the given user, verifies it is active,
    and appends the requested signal dictionary to the action's pending
    signal-input requests list in the database.

    Args:
        action_instance_id: Unique identifier of the action instance (also used
            as the workflow ID).
        user_id: Identifier of the user who owns the action instance.
        signal_fn_name: Name of the workflow signal representing the requested
            signal input.
        title: Optional title to display in the frontend form.
        description: Optional description to display in the frontend form.
        content: Optional markdown content to display in the frontend form.
        timeout: Optional timeout for the signal input activity. Defaults to 1 hour.

    Raises:
        Exception: If the action does not exist for the user or is not active.
    """
    await workflow.execute_activity(
        request_signal_input_activity,
        RequestSignalInputActivityInput(
            action_instance_id=action_instance_id,
            user_id=user_id,
            signal_fn_name=signal_fn_name,
            title=title,
            description=description,
            content=content,
        ),
        start_to_close_timeout=timeout,
        retry_policy=RetryPolicy(maximum_attempts=3),
    )


async def submit_signal_input(
    action_instance_id: str, user_id: str, signal_fn_name: str, data: Any
):
    """
    Submit signal input to a running action by signaling its workflow.

    Validates that the action exists for the user and is active, resolves the
    associated Temporal workflow class from the action registry, and sends the
    specified signal with the provided payload to the workflow execution.

    Args:
        action_instance_id: Unique identifier of the action instance (also used
            as the workflow ID).
        user_id: Identifier of the user who owns the action instance.
        signal_fn_name: Name of the workflow signal to invoke.
        data: Payload to send with the signal. Must be serializable by the
            Temporal payload converter.

    Raises:
        Exception: If the action does not exist for the user or is not active.
        AssertionError: If no workflow entry point is registered for the action.
        temporalio.exceptions.TemporalError: If signaling the workflow fails.
    """
    collection = ActionDocument.get_pymongo_collection()
    action_data = await collection.find_one_and_update(
        {
            'action_instance_id': action_instance_id,
            'user_id': user_id,
            'status': {'$in': ['PENDING', 'RUNNING']},
            'signal_input_requests.signal_fn_name': signal_fn_name,
        },
        {
            '$pull': {'signal_input_requests': {'signal_fn_name': signal_fn_name}},
            '$set': {'updated_at': datetime.now(timezone.utc)},
        },
        return_document=ReturnDocument.BEFORE,
    )
    if not action_data:
        action = await ActionDocument.find_one(
            ActionDocument.action_instance_id == action_instance_id,
            ActionDocument.user_id == user_id,
        )
        if not action:
            raise Exception(
                'The action was not registered in the DB or was registered under a different user.'
            )
        if action.status not in ('PENDING', 'RUNNING'):
            raise Exception('Action is not running.')
        raise Exception(
            f"No pending signal input request found for signal '{signal_fn_name}'."
        )

    pending_requests = action_data.get('signal_input_requests') or []
    matching_request = next(
        (
            req
            for req in pending_requests
            if req.get('signal_fn_name') == signal_fn_name
        ),
        None,
    )
    if not matching_request:
        raise Exception(
            f"No pending signal input request found for signal '{signal_fn_name}'."
        )

    action_id = action_data.get('action_id')
    if action_id is None:
        raise Exception(
            'The action was not registered in the DB or was registered under a different user.'
        )
    action_entry_point = get_actions().get(action_id)
    assert action_entry_point, f'No action data for the given {action_id} ID'
    workflow_cls = action_entry_point.load().workflow
    try:
        await _async_signal_workflow(
            workflow_cls,
            workflow_id=action_instance_id,
            signal_fn_name=signal_fn_name,
            data=data,
        )
    except Exception:
        # Best-effort rollback of the pending request if signaling fails.
        await collection.update_one(
            {
                'action_instance_id': action_instance_id,
                'user_id': user_id,
                'signal_input_requests.signal_fn_name': {'$ne': signal_fn_name},
            },
            {
                '$push': {'signal_input_requests': matching_request},
                '$set': {'updated_at': datetime.now(timezone.utc)},
            },
        )
        raise

    # Save the submitted signal data in db.
    serialized_data: Any
    try:
        serialized_data = _to_dict(data)
    except TypeError:
        serialized_data = (
            str(data)
            if not isinstance(data, (int, float, bool, str, list, dict, type(None)))
            else data
        )

    submitted_entry: dict[str, Any] = {
        'signal_fn_name': signal_fn_name,
        'data': serialized_data,
        'timestamp': datetime.now(timezone.utc).isoformat(),
    }
    if matching_request.get('title') is not None:
        submitted_entry['title'] = matching_request.get('title')
    if matching_request.get('description') is not None:
        submitted_entry['description'] = matching_request.get('description')
    if matching_request.get('content') is not None:
        submitted_entry['content'] = matching_request.get('content')

    await collection.update_one(
        {
            'action_instance_id': action_instance_id,
            'user_id': user_id,
        },
        {
            '$push': {'signal_inputs_submitted': submitted_entry},
            '$set': {'updated_at': datetime.now(timezone.utc)},
        },
    )
