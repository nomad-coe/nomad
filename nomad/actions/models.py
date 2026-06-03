from typing import Any, Literal

from pydantic import BaseModel, Field

from nomad.models.common import UTCDateTime

__all__ = [
    'ActionRecord',
    'ActionSummaryRecord',
    'ActionRecordPage',
    'ActionSchemaInfo',
    'RequestSignalInputActivityInput',
]


class ActionRecord(BaseModel):
    action_id: str
    action_instance_id: str
    user_id: str
    upload_id: str | None = None
    status: str
    input_data: dict[str, Any] = Field(default_factory=dict)
    signal_input_requests: list[dict[str, Any]] = Field(default_factory=list)
    signal_inputs_submitted: list[dict[str, Any]] = Field(default_factory=list)
    results: Any = None
    created_at: UTCDateTime
    updated_at: UTCDateTime
    priority_key: int | None = None
    priority_fairness_key: Literal['user_id'] | None = None


class ActionSummaryRecord(BaseModel):
    action_id: str
    action_instance_id: str
    upload_id: str | None = None
    signal_input_requests: list[dict[str, str | None]] = Field(default_factory=list)
    status: str
    created_at: UTCDateTime
    updated_at: UTCDateTime


class ActionRecordPage(BaseModel):
    """Paginated response for action list queries."""

    items: list[ActionSummaryRecord]
    next_cursor: str | None = None
    total: int


class ActionSchemaInfo(BaseModel):
    action_id: str
    json_schema: dict[str, Any]
    name: str | None = None
    plugin_package: str | None = None
    description: str | None = None
    task_queue: str | None = None
    groups: list[str] | None = None
    users: list[str] | None = None
    signals: list[dict[str, Any]] | None = None


class RequestSignalInputActivityInput(BaseModel):
    """Input parameters for the activity that delegates to request_signal_input."""

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
