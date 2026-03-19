from datetime import datetime, timezone
from typing import Annotated, Any

from beanie import Document, Indexed
from pydantic import Field


class ActionDocument(Document):
    """
    A MongoDB document for storing information about actions.

    Attributes:
        action_id: The ID of the action.
        action_instance_id: The unique ID of the action instance.
        input_data: The input data for the action.
        user_id: The ID of the user who initiated the action.
        upload_id: The ID of the upload associated with the action, if any.
        status: The status of the action.
        results: The results of the action.
        user_input_requests: List of pending user input requests.
        created_at: The timestamp when the action was created.
        updated_at: The timestamp when the action was last updated.
    """

    action_id: Annotated[str, Indexed()]
    action_instance_id: Annotated[str, Indexed(unique=True)]
    input_data: Any
    user_id: Annotated[str, Indexed()]
    upload_id: Annotated[str | None, Indexed()] = None
    status: str
    results: Any = None
    user_input_requests: list[dict] = []
    submitted_user_inputs: list[dict] = []
    created_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))
    updated_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))

    class Settings:
        name = 'action_document'

    async def save(self, *args, **kwargs):
        self.updated_at = datetime.now(timezone.utc)
        return await super().save(*args, **kwargs)
