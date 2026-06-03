from datetime import datetime, timezone
from typing import Annotated, Any

from beanie import Document, Indexed
from pydantic import Field
from pymongo import ASCENDING, DESCENDING, IndexModel


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
        signal_input_requests: List of pending signal input requests.
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
    signal_input_requests: list[dict] = Field(default_factory=list)
    signal_inputs_submitted: list[dict] = Field(default_factory=list)
    priority_key: int | None = None
    priority_fairness_key: str | None = None
    created_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))
    updated_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))

    class Settings:
        name = 'action_document'
        indexes = [
            # Supports cursor-paginated list queries: WHERE user_id = ? ORDER BY created_at DESC
            IndexModel(
                [('user_id', ASCENDING), ('created_at', DESCENDING)],
                name='user_id_created_at_desc',
            ),
        ]

    async def save(self, *args, **kwargs):
        self.updated_at = datetime.now(timezone.utc)
        return await super().save(*args, **kwargs)
