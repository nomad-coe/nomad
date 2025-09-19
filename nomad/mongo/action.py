from datetime import datetime, timezone

from mongoengine import DateTimeField, Document, DynamicField, StringField


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
        created_at: The timestamp when the action was created.
        updated_at: The timestamp when the action was last updated.
    """

    action_id = StringField(required=True)
    action_instance_id = StringField(required=True, unique=True)
    input_data = DynamicField(required=True)
    user_id = StringField(required=True)
    upload_id = StringField()
    status = StringField(required=True)
    results = DynamicField()
    created_at = DateTimeField(default=datetime.now(timezone.utc))
    updated_at = DateTimeField(default=datetime.now(timezone.utc))

    meta = {'indexes': ['action_id', 'action_instance_id', 'user_id', 'upload_id']}

    def save(self, *args, **kwargs):
        if not self.created_at:
            self.created_at = datetime.now(timezone.utc)
        self.updated_at = datetime.now(timezone.utc)
        return super().save(*args, **kwargs)
