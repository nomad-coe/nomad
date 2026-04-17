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

from mongoengine import BooleanField, DateTimeField, Document, ListField, StringField

from nomad.common import now
from nomad.config import config


class PAT(Document):
    """
    A MongoDB document for storing personal access token (PAT).
    """

    # Metadata
    name = StringField(required=True)
    description = StringField()

    # Security
    user_id = StringField(required=True)
    token_digest = StringField(required=True, unique=True, min_length=64, max_length=64)
    scopes = ListField(StringField())

    # Lifecycle
    revoked = BooleanField(default=False)
    revoked_at = DateTimeField()
    expired_at = DateTimeField()
    created_at = DateTimeField(default=now)
    updated_at = DateTimeField(default=now)
    last_used_at = DateTimeField()

    meta = {
        'collection': 'personal_access_tokens',
        'indexes': [
            ('user_id', '-created_at'),
            # Auto-delete expired/revoked tokens after set time
            {
                'fields': ['expired_at'],
                'expireAfterSeconds': config.auth.pat_pruning_time * 86400,
            },
            {
                'fields': ['revoked_at'],
                'expireAfterSeconds': config.auth.pat_pruning_time * 86400,
            },
        ],
    }

    def save(self, *args, **kwargs) -> None:
        self.updated_at = now()
        return super().save(*args, **kwargs)

    @property
    def is_expired(self) -> bool:
        """Checks if the token has passed its expiration date."""
        if self.expired_at is None:
            return False

        current_time = now()

        # If DB timestamp is naive (no timezone), force 'current_time' to be naive too
        if self.expired_at.tzinfo is None:
            current_time = current_time.replace(tzinfo=None)

        return self.expired_at < current_time

    @property
    def is_active(self) -> bool:
        """
        Checks if the token is currently active.
        Returns False if revoked or expired.
        """
        return not self.revoked and not self.is_expired
