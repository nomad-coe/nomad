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

from datetime import datetime

from mongoengine import DateTimeField, Document, StringField

from nomad.common import now


class OwnershipTransferRecord(Document):
    """Ownership transfer records between owners and target users.

    Multiple records can exist per resource for history. At most one active record
    (``pending``, ``responding``, or ``canceling``) is allowed per resource.
    """

    RESOURCE_TYPE_UPLOAD = 'upload'
    RESOURCE_TYPE_GROUP = 'group'

    STATE_PENDING = 'pending'
    STATE_RESPONDING = 'responding'
    STATE_CANCELING = 'canceling'
    STATE_REFUSED = 'refused'
    ACTIVE_STATES = (STATE_PENDING, STATE_RESPONDING, STATE_CANCELING)

    resource_type = StringField(required=True, default=RESOURCE_TYPE_UPLOAD)
    resource_id = StringField(required=True)
    source_user_id = StringField(
        required=True, help_text='current owner who initiated transfer'
    )
    target_user_id = StringField(required=True, help_text='prospective new owner')
    state = StringField(
        required=True,
        default=STATE_PENDING,
        choices=[STATE_PENDING, STATE_RESPONDING, STATE_CANCELING, STATE_REFUSED],
    )
    requested_at = DateTimeField(default=now, required=True)
    updated_at = DateTimeField(default=now, required=True)
    actor_user_id = StringField(help_text='who resolved the transfer (accept/refuse)')

    meta = {
        'collection': 'ownership_transfer',
        'indexes': [
            {
                'fields': ['resource_type', 'resource_id'],
                'unique': True,
                'partialFilterExpression': {'state': STATE_PENDING},
            },
            ('target_user_id', 'state'),
            ('source_user_id', 'state'),
        ],
    }

    @classmethod
    def get_by_transfer_id(cls, transfer_id: str) -> OwnershipTransferRecord | None:
        """Return a transfer record by id, or ``None`` if it does not exist."""
        return cls.objects(id=transfer_id).first()

    @classmethod
    def claim_pending(
        cls,
        transfer_id: str,
        claimed_state: str,
        **filters,
    ) -> OwnershipTransferRecord | None:
        """Atomically claim a pending transfer record by id and optional filters."""
        if claimed_state not in {cls.STATE_RESPONDING, cls.STATE_CANCELING}:
            raise ValueError(f'Invalid claimed state: {claimed_state}')

        return cls.objects(
            id=transfer_id,
            state=cls.STATE_PENDING,
            **filters,
        ).modify(
            new=True,
            set__state=claimed_state,
            set__updated_at=now(),
        )

    @classmethod
    def release_claim(cls, transfer_id: str, claimed_state: str) -> bool:
        """Best-effort release of a claimed transfer record back to pending by id."""
        record = cls.objects(id=transfer_id, state=claimed_state).modify(
            new=True,
            set__state=cls.STATE_PENDING,
            set__updated_at=now(),
        )
        return record is not None

    def is_expired(self, expiry_cutoff: datetime) -> bool:
        """Return whether this record is older than *expiry_cutoff*."""
        requested_at = self.requested_at
        if requested_at.tzinfo is None:
            requested_at = requested_at.replace(tzinfo=expiry_cutoff.tzinfo)
        return requested_at < expiry_cutoff

    @classmethod
    def create_or_replace(
        cls,
        resource_type: str,
        resource_id: str,
        source_user_id: str,
        target_user_id: str,
    ) -> OwnershipTransferRecord:
        """Create or replace the active transfer record for one resource."""
        current_time = now()
        record = cls.objects(
            resource_type=resource_type,
            resource_id=resource_id,
            state__in=list(cls.ACTIVE_STATES),
        ).modify(
            new=True,
            set__resource_type=resource_type,
            set__resource_id=resource_id,
            set__source_user_id=source_user_id,
            set__target_user_id=target_user_id,
            set__state=cls.STATE_PENDING,
            set__requested_at=current_time,
            set__updated_at=current_time,
            unset__actor_user_id=1,
        )
        if record is None:
            record = cls(
                resource_type=resource_type,
                resource_id=resource_id,
                source_user_id=source_user_id,
                target_user_id=target_user_id,
                state=cls.STATE_PENDING,
                requested_at=current_time,
                updated_at=current_time,
            )
            record.save()
        return record

    def refuse(self, actor_user_id: str) -> None:
        """Mark this record as refused and record the actor."""
        self.state = self.STATE_REFUSED
        self.updated_at = now()
        self.actor_user_id = actor_user_id
        self.save()
