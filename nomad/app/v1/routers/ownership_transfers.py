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

from datetime import timedelta
from enum import Enum
from typing import Annotated, Any, Literal

from fastapi import APIRouter, Depends, HTTPException, Path, status
from fastapi import Query as FastApiQuery
from mongoengine.queryset.visitor import Q
from pydantic import BaseModel, Field

from nomad import utils
from nomad.app.v1.routers.auth import get_current_user
from nomad.auth.scopes import Scope
from nomad.common import now
from nomad.config import config
from nomad.datamodel import User as DatamodelUser
from nomad.mongo.users import OwnershipTransferRecord
from nomad.processing import Upload
from nomad.uploads import add_upload_reviewers, remove_upload_reviewers

from ..models import User
from .uploads import (
    UploadProcDataResponse,
    _check_upload_not_processing,
    _get_upload_with_write_access,
    upload_to_pydantic,
)

router = APIRouter()
logger = utils.get_logger(__name__)


class APITag(str, Enum):
    DEFAULT = 'ownership-transfers'


class OwnershipTransferCreateRequest(BaseModel):
    resource_type: str = Field(description='Transfer resource type, e.g. upload.')
    resource_id: str = Field(description='Transfer resource id.')
    target_user: str = Field(
        description='The target user identifier (user_id, username, or email).'
    )
    target_user_type: Literal['user_id', 'username', 'email'] = Field(
        description='Type of identifier provided in target_user.'
    )


class OwnershipTransferRespondRequest(BaseModel):
    action: Literal['accept', 'refuse'] = Field(
        description='Response action by target user: accept or refuse.'
    )


class OwnershipTransferResource(BaseModel):
    transfer_id: str
    resource_type: str
    resource_id: str
    resource_name: str | None = None
    source_user_id: str
    target_user_id: str
    requested_at: str
    updated_at: str
    actor_user_id: str | None = None


class OwnershipTransferResponse(BaseModel):
    transfers: list[OwnershipTransferResource]


class OwnershipTransferActionResponse(BaseModel):
    transfer_id: str
    resource_type: str
    resource_id: str | None = None
    result: dict[str, Any] | None = None


def _ensure_supported_resource_type(resource_type: str) -> None:
    if resource_type != 'upload':
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=(
                f'Unsupported resource_type={resource_type!r}. '
                'Currently only resource_type="upload" is supported.'
            ),
        )


def _map_upload_ownership_transfer_resource(
    record: OwnershipTransferRecord,
    upload_name: str | None,
) -> OwnershipTransferResource:
    return OwnershipTransferResource(
        transfer_id=str(record.id),
        resource_type='upload',
        resource_id=record.resource_id,
        resource_name=upload_name,
        source_user_id=record.source_user_id,
        target_user_id=record.target_user_id,
        requested_at=record.requested_at.isoformat(),
        updated_at=record.updated_at.isoformat(),
        actor_user_id=record.actor_user_id,
    )


def _cleanup_stale_upload_ownership_transfer_records(upload: Upload) -> None:
    """Delete stale pending transfer records and revoke stale reviewer access."""
    expiry_cutoff = now() - timedelta(
        seconds=config.uploads.ownership_transfer_record_ttl_seconds
    )
    stale_reviewer_ids: list[str] = []

    pending_records = OwnershipTransferRecord.objects(
        resource_type=OwnershipTransferRecord.RESOURCE_TYPE_UPLOAD,
        resource_id=upload.upload_id,
        state=OwnershipTransferRecord.STATE_PENDING,
    )
    for record in pending_records:
        stale = upload.main_author != record.source_user_id or record.is_expired(
            expiry_cutoff
        )
        if not stale:
            continue

        stale_reviewer_ids.append(record.target_user_id)
        record.delete()

    removed_count = remove_upload_reviewers(stale_reviewer_ids, upload=upload)
    if removed_count > 0:
        upload.reload()


def _create_upload_ownership_transfer(
    request: OwnershipTransferCreateRequest,
    user: User,
) -> OwnershipTransferResource:
    upload = _get_upload_with_write_access(
        request.resource_id,
        user,
        include_published=True,
        published_requires_admin=False,
        only_main_author=True,
    )
    _check_upload_not_processing(upload)
    _cleanup_stale_upload_ownership_transfer_records(upload)

    new_owner = None
    try:
        new_owner = DatamodelUser.get(
            **{str(request.target_user_type): request.target_user}
        )
    except KeyError:
        pass

    if new_owner is None:
        detail = (
            f'Could not resolve target user by {request.target_user_type}. '
            f'Provide a valid {request.target_user_type}.'
        )
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail=detail,
        )

    if new_owner.user_id == upload.main_author:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail='The specified user is already the owner of this upload.',
        )

    record = OwnershipTransferRecord.create_or_replace(
        resource_type=OwnershipTransferRecord.RESOURCE_TYPE_UPLOAD,
        resource_id=upload.upload_id,
        source_user_id=upload.main_author,
        target_user_id=new_owner.user_id,
    )

    add_upload_reviewers(new_owner.user_id, upload=upload)

    logger.info(
        'upload ownership transfer requested',
        upload_id=upload.upload_id,
        transfer_id=str(record.id),
        actor_user_id=user.user_id,
        actor_is_admin=user.is_admin,
        current_owner_user_id=upload.main_author,
        new_owner_user_id=new_owner.user_id,
    )

    upload.reload()
    return _map_upload_ownership_transfer_resource(record, upload.upload_name)


def _list_upload_ownership_transfers(
    direction: Literal['incoming', 'outgoing', 'all'] | None,
    resource_id: str | None,
    state: str | None,
    user: User,
) -> OwnershipTransferResponse:
    query = OwnershipTransferRecord.objects

    if direction == 'incoming':
        query = query(target_user_id=user.user_id)
    elif direction == 'outgoing':
        query = query(source_user_id=user.user_id)
    else:
        query = query(Q(target_user_id=user.user_id) | Q(source_user_id=user.user_id))

    if resource_id is not None:
        query = query(
            resource_type=OwnershipTransferRecord.RESOURCE_TYPE_UPLOAD,
            resource_id=resource_id,
        )
    else:
        query = query(resource_type=OwnershipTransferRecord.RESOURCE_TYPE_UPLOAD)

    if state is not None:
        query = query(state=state)

    expiry_cutoff = now() - timedelta(
        seconds=config.uploads.ownership_transfer_record_ttl_seconds
    )

    transfers: list[OwnershipTransferResource] = []
    for record in query.order_by('-updated_at'):
        upload = Upload.get(record.resource_id)
        if upload is None:
            continue

        if record.state == OwnershipTransferRecord.STATE_PENDING:
            stale = upload.main_author != record.source_user_id or record.is_expired(
                expiry_cutoff
            )
            if stale:
                continue

        transfers.append(
            _map_upload_ownership_transfer_resource(record, upload.upload_name)
        )

    return OwnershipTransferResponse(transfers=transfers)


def _get_upload_ownership_transfer(
    transfer_id: str, user: User
) -> OwnershipTransferResource:
    record = OwnershipTransferRecord.get_by_transfer_id(transfer_id)
    if record is None:
        raise HTTPException(
            status.HTTP_404_NOT_FOUND, detail='Transfer does not exist.'
        )

    if (
        record.target_user_id != user.user_id
        and record.source_user_id != user.user_id
        and not user.is_admin
    ):
        raise HTTPException(
            status.HTTP_403_FORBIDDEN,
            detail='You are not authorized to access this transfer.',
        )

    upload = Upload.get(record.resource_id)
    if upload is None:
        raise HTTPException(status.HTTP_404_NOT_FOUND, detail='Upload does not exist.')

    if record.state == OwnershipTransferRecord.STATE_PENDING:
        expiry_cutoff = now() - timedelta(
            seconds=config.uploads.ownership_transfer_record_ttl_seconds
        )
        stale = upload.main_author != record.source_user_id or record.is_expired(
            expiry_cutoff
        )
        if stale:
            raise HTTPException(
                status.HTTP_404_NOT_FOUND,
                detail='Transfer does not exist.',
            )

    return _map_upload_ownership_transfer_resource(record, upload.upload_name)


def _respond_to_upload_ownership_transfer(
    transfer_id: str,
    request: OwnershipTransferRespondRequest,
    user: User,
) -> UploadProcDataResponse:
    record = OwnershipTransferRecord.claim_pending(
        transfer_id,
        OwnershipTransferRecord.STATE_RESPONDING,
        target_user_id=user.user_id,
    )
    if record is None:
        existing_record = OwnershipTransferRecord.get_by_transfer_id(transfer_id)
        if existing_record is None:
            raise HTTPException(
                status.HTTP_404_NOT_FOUND,
                detail='Transfer does not exist.',
            )
        if existing_record.target_user_id != user.user_id:
            raise HTTPException(
                status.HTTP_403_FORBIDDEN,
                detail='You are not authorized to respond to this transfer request.',
            )
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail='No pending transfer request for this transfer id. The request may have expired or already been handled.',
        )

    upload = Upload.get(record.resource_id)
    if upload is None:
        record.delete()
        raise HTTPException(status.HTTP_404_NOT_FOUND, detail='Upload does not exist.')

    try:
        _check_upload_not_processing(upload)
        _cleanup_stale_upload_ownership_transfer_records(upload)
    except Exception:
        OwnershipTransferRecord.release_claim(
            transfer_id,
            OwnershipTransferRecord.STATE_RESPONDING,
        )
        raise

    if request.action == 'accept':
        if upload.main_author != record.source_user_id:
            record.delete()
            raise HTTPException(
                status.HTTP_400_BAD_REQUEST,
                detail='The transfer request is no longer valid.',
            )

        if not config.services.admin_user_id:
            OwnershipTransferRecord.release_claim(
                transfer_id,
                OwnershipTransferRecord.STATE_RESPONDING,
            )
            raise HTTPException(
                status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail='No admin user configured for transfer.',
            )

        previous_owner_user_id = upload.main_author
        try:
            upload.transfer_ownership(
                new_owner_user_id=user.user_id,
                previous_owner_user_id=previous_owner_user_id,
            )
        except Exception as e:
            OwnershipTransferRecord.release_claim(
                transfer_id,
                OwnershipTransferRecord.STATE_RESPONDING,
            )
            raise HTTPException(
                status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail=f'Failed to execute transfer workflow: {e}',
            )

        extra_log = dict(
            previous_owner_user_id=previous_owner_user_id,
            new_owner_user_id=user.user_id,
        )
    else:
        remove_upload_reviewers(user.user_id, upload=upload)
        record.refuse(actor_user_id=user.user_id)
        extra_log = dict(owner_user_id=upload.main_author)

    logger.info(
        f'upload transfer {request.action}ed',
        upload_id=upload.upload_id,
        transfer_id=transfer_id,
        actor_user_id=user.user_id,
        actor_is_admin=user.is_admin,
        **extra_log,
    )

    upload.reload()
    return UploadProcDataResponse(
        upload_id=upload.upload_id, data=upload_to_pydantic(upload)
    )


def _cancel_upload_ownership_transfer(
    transfer_id: str,
    user: User,
) -> UploadProcDataResponse:
    record = OwnershipTransferRecord.claim_pending(
        transfer_id,
        OwnershipTransferRecord.STATE_CANCELING,
        source_user_id=user.user_id,
    )
    if record is None:
        existing_record = OwnershipTransferRecord.get_by_transfer_id(transfer_id)
        if existing_record is None:
            raise HTTPException(
                status.HTTP_404_NOT_FOUND,
                detail='Transfer does not exist.',
            )
        if existing_record.source_user_id != user.user_id:
            raise HTTPException(
                status.HTTP_403_FORBIDDEN,
                detail='You are not authorized to cancel this transfer request.',
            )
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail='No pending transfer request for this transfer id. The request may have expired or already been handled.',
        )

    try:
        upload = _get_upload_with_write_access(
            record.resource_id,
            user,
            include_published=True,
            published_requires_admin=False,
            only_main_author=True,
        )
        _check_upload_not_processing(upload)
        _cleanup_stale_upload_ownership_transfer_records(upload)
    except Exception:
        OwnershipTransferRecord.release_claim(
            transfer_id,
            OwnershipTransferRecord.STATE_CANCELING,
        )
        raise

    canceled_user_ids = [record.target_user_id]
    record.delete()
    duplicate_pending_records = list(
        OwnershipTransferRecord.objects(
            resource_type=OwnershipTransferRecord.RESOURCE_TYPE_UPLOAD,
            resource_id=upload.upload_id,
            state=OwnershipTransferRecord.STATE_PENDING,
        )
    )
    if duplicate_pending_records:
        canceled_user_ids.extend(r.target_user_id for r in duplicate_pending_records)
        for duplicate_record in duplicate_pending_records:
            duplicate_record.delete()

    reviewer_access_removed = remove_upload_reviewers(canceled_user_ids, upload=upload)
    if reviewer_access_removed > 0:
        upload.reload()

    logger.info(
        'upload transfer canceled',
        upload_id=upload.upload_id,
        transfer_id=transfer_id,
        actor_user_id=user.user_id,
        actor_is_admin=user.is_admin,
        owner_user_id=upload.main_author,
        reviewer_access_removed=reviewer_access_removed,
    )
    upload.reload()
    return UploadProcDataResponse(
        upload_id=upload.upload_id, data=upload_to_pydantic(upload)
    )


@router.get(
    '',
    tags=[APITag.DEFAULT],
    summary='List transfer resources for the current user.',
    response_model=OwnershipTransferResponse,
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
def list_ownership_transfers(
    direction: Annotated[
        Literal['incoming', 'outgoing', 'all'] | None,
        FastApiQuery(description='Optional transfer direction filter.'),
    ] = 'all',
    resource_type: Annotated[
        str,
        FastApiQuery(description='Transfer resource type.'),
    ] = 'upload',
    resource_id: Annotated[
        str | None,
        FastApiQuery(description='Optional transfer resource id filter.'),
    ] = None,
    state: Annotated[
        str | None,
        FastApiQuery(description='Optional transfer state filter.'),
    ] = None,
    user: Annotated[
        User,
        Depends(get_current_user([Scope.UPLOADS_READ], allow_anonymous=False)),
    ] = None,
):
    if resource_type == 'upload':
        return _list_upload_ownership_transfers(
            direction=direction, resource_id=resource_id, state=state, user=user
        )
    _ensure_supported_resource_type(resource_type)


@router.post(
    '',
    tags=[APITag.DEFAULT],
    summary='Create a transfer resource.',
    response_model=OwnershipTransferResource,
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
def create_ownership_transfer(
    request: OwnershipTransferCreateRequest,
    user: Annotated[
        User,
        Depends(
            get_current_user(
                [Scope.UPLOADS_READ, Scope.UPLOADS_WRITE],
                allow_anonymous=False,
            )
        ),
    ],
):
    if request.resource_type == 'upload':
        return _create_upload_ownership_transfer(request=request, user=user)
    _ensure_supported_resource_type(request.resource_type)


@router.get(
    '/{transfer_id}',
    tags=[APITag.DEFAULT],
    summary='Get one transfer resource by id.',
    response_model=OwnershipTransferResource,
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
def get_ownership_transfer(
    transfer_id: Annotated[str, Path(description='Transfer id.')],
    resource_type: Annotated[
        str,
        FastApiQuery(description='Transfer resource type.'),
    ] = 'upload',
    user: Annotated[
        User,
        Depends(get_current_user([Scope.UPLOADS_READ], allow_anonymous=False)),
    ] = None,
):
    if resource_type == 'upload':
        return _get_upload_ownership_transfer(transfer_id=transfer_id, user=user)
    _ensure_supported_resource_type(resource_type)


@router.post(
    '/{transfer_id}/respond',
    tags=[APITag.DEFAULT],
    summary='Accept or refuse a pending transfer by transfer id.',
    response_model=OwnershipTransferActionResponse,
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
def respond_to_ownership_transfer(
    transfer_id: Annotated[str, Path(description='Transfer id.')],
    request: OwnershipTransferRespondRequest,
    resource_type: Annotated[
        str,
        FastApiQuery(description='Transfer resource type.'),
    ] = 'upload',
    user: Annotated[
        User,
        Depends(
            get_current_user(
                [Scope.UPLOADS_READ, Scope.UPLOADS_WRITE],
                allow_anonymous=False,
            )
        ),
    ] = None,
):
    if resource_type == 'upload':
        result = _respond_to_upload_ownership_transfer(
            transfer_id=transfer_id,
            request=OwnershipTransferRespondRequest(action=request.action),
            user=user,
        )
    else:
        _ensure_supported_resource_type(resource_type)
    return OwnershipTransferActionResponse(
        transfer_id=transfer_id,
        resource_type='upload',
        resource_id=result.upload_id,
        result=result.model_dump(exclude_none=True),
    )


@router.post(
    '/{transfer_id}/cancel',
    tags=[APITag.DEFAULT],
    summary='Cancel a pending transfer by transfer id.',
    response_model=OwnershipTransferActionResponse,
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
def cancel_ownership_transfer(
    transfer_id: Annotated[str, Path(description='Transfer id.')],
    resource_type: Annotated[
        str,
        FastApiQuery(description='Transfer resource type.'),
    ] = 'upload',
    user: Annotated[
        User,
        Depends(
            get_current_user(
                [Scope.UPLOADS_READ, Scope.UPLOADS_WRITE],
                allow_anonymous=False,
            )
        ),
    ] = None,
):
    if resource_type == 'upload':
        result = _cancel_upload_ownership_transfer(transfer_id=transfer_id, user=user)
    else:
        _ensure_supported_resource_type(resource_type)
    return OwnershipTransferActionResponse(
        transfer_id=transfer_id,
        resource_type='upload',
        resource_id=result.upload_id,
        result=result.model_dump(exclude_none=True),
    )
