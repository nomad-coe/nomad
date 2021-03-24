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
from datetime import datetime
from typing import Optional, Tuple, List, Set, Dict, Any
from pydantic import BaseModel, Field, validator
from fastapi import APIRouter, Depends, Query as FastApiQuery, status, HTTPException

from nomad import utils
from nomad.processing import Upload
from nomad.utils import strip

from .auth import get_optional_user
from ..models import (BaseModel, User, Owner, Direction, Pagination,
                      PaginationResponse)
from ..utils import parameter_dependency_from_model

router = APIRouter()
default_tag = 'uploads'

logger = utils.get_logger(__name__)


UPLOAD_FIELDS: Set[str] = set(Upload._fields.keys())
UPLOAD_FIELDS_SORTABLE: Tuple[str, ...] = ('upload_time', 'upload_name', 'current_process',
                                           'process_status')


class UploadPagination(Pagination):
    order_by: Optional[str] = Field(
        None, description=strip('''
        The search results are ordered by the values of this quantity. The response
        either contains the first `size` value or the next `size` values after `after`.
        '''))

    @validator('order_by')
    def validate_order_by(cls, order_by):  # pylint: disable=no-self-argument
        if order_by is None:
            return order_by
        assert order_by in UPLOAD_FIELDS_SORTABLE, 'order_by must be a valid attribute'
        return order_by


upload_pagination_parameters = parameter_dependency_from_model(
    'upload_pagination_parameters', UploadPagination)


class UploadsMetadataResponse(BaseModel):
    pagination: PaginationResponse
    data: List[Dict[str, Any]] = Field(
        None, description=strip('''
        The upload metadata as a list. Each item is a dictionary with the metadata for each
        upload.'''))


@router.get(
    '', tags=[default_tag],
    summary='List uploads of authenticated user.',
    response_model=UploadsMetadataResponse,
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_uploads(
        owner: Owner = FastApiQuery(
            Owner.user,
            description=strip(Owner.__doc__)),
        modified_since: datetime = FastApiQuery(
            datetime(1900, 1, 1),
            description=strip('''
            Specify to retrieve only uploads modified after a specific point in time.
            ''')),
        processing_filters: List[str] = FastApiQuery(
            [],
            description=strip('''
            Specify to retrieve only uploads with specific processing statuses.''')),
        required: List[str] = FastApiQuery(
            [],
            description=strip('''
            Specify to retrieve only selected fields.''')),
        pagination: UploadPagination = Depends(upload_pagination_parameters),
        user: User = Depends(get_optional_user)):
    '''
        TODO: doc me!
    '''
    if owner != Owner.user:
        # In the future, we may support other values
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=strip('''
            Only owner = user is currently supported.
            '''))
    # Check access
    if not user:
        # For now, login is required in all cases
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=strip(f'''
            You need to be logged in to search for uploads with owner = {owner}.
            '''))

    # Check query
    if required:
        for attr in required:
            if attr not in UPLOAD_FIELDS:
                raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=strip(f'''
                    Invalid attribute requested: `{attr}`.
                    '''))
    if processing_filters:
        pass  # TODO: check that all values are valid

    # Build query
    query_kwargs: Dict[str, Any] = {}
    query_kwargs.update(user_id=str(user.user_id))
    # TODO: apply modified_since
    # TODO: apply processing filters

    # Fetch data from DB
    mongodb_query = Upload.objects(**query_kwargs)

    # Create response
    order_by = pagination.order_by or 'upload_time'
    start = 0
    if pagination.after is not None:
        start = int(pagination.after)
    end = start + pagination.size

    order_by_with_sign = order_by if pagination.order == Direction.asc else '-' + order_by
    mongodb_query = mongodb_query.order_by(order_by_with_sign, 'upload_id')

    data = [_upload_to_dict(upload, required) for upload in mongodb_query[start:end]]

    pagination_response = PaginationResponse(
        total=mongodb_query.count(),
        next_after=str(end),
        **pagination.dict())
    pagination_response.order_by = order_by  # Since these may have been defaulted
    pagination_response.after = str(start)

    return UploadsMetadataResponse(
        pagination=pagination_response,
        data=data)


def _upload_to_dict(upload, required):
    ''' Extracts the required fields from the provided upload, and returns as dict. '''
    if required:
        rv = {}
        for attr in required:
            if hasattr(upload, attr):
                rv[attr] = getattr(upload, attr)
            else:
                rv[attr] = None
        return rv
    else:
        rv = dict(upload.to_mongo())
        rv['upload_id'] = rv.pop('_id')
        return rv
