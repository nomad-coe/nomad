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

from typing import Any
from fastapi import APIRouter, Depends, Path, status, HTTPException, Request
from fastapi.exception_handlers import RequestValidationError
from pydantic import BaseModel, Field

from nomad import utils
from nomad.utils import strip
from nomad.search import AuthenticationRequiredError, SearchError
from nomad.search import search, QueryValidationError
from nomad.metainfo.elasticsearch_extension import material_type, material_index

from .auth import create_user_dependency
from ..utils import create_responses
from ..models import (
    User, Owner, WithQuery,
    MetadataResponse, Metadata, MetadataPagination, MetadataRequired,
    metadata_pagination_parameters, metadata_required_parameters, QueryParameters,
    HTTPExceptionModel)


router = APIRouter()

logger = utils.get_logger(__name__)

query_parameters = QueryParameters(doc_type=material_type)

_bad_owner_response = status.HTTP_401_UNAUTHORIZED, {
    'model': HTTPExceptionModel,
    'description': strip('''
        Unauthorized. The given owner requires authorization,
        but no or bad authentication credentials are given.''')}

_bad_id_response = status.HTTP_404_NOT_FOUND, {
    'model': HTTPExceptionModel,
    'description': strip('''
        Material not found. The given id does not match any material.''')}


class MaterialMetadataResponse(BaseModel):
    material_id: str = Field(None)
    required: MetadataRequired = Field(None)
    data: Any = Field(
        None, description=strip('''The material metadata as dictionary.'''))


def perform_search(*args, **kwargs) -> MetadataResponse:
    kwargs.update(index=material_index)
    try:
        search_response = search(*args, **kwargs)
        search_response.es_query = None
        return search_response
    except QueryValidationError as e:
        raise RequestValidationError(errors=e.errors)
    except AuthenticationRequiredError as e:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=str(e))
    except SearchError as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='Elasticsearch could not process your query: %s' % str(e))


@router.post(
    '/query', tags=['materials'],
    summary='Search materials and retrieve their metadata',
    response_model=MetadataResponse,
    responses=create_responses(_bad_owner_response),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def post_entries_metadata_query(
        request: Request,
        data: Metadata,
        user: User = Depends(create_user_dependency())):

    '''
    Executes a *query* and returns a *page* of the results with *required* result data
    as well as *statistics* and *aggregated* data.

    This is the basic search operation to retrieve metadata for entries that match
    certain search criteria (`query` and `owner`). All parameters (including `query`, `owner`)
    are optional. Look at the body schema or parameter documentation for more details.

    By default the *empty* search (that returns everything) is performed. Only a small
    page of the search results are returned at a time; use `pagination` in subsequent
    requests to retrive more data. Each entry has a lot of different *metadata*, use
    `required` to limit the data that is returned.

    The `statistics` and `aggregations` keys will further allow to return statistics
    and aggregated data over all search results.
    '''

    return perform_search(
        owner=data.owner,
        query=data.query,
        pagination=data.pagination,
        required=data.required,
        aggregations=data.aggregations,
        user_id=user.user_id if user is not None else None)


@router.get(
    '', tags=['materials'],
    summary='Search materials and retrieve their metadata',
    response_model=MetadataResponse,
    responses=create_responses(_bad_owner_response),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_entries_metadata(
        request: Request,
        with_query: WithQuery = Depends(query_parameters),
        pagination: MetadataPagination = Depends(metadata_pagination_parameters),
        required: MetadataRequired = Depends(metadata_required_parameters),
        user: User = Depends(create_user_dependency())):
    '''
    Executes a *query* and returns a *page* of the results with *required* result data.
    This is a version of `/entries/query`. Queries work a little different, because
    we cannot put complex queries into URL parameters.

    In addition to the `q` parameter (see parameter documentation for details), you can use all NOMAD
    search quantities as parameters, e.g. `?atoms=H&atoms=O`. Those quantities can be
    used with additional operators attached to their names, e.g. `?n_atoms__gte=3` for
    all entries with more than 3 atoms. Operators are `all`, `any`, `none`, `gte`,
    `gt`, `lt`, `lte`.
    '''

    res = perform_search(
        owner=with_query.owner, query=with_query.query,
        pagination=pagination, required=required,
        user_id=user.user_id if user is not None else None)
    res.pagination.populate_urls(request)
    return res


@router.get(
    '/{material_id}', tags=['materials'],
    summary='Get the metadata of a material by its id',
    response_model=MaterialMetadataResponse,
    responses=create_responses(_bad_id_response),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_material_metadata(
        material_id: str = Path(..., description='The unique material id of the material to retrieve metadata from.'),
        required: MetadataRequired = Depends(metadata_required_parameters),
        user: User = Depends(create_user_dependency())):
    '''
    Retrives the material metadata for the given id.
    '''

    query = {'material_id': material_id}
    response = perform_search(
        owner=Owner.all_, query=query, required=required,
        user_id=user.user_id if user is not None else None)

    if response.pagination.total == 0:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail='The entry with the given id does not exist or is not visible to you.')

    return {
        'material_id': material_id,
        'required': required,
        'data': response.data[0]
    }
