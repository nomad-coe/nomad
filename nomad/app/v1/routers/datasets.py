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

import re
from typing import cast, Optional, List
from fastapi import (
    APIRouter, Request, Depends, Query as FastApiQuery, Path, HTTPException, status)
from pydantic import BaseModel, Field, validator
from datetime import datetime
import enum

from nomad import utils, datamodel, processing, config
from nomad.metainfo.elasticsearch_extension import entry_type
from nomad.utils import strip, create_uuid
from nomad.datamodel import Dataset as DatasetDefinitionCls
from nomad.doi import DOI, DOIException
from nomad.search import search, update_by_query

from .auth import create_user_dependency
from .entries import _do_exaustive_search
from ..utils import create_responses, parameter_dependency_from_model
from ..models import (
    Pagination, PaginationResponse, MetadataPagination, Query, HTTPExceptionModel,
    User, Direction, Owner, Any_)


router = APIRouter()
default_tag = 'datasets'

logger = utils.get_logger(__name__)


_bad_id_response = status.HTTP_404_NOT_FOUND, {
    'model': HTTPExceptionModel,
    'description': strip('''
        Dataset not found. The given id does not match any dataset.''')}

_bad_user_response = status.HTTP_401_UNAUTHORIZED, {
    'model': HTTPExceptionModel,
    'description': strip('''
        The dataset can only be edited by the user who created the dataset.''')}

_bad_owned_dataset_response = status.HTTP_400_BAD_REQUEST, {
    'model': HTTPExceptionModel,
    'description': strip('''
        The requested action cannot be performed for this type of dataset.
        Owned datasets can only have entries that where uploaded by the user that
        creates the dataset.
    ''')}

_existing_name_response = status.HTTP_400_BAD_REQUEST, {
    'model': HTTPExceptionModel,
    'description': strip('''
        The chosen dataset name is already taken. Datesets of the same user must have a
        unique name.
    ''')}

_dataset_is_fixed_response = status.HTTP_400_BAD_REQUEST, {
    'model': HTTPExceptionModel,
    'description': strip('''
        The dataset already has a DOI and cannot be changed anymore.
    ''')}

_dataset_has_unpublished_contents = status.HTTP_400_BAD_REQUEST, {
    'model': HTTPExceptionModel,
    'description': strip('''
        The dataset has unpublished contents. No DOI can be assigned at the moment.
        Publish the dataset contents first.
    ''')}

_dataset_is_empty = status.HTTP_400_BAD_REQUEST, {
    'model': HTTPExceptionModel,
    'description': strip('''
        The dataset is empty. No DOI can be assigned at this moment. Add some published
        contents to the dataset first.
    ''')}

_datacite_did_not_resolve = status.HTTP_500_INTERNAL_SERVER_ERROR, {
    'model': HTTPExceptionModel,
    'description': strip('''
        Datacite server couldn't resolve the request. Please try again later.
    ''')}

_existing_dataset_with_findable_state = status.HTTP_400_BAD_REQUEST, {
    'model': HTTPExceptionModel,
    'description': strip('''
        The dataset was failed to to be submitted previously. It is removed now.
    ''')}

Dataset = datamodel.Dataset.m_def.a_pydantic.model


def _delete_dataset(user: User, dataset_id, dataset):
    es_query = cast(Query, {'datasets.dataset_id': dataset_id})
    entries = _do_exaustive_search(
        owner=Owner.user, query=es_query, user=user,
        include=['entry_id'])
    entry_ids = [entry['entry_id'] for entry in entries]
    mongo_query = {'_id': {'$in': entry_ids}}

    dataset.delete()

    if len(entry_ids) > 0:
        processing.Entry._get_collection().update_many(
            mongo_query, {'$pull': {'datasets': dataset.dataset_id}})
        update_by_query(
            '''
                int index = -1;
                for (int i = 0; i < ctx._source.datasets.length; i++) {
                    if (ctx._source.datasets[i].dataset_id == params.dataset_id) {
                        index = i
                    }
                }
                if (index != -1) {
                    ctx._source.datasets.remove(index);
                }
            ''',
            params=dict(dataset_id=dataset_id),
            query=es_query, user_id=user.user_id, refresh=True)


class DatasetPagination(Pagination):
    @validator('order_by')
    def validate_order_by(cls, order_by):  # pylint: disable=no-self-argument
        if order_by is None:
            return order_by
        assert order_by in ('dataset_create_time', 'dataset_modified_time', 'dataset_name'), 'order_by must be a valid attribute'
        return order_by

    @validator('page_after_value')
    def validate_page_after_value(cls, page_after_value, values):  # pylint: disable=no-self-argument
        # Validation handled elsewhere
        return page_after_value


dataset_pagination_parameters = parameter_dependency_from_model(
    'dataset_pagination_parameters', DatasetPagination)


class DatasetsResponse(BaseModel):
    pagination: PaginationResponse = Field(None)
    data: List[Dataset] = Field(None)  # type: ignore


class DatasetResponse(BaseModel):
    dataset_id: str = Field(..., description=strip('''The unique dataset id. '''))
    data: Dataset = Field()  # type: ignore


class DatasetType(str, enum.Enum):
    owned = 'owned',
    foreign = 'foreign'


class DatasetCreate(BaseModel):  # type: ignore
    dataset_name: Optional[str] = Field(None, description='The new name for the dataset.')
    dataset_type: Optional[DatasetType] = Field(None)
    query: Optional[Query] = Field(None)
    entries: Optional[List[str]] = Field(None)


@router.get(
    '/', tags=[default_tag],
    summary='Get a list of datasets',
    response_model=DatasetsResponse,
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_datasets(
        request: Request,
        dataset_id: str = FastApiQuery(None),
        dataset_name: str = FastApiQuery(None),
        user_id: List[str] = FastApiQuery(None),
        dataset_type: str = FastApiQuery(None),
        doi: str = FastApiQuery(None),
        prefix: str = FastApiQuery(None),
        pagination: DatasetPagination = Depends(dataset_pagination_parameters)):
    '''
    Retrieves all datasets that match the given criteria.
    '''

    mongodb_objects = DatasetDefinitionCls.m_def.a_mongo.objects
    query_params = dict(dataset_id=dataset_id, dataset_name=dataset_name, user_id__in=user_id, dataset_type=dataset_type, doi=doi)
    if prefix is not None and prefix != '':
        query_params.update(dataset_name=re.compile('^%s.*' % prefix, re.IGNORECASE))  # type: ignore
    query_params = {k: v for k, v in query_params.items() if v is not None}

    mongodb_query = mongodb_objects(**query_params)

    order_by = pagination.order_by
    order_by_with_sign = order_by if pagination.order == Direction.asc else '-' + order_by
    if order_by == 'dataset_create_time':
        order_by_args = [order_by_with_sign, 'dataset_id']  # Use upload_id as tie breaker
    else:
        order_by_args = [order_by_with_sign, 'dataset_create_time', 'dataset_id']

    mongodb_query = mongodb_query.order_by(*order_by_args)

    start = pagination.get_simple_index()
    end = start + pagination.page_size

    pagination_response = PaginationResponse(total=mongodb_query.count(), **pagination.dict())
    pagination_response.populate_simple_index_and_urls(request)

    return {
        'pagination': pagination_response,
        'data': list(mongodb_query[start:end])}


@router.get(
    '/{dataset_id}', tags=[default_tag],
    summary='Get a list of datasets',
    response_model=DatasetResponse,
    responses=create_responses(_bad_id_response),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_dataset(
        dataset_id: str = Path(..., description='The unique dataset id of the dataset to retrieve.')):
    '''
    Retrieves the dataset with the given id.
    '''
    mongodb_objects = DatasetDefinitionCls.m_def.a_mongo.objects
    dataset = mongodb_objects(dataset_id=dataset_id).first()

    if dataset is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail='The dataset with the given id does not exist.')

    return {
        'dataset_id': dataset_id,
        'data': dataset}


@router.post(
    '/', tags=[default_tag],
    summary='Create a new dataset',
    response_model=DatasetResponse,
    responses=create_responses(_existing_name_response),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def post_datasets(
        create: DatasetCreate, user: User = Depends(create_user_dependency(required=True))):
    '''
    Create a new dataset.
    '''

    now = datetime.utcnow()
    dataset_type = create.dataset_type if create.dataset_type is not None else DatasetType.owned

    # check if name already exists
    existing_dataset = DatasetDefinitionCls.m_def.a_mongo.objects(
        user_id=user.user_id, dataset_name=create.dataset_name).first()
    if existing_dataset is not None:
        raise HTTPException(
            status_code=_existing_name_response[0],
            detail=_existing_name_response[1]['description'])

    # create dataset
    dataset = DatasetDefinitionCls(
        dataset_id=create_uuid(),
        dataset_name=create.dataset_name,
        user_id=user.user_id,
        dataset_create_time=now,
        dataset_modified_time=now,
        dataset_type=dataset_type)
    dataset.a_mongo.create()

    # add dataset to entries in mongo and elastic
    # TODO this should be part of a new edit API
    if dataset_type != DatasetType.owned:
        dataset.query = create.query
        dataset.entrys = create.entries
        empty = True
    else:
        # add dataset to entries in mongo and elastic
        # TODO this should be part of a new edit API
        if create.entries is not None:
            es_query = cast(Query, {'entry_id': Any_(any=create.entries)})
        elif create.query is not None:
            es_query = create.query
        else:
            es_query = None

        if es_query is None:
            empty = True
        else:
            entries = _do_exaustive_search(
                owner=Owner.user, query=es_query, user=user, include=['entry_id'])
            entry_ids = [entry['entry_id'] for entry in entries]
            mongo_query = {'_id': {'$in': entry_ids}}
            empty = len(entry_ids) == 0

    if not empty:
        processing.Entry._get_collection().update_many(
            mongo_query, {'$push': {'datasets': dataset.dataset_id}})
        update_by_query(
            '''
                if (ctx._source.datasets == null) {
                    ctx._source.datasets = new ArrayList();
                }
                ctx._source.datasets.add(params.dataset);
            ''',
            params=dict(dataset=entry_type.create_index_doc(dataset)),
            query=es_query, user_id=user.user_id, refresh=True)

    return {
        'dataset_id': dataset.dataset_id,
        'data': dataset}


@router.delete(
    '/{dataset_id}', tags=[default_tag],
    summary='Delete a dataset',
    response_model=DatasetResponse,
    responses=create_responses(_bad_id_response, _dataset_is_fixed_response, _bad_user_response),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def delete_dataset(
        dataset_id: str = Path(..., description='The unique dataset id of the dataset to delete.'),
        user: User = Depends(create_user_dependency(required=True))):
    '''
    Delete an dataset.
    '''

    dataset = DatasetDefinitionCls.m_def.a_mongo.objects(dataset_id=dataset_id).first()
    if dataset is None:
        raise HTTPException(
            status_code=_bad_id_response[0],
            detail=_bad_id_response[1]['description'])

    if dataset.doi is not None and not user.is_admin:
        raise HTTPException(
            status_code=_existing_name_response[0],
            detail=_dataset_is_fixed_response[1]['description'])

    if dataset.user_id != user.user_id:
        raise HTTPException(
            status_code=_bad_user_response[0],
            detail=_bad_user_response[1]['description'])

    # delete dataset from entries in mongo and elastic
    # TODO this should be part of a new edit API
    _delete_dataset(user=user, dataset_id=dataset_id, dataset=dataset)

    return {
        'dataset_id': dataset.dataset_id,
        'data': dataset}


@router.post(
    '/{dataset_id}/action/doi', tags=[default_tag],
    summary='Assign a DOI to a dataset',
    response_model=DatasetResponse,
    responses=create_responses(
        _bad_id_response, _dataset_is_fixed_response, _dataset_has_unpublished_contents,
        _bad_user_response, _dataset_is_empty),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def assign_doi(
        dataset_id: str = Path(..., description='The unique dataset id of the dataset to delete.'),
        user: User = Depends(create_user_dependency(required=True))):
    '''
    Assign a DOI to a dataset.
    '''

    dataset = DatasetDefinitionCls.m_def.a_mongo.objects(dataset_id=dataset_id).first()
    if dataset is None:
        raise HTTPException(
            status_code=_bad_id_response[0],
            detail=_bad_id_response[1]['description'])

    if dataset.doi is not None:
        doi = DOI.objects(doi=dataset.doi).first()
        if type(doi) == DOI and not (doi.state == 'findable'):
            _delete_dataset(user=user, dataset_id=dataset_id, dataset=dataset)
            raise HTTPException(
                status_code=_existing_dataset_with_findable_state[0],
                detail=_existing_dataset_with_findable_state[1]['description'])
        raise HTTPException(
            status_code=_existing_name_response[0],
            detail=_dataset_is_fixed_response[1]['description'])

    if dataset.user_id != user.user_id:
        raise HTTPException(
            status_code=_bad_user_response[0],
            detail=_bad_user_response[1]['description'])

    response = search(
        owner='admin',
        query={'datasets.dataset_id': dataset_id},
        pagination=MetadataPagination(page_size=0),
        user_id=config.services.admin_user_id)
    if response.pagination.total == 0:
        raise HTTPException(
            status_code=_dataset_is_empty[0],
            detail=_dataset_is_empty[1]['description'])

    response = search(
        owner='admin',
        query={'datasets.dataset_id': dataset_id, 'published': False},
        pagination=MetadataPagination(page_size=0),
        user_id=config.services.admin_user_id)

    if response.pagination.total > 0:
        raise HTTPException(
            status_code=_dataset_has_unpublished_contents[0],
            detail=_dataset_has_unpublished_contents[1]['description'])

    doi = DOI.create(title='NOMAD dataset: %s' % dataset.dataset_name, user=user)

    try:
        doi.create_draft()
        doi.make_findable()
    except DOIException:
        if doi.doi:
            dataset.doi = doi
            dataset.save()

        raise HTTPException(
            status_code=_datacite_did_not_resolve[0],
            detail=_datacite_did_not_resolve[1]['description'])

    dataset.doi = doi.doi
    dataset.save()

    return {
        'dataset_id': dataset.dataset_id,
        'data': dataset}
