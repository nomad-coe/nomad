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

from typing import Union
from fastapi import APIRouter, Query, Path, HTTPException, status, Depends
from datetime import datetime, date
from elasticsearch_dsl import Q

from nomad import utils
from nomad.utils import strip
from nomad.search import search
from nomad.app.v1.models import MetadataPagination, HTTPExceptionModel
from nomad.app.v1.utils import create_responses

from ..common import rdf_response
from ..mapping import Mapping

router = APIRouter()
default_tag = 'dcat'

logger = utils.get_logger(__name__)


_bad_id_response = status.HTTP_404_NOT_FOUND, {
    'model': HTTPExceptionModel,
    'description': strip('''
        Dataset not found. The given id does not match any dataset.''')}

_raw_response = 200, {
    'content': {'application/octet': {}},
    'description': 'The response. The returned content type depends on the format parameter.'}


@router.get(
    '/datasets/{entry_id}', tags=[default_tag],
    summary='Returns a DCAT dataset for a given NOMAD entry id.',
    responses=create_responses(_bad_id_response, _raw_response)
)
async def get_dataset(
    entry_id: str = Path(..., description='The unique NOMAD entry id.'),
    rdf_respose=Depends(rdf_response)
):
    ''' Returns a DCAT dataset for a given NOMAD entry id. '''

    results = search(owner='public', query=dict(entry_id=entry_id))
    if results.pagination.total == 0:
        raise HTTPException(
            status_code=_bad_id_response[0],
            detail=_bad_id_response[1]['description'])

    entry = results.data[0]

    mapping = Mapping()
    mapping.map_entry(entry)
    return rdf_respose(mapping.g)


@router.get(
    '/catalog/', tags=[default_tag],
    summary='Returns a DCAT dataset for a given NOMAD entry id.',
    responses=create_responses(_raw_response)
)
async def get_catalog(
    after: str = Query(None, description='return entries after the given entry_id'),
    modified_since: Union[datetime, date] = Query(None, description='maximum entry time (e.g. upload time)'),
    rdf_respose=Depends(rdf_response)
):
    ''' Returns a DCAT dataset for a given NOMAD entry id. '''

    search_query = Q()
    if modified_since is not None:
        modified_clause = Q('range', upload_create_time=dict(gte=modified_since))
        modified_clause |= Q('range', last_edit_time=dict(gte=modified_since))
        modified_clause |= Q('range', last_processing_time=dict(gte=modified_since))
        search_query &= modified_clause

    pagination = MetadataPagination(page_after_value=after)

    search_response = search(owner='public', query=search_query, pagination=pagination)

    mapping = Mapping()
    mapping.map_catalog(
        search_response.data,
        total=search_response.pagination.total,
        after=after,
        modified_since=modified_since, slim=False)
    return rdf_respose(mapping.g)
