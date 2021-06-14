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

from typing import List, Dict, Any, Optional
from collections import defaultdict
from pydantic import BaseModel, Field
from fastapi import (
    APIRouter, Depends, status, Request, Path
)
from elasticsearch_dsl import Search
from elasticsearch.exceptions import RequestError

from nomad.utils import strip
from nomad.metainfo.elasticsearch_extension import entry_index

from .auth import create_user_dependency
from ..utils import create_responses
from ..models import User, HTTPExceptionModel, Aggregation, TermsAggregation, MetadataRequired
from .entries import perform_search


router = APIRouter()


class SuggestionError(Exception): pass


class Suggestion(BaseModel):
    value: str = Field(None, description='The returned suggestion.')
    weight: Optional[float] = Field(None, description='The suggestion weight.')


class SuggestionsRequest(BaseModel):
    quantities: List[str] = Field(None, description='List of quantities from which the suggestions are retrieved.')
    input: str = Field(None, description='The input that is used as a basis for returning a suggestion.')


class SuggestionsQuantityRequest(BaseModel):
    query: Any = Field(None, description='The query that is used to restrict the suggestion context.')
    input: str = Field(None, description='The input that is used as a basis for returning a suggestion.')


_bad_quantity_response = status.HTTP_404_NOT_FOUND, {
    'model': HTTPExceptionModel,
    'description': strip('''
        Suggestions for quantity not found. The given quantity name does not have any suggestions associated with it.
    ''')
}


@router.post(
    '',
    tags=['suggestions'],
    summary='Get a list of suggestions for the given quantity names and input.',
    response_model=Dict[str, List[Suggestion]],
    responses=create_responses(_bad_quantity_response),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_suggestions(
        request: Request,
        data: SuggestionsRequest,
        user: User = Depends(create_user_dependency())):

    search = Search(index=entry_index.index_name)
    quantities_es = [x.replace(".", "-") for x in data.quantities]
    for quantity, quantity_es in zip(data.quantities, quantities_es):
        search = search.suggest(quantity_es, data.input, completion={
            'field': '{}.suggestion'.format(quantity),
            'size': 5,
            'skip_duplicates': True,
            # 'fuzziness': {},
        })
    search = search.extra(_source='suggest')

    try:
        # For some reason calling the search.extra()-method messes up the type
        # information for the Search-object. This is why linting is disabled
        # here.
        es_response = search.execute()  # pylint: disable=no-member
    except RequestError as e:
        raise SuggestionError(e)

    response: Dict[str, List[Suggestion]] = defaultdict(list)
    for quantity, quantity_es in zip(data.quantities, quantities_es):
        for option in es_response.suggest[quantity_es][0].options:
            response[quantity].append(Suggestion(value=option.text, weight=option._score))

    return response


@router.post(
    '/{quantity}',
    tags=['suggestions_quantity'],
    summary='''Get a list of suggestions for the given quantity, in the context
    of the given query. In order to give meaningful suggestions, the targeted
    quanitity should be of MEnum type.''',
    response_model=List[Suggestion],
    responses=create_responses(_bad_quantity_response),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_suggestion_quantity(
        data: SuggestionsQuantityRequest,
        request: Request,
        quantity: str = Path(..., description='Name of the the quantity for which suggestions are retrieved for.'),
        user: User = Depends(create_user_dependency())):

    # If there is no input based on which to filter the values, we simply
    # provide a set of options that may not cover all values.
    if not data.input:
        aggs = {}
        aggs[quantity] = Aggregation(terms=TermsAggregation(quantity=quantity, size=100))
        response_es = perform_search(
            owner="visible",
            query=data.query,
            required=MetadataRequired(include=[]),
            aggregations=aggs,
            user_id=user.user_id if user is not None else None)

    # TODO: Perform match_phrase_prefix query if an input is given

    response = [Suggestion(value=x.value) for x in response_es.aggregations[quantity].terms.data]  # pylint: disable=no-member
    return response
