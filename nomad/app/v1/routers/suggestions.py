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

from typing import List
from fastapi import (
    APIRouter, Depends, status, Request
)
from pydantic import BaseModel, Field
from nomad.utils import strip
from .auth import create_user_dependency
from ..utils import create_responses
from ..models import User, HTTPExceptionModel


router = APIRouter()


class Suggestion(BaseModel):
    value: str = Field(None, description='The returned suggestion.')
    weight: float = Field(None, description='The suggestion weight.')


class SuggestionList(BaseModel):
    __root__: List[Suggestion]


class SuggestionRequest(BaseModel):
    quantities: List[str] = Field(None, description='List of quantities from which the suggestions are retrieved.')
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
    response_model=SuggestionList,
    responses=create_responses(_bad_quantity_response),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_suggestions(
        request: Request,
        data: SuggestionRequest,
        user: User = Depends(create_user_dependency())):

    return [
        Suggestion(value="Moi"),
        Suggestion(value="Hei"),
        Suggestion(value="Terve"),
    ]
