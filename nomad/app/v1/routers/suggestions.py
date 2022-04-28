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

from typing import List, Dict, Optional, Set
from collections import defaultdict
from pydantic import BaseModel, Field
from fastapi import APIRouter, Depends, Request, HTTPException, status
from elasticsearch_dsl import Search
from elasticsearch_dsl.utils import AttrList
from elasticsearch.exceptions import RequestError

from nomad.metainfo.elasticsearch_extension import entry_index, entry_type

from .auth import create_user_dependency
from ..models import User


router = APIRouter()

# This is a dynamically create enum class for enumerating all allowed
# quantities. FastAPI uses python enums to validate and document options.
suggestable_quantities: Set[str] = None


class SuggestionError(Exception): pass


class Suggestion(BaseModel):
    value: str = Field(None, description='The returned suggestion.')
    weight: Optional[float] = Field(None, description='The suggestion weight.')


class SuggestionsRequest(BaseModel):
    quantities: List[str] = Field(  # type: ignore
        None,
        description='List of quantities for which the suggestions are retrieved.'
    )
    input: str = Field(
        None,
        description='The input that is used as a basis for returning a suggestion.'
    )


@router.post(
    '',
    tags=['suggestions'],
    summary='Get a list of suggestions for the given quantity names and input.',
    response_model=Dict[str, List[Suggestion]],
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_suggestions(
        request: Request,
        data: SuggestionsRequest,
        user: User = Depends(create_user_dependency())):

    global suggestable_quantities
    if suggestable_quantities is None:
        suggestable_quantities = set(entry_type.suggestions.keys())

    search = Search(index=entry_index.index_name)
    quantities: List[str] = data.quantities
    quantities_es = [x.replace(".", "-") for x in quantities]
    for index, (quantity, quantity_es) in enumerate(zip(quantities, quantities_es)):
        if quantity not in suggestable_quantities:
            raise HTTPException(
                status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
                detail=[dict(
                    msg=(
                        f'The string "{quantity}" does not represent a suggestable quantity. '
                        f'Possible values are {", ".join(suggestable_quantities)}.'),
                    loc=['quantities', index])
                ]
            )

        annotation = entry_type.suggestions[quantity]
        postfix = ".suggestion" if annotation.field else "__suggestion"
        search = search.suggest(quantity_es, data.input, completion={
            'field': f'{quantity}{postfix}',
            'size': 5,
            'skip_duplicates': True,
            # 'fuzzy': {},
        })
    search = search.extra(_source=quantities)

    try:
        # For some reason calling the search.extra()-method messes up the type
        # information for the Search-object. This is why linting is disabled
        # here.
        es_response = search.execute()  # pylint: disable=no-member
    except RequestError as e:
        raise SuggestionError(e)

    # We return the original field in the source document.
    response: Dict[str, List[Suggestion]] = defaultdict(list)
    suggestions: Dict[str, Dict[str, float]] = defaultdict(dict)

    def add_suggestion(quantity, value, weight):
        values = suggestions[quantity]
        if value not in values or values[value] < weight:
            suggestions[quantity][value] = weight

    for quantity, quantity_es in zip(quantities, quantities_es):
        variants = entry_type.suggestions[quantity].variants
        for option in es_response.suggest[quantity_es][0].options:
            weight = option._score

            # We use the original input text to do the matching. This works
            # better than the text returned by the completion suggester
            # (option.text), since it can match several items if there are
            # multiple values per quantity.
            text = data.input.lower().strip()

            # Nested fields use the nested document as _source: we need to
            # modify the path accordingly.
            try:
                nested_field = option._nested.field
            except AttributeError:
                nested_field = None
            quantity_path = quantity[len(nested_field):] if nested_field else quantity

            # Gather all options recursively. The source may contain
            # inner/nested sections and also values may be lists.
            def gather_options(root, parts, options):
                original = root
                for i, part in enumerate(parts):
                    if part != "":
                        original = original[part]
                        if isinstance(original, AttrList):
                            for item in original:
                                gather_options(item, parts[i + 1:], options)
                            return
                original_type = type(original)
                if original_type == str:
                    options.append(original)
                elif original_type == list:
                    for item in original:
                        options.append(item)
            options: List[str] = []
            parts = quantity_path.split('.')
            gather_options(option._source, parts, options)

            # There may be multiple options and we have to look which options
            # were actually matched (the completion suggester does not have a
            # mechanism for saving which item in a list was matched.). If the
            # value has several variants, we have to expand the original source
            # value and return only the best match.
            for option in options:
                if variants:
                    best_match = float("Inf")
                    best_option = option
                    for variant in variants(option):
                        match_start = variant.lower().strip().find(text)
                        if match_start >= 0 and match_start < best_match:
                            best_match = match_start
                            best_option = variant
                    add_suggestion(quantity, best_option, weight)
                elif text in option.lower().strip():
                    add_suggestion(quantity, option, weight)

    for name, suggestion in suggestions.items():
        for value, weight in suggestion.items():
            response[name].append(Suggestion(value=value, weight=weight))

    return response
