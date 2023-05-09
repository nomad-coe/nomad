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

from fastapi import Depends, APIRouter, Body, HTTPException

from nomad.archive.query_reader import MongoReader, ConfigError, GeneralReader, UserReader
from .auth import create_user_dependency
from .entries import EntriesArchive
from ..models import User

router = APIRouter()
default_tag = 'graph'


def normalise_response(query, response):
    if GeneralReader.__CACHE__ in response:
        del response[GeneralReader.__CACHE__]

    if query:
        return {'query': query, 'm_response': response}

    return {'m_response': response}


@router.post(
    '/query',
    tags=[default_tag],
    summary='Query the database with a graph style.',
    description='Query the database with a graph style.')
async def basic_query(query: dict = Body(...), user: User = Depends(create_user_dependency(required=True))):
    try:
        with MongoReader(query, user=user) as reader:
            response: dict = reader.read()
    except ConfigError as e:
        raise HTTPException(400, detail=str(e))
    except Exception as e:
        raise HTTPException(422, detail=str(e))

    return normalise_response(query, response)


@router.post(
    '/archive/query',
    tags=[default_tag],
    summary='Search entries and access their archives',
)
async def archive_query(
        data: EntriesArchive,
        user: User = Depends(create_user_dependency())
):
    graph_dict: dict = {'m_entries': {'m_request': {}}}
    root_request: dict = graph_dict['m_entries']['m_request']

    if data.pagination:
        root_request['pagination'] = data.pagination
    if data.query:
        root_request['query'] = data.query

    if data.required is None:
        root_request['directive'] = 'plain'
    else:
        graph_dict['m_entries']['*'] = {'m_archive': data.required}

    if not root_request:
        del graph_dict['m_entries']['m_request']

    with UserReader(graph_dict, user=user) as reader:
        response: dict = reader.read(user.user_id)

    return normalise_response(data.query, response)
