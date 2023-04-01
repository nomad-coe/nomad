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

from nomad.archive.query_reader import MongoReader, ConfigError, GeneralReader
from .auth import create_user_dependency
from ..models import User

router = APIRouter()
default_tag = 'graph'


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

    if GeneralReader.__CACHE__ in response:
        del response[GeneralReader.__CACHE__]

    return {'query': query, 'm_response': response}
