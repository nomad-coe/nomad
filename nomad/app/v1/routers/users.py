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

from fastapi import Depends, APIRouter, status

from nomad.utils import strip

from .auth import create_user_dependency
from ..models import User, HTTPExceptionModel
from ..utils import create_responses

router = APIRouter()
default_tag = 'users'


_authentication_required_response = status.HTTP_401_UNAUTHORIZED, {
    'model': HTTPExceptionModel,
    'description': strip('''
        Unauthorized. The operation requires authorization,
        but no or bad authentication credentials are given.''')}


@router.get(
    '/me',
    tags=[default_tag],
    summary='Get your account data',
    description='Returnes the account data of the authenticated user.',
    responses=create_responses(_authentication_required_response),
    response_model=User)
async def read_users_me(current_user: User = Depends(create_user_dependency(required=True))):
    return current_user
