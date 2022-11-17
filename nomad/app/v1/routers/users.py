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

from typing import List, Union, Optional
from fastapi import Depends, APIRouter, status, HTTPException, Query
from pydantic.main import BaseModel

from nomad import infrastructure, config, datamodel
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

_bad_invite_response = status.HTTP_400_BAD_REQUEST, {
    'model': HTTPExceptionModel,
    'description': strip('''
        The invite is invalid.''')}


class Users(BaseModel):
    data: List[User]


@router.get(
    '/me',
    tags=[default_tag],
    summary='Get your account data',
    description='Returns the account data of the authenticated user.',
    responses=create_responses(_authentication_required_response),
    response_model=User)
async def read_users_me(current_user: User = Depends(create_user_dependency(required=True))):
    current_user_dict: dict = current_user.m_to_dict(
        with_out_meta=True, include_derived=True)
    additional_info: dict = datamodel.User.get(user_id=current_user.user_id).m_to_dict(
        with_out_meta=True, include_derived=True)
    current_user_dict.update(additional_info)
    return current_user_dict


@router.get(
    '',
    tags=[default_tag],
    summary='Get existing users',
    description='Get existing users for given criteria',
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
    response_model=Users)
async def get_users(
    prefix: Optional[str] = Query(None, description=strip('''
            Search the user with the given prefix.
        ''')),
    user_id: Union[List[str], None] = Query(None, description=strip('''
            To get the user(s) by their user_id(s).
        ''')),
    username: Union[List[str], None] = Query(None, description=strip('''
            To get the user(s) by their username(s).
        ''')),
    email: Union[List[str], None] = Query(None, description=strip('''
            To get the user(s) by their email(s).
        '''))
):
    users: List[User] = []
    for key, values in dict(user_id=user_id, username=username, email=email).items():
        if not values:
            continue

        if isinstance(values, str):
            values = [values]

        for value in values:
            try:
                user = datamodel.User.get(**{key: str(value)}).m_copy()
                user.email = None
                users.append(user)
            except KeyError:
                pass

    if prefix:
        for user in infrastructure.user_management.search_user(prefix):
            user_dict = user.m_to_dict(include_derived=True)
            user_dict['email'] = None
            users.append(user_dict)

    return dict(data=users)


class PublicUserInfo(BaseModel):
    ''' User information that is publicly available. '''
    name: str = None
    first_name: str = None
    last_name: str = None
    affiliation: str = None
    affiliation_address: str = None
    user_id: str = None
    username: str = None


@router.get(
    '/{user_id}',
    tags=[default_tag],
    summary='Get existing users',
    description='Get the user using the given user_id',
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
    response_model=PublicUserInfo)
async def get_user(user_id: str):
    return datamodel.User.get(user_id=str(user_id)).m_to_dict(with_out_meta=True, include_derived=True)


@router.put(
    '/invite',
    tags=[default_tag],
    summary='Invite a new user',
    responses=create_responses(_authentication_required_response, _bad_invite_response),
    response_model=User)
async def invite_user(user: User, current_user: User = Depends(create_user_dependency(required=True))):
    if config.oasis.is_oasis:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='User invite does not work this NOMAD OASIS.')

    json_data = user.dict()
    try:
        user = datamodel.User.m_from_dict(json_data)
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='Invalid user data: %s' % str(e))

    if user.email is None:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='Invalid user data: email is required')

    try:
        error = infrastructure.user_management.add_user(user, invite=True)
    except KeyError as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='Invalid user data: %s' % str(e))

    if error is not None:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='Could not invite user: %s' % str(error))

    return datamodel.User.get(username=user.username), 200
