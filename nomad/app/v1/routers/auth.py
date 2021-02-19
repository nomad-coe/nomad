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

from typing import cast
from fastapi import Depends, APIRouter, HTTPException, status
from fastapi.security import OAuth2PasswordBearer, OAuth2PasswordRequestForm
from pydantic import BaseModel

from nomad import infrastructure, config, datamodel
from nomad.utils import get_logger, strip

from ..common import root_path
from ..models import User, HTTPExceptionModel
from ..utils import create_responses

logger = get_logger(__name__)

router = APIRouter()
default_tag = 'auth'


class Token(BaseModel):
    access_token: str
    token_type: str


oauth2_scheme = OAuth2PasswordBearer(tokenUrl=f'{root_path}/auth/token', auto_error=False)


async def get_optional_user(access_token: str = Depends(oauth2_scheme)) -> User:
    '''
    A dependency that provides the authenticated (if credentials are available) or None.
    '''
    if access_token is None:
        user: datamodel.User = None
    else:
        try:
            user = cast(datamodel.User, infrastructure.keycloak.tokenauth(access_token))
        except infrastructure.KeycloakError as e:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail=str(e), headers={'WWW-Authenticate': 'Bearer'})

    if config.oasis.allowed_users is not None:
        if user is None:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail='Authentication is required for this Oasis',
                headers={'WWW-Authenticate': 'Bearer'})

        if user.email not in config.oasis.allowed_users:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail='You are not authorized to access this Oasis',
                headers={'WWW-Authenticate': 'Bearer'})

    return user


async def get_required_user(user: User = Depends(get_optional_user)) -> User:
    '''
    A dependency that provides the authenticated user or raises 401 if no user is
    authenticated.
    '''
    if user is None:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail='Authentication required',
            headers={'WWW-Authenticate': 'Bearer'})

    return user


_bad_credentials_response = status.HTTP_401_UNAUTHORIZED, {
    'model': HTTPExceptionModel,
    'description': strip('''
        Unauthorized. The provided credentials were not recognized.''')}


@router.post(
    '/token',
    tags=[default_tag],
    summary='Get an access token',
    responses=create_responses(_bad_credentials_response),
    response_model=Token)
async def get_token(form_data: OAuth2PasswordRequestForm = Depends()):
    '''
    This API uses OAuth as an authentication mechanism. This operation allows you to
    retrieve an *access token* by posting username and password as form data.

    This token can be used on subsequent API calls to authenticate
    you. Operations that support or require authentication will expect the *access token*
    in an HTTP Authorization header like this: `Authorization: Bearer <access token>`.

    On the OpenAPI dashboard, you can use the *Authorize* button at the top.

    You only need to provide `username` and `password` values. You can ignore the other
    parameters.
    '''

    try:
        access_token = infrastructure.keycloak.basicauth(
            form_data.username, form_data.password)
    except infrastructure.KeycloakError:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail='Incorrect username or password',
            headers={'WWW-Authenticate': 'Bearer'})

    return {'access_token': access_token, 'token_type': 'bearer'}


@router.get(
    '/token',
    tags=[default_tag],
    summary='Get an access token',
    responses=create_responses(_bad_credentials_response),
    response_model=Token)
async def get_token_via_query(username: str, password: str):
    '''
    This is an convenience alternative to the **POST** version of this operation.
    It allows you to retrieve an *access token* by providing username and password.
    '''

    try:
        access_token = infrastructure.keycloak.basicauth(username, password)
    except infrastructure.KeycloakError:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail='Incorrect username or password',
            headers={'WWW-Authenticate': 'Bearer'})

    return {'access_token': access_token, 'token_type': 'bearer'}
