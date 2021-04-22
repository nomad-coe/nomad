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

import hmac
import hashlib
import uuid
from typing import Callable, cast
from inspect import Parameter, signature
from functools import wraps
from fastapi import APIRouter, Depends, Query as FastApiQuery, HTTPException, status
from fastapi.security import OAuth2PasswordBearer, OAuth2PasswordRequestForm
from pydantic import BaseModel

from nomad import utils, infrastructure, config, datamodel
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


def create_user_dependency(
        required: bool = False,
        basic_auth_allowed: bool = False,
        bearer_token_auth_allowed: bool = True,
        upload_token_auth_allowed: bool = False) -> Callable:
    '''
    Creates a dependency for getting the authenticated user. The parameters define if
    the authentication is required or not, and which authentication methods are allowed.
    '''

    def user_dependency(**kwargs) -> User:
        user = None
        if basic_auth_allowed:
            user = _get_user_basic_auth(kwargs.get('form_data'))
        if not user and bearer_token_auth_allowed:
            user = _get_user_bearer_token_auth(kwargs.get('bearer_token'))
        if not user and upload_token_auth_allowed:
            user = _get_user_upload_token_auth(kwargs.get('token'))

        if required and not user:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail='Authorization required.')

        if config.oasis.allowed_users is not None:
            # We're an oasis, and have allowed_users set
            if not user:
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

    # Create the desired function signature (as it depends on which auth options are allowed)
    new_parameters = []
    if basic_auth_allowed:
        new_parameters.append(
            Parameter(
                name='form_data',
                annotation=OAuth2PasswordRequestForm,
                default=Depends(),
                kind=Parameter.KEYWORD_ONLY))
    if bearer_token_auth_allowed:
        new_parameters.append(
            Parameter(
                name='bearer_token',
                annotation=str,
                default=Depends(oauth2_scheme),
                kind=Parameter.KEYWORD_ONLY))
    if upload_token_auth_allowed:
        new_parameters.append(
            Parameter(
                name='token',
                annotation=str,
                default=FastApiQuery(
                    None,
                    description='Token for simplified authorization for uploading.'),
                kind=Parameter.KEYWORD_ONLY))

    # Create a wrapper around user_dependency, and set the signature on it
    @wraps(user_dependency)
    def wrapper(**kwargs) -> Callable:
        return user_dependency(**kwargs)

    sig = signature(user_dependency)
    sig = sig.replace(parameters=tuple(new_parameters))
    wrapper.__signature__ = sig  # type: ignore
    return wrapper


def _get_user_basic_auth(form_data: OAuth2PasswordRequestForm) -> User:
    '''
    Verifies basic auth (username and password), throwing an exception if illegal credentials
    are provided, and returns the corresponding user object if successful, None if no
    credentials provided.
    '''
    if form_data and form_data.username and form_data.password:
        try:
            infrastructure.keycloak.basicauth(form_data.username, form_data.password)
            user = cast(datamodel.User, infrastructure.keycloak.get_user(form_data.username))
            return user
        except infrastructure.KeycloakError:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail='Incorrect username or password',
                headers={'WWW-Authenticate': 'Bearer'})
    return None


def _get_user_bearer_token_auth(bearer_token: str) -> User:
    '''
    Verifies bearer_token (throwing exception if illegal value provided) and returns the
    corresponding user object, or None, if no bearer_token provided.
    '''
    if bearer_token:
        try:
            user = cast(datamodel.User, infrastructure.keycloak.tokenauth(bearer_token))
            return user
        except infrastructure.KeycloakError as e:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail=str(e), headers={'WWW-Authenticate': 'Bearer'})
    return None


def _get_user_upload_token_auth(upload_token: str) -> User:
    '''
    Verifies the upload token (throwing exception if illegal value provided) and returns the
    corresponding user object, or None, if no upload_token provided.
    '''
    if upload_token:
        try:
            payload, signature = upload_token.split('.')
            payload_bytes = utils.base64_decode(payload)
            signature_bytes = utils.base64_decode(signature)

            compare = hmac.new(
                bytes(config.services.api_secret, 'utf-8'),
                msg=payload_bytes,
                digestmod=hashlib.sha1)

            if signature_bytes == compare.digest():
                user_id = str(uuid.UUID(bytes=payload_bytes))
                user = cast(datamodel.User, infrastructure.keycloak.get_user(user_id))
                return user
        except Exception:
            # Decode error, format error, user not found, etc.
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail='A invalid upload token was supplied.')
    return None


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


def generate_upload_token(user):
    payload = uuid.UUID(user.user_id).bytes
    signature = hmac.new(
        bytes(config.services.api_secret, 'utf-8'),
        msg=payload,
        digestmod=hashlib.sha1)

    return '%s.%s' % (
        utils.base64_encode(payload),
        utils.base64_encode(signature.digest()))
