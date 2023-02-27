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

import os
import requests

from typing import List, Dict, cast, Optional
from enum import Enum
from pydantic import BaseModel
from fastapi import APIRouter, Depends, status, HTTPException
from mongoengine.queryset.visitor import Q

from nomad import config
from nomad.config.models import NORTHTool
from nomad.utils import strip, get_logger, slugify
from nomad.processing import Upload
from .auth import create_user_dependency, oauth2_scheme
from ..models import User, HTTPExceptionModel
from ..utils import create_responses


default_tag = 'north'
router = APIRouter()

hub_api_headers = {'Authorization': f'Bearer {config.north.hub_service_api_token}'}
logger = get_logger(__name__)


class ToolStateEnum(str, Enum):
    running = 'running',
    starting = 'starting',
    stopping = 'stopping',
    stopped = 'stopped'


class ToolModel(NORTHTool):
    name: str
    state: Optional[ToolStateEnum]


class ToolResponseModel(BaseModel):
    tool: str
    username: str
    data: ToolModel


class ToolsResponseModel(BaseModel):
    data: List[ToolModel] = []


_bad_tool_response = status.HTTP_404_NOT_FOUND, {
    'model': HTTPExceptionModel,
    'description': strip('''
        The tool does not exist.''')}


def _get_status(tool: ToolModel, user: User) -> ToolModel:
    if not user:
        return tool

    url = f'{config.hub_url()}/api/users/{user.username}/servers/{tool.name}/progress'
    response = requests.get(url, headers=hub_api_headers)

    if response.status_code == 404:
        # The user or the tool does not yet exist
        tool.state = ToolStateEnum.stopped
    elif response.status_code == 200:
        if '"ready": true' in response.text:
            tool.state = ToolStateEnum.running
        else:
            tool.state = ToolStateEnum.starting
    else:
        logger.error(
            'unexpected jupyterhub response', data=dict(status_code=response.status_code),
            text=response.text)
        tool.state = ToolStateEnum.stopped

    return tool


@router.get(
    '/', tags=[default_tag],
    response_model=ToolsResponseModel,
    summary='Get a list of all configured tools and their current state.',
    response_model_exclude_unset=True,
    response_model_exclude_none=True
)
async def get_tools(user: User = Depends(create_user_dependency())):
    return ToolsResponseModel(
        data=[
            _get_status(ToolModel(name=name, **tool.dict()), user)
            for name, tool in cast(Dict[str, NORTHTool], config.north.tools).items()
        ]
    )


async def tool(name: str) -> ToolModel:
    if name not in config.north.tools:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail='The tools does not exist.')

    tool = cast(Dict[str, NORTHTool], config.north.tools)[name]
    return ToolModel(name=name, **tool.dict())


@router.get(
    '/{name}', tags=[default_tag],
    summary='Get information for a specific tool.',
    response_model=ToolResponseModel,
    responses=create_responses(_bad_tool_response),
    response_model_exclude_unset=True,
    response_model_exclude_none=True
)
async def get_tool(
    tool: ToolModel = Depends(tool),
    user: User = Depends(create_user_dependency(required=True))
):
    return ToolResponseModel(
        tool=tool.name,
        username=user.username,
        data=_get_status(tool, user))


@router.post(
    '/{name}', tags=[default_tag],
    response_model=ToolResponseModel,
    summary='Start a tool.',
    response_model_exclude_unset=True,
    response_model_exclude_none=True
)
async def start_tool(
    tool: ToolModel = Depends(tool),
    access_token: str = Depends(oauth2_scheme),
    user: User = Depends(create_user_dependency(required=True))
):
    tool.state = ToolStateEnum.stopped

    # Make sure the user exists
    url = f'{config.hub_url()}/api/users/{user.username}'
    response = requests.get(url, headers=hub_api_headers)
    if response.status_code == 404:
        response = requests.post(url, headers=hub_api_headers)
        if response.status_code == 200:
            logger.info('created north user', user_id=user.user_id)
        else:
            # TODO
            logger.error('could not create north user', user_id=user.user_id)

    # Make sure that the home folder of the user exists
    user_home = os.path.join(config.fs.north_home, user.user_id)
    if not os.path.exists(user_home):
        os.makedirs(user_home)

    def truncate(path_name):
        # On Linux: The maximum length for a file name is 255 bytes
        return path_name[:230]

    user_id = str(user.user_id)
    upload_query = Q()
    upload_query &= Q(main_author=user_id) | Q(coauthors=user_id)
    upload_query &= Q(publish_time=None)

    uploads: List[Dict] = []
    for upload in Upload.objects.filter(upload_query):

        if not hasattr(upload.upload_files, 'external_os_path'):
            # In case the files are missing for one reason or another
            logger.info('upload: the external path is missing for one reason or another')
            continue

        if upload.upload_name:
            upload_dir = f'uploads/{truncate(slugify(upload.upload_name))}-{upload.upload_id}'
        else:
            upload_dir = f'uploads/{upload.upload_id}'

        uploads.append(
            {
                'host_path': os.path.join(upload.upload_files.external_os_path, 'raw'),
                'mount_path': os.path.join(tool.mount_path, upload_dir)
            }
        )

    # Check if the tool/named server already exists
    _get_status(tool, user)
    if tool.state != ToolStateEnum.stopped:
        return ToolResponseModel(
            tool=tool.name,
            username=user.username,
            data=_get_status(tool, user))

    url = f'{config.hub_url()}/api/users/{user.username}/servers/{tool.name}'
    body = {
        'tool': {
            'image': tool.image,
            'cmd': tool.cmd,
            'privileged': tool.privileged
        },
        'environment': {
            'SUBFOLDER': f'{config.services.api_base_path.rstrip("/")}/north/user/{user.username}/',
            'JUPYTERHUB_CLIENT_API_URL': f'{config.north_url()}/hub/api'
        },
        'user_home': {
            'host_path': os.path.join(config.fs.north_home_external, user.user_id),
            'mount_path': os.path.join(tool.mount_path, 'work')
        },
        'uploads': uploads
    }

    logger.info('body of the post call', body=body)

    response = requests.post(url, json=body, headers=hub_api_headers)

    if response.status_code == 400 and 'is already running' in response.json()['message']:
        tool.state = ToolStateEnum.running
    elif response.status_code == 201:
        tool.state = ToolStateEnum.running
    elif response.status_code == 202:
        tool.state = ToolStateEnum.starting
    else:
        logger.error(
            'unexpected jupyterhub response', data=dict(status_code=response.status_code),
            text=response.text)
        tool.state = ToolStateEnum.stopped

    return ToolResponseModel(
        tool=tool.name,
        username=user.username,
        data=_get_status(tool, user))


@router.delete(
    '/{name}', tags=[default_tag],
    response_model=ToolResponseModel,
    summary='Stop a tool.',
    response_model_exclude_unset=True,
    response_model_exclude_none=True
)
async def stop_tool(
    tool: ToolModel = Depends(tool),
    user: User = Depends(create_user_dependency(required=True))
):
    url = f'{config.hub_url()}/api/users/{user.username}/servers/{tool.name}'
    response = requests.delete(url, json={'remove': True}, headers=hub_api_headers)

    if response.status_code == 404:
        tool.state = ToolStateEnum.stopped
    elif response.status_code == 204:
        tool.state = ToolStateEnum.stopped
    elif response.status_code == 202:
        tool.state = ToolStateEnum.stopping
    else:
        logger.error(
            'unexpected jupyterhub response', data=dict(status_code=response.status_code),
            text=response.text)
        tool.state = ToolStateEnum.stopped

    return ToolResponseModel(
        tool=tool.name,
        username=user.username,
        data=_get_status(tool, user))
