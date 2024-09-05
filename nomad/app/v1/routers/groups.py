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

from typing import List, Optional, Set

from fastapi import APIRouter, Depends, HTTPException, Query, status
from pydantic import BaseModel, Field

from nomad.datamodel import User as UserDataModel
from nomad.groups import (
    UserGroup as MongoUserGroup,
)
from nomad.groups import (
    create_user_group as create_mongo_user_group,
)
from nomad.utils import strip

from ..models import User
from .auth import create_user_dependency

router = APIRouter()
default_tag = 'groups'

group_name_description = 'Name of the group.'
group_members_description = 'User ids of the group members.'


class UserGroupEdit(BaseModel):
    group_name: Optional[str] = Field(
        default=None,
        description=group_name_description,
        min_length=3,
        max_length=32,
        regex=r'^[a-zA-Z0-9][a-zA-Z0-9 \.\_\-]+[a-zA-Z0-9]$',
    )
    members: Optional[Set[str]] = Field(
        default=None, description=group_members_description
    )


class UserGroup(BaseModel):
    group_id: str = Field(description='Unique id of the group.')
    group_name: str = Field(
        default='Default Group Name', description=group_name_description
    )
    owner: str = Field(description='User id of the group owner.')
    members: List[str] = Field(
        default_factory=list, description=group_members_description
    )

    class Config:
        orm_mode = True


class UserGroups(BaseModel):
    data: List[UserGroup]


def get_mongo_user_group(group_id: str) -> MongoUserGroup:
    user_group = MongoUserGroup.objects(group_id=group_id).first()
    if user_group is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"User group '{group_id}' was not found.",
        )

    return user_group


def check_user_ids(user_ids):
    """Raise 404 NOT FOUND if user id is not in database."""
    for user_id in user_ids:
        try:
            UserDataModel.get(user_id=user_id)
        except KeyError:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"User '{user_id}' was not found.",
            )


def check_user_may_edit_user_group(user: User, user_group: MongoUserGroup):
    if user.is_admin or user.user_id == user_group.owner:
        return

    raise HTTPException(
        status_code=status.HTTP_401_UNAUTHORIZED,
        detail=strip(
            f"Not authorized to edit user group '{user_group.group_id}'."
            ' Only group owners and admins are allowed to edit a group.'
        ),
    )


@router.get(
    '',
    tags=[default_tag],
    summary='List user groups. Use at most one filter.',
    response_model=UserGroups,
)
async def get_user_groups(
    group_id: Optional[List[str]] = Query(
        None, description='Search groups by their full id.'
    ),
    user_id: Optional[str] = Query(
        None, description='Search groups by their owner or members ids.'
    ),
    search_terms: Optional[str] = Query(
        None, description='Search groups by parts of their name.'
    ),
):
    """Get data about user groups."""
    if sum(param is not None for param in (group_id, user_id, search_terms)) > 1:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='Only one of (group_id, user_id, search_terms) may be used at a time.',
        )

    if group_id is not None:
        user_groups = MongoUserGroup.get_by_ids(group_id)
    elif user_id is not None:
        user_groups = MongoUserGroup.get_by_user_id(user_id)
    elif search_terms is not None:
        user_groups = MongoUserGroup.get_by_search_terms(search_terms)
    else:
        user_groups = MongoUserGroup.objects
    data = [UserGroup.from_orm(user_group) for user_group in user_groups]
    return {'data': data}


@router.get(
    '/{group_id}',
    tags=[default_tag],
    summary='Get data about user group.',
    response_model=UserGroup,
)
async def get_user_group(group_id: str):
    """Get data about user group."""
    user_group = get_mongo_user_group(group_id)

    return user_group


@router.post(
    '',
    tags=[default_tag],
    status_code=status.HTTP_201_CREATED,
    summary='Create user group.',
    response_model=UserGroup,
)
async def create_user_group(
    user_group_edit: UserGroupEdit,
    user: User = Depends(create_user_dependency(required=True)),
):
    """Create user group."""
    user_group_dict = user_group_edit.dict(exclude_none=True)
    members = user_group_dict.get('members')
    if members is not None:
        check_user_ids(members)
    user_group_dict['owner'] = user.user_id

    user_group = create_mongo_user_group(**user_group_dict)
    return user_group


@router.post(
    '/{group_id}/edit',
    tags=[default_tag],
    summary='Update user group.',
    response_model=UserGroup,
)
async def update_user_group(
    group_id: str,
    user_group_edit: UserGroupEdit,
    user: User = Depends(create_user_dependency(required=True)),
):
    """Update user group."""
    user_group = get_mongo_user_group(group_id)
    check_user_may_edit_user_group(user, user_group)

    user_group_dict = user_group_edit.dict(exclude_none=True)
    members = user_group_dict.get('members')
    if members is not None:
        check_user_ids(members)

    user_group.update(**user_group_dict)
    user_group.save()
    user_group.reload()
    return user_group


@router.delete(
    '/{group_id}',
    tags=[default_tag],
    status_code=status.HTTP_204_NO_CONTENT,
    summary='Delete user group.',
)
async def delete_user_group(
    group_id: str, user: User = Depends(create_user_dependency(required=True))
):
    """Delete user group."""
    user_group = get_mongo_user_group(group_id)
    check_user_may_edit_user_group(user, user_group)

    user_group.delete()
