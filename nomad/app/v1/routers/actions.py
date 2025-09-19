import asyncio
from typing import Final

from fastapi import APIRouter, Depends, HTTPException, status
from fastapi_cache.decorator import cache
from pydantic import BaseModel

from nomad.actions.utils import (
    ActionModel,
    ActionModelSummary,
    ActionSchemaInfo,
    get_action_result,
    get_action_status,
    get_all_action_schemas,
    get_all_user_actions,
    get_user_action,
    start_action,
    stop_action,
    validate_action_arg,
)
from nomad.app.v1.models import User
from nomad.app.v1.routers.auth import create_user_dependency
from nomad.utils import strip

from ..models import HTTPExceptionModel, User
from ..utils import create_responses

router = APIRouter()


class ActionStart(BaseModel):
    data: dict


SCHEMA_CACHE_TTL: Final[int] = 1 * 24 * 60 * 60  # 1 day in seconds


@router.post('/{action_id}/start')
async def action_start(
    action_id: str,
    start_data: ActionStart,
    user: User = Depends(create_user_dependency(required=True)),
):
    """
    Starts a new action.

    Args:
        action_id: The ID of the action to start.
        start_data: The input data for the action.
        user: The authenticated user.

    Returns:
        The ID of the started action instance.
    """
    start_data.data['user_id'] = user.user_id
    try:
        input_data = validate_action_arg(action_id, start_data.data)
        action_instance_id = await asyncio.to_thread(
            lambda: start_action(action_id=action_id, data=input_data)
        )
        return {'action_instance_id': action_instance_id}
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post('/{action_instance_id}/stop')
async def action_stop(
    action_instance_id: str,
    user: User = Depends(create_user_dependency(required=True)),
):
    """
    Stops an action.

    Args:
        action_instance_id: The ID of the action instance to stop.
        user: The authenticated user.
    """
    try:
        await asyncio.to_thread(
            lambda: stop_action(
                action_instance_id=action_instance_id, user_id=user.user_id
            )
        )
        return {'status': 'stopped'}
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get('/{action_instance_id}/status')
async def action_status(
    action_instance_id: str, user: User = Depends(create_user_dependency(required=True))
):
    """
    Gets the status of an action.

    Args:
        action_instance_id: The ID of the action instance.
        user: The authenticated user.

    Returns:
        The status of the action.
    """
    try:
        status = await asyncio.to_thread(
            lambda: get_action_status(
                action_instance_id=action_instance_id, user_id=user.user_id
            )
        )
        return {'status': status.name}
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get('/{action_instance_id}/result')
async def action_result(
    action_instance_id: str, user: User = Depends(create_user_dependency(required=True))
):
    """
    Gets the result of an action.

    Args:
        action_instance_id: The ID of the action instance.
        user: The authenticated user.

    Returns:
        The result of the action.
    """
    try:
        result = await asyncio.to_thread(
            lambda: get_action_result(
                action_instance_id=action_instance_id, user_id=user.user_id
            )
        )
        return result
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get(
    '/schemas',
    response_model=list[ActionSchemaInfo],
)
@cache(expire=SCHEMA_CACHE_TTL)
async def action_input_schemas(
    user: User = Depends(create_user_dependency(required=True)),
):
    """
    Gets the input schemas for all available actions.

    Args:
        user: The authenticated user.

    Returns:
        A list of action input schemas.
    """
    try:
        result = await asyncio.to_thread(lambda: get_all_action_schemas())
        return result
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


_not_authorized = (
    status.HTTP_401_UNAUTHORIZED,
    {
        'model': HTTPExceptionModel,
        'description': strip(
            """
        Unauthorized. Authorization is required, but no or bad authentication credentials provided."""
        ),
    },
)


@router.get(
    '/{action_instance_id}',
    summary='Get a specific action of the authenticated user.',
    response_model=ActionModel,
    responses=create_responses(_not_authorized),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
async def action(
    action_instance_id: str,
    user: User = Depends(create_user_dependency(required=True)),
):
    """
    Gets a specific action for the authenticated user.

    Args:
        action_instance_id: The ID of the action instance.
        user: The authenticated user.

    Returns:
        The action.
    """
    try:
        result = await asyncio.to_thread(
            lambda: get_user_action(
                action_instance_id=action_instance_id, user_id=user.user_id
            )
        )
        if result is None:
            raise HTTPException(status_code=404, detail='Action not found.')
        return result
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get(
    '',
    summary='List uploads of authenticated user.',
    response_model=list[ActionModelSummary],
    responses=create_responses(_not_authorized),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
async def actions(
    user: User = Depends(create_user_dependency(required=True)),
):
    """
    Lists all actions for the authenticated user.

    Args:
        user: The authenticated user.

    Returns:
        A list of actions.
    """
    try:
        result = await asyncio.to_thread(
            lambda: get_all_user_actions(user_id=user.user_id)
        )
        return result
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
