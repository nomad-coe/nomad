from fastapi import Depends, APIRouter, status

from nomad.app_fastapi.routers.auth import get_required_user
from nomad.app_fastapi.models import User, HTTPExceptionModel
from nomad.app_fastapi.utils import create_responses
from nomad.utils import strip

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
async def read_users_me(current_user: User = Depends(get_required_user)):
    return current_user
