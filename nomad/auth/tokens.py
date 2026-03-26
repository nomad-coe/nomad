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

from __future__ import annotations

import datetime
import hashlib
import hmac
import secrets
import uuid
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal, NamedTuple, cast

from mongoengine import DoesNotExist, Q
from mongoengine.errors import ValidationError
from pydantic import BaseModel, ConfigDict, Field

from nomad import datamodel, utils
from nomad.auth import keycloak, user_management
from nomad.auth.keycloak import KeycloakError
from nomad.auth.scopes import Scope, _resolve_scopes
from nomad.common import now
from nomad.config import config
from nomad.config.models.config import _DEFAULT_API_KEY, ModeEnum
from nomad.datamodel import User
from nomad.mongo.pat import PAT

if TYPE_CHECKING:
    from typing import Final


@dataclass(frozen=True)
class AuthResult:
    user: User
    scopes: set[str]


# Keycloak token


def get_user_from_keycloak_token(keycloak_token: str | None) -> AuthResult | None:
    """
    Verifies keycloak bearer token.

    Returns:
        The corresponding AuthResult object,
        or None if cannot resolve.
    """
    from fastapi import HTTPException, status

    if keycloak_token is None or keycloak_token.startswith(PAT_PREFIX):
        return None

    try:
        user = cast(datamodel.User, keycloak.keycloak.tokenauth(keycloak_token))

    except KeycloakError as e:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail=str(e),
            headers={'WWW-Authenticate': 'Bearer'},
        )

    return AuthResult(user, _resolve_scopes(['*:*']))


# Personal access token (PAT) endpoint models

PAT_PREFIX: Final[str] = 'nomad_pat_'


class PATCreationResult(NamedTuple):
    """PAT creation return type, added to avoid using a regular tuple
    where one could easily get confused by accessing via index.
    """

    pat: PAT
    raw_token: str


class PATMetadata(BaseModel):
    """
    Holds all user-defined, copyable properties of a token.
    """

    name: str
    scopes: list[str]
    description: str | None = None

    class Config:
        from_attributes = True


class PATQuery(BaseModel):
    """Query parameters for filtering PATs."""

    search: str | None = Field(
        None, description='Search by token name (case-insensitive)'
    )
    revoked: bool | None = Field(
        None, description='Filter by explicitly revoked status'
    )
    state: Literal['active', 'inactive'] | None = Field(
        None, description='Filter by active/inactive state'
    )

    created_after: datetime.datetime | None = None
    created_before: datetime.datetime | None = None

    last_used_after: datetime.datetime | None = None
    last_used_before: datetime.datetime | None = None

    expires_after: datetime.datetime | None = None
    expires_before: datetime.datetime | None = None


class PATQueryResult(BaseModel):
    data: list[PAT]
    total: int = Field(ge=0, description='Total number of matching tokens')

    model_config = ConfigDict(arbitrary_types_allowed=True)


# PAT Service layer


def _hash_token(raw_token: str) -> str:
    """Hashes the token using SHA256 for storage."""
    return hashlib.sha256(raw_token.encode('utf-8')).hexdigest()


def create_pat(
    *,
    user_id: str,
    metadata: PATMetadata,
    expires_in_days: int | None,
) -> PATCreationResult:
    """
    Creates a new PAT.
    WARNING: The raw secret is only available once!

    Returns: (The saved PAT object, The RAW secret string)
    """
    # Check lifetime
    if expires_in_days is not None and expires_in_days <= 0:
        raise ValueError('Cannot create an already expired token.')

    if config.auth.pat_max_lifetime is not None:
        if expires_in_days is None:
            raise ValueError(
                f'Infinite tokens are disabled. Max lifetime is {config.auth.pat_max_lifetime} days.'
            )

        if expires_in_days > config.auth.pat_max_lifetime:
            raise ValueError(
                f'Requested lifetime ({expires_in_days} days) exceeds the maximum allowed ({config.auth.pat_max_lifetime} days).'
            )

    # PAT shouldn't be able to operate on PATs
    if (
        Scope.TOKENS_CREATE in metadata.scopes
        or Scope.TOKENS_DELETE in metadata.scopes
        or Scope.TOKENS_READ in metadata.scopes
    ):
        raise ValueError('Personal access tokens are not allowed to operate on PATs.')

    # Check allowed number of active tokens per user
    if (
        list_pat(user_id=user_id, query=PATQuery(state='active')).total
        >= config.auth.pat_max_active_per_user
    ):
        raise ValueError(
            f'Maximum number of active PAT ({config.auth.pat_max_active_per_user}) reached.'
        )

    # Generates a random 32-byte URL-safe string (the secret)
    raw_token = f'{PAT_PREFIX}{secrets.token_urlsafe(32)}'
    token_digest = _hash_token(raw_token)

    current_time = now()

    expires_at: datetime.datetime | None = (
        current_time + datetime.timedelta(days=expires_in_days)
        if expires_in_days is not None
        else None
    )

    # Create the token
    pat = PAT(
        user_id=user_id,
        token_digest=token_digest,
        expired_at=expires_at,
        created_at=current_time,
        updated_at=current_time,
        **metadata.model_dump(),
    )
    pat.save()

    return PATCreationResult(pat=pat, raw_token=raw_token)


def rotate_pat(*, user_id: str, pat_id: str) -> PATCreationResult | None:
    """
    Rotate (regenerate) a token (not expired nor revoked).
    1. Finds the old token.
    2. Copies its Name and Scopes.
    3. Revokes the old token.
    4. Creates a new token with a fresh secret.
    """
    # Find the old token (MUST belong to user)
    try:
        old_pat = PAT.objects(id=pat_id, user_id=user_id).first()
    except ValidationError:
        return None

    if not old_pat:
        return None

    # Expired/revoked (inactive) tokens are not allowed to be rotated
    if not old_pat.is_active:
        raise ValueError('Cannot rotate an expired/revoked token.')

    # Calculate the original token lifetime
    if old_pat.expired_at is None:
        expires_in_days = None
    else:
        duration = old_pat.expired_at - old_pat.created_at
        expires_in_days = duration.days

    # Revoke the old token
    revoke_pat(user_id=user_id, pat_id=pat_id)

    extracted_metadata = PATMetadata.model_validate(old_pat)

    return create_pat(
        user_id=user_id,
        metadata=extracted_metadata,
        expires_in_days=expires_in_days,
    )


PATSortOrder = Literal[
    'created_asc',
    'created_desc',
    'expires_asc',
    'expires_desc',
    'last_used_asc',
    'last_used_desc',
    'name_asc',
    'name_desc',
]


def list_pat(
    *,
    user_id: str,
    query: PATQuery | None = None,
    start: int = 0,
    limit: int | None = None,
    order_by: PATSortOrder = 'created_desc',
) -> PATQueryResult:
    """
    Lists tokens for a user based on filters, with pagination.
    Returns a tuple containing the page of data and the total count of matches.
    """
    # Security: only get tokens for current user
    mongo_query = Q(user_id=user_id)

    if query is not None:
        # Search (partial string match on the token name)
        if query.search:
            mongo_query &= Q(name__icontains=query.search)

        # Revoked
        if query.revoked is not None:
            mongo_query &= Q(revoked=query.revoked)

        # State (active/inactive)
        if query.state == 'active':
            current_time = now()
            mongo_query &= Q(revoked=False)
            mongo_query &= Q(expired_at=None) | Q(expired_at__gt=current_time)
        elif query.state == 'inactive':  # revoked OR expired
            current_time = now()
            mongo_query &= Q(revoked=True) | Q(expired_at__lte=current_time)

        # Created
        if query.created_after:
            mongo_query &= Q(created_at__gte=query.created_after)
        if query.created_before:
            mongo_query &= Q(created_at__lte=query.created_before)

        # Last Used
        if query.last_used_after:
            mongo_query &= Q(last_used_at__gte=query.last_used_after)
        if query.last_used_before:
            mongo_query &= Q(last_used_at__lte=query.last_used_before)

        # Expiration
        if query.expires_after:
            mongo_query &= Q(expired_at__gte=query.expires_after)
        if query.expires_before:
            mongo_query &= Q(expired_at__lte=query.expires_before)

    # The base QuerySet
    base_qs = PAT.objects.filter(mongo_query).exclude('token_digest')

    # Apply sorting
    sort_mapping: dict[str, str] = {
        'created_asc': 'created_at',
        'created_desc': '-created_at',
        'expires_asc': 'expired_at',
        'expires_desc': '-expired_at',
        'last_used_asc': 'last_used_at',
        'last_used_desc': '-last_used_at',
        'name_asc': 'name',
        'name_desc': '-name',
    }
    base_qs = base_qs.order_by(sort_mapping[order_by])

    # Apply pagination slicing
    if limit is None:
        results = list(base_qs[start:])
    else:
        results = list(base_qs[start : start + limit])

    return PATQueryResult(data=results, total=base_qs.count())


def get_pat(*, pat_id: str, user_id: str) -> PAT | None:
    """
    Retrieves a token by ID.
    """
    try:
        return PAT.objects(id=pat_id, user_id=user_id).exclude('token_digest').first()
    except ValidationError:  # Invalid PAT ID format
        return None


def revoke_pat(*, user_id: str, pat_id: str) -> bool:
    """Revoke a token."""
    try:
        # Ensure user_id matches to prevent users revoking others' tokens
        pat = PAT.objects.get(id=pat_id, user_id=user_id)
    except (DoesNotExist, ValidationError):
        return False

    if pat.revoked:
        return True

    pat.revoked_at = now()
    pat.revoked = True
    pat.save()

    return True


def authenticate_pat(raw_token: str) -> PAT | None:
    """
    Validates a raw token string.
    1. Find its hash in database
    2. Check validity (revocation and expiry)
    3. Update usage (last_used_at)
    """
    if not raw_token.startswith(PAT_PREFIX):
        return None

    # Find the token
    pat = PAT.objects(token_digest=_hash_token(raw_token)).first()

    if not pat or not pat.is_active:
        return None

    # Update usage
    pat.update(set__last_used_at=now())
    return pat


# [Deprecated] Other NOMAD custom tokens (signature/simple/upload token)


JWT_ALGORITHM = 'HS256'
HMAC_DIGESTMOD = hashlib.sha256


def check_api_secret() -> None:
    if (
        config.services.mode == ModeEnum.PRODUCTION
        and config.services.api_secret == _DEFAULT_API_KEY
    ):
        raise ValueError(
            'When running NOMAD in production mode, value for config.services.api_secret must be set to a minimum 32 character string through the environment variable NOMAD_SERVICES_API_SECRET. '
            'Alternatively you can run NOMAD in an insecure development mode by setting config.services.mode to development.'
        )


def generate_simple_token(user_id: str, expires_in: float) -> str:
    """
    Generate a simple token: JWT encoded user_id and expiration time,
    signed with the API secret.
    """
    import jwt

    check_api_secret()
    expires_at = now() + datetime.timedelta(seconds=expires_in)
    payload = dict(user=user_id, exp=expires_at)
    return jwt.encode(
        payload=payload, key=config.services.api_secret, algorithm=JWT_ALGORITHM
    )


def get_user_from_simple_token(simple_token: str | None) -> AuthResult | None:
    """
    Verifies a simple token (throwing HTTPException if illegal value provided).

    Returns:
        The corresponding AuthResult object,
        or None if cannot resolve.
    """
    import jwt
    from fastapi import HTTPException, status

    if simple_token is None:
        return None

    check_api_secret()

    try:
        decoded = jwt.decode(
            simple_token, config.services.api_secret, algorithms=[JWT_ALGORITHM]
        )
        user = User.get(user_id=decoded['user'])
        scopes = _resolve_scopes(['*:*']) - _resolve_scopes(['tokens:*'])
        return AuthResult(user, scopes)

    except KeyError:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail='Token with invalid/unexpected payload.',
            headers={'WWW-Authenticate': 'Bearer'},
        )
    except jwt.ExpiredSignatureError:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail='Expired token.',
            headers={'WWW-Authenticate': 'Bearer'},
        )
    except jwt.InvalidTokenError:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail='Invalid token.',
            headers={'WWW-Authenticate': 'Bearer'},
        )


def generate_upload_token(user: User) -> str:
    """Generate an upload token for user."""
    check_api_secret()
    payload = uuid.UUID(user.user_id).bytes
    signature = hmac.new(
        config.services.api_secret.encode('utf-8'),
        msg=payload,
        digestmod=HMAC_DIGESTMOD,
    )

    return f'{utils.base64_encode(payload)}.{utils.base64_encode(signature.digest())}'


def get_user_from_upload_token(upload_token: str | None) -> AuthResult | None:
    """
    Verifies the upload token (throwing HTTPException if illegal value provided).

    Returns:
        The corresponding AuthResult object,
        or None if cannot resolve.
    """
    from fastapi import HTTPException, status

    if upload_token is None:
        return None

    check_api_secret()

    try:
        payload, signature = upload_token.split('.', 1)
        payload_bytes = utils.base64_decode(payload)
        signature_bytes = utils.base64_decode(signature)

        expected = hmac.new(
            config.services.api_secret.encode('utf-8'),
            msg=payload_bytes,
            digestmod=HMAC_DIGESTMOD,
        )

        if not hmac.compare_digest(signature_bytes, expected.digest()):
            raise ValueError('Invalid HMAC signature')

        user_id = str(uuid.UUID(bytes=payload_bytes))
        user = cast(datamodel.User, user_management.user_management.get_user(user_id))
        return AuthResult(user, _resolve_scopes(['uploads:*']))

    except Exception:
        # Decode error, format error, user not found, etc.
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail='An invalid upload token was supplied.',
            headers={'WWW-Authenticate': 'Bearer'},
        )
