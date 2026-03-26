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


import datetime
import time
from unittest.mock import MagicMock

import pytest
from bson.objectid import ObjectId

from nomad.auth.scopes import Scope
from nomad.auth.tokens import (
    PAT_PREFIX,
    PATMetadata,
    PATQuery,
    _hash_token,
    authenticate_pat,
    create_pat,
    get_pat,
    list_pat,
    revoke_pat,
    rotate_pat,
)
from nomad.common import now
from nomad.config import config
from nomad.mongo.pat import PAT

# Test `create`


def test_create_returns_correct_structure(mongo_function):
    """
    Verifies that create() returns the NamedTuple with a raw token
    and a saved DB object.
    """
    result = create_pat(
        user_id='user_123',
        expires_in_days=30,
        metadata=PATMetadata(
            name='CI/CD Token',
            scopes=['read', 'write'],
        ),
    )

    # Check return type
    assert result.raw_token is not None
    assert isinstance(result.pat, PAT)

    # Check raw token format
    assert result.raw_token.startswith(PAT_PREFIX)
    assert len(result.raw_token) > 40  # Prefix + 32 bytes (b64 encoded)

    # Check DB object
    saved_pat = PAT.objects.get(id=result.pat.id)
    assert saved_pat.user_id == 'user_123'
    assert saved_pat.name == 'CI/CD Token'
    assert saved_pat.scopes == ['read', 'write']


def test_create_hashing_security(mongo_function):
    """
    Verify that we NEVER store the raw token, only the hash.
    """
    result = create_pat(
        user_id='u1',
        metadata=PATMetadata(name='Security Test', scopes=[]),
        expires_in_days=1,
    )

    raw = result.raw_token
    stored_digest = result.pat.token_digest

    # Raw token should NOT be equal to stored digest
    assert raw != stored_digest

    # Re-hashing the raw token should match the stored digest
    assert _hash_token(raw) == stored_digest


def test_create_expiration_logic(mongo_function):
    """Test standard expiration math."""
    days = 10
    result = create_pat(
        user_id='u1',
        metadata=PATMetadata(name='Exp Test', scopes=[]),
        expires_in_days=days,
    )

    expected_date = now() + datetime.timedelta(days=days)
    # Check if dates are close (within 10 seconds)
    delta = abs((result.pat.expired_at - expected_date).total_seconds())
    assert delta < 10


@pytest.mark.parametrize('expires_in_days', [0, -1])
def test_create_invalid_lifespan(mongo_function, expires_in_days):
    """
    Verifies that the service rejects requests to create tokens
    with a lifespan of 0 or negative days.
    """
    with pytest.raises(ValueError, match='already expired'):
        create_pat(
            user_id='u_invalid_life',
            metadata=PATMetadata(name='incorrect expiration', scopes=[]),
            expires_in_days=expires_in_days,
        )


@pytest.mark.parametrize(
    'requested_days, expected_match',
    [
        (None, 'Infinite tokens are disabled'),
        (60, 'exceeds the maximum allowed'),
    ],
)
def test_create_exceeds_configurable_lifetime(
    mongo_function, monkeypatch, requested_days, expected_match
):
    """
    Verifies that the service rejects requests that violate the
    configured maximum token lifetime (e.g., requesting infinite or exceeding max).
    """
    monkeypatch.setattr(config.auth, 'pat_max_lifetime', 30)

    with pytest.raises(ValueError, match=expected_match):
        create_pat(
            user_id='u_exceeds_life',
            metadata=PATMetadata(name='configurable expiration', scopes=[]),
            expires_in_days=requested_days,
        )


def test_create_within_configurable_lifetime(mongo_function, monkeypatch):
    """
    Verifies that the service successfully creates a token when
    the requested lifetime is within the configured limits.
    """
    monkeypatch.setattr(config.auth, 'pat_max_lifetime', 30)

    result = create_pat(
        user_id='u_valid_life',
        metadata=PATMetadata(name='valid expiration', scopes=[]),
        expires_in_days=15,
    )

    pat_obj, raw_secret = result
    assert pat_obj is not None
    assert raw_secret is not None

    # Test that we can create non-expiring tokens
    monkeypatch.setattr(config.auth, 'pat_max_lifetime', None)
    result = create_pat(
        user_id='u1',
        metadata=PATMetadata(name='Forever Token', scopes=[]),
        expires_in_days=None,
    )

    assert result.pat.expired_at is None


def test_create_exceeds_max_per_user(mongo_function, monkeypatch):
    """
    Verifies that the service raises a ValueError when a user attempts
    to create more tokens than the configured maximum limit.
    """
    # Restrict the maximum number of ative tokens per user to 2
    monkeypatch.setattr(config.auth, 'pat_max_active_per_user', 2)
    user_id = 'u_limit_test'

    # Create the first token (should succeed)
    result_1 = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Token 1', scopes=[]),
        expires_in_days=30,
    )
    assert result_1.pat is not None

    # Create the second token (should succeed - now at the limit)
    result_2 = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Token 2', scopes=[]),
        expires_in_days=30,
    )
    assert result_2.pat is not None

    # Attempt to create a third token
    with pytest.raises(ValueError, match='Maximum number of active'):
        create_pat(
            user_id=user_id,
            metadata=PATMetadata(name='Token 3', scopes=[]),
            expires_in_days=30,
        )


@pytest.mark.parametrize(
    'scope', [Scope.TOKENS_CREATE, Scope.TOKENS_DELETE, Scope.TOKENS_READ]
)
def test_create_rejects_token_operating_scope(mongo_function, scope):
    """
    Verifies that requesting any tokens scope immediately throws a ValueError.
    """
    with pytest.raises(
        ValueError, match='Personal access tokens are not allowed to operate on PATs'
    ):
        create_pat(
            user_id='user_123',
            expires_in_days=30,
            metadata=PATMetadata(
                name='Malicious Token',
                scopes=[
                    scope,
                ],
            ),
        )


# Test `rotate`


def test_rotate_basic_success(mongo_function):
    """
    Test the sucess path:
    1. Old token is revoked.
    2. New token is created with a NEW secret.
    3. Metadata (name, scopes) is copied.
    """
    user_id = 'u_rotate_happy'

    # Create original
    original_res = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='CI Token', scopes=['read']),
        expires_in_days=30,
    )
    original_id = original_res.pat.id
    original_secret = original_res.raw_token

    # Rotate
    new_res = rotate_pat(user_id=user_id, pat_id=str(original_id))

    # Verify old token is revoked
    old_db = PAT.objects.get(id=original_id)
    assert old_db.revoked is True

    # Verify new token
    assert new_res.pat.id != original_id
    assert new_res.raw_token != original_secret
    assert new_res.pat.name == 'CI Token'
    assert new_res.pat.scopes == ['read']
    assert new_res.pat.is_active is True


def test_rotate_preserves_original_lifespan(mongo_function):
    """
    If one creates a 1-Year token and rotate it 6 months later,
    the new token should be valid for 1 Year.
    """
    user_id = 'u_lifespan'
    lifespan_days = 365

    # Create the token
    original_res = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Yearly', scopes=[]),
        expires_in_days=lifespan_days,
    )
    original = original_res.pat

    # Time Travel: Move BOTH dates back by 180 days
    # This simulates a token that was created 6 months ago
    # and was intended to expire 6 months from now.
    time_shift = datetime.timedelta(days=180)

    original.created_at -= time_shift
    original.expired_at -= time_shift
    original.save()

    # Ensure the gap is still 365 days before we rotate
    initial_duration = (original.expired_at - original.created_at).days
    assert abs(initial_duration - lifespan_days) <= 1

    # Rotate
    new_res = rotate_pat(user_id=user_id, pat_id=str(original.id))

    # Check duration of NEW token
    # It should be a fresh 365 days from now
    new_duration = new_res.pat.expired_at - new_res.pat.created_at
    assert abs(new_duration.days - lifespan_days) <= 1


def test_rotate_infinite_token(mongo_function, monkeypatch):
    """
    Test that rotating a token with NO expiration (None)
    creates a new token with NO expiration.
    """
    monkeypatch.setattr(config.auth, 'pat_max_lifetime', None)

    user_id = 'u_infinite'
    original = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Forever', scopes=[]),
        expires_in_days=None,
    ).pat

    new_res = rotate_pat(user_id=user_id, pat_id=str(original.id))

    original.reload()
    assert original.expired_at is None
    assert new_res.pat.expired_at is None


def test_rotate_wrong_user(mongo_function):
    """
    User A cannot rotate User B's token.
    """
    victim_id = 'u_victim'
    attacker_id = 'u_attacker'

    original = create_pat(
        user_id=victim_id,
        metadata=PATMetadata(name='Secret', scopes=[]),
        expires_in_days=30,
    ).pat

    # Attacker tries to rotate
    result = rotate_pat(user_id=attacker_id, pat_id=str(original.id))

    # Should fail
    assert result is None

    # Original should NOT be revoked
    original.reload()
    assert original.revoked is False


def test_rotate_non_existent_invalid(mongo_function):
    """Sanity check for non-existent or invalid PAT ID."""
    assert rotate_pat(user_id='user_1', pat_id=str(ObjectId())) is None
    assert rotate_pat(user_id='user_1', pat_id='invalid') is None


def test_rotate_revoked(mongo_function):
    """Shouldn't be able to rotate revoked token."""
    user_id = 'u_owner'
    # Create a token
    result = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='To Revoke', scopes=[]),
        expires_in_days=30,
    )
    token_id = result.pat.id

    # Revoke it
    success = revoke_pat(user_id=user_id, pat_id=str(token_id))
    assert PAT.objects.get(id=token_id).revoked is True

    with pytest.raises(ValueError, match='Cannot rotate an expired/revoked token'):
        rotate_pat(user_id=user_id, pat_id=token_id)


def test_rotate_expired(mongo_function):
    """Shouldn't be able to rotate expired token."""
    user_id = 'u_owner'

    result = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='To Expire', scopes=[]),
        expires_in_days=30,
    )
    token_id = result.pat.id

    # Manually expire the token
    expired_time = now() - datetime.timedelta(minutes=1)
    PAT.objects(id=token_id).update(set__expired_at=expired_time)

    token = PAT.objects.get(id=token_id)
    assert token.is_active is False

    current_time = now()
    if token.expired_at.tzinfo is None and current_time.tzinfo is not None:
        current_time = current_time.replace(tzinfo=None)
    assert token.expired_at < current_time

    # Attempt rotation
    with pytest.raises(ValueError, match='Cannot rotate an expired/revoked token'):
        rotate_pat(user_id=user_id, pat_id=token_id)


def test_rotate_at_max_limit(mongo_function, monkeypatch):
    """
    Verifies that rotating a PAT succeeds when the user is at their
    maximum active token limit.
    """
    monkeypatch.setattr(config.auth, 'pat_max_active_per_user', 1)
    user_id = 'u_rotate_limit_test'

    create_result = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Token to Rotate', scopes=[]),
        expires_in_days=30,
    )
    assert create_result.pat is not None

    rotate_pat(user_id=user_id, pat_id=create_result.pat.id)


# Test `list`


def test_list_isolation(mongo_function):
    """
    Test that we only see tokens for the requested user.
    """
    # Create tokens for User A
    create_pat(
        user_id='user_A', metadata=PATMetadata(name='A1', scopes=[]), expires_in_days=30
    )
    create_pat(
        user_id='user_A', metadata=PATMetadata(name='A2', scopes=[]), expires_in_days=30
    )

    # Create token for User B
    create_pat(
        user_id='user_B', metadata=PATMetadata(name='A1', scopes=[]), expires_in_days=30
    )

    # List user A
    tokens_A = list_pat(user_id='user_A')
    assert tokens_A.total == 2
    assert all(t.user_id == 'user_A' for t in tokens_A.data)
    # Ensure digest is dropped by the DB query
    assert all(t.token_digest is None for t in tokens_A.data)

    # List user B
    tokens_B = list_pat(user_id='user_B')
    assert tokens_B.total == 1
    assert tokens_B.data[0].name == 'A1'
    assert tokens_B.data[0].token_digest is None


def test_list_pagination(mongo_function):
    """
    Test that start and limit correctly slice the results while maintaining the total count.
    """
    user_id = 'u_page_test'

    # Create 5 tokens with predictable names
    for i in range(5):
        create_pat(
            user_id=user_id,
            metadata=PATMetadata(name=f'Token_{i}', scopes=[]),
            expires_in_days=30,
        )

    # Sort by name_asc so the order is perfectly predictable: 0, 1, 2, 3, 4

    # Page 1: Items 0, 1
    page_1 = list_pat(user_id=user_id, start=0, limit=2, order_by='name_asc')
    assert page_1.total == 5
    assert len(page_1.data) == 2
    assert [t.name for t in page_1.data] == ['Token_0', 'Token_1']

    # Page 2: Items 2, 3
    page_2 = list_pat(user_id=user_id, start=2, limit=2, order_by='name_asc')
    assert page_2.total == 5
    assert len(page_2.data) == 2
    assert [t.name for t in page_2.data] == ['Token_2', 'Token_3']

    # Page 3: Item 4 (Partial page)
    page_3 = list_pat(user_id=user_id, start=4, limit=2, order_by='name_asc')
    assert page_3.total == 5
    assert len(page_3.data) == 1
    assert page_3.data[0].name == 'Token_4'

    # Page 4: Out of bounds
    page_empty = list_pat(user_id=user_id, start=10, limit=2)
    assert page_empty.total == 5
    assert len(page_empty.data) == 0


def test_list_pat_filter_search(mongo_function):
    """
    Test filtering by token name search.
    """
    user_id = 'u_test'

    create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Active Token', scopes=[]),
        expires_in_days=30,
    )

    create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Other', scopes=[]),
        expires_in_days=30,
    )

    results = list_pat(user_id=user_id, query=PATQuery(search='active'))

    assert results.total == 1
    assert results.data[0].name == 'Active Token'


def test_list_pat_filter_revoked(mongo_function):
    """
    Test filtering by revoked status.
    """
    user_id = 'u_test'

    t_revoked = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Revoked Token', scopes=[]),
        expires_in_days=30,
    ).pat
    t_revoked.revoked = True
    t_revoked.save()

    create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Valid Token', scopes=[]),
        expires_in_days=30,
    )

    results = list_pat(user_id=user_id, query=PATQuery(revoked=True))

    assert results.total == 1
    assert results.data[0].name == 'Revoked Token'


def test_list_pat_filter_state_active(mongo_function):
    """
    Test filtering by active state.
    """
    user_id = 'u_test'
    current_time = now()

    create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Active Token', scopes=[]),
        expires_in_days=30,
    )

    t_expired = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Expired Token', scopes=[]),
        expires_in_days=30,
    ).pat
    t_expired.expired_at = current_time - datetime.timedelta(days=1)
    t_expired.save()

    results = list_pat(user_id=user_id, query=PATQuery(state='active'))

    assert results.total == 1
    assert results.data[0].name == 'Active Token'


def test_list_pat_filter_state_inactive_expired(mongo_function):
    """
    Test that inactive state matches expired tokens.
    """
    user_id = 'u_test'
    current_time = now()

    t_expired = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Expired Token', scopes=[]),
        expires_in_days=30,
    ).pat
    t_expired.expired_at = current_time - datetime.timedelta(days=1)
    t_expired.save()

    create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Active Token', scopes=[]),
        expires_in_days=30,
    )

    results = list_pat(user_id=user_id, query=PATQuery(state='inactive'))

    assert results.total == 1
    assert results.data[0].name == 'Expired Token'


def test_list_pat_filter_state_inactive_revoked(mongo_function):
    """
    Test that inactive state matches revoked tokens.
    """
    user_id = 'u_test'

    t_revoked = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Revoked Token', scopes=[]),
        expires_in_days=30,
    ).pat
    t_revoked.revoked = True
    t_revoked.save()

    create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Active Token', scopes=[]),
        expires_in_days=30,
    )

    results = list_pat(user_id=user_id, query=PATQuery(state='inactive'))

    assert results.total == 1
    assert results.data[0].name == 'Revoked Token'


@pytest.mark.parametrize(
    ('query_field', 'matching_name'),
    [
        pytest.param('created_before', 'Old Token', id='created_before'),
        pytest.param('created_after', 'New Token', id='created_after'),
    ],
)
def test_list_pat_filter_created_bounds(
    mongo_function, query_field: str, matching_name: str
):
    """
    Test filtering by created_at lower/upper bounds.
    """
    user_id = 'u_test'
    current_time = now()

    t_old = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Old Token', scopes=[]),
        expires_in_days=30,
    ).pat
    t_old.created_at = current_time - datetime.timedelta(days=10)
    t_old.save()

    t_new = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='New Token', scopes=[]),
        expires_in_days=30,
    ).pat
    t_new.created_at = current_time - datetime.timedelta(days=1)
    t_new.save()

    cutoff = current_time - datetime.timedelta(days=5)
    results = list_pat(
        user_id=user_id,
        query=PATQuery(**{query_field: cutoff}),
    )

    assert results.total == 1
    assert results.data[0].name == matching_name


@pytest.mark.parametrize(
    ('query_field', 'matching_name'),
    [
        pytest.param('last_used_before', 'Old Used Token', id='last_used_before'),
        pytest.param('last_used_after', 'Recent Used Token', id='last_used_after'),
    ],
)
def test_list_pat_filter_last_used_bounds(
    mongo_function, query_field: str, matching_name: str
):
    """
    Test filtering by last_used_at lower/upper bounds.
    """
    user_id = 'u_test'
    current_time = now()

    t_old_used = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Old Used Token', scopes=[]),
        expires_in_days=30,
    ).pat
    t_old_used.last_used_at = current_time - datetime.timedelta(days=10)
    t_old_used.save()

    t_recent_used = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Recent Used Token', scopes=[]),
        expires_in_days=30,
    ).pat
    t_recent_used.last_used_at = current_time - datetime.timedelta(days=1)
    t_recent_used.save()

    cutoff = current_time - datetime.timedelta(days=5)
    results = list_pat(
        user_id=user_id,
        query=PATQuery(**{query_field: cutoff}),
    )

    assert results.total == 1
    assert results.data[0].name == matching_name


@pytest.mark.parametrize(
    ('query_field', 'matching_name'),
    [
        pytest.param('expires_before', 'Expired Token', id='expires_before'),
        pytest.param('expires_after', 'Valid Token', id='expires_after'),
    ],
)
def test_list_pat_filter_expires_bounds(
    mongo_function, query_field: str, matching_name: str
):
    """
    Test filtering by expired_at lower/upper bounds.
    """
    user_id = 'u_test'
    current_time = now()

    t_expired = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Expired Token', scopes=[]),
        expires_in_days=30,
    ).pat
    t_expired.expired_at = current_time - datetime.timedelta(days=1)
    t_expired.save()

    t_valid = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Valid Token', scopes=[]),
        expires_in_days=30,
    ).pat
    t_valid.expired_at = current_time + datetime.timedelta(days=10)
    t_valid.save()

    results = list_pat(
        user_id=user_id,
        query=PATQuery(**{query_field: current_time}),
    )

    assert results.total == 1
    assert results.data[0].name == matching_name


@pytest.mark.parametrize(
    ('order_by', 'expected_names'),
    [
        pytest.param('created_asc', ['Old', 'New'], id='created_asc'),
        pytest.param('created_desc', ['New', 'Old'], id='created_desc'),
    ],
)
def test_list_pat_order_created(mongo_function, order_by, expected_names: list[str]):
    """
    Test ordering by created_at.
    """
    user_id = 'u_order_created'
    current_time = now()

    t_old = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Old', scopes=[]),
        expires_in_days=30,
    ).pat
    t_old.created_at = current_time - datetime.timedelta(days=1)
    t_old.save()

    t_new = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='New', scopes=[]),
        expires_in_days=30,
    ).pat
    t_new.created_at = current_time
    t_new.save()

    results = list_pat(user_id=user_id, order_by=order_by)

    assert results.total == 2
    assert [token.name for token in results.data] == expected_names


@pytest.mark.parametrize(
    ('order_by', 'expected_names'),
    [
        pytest.param('name_asc', ['Apple', 'Zebra'], id='name_asc'),
        pytest.param('name_desc', ['Zebra', 'Apple'], id='name_desc'),
    ],
)
def test_list_pat_order_name(mongo_function, order_by, expected_names: list[str]):
    """
    Test ordering by name.
    """
    user_id = 'u_order_name'

    create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Zebra', scopes=[]),
        expires_in_days=30,
    )

    create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Apple', scopes=[]),
        expires_in_days=30,
    )

    results = list_pat(user_id=user_id, order_by=order_by)

    assert results.total == 2
    assert [token.name for token in results.data] == expected_names


@pytest.mark.parametrize(
    ('order_by', 'expected_names'),
    [
        pytest.param('last_used_asc', ['Old', 'New'], id='last_used_asc'),
        pytest.param('last_used_desc', ['New', 'Old'], id='last_used_desc'),
    ],
)
def test_list_pat_order_last_used(mongo_function, order_by, expected_names: list[str]):
    """
    Test ordering by last_used_at.
    """
    user_id = 'u_order_last_used'
    current_time = now()

    t_old = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Old', scopes=[]),
        expires_in_days=30,
    ).pat
    t_old.last_used_at = current_time - datetime.timedelta(days=1)
    t_old.save()

    t_new = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='New', scopes=[]),
        expires_in_days=30,
    ).pat
    t_new.last_used_at = current_time
    t_new.save()

    results = list_pat(user_id=user_id, order_by=order_by)

    assert results.total == 2
    assert [token.name for token in results.data] == expected_names


@pytest.mark.parametrize(
    ('order_by', 'expected_names'),
    [
        pytest.param('expires_asc', ['Soon', 'Later'], id='expires_asc'),
        pytest.param('expires_desc', ['Later', 'Soon'], id='expires_desc'),
    ],
)
def test_list_pat_order_expires(mongo_function, order_by, expected_names: list[str]):
    """
    Test ordering by expired_at.
    """
    user_id = 'u_order_expires'
    current_time = now()

    t_soon = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Soon', scopes=[]),
        expires_in_days=30,
    ).pat
    t_soon.expired_at = current_time + datetime.timedelta(days=1)
    t_soon.save()

    t_later = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Later', scopes=[]),
        expires_in_days=30,
    ).pat
    t_later.expired_at = current_time + datetime.timedelta(days=10)
    t_later.save()

    results = list_pat(user_id=user_id, order_by=order_by)

    assert results.total == 2
    assert [token.name for token in results.data] == expected_names


def test_list_empty(mongo_function):
    """Test that a user with no tokens gets an empty list, not None."""
    tokens = list_pat(user_id='ghost_user')
    assert isinstance(tokens.data, list)
    assert tokens.total == 0


# Test `get`


def test_get_own_token_success(mongo_function):
    """
    Test that a user can successfully retrieve their own token.
    """
    user_id = 'u_owner'

    # Create token
    created = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='My Token', scopes=[]),
        expires_in_days=30,
    )
    pat_id = str(created.pat.id)

    # Retrieve it
    retrieved = get_pat(user_id=user_id, pat_id=pat_id)

    # Verify
    assert retrieved is not None
    assert retrieved.id == created.pat.id
    assert retrieved.name == 'My Token'
    assert retrieved.user_id == user_id

    # Digest should be dropped by the DB query
    assert retrieved.token_digest is None


def test_get_others_token_fails(mongo_function):
    """
    Ensure User A cannot fetch User B's token.
    """
    victim_id = 'u_victim'
    attacker_id = 'u_attacker'

    # Victim creates a token
    victim_token = create_pat(
        user_id=victim_id,
        metadata=PATMetadata(name='Secret Token', scopes=[]),
        expires_in_days=30,
    ).pat
    pat_id = str(victim_token.id)

    # Attacker tries to get it
    result = get_pat(user_id=attacker_id, pat_id=pat_id)

    assert result is None


def test_get_expired_token(mongo_function):
    """
    Test that an expired token can still be retrieved directly.
    """
    user_id = 'u_expired_get_test'

    created = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Expired Token', scopes=[]),
        expires_in_days=30,
    )
    pat = created.pat

    # Backdate the expiration date to simulate natural expiration
    pat.expired_at = pat.created_at - datetime.timedelta(days=5)
    pat.save()

    pat_id = str(pat.id)

    retrieved = get_pat(user_id=user_id, pat_id=pat_id)

    assert retrieved is not None
    assert retrieved.id == pat.id
    assert retrieved.name == 'Expired Token'
    assert retrieved.user_id == user_id


def test_get_revoked_token(mongo_function):
    """
    Test that can get metadata for a revoked token.
    """
    user_id = 'u_revoked_get_test'

    created = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Revoked Token', scopes=[]),
        expires_in_days=30,
    )
    pat_id = str(created.pat.id)

    # Revoke it
    revoke_pat(user_id=user_id, pat_id=pat_id)

    # Retrieve it
    retrieved = get_pat(user_id=user_id, pat_id=pat_id)

    assert retrieved is not None
    assert str(retrieved.id) == pat_id
    assert retrieved.name == 'Revoked Token'
    assert retrieved.revoked is True
    assert retrieved.is_active is False


def test_get_non_existent_invalid_token(mongo_function):
    """
    Test retrieving a nonexistent or invalid ID returns None.
    """
    assert get_pat(user_id='any_user', pat_id=str(ObjectId())) is None
    assert get_pat(user_id='any_user', pat_id='123') is None


# Test `revoke`


def test_revoke_success(mongo_function):
    """
    Test that a user can successfully revoke their own token.

    Revocation should mark the token as revoked and record `revoked_at`.
    """
    user_id = 'u_owner'

    # Create a token
    result = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='To Revoke', scopes=[]),
        expires_in_days=30,
    )
    token_id = result.pat.id
    original_expires_at = result.pat.expired_at

    # Revoke it
    success = revoke_pat(user_id=user_id, pat_id=str(token_id))
    assert success is True

    # Verify database state
    updated_pat = PAT.objects.get(id=token_id)
    assert updated_pat.revoked is True
    assert updated_pat.revoked_at is not None

    # Strip timezones for comparison
    revoked_at_naive = updated_pat.revoked_at.replace(tzinfo=None)
    updated_at_naive = updated_pat.updated_at.replace(tzinfo=None)

    original_expired_naive = original_expires_at.replace(tzinfo=None)
    updated_expired_naive = updated_pat.expired_at.replace(tzinfo=None)
    expire_diff = abs((updated_expired_naive - original_expired_naive).total_seconds())
    assert expire_diff < 0.1  # allow 100 ms drift

    # The revocation time should be roughly equal to `updated_at`
    time_diff = abs((revoked_at_naive - updated_at_naive).total_seconds())
    assert time_diff < 5  # seconds


def test_revoke_wrong_user_security(mongo_function):
    """
    Test that User A cannot revoke User B's token.
    """
    owner_id = 'u_victim'
    attacker_id = 'u_attacker'

    # Create token for victim
    result = create_pat(
        user_id=owner_id,
        metadata=PATMetadata(name='Victim Token', scopes=[]),
        expires_in_days=30,
    )
    token_id = result.pat.id

    # Attacker tries to revoke it
    success = revoke_pat(user_id=attacker_id, pat_id=str(token_id))

    # Should fail
    assert success is False

    # Token should still be active
    updated_pat = PAT.objects.get(id=token_id)
    assert updated_pat.revoked is False


def test_revoke_non_existent_invalid(mongo_function):
    """
    Test revoking a non-existent or invalid ID returns False.
    """
    assert revoke_pat(user_id='user_1', pat_id=str(ObjectId())) is False
    assert revoke_pat(user_id='user_1', pat_id='invalid') is False


def test_revoke_already_revoked(mongo_function, monkeypatch):
    """
    Test that revoking an already revoked token is safe.
    It should return True (because the end state is 'revoked=True')
    and short-circuit without touching the database.
    """
    user_id = 'u_double'
    result = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Double Tap', scopes=[]),
        expires_in_days=30,
    )
    token_id = str(result.pat.id)

    # First revoke
    assert revoke_pat(user_id=user_id, pat_id=token_id) is True

    result.pat.reload()
    first_updated_at = result.pat.updated_at
    first_revoked_at = result.pat.revoked_at
    first_expired_at = result.pat.expired_at

    assert result.pat.revoked is True
    assert first_revoked_at is not None

    mock_save = MagicMock()
    monkeypatch.setattr('nomad.mongo.pat.PAT.save', mock_save)

    # Second revoke (already revoked, should be a no-op)
    assert revoke_pat(user_id=user_id, pat_id=token_id) is True

    # Make sure DB state is not changed
    mock_save.assert_not_called()
    result.pat.reload()

    assert result.pat.revoked is True
    assert result.pat.updated_at == first_updated_at
    assert result.pat.revoked_at == first_revoked_at
    assert result.pat.expired_at == first_expired_at


# Test `authenticate`


def test_authenticate_success(mongo_function):
    """
    A valid raw token should return the PAT object.
    """
    user_id = 'u_auth_success'

    # Create a valid token
    result = create_pat(
        user_id=user_id,
        metadata=PATMetadata(name='Login Key', scopes=[]),
        expires_in_days=30,
    )
    raw_token = result.raw_token

    # Authenticate
    authenticated_pat = authenticate_pat(raw_token)

    # Verify
    assert authenticated_pat is not None
    assert authenticated_pat.id == result.pat.id
    assert authenticated_pat.user_id == user_id


def test_authenticate_updates_last_used(mongo_function):
    """
    Verify that authenticating successfully updates the 'last_used_at' timestamp.
    """
    # Create token
    result = create_pat(
        user_id='u_usage',
        metadata=PATMetadata(name='Usage Test', scopes=[]),
        expires_in_days=30,
    )

    # Ensure usage is initially None
    assert result.pat.last_used_at is None

    # Sleep briefly to ensure the timestamp will be different
    time.sleep(0.01)

    # Authenticate
    authenticate_pat(result.raw_token)

    # Reload from DB
    updated_pat = result.pat.reload()
    assert updated_pat.last_used_at is not None

    # mongo doesn't store TZ
    current_time = now().replace(tzinfo=None)

    # Verify it happened "just now"
    assert (current_time - updated_pat.last_used_at).total_seconds() < 1


def test_authenticate_fails_wrong_prefix(mongo_function):
    """
    Strings without prefix should fail immediately.
    """
    assert authenticate_pat('wrong-prefix-token') is None


def test_authenticate_fails_wrong_secret(mongo_function):
    """
    A token with correct prefix but wrong random bytes fails.
    """
    # Create a real token
    real_result = create_pat(
        user_id='u_hacker',
        metadata=PATMetadata(name='Real', scopes=[]),
        expires_in_days=30,
    )

    # Modify the secret slightly
    fake_token = real_result.raw_token[:-1] + 'X'

    # Attempt auth
    assert authenticate_pat(fake_token) is None


def test_authenticate_fails_revoked(mongo_function):
    """
    A revoked token is rejected.
    """
    # Create & Revoke
    result = create_pat(
        user_id='u_revoked',
        metadata=PATMetadata(name='Revoked', scopes=[]),
        expires_in_days=30,
    )
    revoke_pat(user_id='u_revoked', pat_id=str(result.pat.id))

    pat = authenticate_pat(result.raw_token)
    assert pat is None


def test_authenticate_fails_expired(mongo_function):
    """
    An expired token is rejected.
    """
    # Create token
    result = create_pat(
        user_id='u_expired',
        metadata=PATMetadata(name='Expired', scopes=[]),
        expires_in_days=30,
    )
    pat = result.pat

    # Manually expire it
    pat.expired_at = now() - datetime.timedelta(days=1)
    pat.save()

    assert authenticate_pat(result.raw_token) is None
