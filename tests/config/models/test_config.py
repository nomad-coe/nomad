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

import pytest

from nomad.auth.scopes import _resolve_scopes
from nomad.config.models.config import Auth

test_cases = [
    pytest.param(
        {},
        _resolve_scopes({'*:*'}),
        id='empty-dict-gives-all',
    ),
    pytest.param(
        {'include': ['*:read']},
        _resolve_scopes({'*:read'}),
        id='include-only',
    ),
    pytest.param(
        {'exclude': ['*:read']},
        _resolve_scopes({'*:*'}) - _resolve_scopes({'*:read'}),
        id='exclude-only',
    ),
    pytest.param(
        {'include': ['*:*'], 'exclude': ['tokens:*']},
        _resolve_scopes({'*:*'}) - _resolve_scopes({'tokens:*'}),
        id='include-and-exclude',
    ),
    pytest.param(
        {'include': ['*:*'], 'exclude': ['*:*']},
        set(),
        id='equal-include-exclude-all',
    ),
    pytest.param(
        {'include': ['tokens:*'], 'exclude': ['tokens:*']},
        set(),
        id='equal-include-exclude_specific',
    ),
    pytest.param(
        {'include': [], 'exclude': ['tokens:*']},
        set(),
        id='empty-include',
    ),
]


@pytest.mark.parametrize('scopes, expected', test_cases)
def test_unauthenticated_user_scopes_resolved(scopes, expected):
    """Test that includes and excludes are correctly resolved."""
    auth = Auth.model_validate({'unauthenticated_user_scopes': scopes})
    assert auth.unauthenticated_user_scopes_resolved == expected


@pytest.mark.parametrize('scopes, expected', test_cases)
def test_unauthorized_user_scopes_resolved(scopes, expected):
    """Test that includes and excludes are correctly resolved for unauthorized users."""
    auth = Auth.model_validate({'unauthorized_user_scopes': scopes})
    assert auth.unauthorized_user_scopes_resolved == expected
