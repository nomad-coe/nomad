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

from enum import Enum
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable


class Scope(str, Enum):
    """
    Backend authorization scopes.

    Each scope is represented as a colon-separated string of the form:
        `resource:action`

    where:
    - `resource` identifies a protected backend domain or API surface
      (e.g. `datasets`, `uploads`, `users`).
    - `action` identifies the permitted operation on that resource
      (e.g. `read`, `write`, `delete`, `run`).
    """

    # actions
    ACTIONS_READ = 'actions:read'
    ACTIONS_RUN = 'actions:run'

    # apps
    APPS_READ = 'apps:read'

    # Personal access tokens (PAT) and custom NOMAD tokens
    TOKENS_CREATE = 'tokens:create'
    TOKENS_READ = 'tokens:read'
    TOKENS_DELETE = 'tokens:delete'

    # datasets
    DATASETS_READ = 'datasets:read'
    DATASETS_WRITE = 'datasets:write'
    DATASETS_DELETE = 'datasets:delete'
    DATASETS_ASSIGN_DOI = 'datasets:assign_doi'

    # entries
    ENTRIES_READ = 'entries:read'
    ENTRIES_WRITE = 'entries:write'

    # federation
    FEDERATION_WRITE = 'federation:write'

    # graph
    GRAPH_READ = 'graph:read'

    # groups
    GROUPS_READ = 'groups:read'
    GROUPS_WRITE = 'groups:write'
    GROUPS_DELETE = 'groups:delete'

    # info
    INFO_READ = 'info:read'

    # materials
    MATERIALS_READ = 'materials:read'

    # metainfo
    METAINFO_READ = 'metainfo:read'

    # north
    NORTH_READ = 'north:read'
    NORTH_RUN = 'north:run'

    # schemas
    SCHEMAS_READ = 'schemas:read'

    # suggestions
    SUGGESTIONS_READ = 'suggestions:read'

    # systems
    SYSTEMS_READ = 'systems:read'

    # uploads
    UPLOADS_READ = 'uploads:read'
    UPLOADS_WRITE = 'uploads:write'
    UPLOADS_PUBLISH = 'uploads:publish'
    UPLOADS_PROCESS = 'uploads:process'
    UPLOADS_ASSIGN_DOI = 'uploads:assign_doi'

    # uploads bundle
    UPLOADS_BUNDLE_READ = 'uploads_bundle:read'
    UPLOADS_BUNDLE_WRITE = 'uploads_bundle:write'

    # users
    USERS_READ = 'users:read'
    USERS_INVITE = 'users:invite'

    # external apps
    EXTERNAL_OPTIMADE_READ = 'external_optimade:read'
    EXTERNAL_DCAT_READ = 'external_dcat:read'
    EXTERNAL_H5GROVE_READ = 'external_h5grove:read'

    @classmethod
    def all_values(cls) -> set[str]:
        """Return all concrete scope strings."""
        return {scope.value for scope in cls}


def _resolve_scopes(
    scopes: str | Iterable[str],
) -> set[str]:
    """Resolve and validate configured scopes, supporting '*' as wildcards.

    Supported wildcard forms:
      - '*:*'         -> all known scopes
      - 'resource:*'  -> all actions for a resource
      - '*:action'    -> all resources for an action

    Unsupported:
      - partial wildcards like 'upl*', '*load'
      - glob syntax beyond '*'

    Args:
        scopes: Scope strings from configuration.
    """
    if isinstance(scopes, str):
        scopes = {scopes}
    else:
        scopes = set(scopes)
    known_scopes: set[str] = Scope.all_values()
    resolved_scopes: set[str] = set()

    unknown_concrete: set[str] = set()
    unmatched_wildcards: set[str] = set()

    for raw in scopes:
        scope = raw.strip()
        if not scope:
            continue

        scope_parts: list[str] = scope.split(':')

        if len(scope_parts) != 2:
            raise ValueError(
                f"Illegal scope {scope}; expected format is 'resource:action'."
            )

        resource, action = scope_parts

        # Reject partial wildcard like "up*:read" or "uploads:r*".
        for seg in (resource, action):
            if seg != '*' and '*' in seg:
                raise ValueError(f'Partial wildcard in {scope} is not allowed.')

        # Concrete scope
        if '*' not in scope:
            if scope in known_scopes:
                resolved_scopes.add(scope)
            else:
                unknown_concrete.add(scope)
            continue

        # Wildcards
        if scope == '*:*':
            resolved_scopes |= known_scopes
            continue

        if resource == '*' and action != '*':  # *:action
            matches = {s for s in known_scopes if s.endswith(f':{action}')}
        else:  # resource:*
            matches = {s for s in known_scopes if s.startswith(f'{resource}:')}

        if matches:
            resolved_scopes |= matches
        else:
            unmatched_wildcards.add(scope)
        continue

    if unknown_concrete or unmatched_wildcards:
        err_msg: list[str] = ['Invalid scope configuration.']
        if unknown_concrete:
            err_msg.append(f'Unknown concrete scopes: {sorted(unknown_concrete)}')
        if unmatched_wildcards:
            err_msg.append(f'Wildcards matched nothing: {sorted(unmatched_wildcards)}')
        raise ValueError(' '.join(err_msg))

    return resolved_scopes
