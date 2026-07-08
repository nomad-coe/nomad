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

import copy
import json
import os
from typing import Any

import jmespath

from nomad.layouts.constants import BUILTIN_NODE_TYPE_DEFAULTS
from nomad.utils import get_logger

logger = get_logger(__name__)


def _definitions_path() -> str:
    """Return the directory containing bundled layout JSON definitions."""
    return os.path.join(os.path.dirname(__file__), 'definitions')


def get_builtin_layout_definitions() -> list[dict[str, Any]]:
    """Load all bundled layout definitions in deterministic filename order."""
    definitions_path = _definitions_path()
    definitions: list[dict[str, Any]] = []
    if not os.path.exists(definitions_path):
        return definitions

    for filename in sorted(os.listdir(definitions_path)):
        if not filename.endswith('.json') or filename == 'defaults.json':
            continue
        with open(os.path.join(definitions_path, filename)) as handle:
            definitions.append(json.load(handle))
    return definitions


def get_builtin_node_type_defaults() -> dict[str, Any]:
    """Return built-in widget requests extended by optional JSON defaults."""
    defaults = copy.deepcopy(BUILTIN_NODE_TYPE_DEFAULTS)
    defaults_path = os.path.join(_definitions_path(), 'defaults.json')
    if not os.path.exists(defaults_path):
        return defaults
    with open(defaults_path) as handle:
        defaults.update(json.load(handle).get('node_types', {}))
    return defaults


def evaluate_query(query: dict[str, Any], data: dict[str, Any]) -> bool:
    """Evaluate layout match conditions with frontend ``evaluateQuery`` parity."""
    if not isinstance(query, dict):
        raise TypeError('conditions must be a plain object')

    for key, expected in query.items():
        sep_idx = key.rfind(':')
        if sep_idx == -1:
            raise ValueError(
                f'Condition key "{key}" must contain a ":" separating the JMESPath and operator'
            )

        path = key[:sep_idx]
        operator = key[sep_idx + 1 :].strip().lower()

        # JMESPath search
        actual = jmespath.search(path, data)

        expected_list = expected if isinstance(expected, list) else [expected]
        actual_list = actual if isinstance(actual, list) else [actual]

        # Use the same matching logic as the frontend
        matches = [v in actual_list for v in expected_list]

        if operator == 'any':
            if not any(matches):
                return False
        elif operator == 'all':
            if not all(matches):
                return False
        elif operator == 'none':
            if any(matches):
                return False
        else:
            raise ValueError(
                f'Unsupported operator "{operator}" in condition key "{key}" (use any | all | none)'
            )

    return True


class LayoutRegistry:
    """Load, validate, and match the layouts available to entry pages.

    Definitions and widget request defaults are loaded lazily from bundled JSON.
    Matching layouts are ordered by priority. A layout marked ``is_fallback`` is
    selected only when no specialized layout matches the entry; ``default`` is the
    final safety net when no definition carries that marker.

    The registry is deliberately internal today. It also provides the natural
    registration boundary for external plugins once a public custom-layout API is
    designed.
    """

    def __init__(self):
        """Create an unloaded registry."""
        self._registry: dict[str, Any] | None = None
        self._node_type_defaults: dict[str, Any] = {}

    def reset(self) -> None:
        """Discard loaded state so definitions are re-read on next access."""
        self._registry = None
        self._node_type_defaults = {}

    def _ensure_loaded(self) -> None:
        """Load registry contents on first access."""
        if self._registry is None:
            self._load_registry()

    def _load_registry(self) -> None:
        """Load bundled definitions and reject missing or duplicate IDs."""
        items: dict[str, dict[str, Any]] = {}
        for definition in get_builtin_layout_definitions():
            layout_id = definition.get('id')
            if not layout_id:
                raise ValueError('Layout definition must define an "id".')
            if layout_id in items:
                raise ValueError(f'Duplicate layout id registered: {layout_id}')
            items[layout_id] = copy.deepcopy(definition)

        self._node_type_defaults = get_builtin_node_type_defaults()
        self._registry = {'items': items}

    @property
    def node_type_defaults(self) -> dict[str, Any]:
        """Return default graph requests keyed by layout node type."""
        self._ensure_loaded()
        return self._node_type_defaults

    @property
    def registry(self) -> dict[str, Any]:
        """Return the loaded layout catalog keyed by layout ID."""
        self._ensure_loaded()
        return self._registry

    @property
    def fallback_layout_id(self) -> str:
        """Return the layout used when no specialized match rule succeeds."""
        for item in sorted(
            self.registry['items'].values(),
            key=lambda layout: (layout.get('priority', 0), layout.get('id', '')),
        ):
            if item.get('is_fallback'):
                return item['id']
        return 'default'

    def get_matching_layout_ids(self, data: dict[str, Any]) -> list[str]:
        """Return enabled layouts whose match rules evaluate against the context."""
        matching_ids: list[str] = []
        fallback_id: str | None = None

        sorted_items = sorted(
            (
                item
                for item in self.registry['items'].values()
                if item.get('enabled', True)
            ),
            key=lambda item: (item.get('priority', 0), item.get('id', '')),
        )

        for item in sorted_items:
            if item.get('is_fallback'):
                fallback_id = item['id']
                continue

            query = item.get('query')
            if not query:
                continue

            try:
                if evaluate_query(query, data):
                    matching_ids.append(item['id'])
            except Exception:
                logger.exception(
                    'Failed to evaluate layout match rule.', layout_id=item.get('id')
                )

        if not matching_ids and fallback_id:
            matching_ids.append(fallback_id)

        return matching_ids


registry = LayoutRegistry()
