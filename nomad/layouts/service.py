from __future__ import annotations

import copy
from dataclasses import dataclass
from typing import Any

from nomad.graph.model import DirectiveType, RequestConfig
from nomad.layouts.constants import (
    _LAYOUT_CONTEXT_OVERRIDE_KEYS,
    _LAYOUT_RESOLUTION_KEYS,
    _LAYOUT_SEARCH_METADATA_KEYS,
    _PROPERTY_VALUE_NODE_TYPES,
)
from nomad.layouts.registry import evaluate_query, registry
from nomad.utils import get_logger

logger = get_logger(__name__)


@dataclass(frozen=True)
class LayoutResolution:
    """Matching, default, and explicitly resolved IDs for one entry context."""

    matching_layout_ids: list[str]
    default_layout_id: str
    resolved_layout_id: str


def build_layout_context(
    metadata: dict[str, Any], extra_context: dict[str, Any] | None = None
) -> dict[str, Any]:
    """Overlay archive metadata onto entry data for matching and compilation."""
    # Only a top-level copy is needed for overlaying extra context.
    # Downstream layout matching/request derivation treats nested values as read-only.
    context = dict(metadata)
    for key, value in (extra_context or {}).items():
        if key not in context or key in _LAYOUT_CONTEXT_OVERRIDE_KEYS:
            context[key] = value
    return context


def resolve_layouts(
    context: dict[str, Any], *, requested_layout_id: str | None = None
) -> LayoutResolution:
    """Choose and validate matching, default, and requested layout IDs."""
    matching_layout_ids = registry.get_matching_layout_ids(context)
    default_layout_id = (
        matching_layout_ids[0] if matching_layout_ids else registry.fallback_layout_id
    )
    if requested_layout_id is not None:
        requested_layout = registry.registry['items'].get(requested_layout_id)
        if requested_layout is None:
            raise ValueError(f'Unknown layout id: {requested_layout_id}')
        if not requested_layout.get('enabled', True):
            raise ValueError(f'Layout is disabled: {requested_layout_id}')
        if requested_layout_id not in matching_layout_ids:
            raise ValueError(f'Layout {requested_layout_id} does not match this entry.')
    resolved_layout_id = requested_layout_id or default_layout_id
    return LayoutResolution(
        matching_layout_ids=matching_layout_ids,
        default_layout_id=default_layout_id,
        resolved_layout_id=resolved_layout_id,
    )


def get_layout_definition(layout_id: str) -> dict[str, Any]:
    """Return a registered layout definition or raise a user-facing error."""
    layout_definition = registry.registry['items'].get(layout_id)
    if layout_definition is None:
        raise ValueError(f'Unknown layout id: {layout_id}')
    return layout_definition


# ==============================================================================
# Compile & Request Derivation
# ==============================================================================


def merge_requests(request: dict[str, Any], additional_request: dict[str, Any]) -> None:
    """Merge graph requests, preserving wildcard requests as the broadest request."""
    for key, value in additional_request.items():
        if key not in request:
            request[key] = copy.deepcopy(value)
            continue

        existing_value = request[key]
        if isinstance(existing_value, dict) and isinstance(value, dict):
            merge_requests(existing_value, value)
        elif existing_value == value:
            continue
        elif existing_value == '*' or value == '*':
            request[key] = '*'
        else:
            raise ValueError(
                f'Conflicting layout data requests for property {key}: {existing_value} vs {value}'
            )


def _compile_children(
    children: list[dict[str, Any]] | None,
    data: dict[str, Any],
    parent_editable: bool | None,
    _depth: int = 0,
) -> list[dict[str, Any]]:
    """Compile children and flatten children emitted by condition nodes."""
    compiled_children: list[dict[str, Any]] = []
    for child in children or []:
        result = compile_layout(child, data, parent_editable, _depth=_depth)
        if result is None:
            continue
        if isinstance(result, list):
            compiled_children.extend(result)
        else:
            compiled_children.append(result)
    return compiled_children


def compile_layout(
    node: dict[str, Any],
    data: dict[str, Any],
    parent_editable: bool | None = None,
    _depth: int = 0,
) -> dict[str, Any] | list[dict[str, Any]] | None:
    """
    Compile a declarative layout tree for a specific archive/search context.

    Conditional nodes disappear when their query does not match. Matching
    conditions are replaced by their compiled children so the GUI only receives
    concrete render nodes.
    """
    if _depth > 50:
        raise ValueError('Max layout recursion depth exceeded (circular reference?)')

    if node.get('type') == 'condition':
        if not node.get('query') or not evaluate_query(node['query'], data):
            return None
        return _compile_children(
            node.get('children'), data, parent_editable, _depth=_depth + 1
        )

    editable = node.get('editable')
    effective_editable = editable if editable is not None else parent_editable
    new_node = copy.deepcopy(node)
    if effective_editable is not None:
        new_node['editable'] = effective_editable

    # Inline apply built-in request defaults before node-local request overrides
    defaults = registry.node_type_defaults.get(new_node.get('type'), {})
    default_request = defaults.get('request') if isinstance(defaults, dict) else None
    if default_request:
        if new_node.get('request') is None:
            new_node['request'] = copy.deepcopy(default_request)
        else:
            merged_request = copy.deepcopy(default_request)
            merge_requests(merged_request, new_node['request'])
            new_node['request'] = merged_request

    if new_node.get('children'):
        new_node['children'] = _compile_children(
            new_node['children'], data, effective_editable, _depth=_depth + 1
        )

    return new_node


def _build_section_request(
    *,
    depth: int | None = None,
    exclude: list[str] | None = None,
    export_whole_package: bool = False,
) -> dict[str, Any]:
    """Build the graph request shared by archive section layout nodes."""
    m_request: dict[str, Any] = {
        'directive': 'plain',
        'include_definition': 'both',
        'm_def_format': 'short',
    }
    if export_whole_package:
        m_request['export_whole_package'] = True
    if depth is not None:
        m_request['depth'] = depth
    if exclude is not None:
        m_request['exclude'] = exclude

    return {
        'm_request': m_request,
        'm_def': {
            'm_request': {
                'directive': 'plain',
                'm_def_format': 'short',
                **({'export_whole_package': True} if export_whole_package else {}),
            }
        },
    }


def _add_prop_request(
    prop_request: Any, property_path: str, parent: dict[str, Any]
) -> None:
    """Translate slash-separated layout property paths into graph request trees."""
    new_request: dict[str, Any] = {}
    current = new_request
    parts = [part for part in property_path.split('/') if not part.isdigit()]
    for index, part in enumerate(parts):
        if index == len(parts) - 1:
            current[part] = prop_request
        else:
            is_top_level_data_section = index == 0 and part == 'data'
            if part not in current:
                current[part] = (
                    _build_section_request(depth=2, export_whole_package=True)
                    if is_top_level_data_section
                    else _build_section_request(exclude=['*'])
                )
            current = current[part]
    merge_requests(parent, new_request)


def calculate_request_from_layout(
    node: dict[str, Any] | list[dict[str, Any]] | None, request: dict[str, Any]
) -> None:
    """Derive the archive graph request needed to render a compiled layout tree."""
    if node is None:
        return
    if isinstance(node, list):
        for child in node:
            calculate_request_from_layout(child, request)
        return

    if not request:
        request.update(_build_section_request(exclude=['*']))

    for child in node.get('children') or []:
        calculate_request_from_layout(child, request)

    # Translate built-in property-bearing nodes into archive graph requests.
    property_path = node.get('property')
    if property_path:
        if node.get('type') in _PROPERTY_VALUE_NODE_TYPES:
            _add_prop_request('*', property_path, request)
        elif node.get('type') == 'sub_section':
            _add_prop_request(_build_section_request(depth=2), property_path, request)
        elif node.get('type') == 'table':
            sub_request = _build_section_request(exclude=['*'])
            for column in node.get('columns') or []:
                _add_prop_request('*', column['property'], sub_request)
            _add_prop_request(sub_request, property_path, request)

    if node.get('request'):
        merge_requests(request, node['request'])


def derive_request_from_compiled_layout(
    compiled: dict[str, Any] | list[dict[str, Any]] | None,
) -> dict[str, Any]:
    """Return the archive graph request needed by a compiled layout tree."""
    derived_request: dict[str, Any] = {}
    if compiled is not None:
        calculate_request_from_layout(compiled, derived_request)
    return derived_request


def derive_request_from_layout(
    layout_id: str,
    context: dict[str, Any],
) -> dict[str, Any]:
    """Compile a layout and return the archive graph request required by it."""
    layout_definition = get_layout_definition(layout_id)
    compiled = compile_layout(layout_definition['overview'], context)
    return derive_request_from_compiled_layout(compiled)


@dataclass(frozen=True)
class LayoutPlan:
    """Compiled layouts and the selected archive request for one entry."""

    matching_layouts: list[dict[str, Any]]
    default_layout_id: str
    resolved_layout_id: str
    archive_request: dict[str, Any]


def _strip_layout_requests(
    node: dict[str, Any] | list[dict[str, Any]] | None,
) -> dict[str, Any] | list[dict[str, Any]] | None:
    """Remove server-only request trees from a compiled layout response."""
    if isinstance(node, list):
        for child in node:
            _strip_layout_requests(child)
        return node
    if not isinstance(node, dict):
        return node
    node.pop('request', None)
    _strip_layout_requests(node.get('children'))
    return node


def create_layout_plan(
    context: dict[str, Any], *, requested_layout_id: str | None = None
) -> LayoutPlan:
    """Resolve and compile layouts once, then derive the selected archive request."""
    resolution = resolve_layouts(context, requested_layout_id=requested_layout_id)
    matching_layouts: list[dict[str, Any]] = []
    archive_request: dict[str, Any] = {}

    for layout_id in resolution.matching_layout_ids:
        definition = copy.deepcopy(get_layout_definition(layout_id))
        compiled_overview = compile_layout(definition['overview'], context)
        if layout_id == resolution.resolved_layout_id:
            archive_request = derive_request_from_compiled_layout(compiled_overview)
        _strip_layout_requests(compiled_overview)
        definition['overview'] = compiled_overview or {
            'type': 'container',
            'children': [],
        }
        matching_layouts.append(definition)

    return LayoutPlan(
        matching_layouts=matching_layouts,
        default_layout_id=resolution.default_layout_id,
        resolved_layout_id=resolution.resolved_layout_id,
        archive_request=archive_request,
    )


# ==============================================================================
# Graph Intent
# ==============================================================================


def _get_request_config_value(config: Any, key: str) -> Any:
    """Read a request option from either a dict or ``RequestConfig``."""
    if isinstance(config, dict):
        return config.get(key)
    return getattr(config, key, None)


@dataclass(frozen=True)
class LayoutQueryIntent:
    """What layout-related work an entry graph request needs before walking."""

    requires_search_metadata: bool
    requires_layout_resolution: bool
    auto_from_layout: bool
    requested_layout_id: str | None


def get_layout_query_intent(required_query: Any) -> LayoutQueryIntent:
    """Inspect a graph query for layout metadata/resolution requirements."""
    if not isinstance(required_query, dict):
        return LayoutQueryIntent(False, False, False, None)

    # Locate the archive request config without introducing a one-use wrapper.
    archive_query = required_query.get('archive')
    if isinstance(archive_query, RequestConfig):
        archive_config: Any | None = archive_query
    elif isinstance(archive_query, dict):
        archive_config = archive_query.get('m_request')
    else:
        archive_config = None

    directive = _get_request_config_value(archive_config, 'directive')
    auto_from_layout = directive in (
        DirectiveType.auto_from_layout,
        DirectiveType.auto_from_layout.value,
    )
    requested_layout_id = (
        _get_request_config_value(archive_config, 'layout_id')
        if auto_from_layout
        else None
    )

    requires_layout_resolution = auto_from_layout or any(
        key in required_query for key in _LAYOUT_RESOLUTION_KEYS
    )
    requires_search_metadata = requires_layout_resolution or any(
        key in required_query for key in _LAYOUT_SEARCH_METADATA_KEYS
    )

    return LayoutQueryIntent(
        requires_search_metadata=requires_search_metadata,
        requires_layout_resolution=requires_layout_resolution,
        auto_from_layout=auto_from_layout,
        requested_layout_id=requested_layout_id,
    )
