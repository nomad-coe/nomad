"""
Server-side layout support for entry overview rendering.
"""

from nomad.layouts.registry import registry
from nomad.layouts.service import (
    build_layout_context,
    resolve_layouts,
    derive_request_from_layout,
    derive_request_from_compiled_layout,
    create_layout_plan,
    compile_layout,
    calculate_request_from_layout,
    LayoutQueryIntent,
    get_layout_query_intent,
)

__all__ = [
    'registry',
    'build_layout_context',
    'resolve_layouts',
    'derive_request_from_layout',
    'derive_request_from_compiled_layout',
    'create_layout_plan',
    'compile_layout',
    'calculate_request_from_layout',
    'LayoutQueryIntent',
    'get_layout_query_intent',
]
