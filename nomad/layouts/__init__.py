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
