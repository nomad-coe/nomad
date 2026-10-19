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

import json
import os


def merge_dicts(base, update):
    result = base.copy()
    for key, value in update.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = merge_dicts(result[key], value)
        else:
            result[key] = value
    return result


def update_plotly_layout(figure, theme='light', spacing=8):
    """
    Update plotly template layout based on GUI logic

    Args:
        figure: figure object
        theme: 'light' or 'dark'
        spacing: Spacing value (default 8 based on GUI)

    Returns:
        Updated layout dict
    """
    layout = figure['layout']

    with_title = False
    if layout is not None:
        with_title = (
            layout.get('title', {}).get('text')
            or layout.get('template', {}).get('title', {}).get('text')
            or any(item.get('text') for item in layout.get('annotations', []))
        )

    default_layout = {
        'paper_bgcolor': '#1A1A1A' if theme == 'dark' else '#FFFFFF',
        'margin': {
            'l': 4 * spacing,
            'r': 1.5 * spacing,
            't': (5 if with_title else 1) * spacing,
            'b': 6 * spacing,
        },
    }

    updated_layout = merge_dicts(layout, default_layout)

    return updated_layout


def get_plotly_template(theme):
    base_dir = os.path.dirname(__file__)
    if theme == 'dark':
        file_path = os.path.join(base_dir, 'dark.json')
    else:
        file_path = os.path.join(base_dir, 'light.json')
    with open(file_path) as file:
        template = json.load(file)

    return template
