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

import json
from pathlib import Path


def _definitions_dir() -> Path:
    return Path(__file__).resolve().parents[2] / 'nomad' / 'layouts' / 'definitions'


def test_builtin_layout_definitions_have_required_fields():
    for definition_file in sorted(_definitions_dir().glob('*.json')):
        if definition_file.name == 'defaults.json':
            continue

        with definition_file.open() as f:
            payload = json.load(f)

        assert payload['id']
        assert payload['label']
        assert payload['overview']['type']
