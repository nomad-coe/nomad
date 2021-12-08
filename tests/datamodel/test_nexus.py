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

from typing import cast, Any
import pytest

from nomad.metainfo import Definition
from nomad.datamodel.metainfo import nexus


@pytest.mark.parametrize('path,value', [
    pytest.param('base_classes.name', 'nexus_base_classes'),
    pytest.param('base_classes.NXobject.name', 'NXobject'),
    pytest.param('base_classes.NXentry.nx_kind', 'group'),
    pytest.param('base_classes.NXentry.default.type', str),
    pytest.param('base_classes.NXentry.NXdataGroup', '*'),
    pytest.param('base_classes.NXentry.NXdataGroup.nx_optional', True),
    pytest.param('base_classes.NXentry.nx_group_data.section_def.nx_kind', 'group'),
    pytest.param('base_classes.NXentry.nx_group_data.section_def.nx_optional', True),
    pytest.param('base_classes.NXentry.nx_group_data.section_def.nx_name.type', str),
    pytest.param('applications.NXarpes.NXentryGroup.NXdataGroup.nx_optional', False)
])
def test_assert_nexus_metainfo(path: str, value: Any):
    segments = path.split('.')
    package, definition_names, property = segments[0], segments[1:-1], segments[-1]

    current: Definition = getattr(nexus, package)
    for name in definition_names:
        assert current or name == definition_names[-1], f'{path} has to exist in nexus_{package}'
        for content in current.m_contents():
            if getattr(content, 'name', None) == name:
                current = cast(Definition, content)
                break

        else:
            current = getattr(current, name, None)

    if value is '*':
        assert current is not None, f'{path} does not exist'
    elif value is None:
        assert current is None, f'{path} does exist'
    else:
        assert getattr(current, property) == value, f'{path} has wrong value'
