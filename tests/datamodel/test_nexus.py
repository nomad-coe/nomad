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

from nomad.metainfo import Definition, MSection, Section
from nomad.datamodel.metainfo import nexus
from nomad.datamodel import EntryArchive


@pytest.mark.parametrize('path,value', [
    pytest.param('base_classes.name', 'nexus_base_classes'),
    pytest.param('base_classes.NXobject.name', 'NXobject'),
    pytest.param('base_classes.NXentry.nx_kind', 'group'),
    pytest.param('base_classes.NXentry.DefaultAttribute.nx_value.type', str),
    pytest.param('base_classes.NXentry.nx_attribute_default', '*'),
    pytest.param('base_classes.NXentry.NXdataGroup', '*'),
    pytest.param('base_classes.NXdetector.RealTimeField', '*'),
    pytest.param('base_classes.NXentry.NXdataGroup.nx_optional', True),
    pytest.param('base_classes.NXentry.nx_group_data.section_def.nx_kind', 'group'),
    pytest.param('base_classes.NXentry.nx_group_data.section_def.nx_optional', True),
    pytest.param('base_classes.NXentry.nx_group_data.section_def.nx_name.type', str),
    pytest.param('base_classes.NXdetector.RealTimeField.name', 'RealTimeField'),
    pytest.param('base_classes.NXdetector.RealTimeField.nx_type', 'NX_NUMBER'),
    pytest.param('base_classes.NXdetector.RealTimeField.nx_units', 'NX_TIME'),
    pytest.param('base_classes.NXdetector.RealTimeField.nx_unit', '*'),
    pytest.param('base_classes.NXdetector.RealTimeField.nx_value', '*'),
    pytest.param('applications.NXarpes.NXentryGroup.NXdataGroup.nx_optional', False)
])
def test_assert_nexus_metainfo(path: str, value: Any):
    segments = path.split('.')
    package, definition_names = segments[0], segments[1:]

    current: Definition = getattr(nexus, package)
    for name in definition_names:
        for content in current.m_contents():
            if getattr(content, 'name', None) == name:
                current = cast(Definition, content)
                break

        else:
            current = getattr(current, name, None)

        if current is None:
            assert False, f'{path} does not exist'

    if value is '*':
        assert current is not None, f'{path} does not exist'
    elif value is None:
        assert current is None, f'{path} does exist'
    else:
        assert current == value, f'{path} has wrong value'

    if isinstance(current, Section):
        assert current.nx_kind is not None
        for base_section in current.all_base_sections:
            assert base_section.nx_kind == current.nx_kind


def test_use_nexus_metainfo():
    # pylint: disable=no-member
    archive = EntryArchive()
    archive.nexus = nexus.Nexus()
    archive.nexus.nx_application_arpes = nexus.NXarpes()
    archive.nexus.nx_application_arpes.nx_group_entry = nexus.NXarpes.NXentryGroup()
    archive.nexus.nx_application_arpes.nx_group_entry.nx_field_title = nexus.NXarpes.NXentryGroup.TitleField()
    archive.nexus.nx_application_arpes.nx_group_entry.nx_field_title.nx_value = 'my title'

    # Entry/default is not overwritten in NXarpes. Therefore technically, there is no attribute section
    # nexus.NXarpes.NXentryGroup.DefaultAttribute. We artifically extented inheritence to
    # include inner section/classes. So both options work:
    # archive.nexus.nx_application_arpes.nx_group_entry.nx_attribute_default = nexus.NXentry.DefaultAttribute()
    archive.nexus.nx_application_arpes.nx_group_entry.nx_attribute_default = nexus.NXarpes.NXentryGroup.DefaultAttribute()
    archive.nexus.nx_application_arpes.nx_group_entry.nx_attribute_default.nx_value = 'my default'
    # pylint: enable=no-member

    archive = EntryArchive.m_from_dict(archive.m_to_dict())
    assert archive.nexus.nx_application_arpes.nx_group_entry.nx_attribute_default.nx_value == 'my default'
    assert archive.nexus.nx_application_arpes.nx_group_entry.nx_field_title.nx_value == 'my title'

    # TODO remove
    # print(json.dumps(archive.m_to_dict(), indent=2))


@pytest.mark.parametrize('path', [
    pytest.param('NXarpes:app/NXentry:group/title:field/my title:value', id='field'),
    pytest.param('NXarpes:app/NXentry:group/default:attribute/my default:value', id='attribute')
])
def test_use_nexus_metainfo_reflectivly(path):
    archive = EntryArchive()
    archive.nexus = nexus.Nexus()  # pylint: disable=no-member
    parent_object: MSection = archive.nexus
    parent_definition: Section = nexus.Nexus.m_def  # pylint: disable=no-member

    segments = path.split('/')
    for segment in segments:
        name_or_value, kind = segment.split(':')
        if kind in ['app', 'group', 'field', 'attribute']:
            if kind == 'app':
                section_definition = nexus.applications.all_definitions[name_or_value]
                sub_section_definition = parent_definition.all_sub_sections[name_or_value.replace('NX', 'nx_application_')]

            if kind == 'group':
                section_definition = parent_definition.all_inner_section_definitions[f'{name_or_value}Group']
                sub_section_definition = parent_definition.all_sub_sections[f'nx_group_{name_or_value.replace("NX", "")}']

            if kind == 'field':
                section_definition = parent_definition.all_inner_section_definitions[f'{name_or_value.title()}Field']
                sub_section_definition = parent_definition.all_sub_sections[f'nx_field_{name_or_value}']

            if kind == 'attribute':
                section_definition = parent_definition.all_inner_section_definitions[f'{name_or_value.title()}Attribute']
                sub_section_definition = parent_definition.all_sub_sections[f'nx_attribute_{name_or_value}']

            new_object = section_definition.section_cls()
            parent_object.m_add_sub_section(sub_section_definition, new_object)

        elif kind == 'value':
            parent_object.m_set(parent_definition.all_quantities['nx_value'], name_or_value)

        else:
            assert False, 'kind does not exist'

        parent_object = new_object
        parent_definition = section_definition
