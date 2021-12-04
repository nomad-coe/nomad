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

from typing import Dict, List, Any
import xml.etree.ElementTree as ET
import os.path
import os
import sys
import numpy as np
import re

from nomad.utils import strip
from nomad.metainfo import (
    Section, Package, SubSection, Definition, Datetime, Bytes, MEnum, Quantity)
from nomad.datamodel import EntryArchive


# url_regexp from https://stackoverflow.com/questions/3809401/what-is-a-good-regular-expression-to-match-a-url
url_regexp = re.compile(r'(https?:\/\/(www\.)?[-a-zA-Z0-9@:%._\+~#=]{1,256}\.[a-zA-Z0-9()]{1,6}\b([-a-zA-Z0-9()@:%_\+.~#?&//=]*))')
xml_namespaces = {'nx': 'http://definition.nexusformat.org/nxdl/3.1'}

# TODO the validation still show some problems. Most notably there are a few higher
# dimensional fields with non number types, which the metainfo does not support
validate = False
current_package: Package = None

_definition_sections: Dict[str, Section] = dict()
_xml_parent_map: Dict[ET.Element, ET.Element] = None
_nx_doc_base = 'https://manual.nexusformat.org/classes'


def to_camel_case(snake_str: str, upper: bool = False) -> str:
    components = snake_str.split('_')
    if upper:
        return ''.join(f'{x[0].upper()}{x[1:]}' for x in components)

    return components[0] + ''.join(f'{x[0].upper()}{x[1:]}' for x in components[1:])


def nx_documenation_url(xml_node: ET.Element, nx_type: str):
    anchor_segments = []
    if nx_type != 'class':
        anchor_segments.append(nx_type)

    while xml_node is not None:
        if 'name' in xml_node.attrib:
            anchor_segments.append(xml_node.attrib['name'].replace('_', '-'))

        xml_node = _xml_parent_map.get(xml_node)

    anchor = "-".join([name.lower() for name in reversed(anchor_segments)])
    nx_package = current_package.name.replace('nexus_', '')
    doc_url = f'{_nx_doc_base}/{nx_package}/{anchor_segments[-1]}.html#{anchor}'
    return doc_url


def get_section(name: str, **kwargs) -> Section:
    '''
    Returns the 'existing' metainfo section for a given top-level nexus base-class name.

    This function ensures that sections for these base-classes are only created one.
    This allows to access the metainfo section even before it is generated from the base-class
    nexus definition.
    '''
    if name in _definition_sections:
        section = _definition_sections[name]
        section.more.update(**kwargs)
        return section

    section = Section(validate=validate, name=name, more=kwargs)
    current_package.section_definitions.append(section)
    _definition_sections[section.name] = section

    return section


def add_definition_properties(xml_node: ET.Element, definition: Definition, nx_type: str = None):
    '''
    Adds general metainfo definition properties (e.g. name, deprecated, description)
    from the given nexus XML node to the given metainfo definition.
    '''
    xml_attrs = xml_node.attrib

    if 'name' in xml_attrs:
        name = xml_attrs['name']
        if nx_type and nx_type != 'class':
            name = f'nx_{nx_type}_{name}'
        definition.name = name

    # Add a link to the nexus classes documentation
    links = []
    if nx_type is not None:
        doc_url = nx_documenation_url(xml_node, nx_type)
        if doc_url:
            links.append(doc_url)

    # Read additional urls from the XML doc element
    doc = xml_node.find('nx:doc', xml_namespaces)
    if doc is not None and doc.text is not None:
        definition.description = strip(doc.text)
        for match in url_regexp.findall(definition.description):
            links.append(match[0])

    if len(links) > 0:
        definition.links = links

    if 'deprecated' in xml_attrs:
        definition.deprecated = xml_attrs['deprecated']

    definition.more['nx_optional'] = xml_attrs.get(
        'optional', current_package.name == 'nexus_base_classes')

    if 'minOccurs' in xml_attrs:
        definition.more['nx_min_occurs'] = xml_attrs['minOccurs']
    if 'maxOccurs' in xml_attrs:
        definition.more['nx_max_occurs'] = xml_attrs['maxOccurs']
    if 'required' in xml_attrs:
        definition.more['nx_required'] = xml_attrs['required']
    if 'recommended' in xml_attrs:
        definition.more['nx_recommended'] = xml_attrs['recommended']

    # TODO there are probably even more nxdl attributes


def get_enum(xml_node: ET.Element):
    enumeration = xml_node.find('nx:enumeration', xml_namespaces)
    if enumeration is not None:
        enum_values = []
        for enum_value in enumeration.findall('nx:item', xml_namespaces):
            enum_values.append(enum_value.attrib['value'])
        return MEnum(*enum_values)
    return None


def add_attributes(xml_node: ET.Element, section: Section):
    '''
    Adds quantities for all attributes in the given nexus XML node to the given
    section.
    '''
    for attribute in xml_node.findall('nx:attribute', xml_namespaces):
        type: Any = get_enum(xml_node)
        if type is None:
            type = str
        quantity = Quantity(type=type, nx_kind='attribute')
        add_definition_properties(attribute, quantity)
        section.quantities.append(quantity)


# TODO There are more types in nxdl, but they are not used by the current base classes and
# application definitions.
_nx_types = {
    'NX_FLOAT': np.dtype(np.float64),
    'NX_CHAR': str,
    'NX_BOOLEAN': bool,
    'NX_INT': np.dtype(np.int64),
    'NX_NUMBER': np.dtype(np.number),
    'NX_POSINT': np.dtype(np.uint64),
    'NX_BINARY': Bytes,
    'NX_DATE_TIME': Datetime
}


def section_from_field(xml_node: ET.Element) -> Section:
    '''
    Generates a metainfo section for the given nexus field XML node.
    '''
    xml_attrs = xml_node.attrib
    nx_name = xml_attrs['name']
    name = to_camel_case(nx_name, True) + 'Field'
    section = Section(validate=validate, name=name, more=dict(nx_kind='field'))
    if nx_name.lower() != nx_name:
        section.quantities.append(Quantity(
            name='nx_name', type=str, description='''
                This is a nexus template property. This quantity holds the actual name used
                in the nexus data.'''))
    value = Quantity(name='nx_value', description='The value for this nexus field')
    section.quantities.append(value)

    if 'type' in xml_attrs:
        nx_type = xml_attrs['type']
        if nx_type not in _nx_types:
            raise NotImplementedError(f'type {nx_type} is not supported')
        value.type = _nx_types[nx_type]
    else:
        value.type = get_enum(xml_node)
        if value.type is None:
            value.type = Any

    if 'unit' in xml_attrs:
        # TODO map unit
        pass

    dimensions = xml_node.find('nx:dimensions', xml_namespaces)
    if dimensions is not None:
        shape = []
        for dimension in dimensions.findall('nx:dim', xml_namespaces):
            dimension_value: Any = dimension.attrib.get('value', '*')
            try:
                dimension_value = int(dimension_value)
            except ValueError:
                pass
            shape.append(dimension_value)
        value.shape = shape

    add_attributes(xml_node, section)

    return section


def add_group_properties(xml_node: ET.Element, definition_section: Section):
    '''
    Adds all properties that can be generated from the given nexus group XML node to
    the given (empty) metainfo section definition.
    '''
    add_attributes(xml_node, definition_section)

    for group in xml_node.findall('nx:group', xml_namespaces):
        assert 'type' in group.attrib, 'group has not type'

        type = group.attrib['type']
        base_section = get_section(type)
        empty_definition = len(group) == 0 or (
            len(group) == 1 and group.find('nx:doc', xml_namespaces) is not None)
        if empty_definition:
            # The group does not define anything new, we can directly use the base definition
            group_section = base_section
        else:
            if 'name' in group.attrib:
                name = to_camel_case(group.attrib['name'], True) + 'Group'
            else:
                name = to_camel_case(type, True) + 'Group'

            group_section = Section(validate=validate, name=name)
            group_section.base_sections = [base_section]
            definition_section.inner_section_definitions.append(group_section)
            add_group_properties(group, group_section)

        sub_section = SubSection(section_def=group_section, nx_kind='group')
        add_definition_properties(group, sub_section, nx_type='group')
        if sub_section.name is None:
            sub_section.name = type.replace('NX', 'nx_group_')
        definition_section.sub_sections.append(sub_section)

    for field in xml_node.findall('nx:field', xml_namespaces):
        assert 'name' in field.attrib, 'field has not name'

        field_section = section_from_field(field)
        definition_section.inner_section_definitions.append(field_section)

        more = dict(nx_kind='field')
        more.update(**{
            f'nx_{key}': field.attrib[key]
            for key in ['type', 'units'] if key in field.attrib})
        field_sub_section = SubSection(section_def=field_section, **more)
        add_definition_properties(field, field_sub_section, nx_type='field')
        definition_section.sub_sections.append(field_sub_section)

    return definition_section


def section_from_class(nxdl_file: str):
    '''
    Creates a metainfo section from the top-level nexus group definition in the given
    nxdl file and adds it to the current metainfo package.
    '''
    xml_tree = ET.parse(nxdl_file)
    definition = xml_tree.getroot()
    xml_attrs = definition.attrib
    global _xml_parent_map
    _xml_parent_map = {child: parent for parent in xml_tree.iter() for child in parent}

    assert xml_attrs.get('type') == 'group', 'definition is not a group'
    assert 'name' in xml_attrs

    definition_section = get_section(xml_attrs['name'], nx_kind=xml_attrs['type'])
    if 'extends' in xml_attrs:
        base_section = get_section(xml_attrs['extends'])
        definition_section.base_sections = [base_section]
    add_group_properties(definition, definition_section)
    add_definition_properties(definition, definition_section, nx_type='class')


def package_from_directory(path: str, package_name: str) -> Package:
    '''
    Creates a metainfo package from the given nexus directory. Will the respective
    metainfo definitions generated from all the nxdl files in that directory.
    '''
    global current_package
    current_package = Package(name=package_name)

    for nxdl_file in sorted(os.listdir(path)):
        if not nxdl_file.endswith('.nxdl.xml'):
            continue

        try:
            section_from_class(
                os.path.join(path, nxdl_file))
        except Exception as e:
            print(f'Exception while mapping {nxdl_file}', file=sys.stderr)
            raise e

    return current_package


nx_definitions_path = os.path.join(
    os.path.dirname(__file__),
    '../../../dependencies/nexus_definitions')


# We generate separated metainfo package for the nexus base classes and application
# definitions.
packages: List[Package] = []
for nx_package in ['base_classes', 'applications']:
    path = os.path.join(nx_definitions_path, nx_package)
    packages.append(package_from_directory(path, f'nexus_{nx_package}'))
# TODO there are problems generating with nx_package='contributed_definitions'


# We take the application definitions and create a common parent section that allows to
# include nexus in an EntryArchive.
application_package = packages[1]

nexus_section = Section(validate=validate, name='nexus')

for application_section in application_package.section_definitions:  # pylint: disable=not-an-iterable
    sub_section = SubSection(
        section_def=application_section,
        name=application_section.name.replace('NX', 'nx_application_'))
    nexus_section.sub_sections.append(sub_section)

application_package.section_definitions.append(nexus_section)

entry_archive_nexus_sub_section = SubSection(name='nexus', section_def=nexus_section)
EntryArchive.nexus = entry_archive_nexus_sub_section  # type: ignore
EntryArchive.m_def.sub_sections.append(entry_archive_nexus_sub_section)


# We need to initialize the metainfo definitions. This is usually done automatically,
# when the metainfo schema is defined though MSection Python classes.
for package in packages:
    package.init_metainfo()


# We skip the Python code generation for now and offer Python classes as variables
# TODO not necessary right now, could also be done case-by-case by the nexus parser
# python_module = sys.modules[__name__]
# for package in packages:
#     for section in package.section_definitions:
#         setattr(python_module, section.name, section.section_cls)
