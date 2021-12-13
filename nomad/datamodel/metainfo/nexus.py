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

from typing import Dict, Any
import xml.etree.ElementTree as ET
import os.path
import os
import sys
import numpy as np
import re

from nomad.utils import strip
from nomad.metainfo import (
    Section, Package, SubSection, Definition, Datetime, Bytes, Unit, MEnum, Quantity)
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


def get_or_create_section(name: str, **kwargs) -> Section:
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

    section = Section(validate=validate, name=name, **kwargs)
    current_package.section_definitions.append(section)
    _definition_sections[section.name] = section

    return section


def get_enum(xml_node: ET.Element):
    enumeration = xml_node.find('nx:enumeration', xml_namespaces)
    if enumeration is not None:
        enum_values = []
        for enum_value in enumeration.findall('nx:item', xml_namespaces):
            enum_values.append(enum_value.attrib['value'])
        return MEnum(*enum_values)
    return None


def add_common_properties(xml_node: ET.Element, definition: Definition):
    '''
    Adds general metainfo definition properties (e.g. deprecated, docs, optional, ...)
    from the given nexus XML node to the given metainfo definition.
    '''
    nx_kind = definition.more.get('nx_kind')
    xml_attrs = xml_node.attrib

    # Read properties from potential base section. Those are not inherited, but we
    # duplicate them for a nicer presentation
    if isinstance(definition, Section) and len(definition.base_sections) > 0:
        base_section = definition.base_sections[0]
        if base_section.description:
            definition.description = base_section.description
        if base_section.deprecated:
            definition.deprecated = base_section.deprecated
        if len(base_section.more) > 0:
            definition.more.update(**base_section.more)

    links = []
    if nx_kind is not None:
        doc_url = nx_documenation_url(xml_node, nx_kind)
        if doc_url:
            links.append(doc_url)

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

    # TODO there are probably even more nxdl attributes?


def add_attributes(xml_node: ET.Element, section: Section):
    '''
    Adds quantities for all attributes in the given nexus XML node to the given
    section.
    '''
    for attribute in xml_node.findall('nx:attribute', xml_namespaces):
        attribute_section = create_attribute_section(attribute, section)

        section.sub_sections.append(SubSection(
            section_def=attribute_section, nx_kind='attribute',
            name=f'nx_attribute_{attribute.attrib["name"]}'))


def add_group_properties(xml_node: ET.Element, section: Section):
    '''
    Adds all properties that can be generated from the given nexus group XML node to
    the given (empty) metainfo section definition.
    '''
    for group in xml_node.findall('nx:group', xml_namespaces):
        group_section = create_group_section(group, section)
        section.inner_section_definitions.append(group_section)

        if 'name' in group.attrib:
            name = f'nx_group_{group.attrib["name"]}'
        else:
            name = group.attrib['type'].replace('NX', 'nx_group_')

        max_occurs = group.attrib.get('maxOccurs', '0')
        repeats = max_occurs == 'unbounded' or int(max_occurs) > 1
        section.sub_sections.append(SubSection(
            section_def=group_section, nx_kind='group', name=name, repeats=repeats))

    for field in xml_node.findall('nx:field', xml_namespaces):
        field_section = create_field_section(field, section)

        section.sub_sections.append(SubSection(
            section_def=field_section, nx_kind='field',
            name=f'nx_field_{field.attrib["name"]}'))

    add_attributes(xml_node, section)


def add_template_properties(xml_node: ET.Element, section: Section):
    '''
    Adds potential abilities of a group or field section to act as a TEMPLATE or
    nameType="any" definition.
    '''
    is_template = section.name.lower() != section.name or 'name' not in xml_node.attrib
    if is_template:
        section.quantities.append(Quantity(
            name='nx_name', type=str, default=xml_node.attrib.get('name'), description='''
                This is a nexus template property. This quantity holds the actual name used
                in the nexus data.'''))


def add_base_section(section: Section, container: Section, default_base_section: Section = None):
    '''
    Potentially adds a base section to the given section, if the given container has
    a base-section with a suitable base.
    '''
    base_section = container.all_inner_section_definitions.get(section.name, None)
    if base_section:
        assert base_section.nx_kind == section.nx_kind, 'base section has wrong nexus kind'
    else:
        base_section = default_base_section

    if base_section:
        section.base_sections = [base_section]


def create_attribute_section(xml_node: ET.Element, container: Section) -> Section:
    '''
    Creates a metainfo section from the nexus attribute given as xml node.
    '''
    xml_attrs = xml_node.attrib
    assert 'name' in xml_attrs, 'attribute has not name'

    attribute_section = Section(
        validate=validate, nx_kind='attribute',
        name=to_camel_case(xml_attrs['name'], True) + 'Attribute')
    add_base_section(attribute_section, container)
    container.inner_section_definitions.append(attribute_section)

    base_value_quantity = attribute_section.all_quantities.get('nx_value')
    if base_value_quantity is None:
        value_quantity = Quantity(
            name='nx_value', description='The value for this nexus attribute')
    else:
        value_quantity = base_value_quantity.m_copy()
    attribute_section.quantities.append(value_quantity)

    enum_type = get_enum(xml_node)
    if enum_type is not None:
        value_quantity.tyoe = enum_type

    if value_quantity.type is None:
        value_quantity.type = str

    add_common_properties(xml_node, attribute_section)
    add_template_properties(xml_node, attribute_section)

    return attribute_section


def create_field_section(xml_node: ET.Element, container: Section):
    '''
    Creates a metainfo section from the nexus field given as xml node.
    '''
    xml_attrs = xml_node.attrib

    assert 'name' in xml_attrs, 'field has not name'
    name = to_camel_case(xml_attrs['name'], True) + 'Field'
    field_section = Section(validate=validate, nx_kind='field', name=name)
    add_base_section(field_section, container)
    container.inner_section_definitions.append(field_section)

    add_template_properties(xml_node, field_section)

    base_value_quantity = field_section.all_quantities.get('nx_value')
    if base_value_quantity:
        value_quantity = base_value_quantity.m_copy()
    else:
        value_quantity = Quantity(name='nx_value', description='The value for this nexus field')
    field_section.quantities.append(value_quantity)

    if 'type' in xml_attrs:
        nx_type = xml_attrs['type']
        if nx_type not in _nx_types:
            raise NotImplementedError(f'type {nx_type} is not supported')
        field_section.more['nx_type'] = nx_type

        if value_quantity.type is None or value_quantity.type is Any or nx_type != 'NX_CHAR':
            value_quantity.type = _nx_types[nx_type]

    enum_type = get_enum(xml_node)
    if enum_type:
        value_quantity.type = enum_type

    if value_quantity.type is None:
        value_quantity.type = Any

    if 'units' in xml_attrs:
        field_section.more['nx_units'] = xml_attrs['units']
        if xml_attrs['units'] != 'NX_UNITLESS':
            # TODO a default could be created from the nx_units value
            field_section.quantities.append(Quantity(
                name='nx_unit', type=Unit,
                description='The specific unit for that this fields data has.'))

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
        value_quantity.shape = shape

    add_common_properties(xml_node, field_section)
    add_attributes(xml_node, field_section)

    return field_section


def create_group_section(xml_node: ET.Element, container: Section) -> Section:
    '''
    Creates a metainfo section from the nexus group given as xml node.
    '''
    xml_attrs = xml_node.attrib
    type = xml_attrs['type']

    if 'name' in xml_attrs:
        name = to_camel_case(xml_attrs['name'], True) + 'Group'
    else:
        name = to_camel_case(type, True) + 'Group'

    group_section = Section(validate=validate, nx_kind='group', name=name)
    add_base_section(group_section, container, get_or_create_section(type))

    add_common_properties(xml_node, group_section)
    add_template_properties(xml_node, group_section)
    add_group_properties(xml_node, group_section)

    return group_section


def create_class_section(xml_node: ET.Element) -> Section:
    '''
    Creates a metainfo section from the top-level nexus definition given as xml node.
    '''
    xml_attrs = xml_node.attrib
    assert 'name' in xml_attrs

    class_section = get_or_create_section(xml_attrs['name'], nx_kind=xml_attrs['type'])

    if 'extends' in xml_attrs:
        base_section = get_or_create_section(xml_attrs['extends'])
        class_section.base_sections = [base_section]

    add_common_properties(xml_node, class_section)
    add_group_properties(xml_node, class_section)

    return class_section


def create_package_from_nxdl_directory(path: str) -> Package:
    '''
    Creates a metainfo package from the given nexus directory. Will generate the respective
    metainfo definitions from all the nxdl files in that directory.
    '''
    global current_package
    current_package = Package(name=f'nexus_{os.path.basename(path)}')

    for nxdl_file in sorted(os.listdir(path)):
        if not nxdl_file.endswith('.nxdl.xml'):
            continue

        try:
            nxdl_path = os.path.join(path, nxdl_file)
            xml_tree = ET.parse(nxdl_path)
            xml_node = xml_tree.getroot()

            global _xml_parent_map
            _xml_parent_map = {child: parent for parent in xml_tree.iter() for child in parent}

            assert xml_node.attrib.get('type') == 'group', 'definition is not a group'

            # The section gets already implicitly added to current_package by
            # get_or_create_section
            create_class_section(xml_node)

        except Exception as e:
            print(f'Exception while mapping {nxdl_file}', file=sys.stderr)
            raise e

    return current_package


nx_definitions_path = os.path.join(
    os.path.dirname(__file__),
    '../../../dependencies/nexus_definitions')


# We generate separated metainfo package for the nexus base classes and application
# definitions.
base_classes = create_package_from_nxdl_directory(os.path.join(nx_definitions_path, 'base_classes'))
applications = create_package_from_nxdl_directory(os.path.join(nx_definitions_path, 'applications'))
# TODO there are problems generating with nx_package='contributed_definitions'
packages = (base_classes, applications)

# We take the application definitions and create a common parent section that allows to
# include nexus in an EntryArchive.
nexus_section = Section(validate=validate, name='Nexus')

for application_section in applications.section_definitions:  # pylint: disable=not-an-iterable
    sub_section = SubSection(
        section_def=application_section,
        name=application_section.name.replace('NX', 'nx_application_'))
    nexus_section.sub_sections.append(sub_section)

applications.section_definitions.append(nexus_section)

entry_archive_nexus_sub_section = SubSection(name='nexus', section_def=nexus_section)
EntryArchive.nexus = entry_archive_nexus_sub_section  # type: ignore
EntryArchive.m_def.sub_sections.append(entry_archive_nexus_sub_section)


# We need to initialize the metainfo definitions. This is usually done automatically,
# when the metainfo schema is defined though MSection Python classes.
for package in packages:
    package.init_metainfo()


# We skip the Python code generation for now and offer Python classes as variables
# TODO not necessary right now, could also be done case-by-case by the nexus parser
python_module = sys.modules[__name__]
for package in packages:
    for section in package.section_definitions:
        setattr(python_module, section.name, section.section_cls)
