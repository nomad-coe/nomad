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

import os
import os.path
import re
import sys
# noinspection PyPep8Naming
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional, Union

import numpy as np
from toposort import toposort_flatten

from nexusutils.nexus import nexus  # pylint: disable=import-error
from nomad.datamodel import EntryArchive
from nomad.metainfo import (
    Attribute, Bytes, Datetime, Definition, MEnum, Package, Property, Quantity, Section, SubSection)
from nomad.utils import get_logger, strip

# __URL_REGEXP from
# https://stackoverflow.com/questions/3809401/what-is-a-good-regular-expression-to-match-a-url
__URL_REGEXP = re.compile(
    r'(?i)\b((?:https?://|www\d{0,3}[.]|[a-z0-9.\-]+[.][a-z]{2,4}/)'
    r'(?:[^\s()<>]+|\(([^\s()<>]+|(\([^\s()<>]+\)))*\))+'
    r'(?:\(([^\s()<>]+|(\([^\s()<>]+\)))*\)|[^\s`!()\[\]{};:\'".,<>?«»“”‘’]))')
# noinspection HttpUrlsUsage
__XML_NAMESPACES = {'nx': 'http://definition.nexusformat.org/nxdl/3.1'}

# TO DO the validation still show some problems. Most notably there are a few higher
# dimensional fields with non number types, which the metainfo does not support

__section_definitions: Dict[str, Section] = dict()

__logger = get_logger(__name__)

VALIDATE = False

__XML_PARENT_MAP: Dict[ET.Element, ET.Element]
__NX_DOC_BASE = 'https://manual.nexusformat.org/classes'
__NX_TYPES = {  # Primitive Types,  'ISO8601' is the only type not defined here
    'NX_COMPLEX': np.float64,
    'NX_FLOAT': np.float64,
    'NX_CHAR': str,
    'NX_BOOLEAN': bool,
    'NX_INT': np.int64,
    'NX_UINT': np.uint64,
    'NX_NUMBER': np.float64,
    'NX_POSINT': np.uint64,
    'NX_BINARY': Bytes,
    'NX_DATE_TIME': Datetime
}


class NXUnitSet:
    '''
    maps from `NX_` token to dimensionality
    None -> disable dimensionality check
    '1' -> dimensionless quantities
    'transformation' -> Specially handled in metainfo
    '''
    mapping: dict = {
        'NX_ANGLE': '[angle]',
        'NX_ANY': None,
        'NX_AREA': '[area]',
        'NX_CHARGE': '[charge]',
        'NX_COUNT': '1',
        'NX_CROSS_SECTION': '[area]',
        'NX_CURRENT': '[current]',
        'NX_DIMENSIONLESS': '1',
        'NX_EMITTANCE': '[length] * [angle]',
        'NX_ENERGY': '[energy]',
        'NX_FLUX': '1 / [time] / [area]',
        'NX_FREQUENCY': '[frequency]',
        'NX_LENGTH': '[length]',
        'NX_MASS': '[mass]',
        'NX_MASS_DENSITY': '[mass] / [volume]',
        'NX_MOLECULAR_WEIGHT': '[mass] / [substance]',
        'NX_PERIOD': '[time]',
        'NX_PER_AREA': '1 / [area]',
        'NX_PER_LENGTH': '1 / [length]',
        'NX_POWER': '[power]',
        'NX_PRESSURE': '[pressure]',
        'NX_PULSES': '1',
        'NX_SCATTERING_LENGTH_DENSITY': '1 / [area]',
        'NX_SOLID_ANGLE': '[angle] * [angle]',
        'NX_TEMPERATURE': '[temperature]',
        'NX_TIME': '[time]',
        'NX_TIME_OF_FLIGHT': '[time]',
        'NX_TRANSFORMATION': 'transformation',
        'NX_UNITLESS': '1',
        'NX_VOLTAGE': '[energy] / [current] / [time]',
        'NX_VOLUME': '[volume]',
        'NX_WAVELENGTH': '[length]',
        'NX_WAVENUMBER': '1 / [length]'
    }

    @staticmethod
    def normalise(value: str) -> str:
        '''
        Normalise the given token
        '''
        value = value.upper()
        if not value.startswith('NX_'):
            value = 'NX_' + value
        return value

    @staticmethod
    def is_nx_token(value: str) -> bool:
        '''
        Check if a given token is one of NX tokens
        '''
        return NXUnitSet.normalise(value) in NXUnitSet.mapping.keys()


def __to_camel_case(snake_str: str, upper: bool = False) -> str:
    '''
    Take as input a snake case variable and return a camel case one
    '''
    components = snake_str.split('_')

    if upper:
        return ''.join(x.capitalize() for x in components)

    return components[0] + ''.join(x.capitalize() for x in components[1:])


def __to_root(xml_node: ET.Element) -> ET.Element:
    '''
    get the root element
    '''
    elem = xml_node
    while True:
        parent = __XML_PARENT_MAP.get(elem)
        if parent is None:
            break
        elem = parent

    return elem


def __if_base(xml_node: ET.Element) -> bool:
    '''
    retrieves the category from the root element
    '''
    return __to_root(xml_node).get('category') == 'base'


def __if_repeats(name: str, max_occurs: str) -> bool:
    repeats = any(char.isupper() for char in name) or max_occurs == 'unbounded'

    if max_occurs.isdigit():
        repeats = repeats or int(max_occurs) > 1

    return repeats


def __if_template(name: Optional[str]) -> bool:
    return name is None or name.lower() != name


def __get_documentation_url(
        xml_node: ET.Element, nx_type: Optional[str]) -> Optional[str]:
    '''
    Get documentation url
    '''
    if nx_type is None:
        return None

    anchor_segments = []
    if nx_type != 'class':
        anchor_segments.append(nx_type)

    while True:
        nx_type = xml_node.get('type')
        if nx_type:
            nx_type = nx_type.replace('NX', '')
        segment = xml_node.get('name', nx_type)  # type: ignore
        anchor_segments.append(segment.replace('_', '-'))

        xml_parent = xml_node
        xml_node = __XML_PARENT_MAP.get(xml_node)
        if xml_node is None:
            break

    nx_package = xml_parent.get('nxdl_base').split('/')[-1]
    anchor = "-".join([name.lower() for name in reversed(anchor_segments)])
    return f'{__NX_DOC_BASE}/{nx_package}/{anchor_segments[-1]}.html#{anchor}'


def __to_section(name: str, **kwargs) -> Section:
    '''
    Returns the 'existing' metainfo section for a given top-level nexus base-class name.

    This function ensures that sections for these base-classes are only created one.
    This allows to access the metainfo section even before it is generated from the base
    class nexus definition.
    '''
    if name in __section_definitions:
        section = __section_definitions[name]
        section.more.update(**kwargs)
        return section

    section = Section(validate=VALIDATE, name=name, **kwargs)

    __section_definitions[name] = section

    return section


def __get_enumeration(xml_node: ET.Element) -> Optional[MEnum]:
    '''
    Get the enumeration field from xml node
    '''
    enumeration = xml_node.find('nx:enumeration', __XML_NAMESPACES)
    if enumeration is None:
        return None

    items = enumeration.findall('nx:item', __XML_NAMESPACES)

    return MEnum([value.attrib['value'] for value in items])


def __add_common_properties(xml_node: ET.Element, definition: Definition):
    '''
    Adds general metainfo definition properties (e.g., deprecated, docs, optional, ...)
    from the given nexus XML node to the given metainfo definition.
    '''
    xml_attrs = xml_node.attrib

    # Read properties from potential base section. Those are not inherited, but we
    # duplicate them for a nicer presentation
    if isinstance(definition, Section) and definition.base_sections:
        base_section = definition.base_sections[0]
        if base_section.description:
            definition.description = base_section.description
        if base_section.deprecated:
            definition.deprecated = base_section.deprecated
        if base_section.more:
            definition.more.update(**base_section.more)

    links = []
    doc_url = __get_documentation_url(xml_node, definition.more.get('nx_kind'))
    if doc_url:
        links.append(doc_url)

    doc = xml_node.find('nx:doc', __XML_NAMESPACES)
    if doc is not None and doc.text is not None:
        definition.description = strip(doc.text)
        links.extend([match[0] for match in __URL_REGEXP.findall(definition.description)])

    if links:
        definition.links = links

    for key, value in xml_attrs.items():
        if key == 'deprecated':
            definition.deprecated = value
            continue
        if 'nxdl_base' in key or 'schemaLocation' in key:
            continue
        definition.more['nx_' + key] = value

    if 'optional' not in xml_attrs:
        definition.more['nx_optional'] = __if_base(xml_node)


def __create_attributes(xml_node: ET.Element, definition: Union[Section, Property]):
    '''
    Add all attributes in the given nexus XML node to the given
    Quantity or SubSection using the Attribute class (new mechanism).

    todo: account for more attributes of attribute, e.g., default, minOccurs
    '''
    for attribute in xml_node.findall('nx:attribute', __XML_NAMESPACES):
        name = attribute.get('name') + "__attribute"

        nx_enum = __get_enumeration(attribute)
        if nx_enum:
            nx_type = nx_enum
            nx_shape: List[str] = []
        else:
            nx_type = __NX_TYPES[attribute.get('type', 'NX_CHAR')]  # type: ignore
            has_bound = False
            has_bound |= 'minOccurs' in attribute.attrib
            has_bound |= 'maxOccurs' in attribute.attrib
            if has_bound:
                nx_min_occurs = attribute.get('minOccurs', '0')  # type: ignore
                nx_max_occurs = attribute.get('maxOccurs', '*')  # type: ignore
                if nx_max_occurs == 'unbounded':
                    nx_max_occurs = '*'
                nx_shape = [f'{nx_min_occurs}..{nx_max_occurs}']
            else:
                nx_shape = []

        # check if the attribute exist
        # if yes then modify directly
        # if not create a new one and append to the list
        for m_attribute in definition.attributes:
            if m_attribute.name == name:
                m_attribute.shape = nx_shape
                m_attribute.type = nx_type

                for name, value in attribute.items():
                    m_attribute.more[f'nx_{name}'] = value

                break
        else:
            m_attribute = Attribute(
                name=name, variable=__if_template(name), shape=nx_shape, type=nx_type)

            for name, value in attribute.items():
                m_attribute.more[f'nx_{name}'] = value

            definition.attributes.append(m_attribute)

    m_nx_data_path_exists = False
    for attrib in definition.attributes:
        if attrib.name == "m_nx_data_path":
            m_nx_data_path_exists = True
            break

        # TODO: Add call to __add_common for description of attribs

    if not m_nx_data_path_exists:
        definition.attributes.append(Attribute(
            name="m_nx_data_path", variable=False, shape=[], type=str,
            description="This is a nexus template property. "
            "This attribute holds the actual path of the value in the nexus data."))

    if isinstance(definition, Quantity):
        for nx_array_attr in ['nx_data_mean', 'nx_data_var', 'nx_data_min', 'nx_data_max']:
            for attrib in definition.attributes:
                if attrib.name == nx_array_attr:
                    break
            else:
                definition.attributes.append(Attribute(
                    name=nx_array_attr, variable=False, shape=[], type=np.float64,
                    description="This is a nexus template property. "
                    "This attribute holds specific  statistics on the nexus data array."))


def __create_field(xml_node: ET.Element, container: Section) -> Quantity:
    '''
    Creates a metainfo quantity from the nexus field given as xml node.
    '''
    xml_attrs = xml_node.attrib

    # name
    assert 'name' in xml_attrs, 'Expecting name to be present'
    name = xml_attrs['name'] + '__field'

    # type
    nx_type = xml_attrs.get('type', 'NX_CHAR')
    if nx_type not in __NX_TYPES:
        raise NotImplementedError(f'Type {nx_type} is not supported for the moment for {name}.')

    # enumeration
    enum_type = __get_enumeration(xml_node)

    # dimensionality
    nx_dimensionality = xml_attrs.get('units', None)
    if nx_dimensionality:
        if nx_dimensionality not in NXUnitSet.mapping:
            raise NotImplementedError(f'Unit {nx_dimensionality} is not supported for {name}.')
        dimensionality = NXUnitSet.mapping[nx_dimensionality]
    else:
        dimensionality = None

    # shape
    shape: list = []
    nx_shape: list = []
    dimensions = xml_node.find('nx:dimensions', __XML_NAMESPACES)
    if dimensions is not None:
        for dimension in dimensions.findall('nx:dim', __XML_NAMESPACES):
            dimension_value: str = dimension.attrib.get('value', '0..*')
            nx_shape.append(dimension_value)

    value_quantity: Quantity = None  # type: ignore

    # copy from base to inherit from it
    if container.base_sections is not None:
        base_quantity: Quantity = container.base_sections[0].all_quantities.get(name)
        if base_quantity:
            value_quantity = base_quantity.m_copy(deep=True)

    # create quantity
    if value_quantity is None:
        value_quantity = Quantity(name=name, flexible_unit=True)

    value_quantity.variable = __if_template(name)

    # check parent type compatibility
    parent_type = getattr(value_quantity, 'type', None)
    if not isinstance(parent_type, MEnum):
        # if parent type is not MEnum then overwrite whatever given
        value_quantity.type = enum_type if enum_type else __NX_TYPES[nx_type]
    elif enum_type:
        # only when derived type is also MEnum to allow overwriting
        value_quantity.type = enum_type

    value_quantity.dimensionality = dimensionality
    value_quantity.shape = shape
    value_quantity.more.update(dict(nx_kind='field', nx_type=nx_type, nx_shape=nx_shape))

    __add_common_properties(xml_node, value_quantity)
    __create_attributes(xml_node, value_quantity)

    container.quantities.append(value_quantity)

    return value_quantity


def __create_group(xml_node: ET.Element, root_section: Section):
    '''
    Adds all properties that can be generated from the given nexus group XML node to
    the given (empty) metainfo section definition.
    '''
    __create_attributes(xml_node, root_section)

    for group in xml_node.findall('nx:group', __XML_NAMESPACES):
        xml_attrs = group.attrib

        assert 'type' in xml_attrs, 'Expecting type to be present'
        nx_type = xml_attrs['type']

        nx_name = xml_attrs.get('name', nx_type)
        group_section = Section(validate=VALIDATE, nx_kind='group', name=nx_name)

        __attach_base_section(group_section, root_section, __to_section(nx_type))
        __copy_base_attributes(group_section, nx_type)
        __add_common_properties(group, group_section)

        nx_name = xml_attrs.get('name', nx_type.replace('NX', '').upper())
        group_subsection = SubSection(
            section_def=group_section,
            nx_kind='group',
            name=nx_name,
            repeats=__if_repeats(nx_name, xml_attrs.get('maxOccurs', '0')),
            variable=__if_template(nx_name))

        root_section.inner_section_definitions.append(group_section)

        root_section.sub_sections.append(group_subsection)

        __create_group(group, group_section)

    for field in xml_node.findall('nx:field', __XML_NAMESPACES):
        __create_field(field, root_section)


def __attach_base_section(section: Section, container: Section, default: Section):
    '''
    Potentially adds a base section to the given section, if the given container has
    a base-section with a suitable base.
    '''
    base_section = container.all_inner_section_definitions.get(section.name)
    if base_section:
        assert base_section.nx_kind == section.nx_kind, 'Base section has wrong kind'
    else:
        base_section = default

    section.base_sections = [base_section]


def __copy_base_attributes(destination: Section, source_name: str):
    '''
    Copy attributes from base subsection to derived subsection.

    Attributes are stored in `SubSection.attributes` list. They are not inherited
    thus need to be manually copied.
    '''
    source: Section = __section_definitions.get(source_name)

    if not source or not destination or source is destination:
        return

    for m_attribute in source.attributes:
        destination.attributes.append(Attribute(
            name=m_attribute.name,
            type=m_attribute.type,
            shape=m_attribute.shape,
            variable=m_attribute.variable,
            **m_attribute.more))


def __create_class_section(xml_node: ET.Element) -> Section:
    '''
    Creates a metainfo section from the top-level nexus definition given as xml node.
    '''
    xml_attrs = xml_node.attrib
    assert 'name' in xml_attrs, 'Expecting name to be present'
    assert 'type' in xml_attrs, 'Expecting type to be present'
    assert 'category' in xml_attrs, 'Expecting category to be present'

    nx_name = xml_attrs['name']
    nx_type = xml_attrs['type']
    nx_category = xml_attrs['category']

    class_section: Section = __to_section(
        nx_name, nx_kind=nx_type, nx_category=nx_category)

    if 'extends' in xml_attrs:
        base_section = __to_section(xml_attrs['extends'])
        class_section.base_sections = [base_section]
        __copy_base_attributes(class_section, base_section.name)

    __add_common_properties(xml_node, class_section)

    __create_group(xml_node, class_section)

    return class_section


def __sort_nxdl_files(paths):
    '''
    Sort all definitions based on dependencies
    '''

    name_node_map = {}
    name_dependency_map = {}
    for path in paths:
        for nxdl_file in os.listdir(path):
            if not nxdl_file.endswith('.nxdl.xml'):
                continue
            xml_node = ET.parse(os.path.join(path, nxdl_file)).getroot()
            xml_node.set('nxdl_base', path)
            assert xml_node.get('type') == 'group', 'definition is not a group'
            xml_name = xml_node.get('name')
            name_node_map[xml_name] = xml_node
            dependency_list = []
            if 'extends' in xml_node.attrib:
                dependency_list.append(xml_node.get('extends'))
            for child in xml_node.iter():
                if child.tag.endswith('group') and child.get('type') != xml_name:
                    dependency_list.append(child.get('type'))
            name_dependency_map[xml_name] = set(dependency_list)

    # manually remove deprecated circular dependency
    name_dependency_map['NXgeometry'].remove('NXorientation')
    name_dependency_map['NXgeometry'].remove('NXtranslation')

    sorted_nodes = toposort_flatten(name_dependency_map)
    validated_names = []
    for node in sorted_nodes:
        if node in name_node_map:
            validated_names.append(name_node_map[node])
        else:
            parent_nodes = []
            for name, dependencies in name_dependency_map.items():
                if node in dependencies:
                    parent_nodes.append(name)
            __logger.error('Missing dependency (incorrect group type).', target_name=node, used_by=parent_nodes)

    return validated_names


def __add_section_from_nxdl(xml_node: ET.Element) -> Optional[Section]:
    '''
    Creates a metainfo section from a nxdl file.
    '''
    try:
        global __XML_PARENT_MAP  # pylint: disable=global-statement
        __XML_PARENT_MAP = {
            child: parent for parent in xml_node.iter() for child in parent}

        return __create_class_section(xml_node)

    except NotImplementedError as err:
        __logger.error('Fail to generate metainfo.', target_name=xml_node.attrib['name'], exe_info=str(err))
        return None


def __create_package_from_nxdl_directories(nexus_section: Section) -> Package:
    '''
    Creates a metainfo package from the given nexus directory. Will generate the
    respective metainfo definitions from all the nxdl files in that directory.
    '''
    package = Package(name='nexus')

    folder_list = ('base_classes', 'contributed_definitions', 'applications')
    paths = [os.path.join(
        nexus.get_nexus_definitions_path(), folder) for folder in folder_list]

    sections = []
    for nxdl_file in __sort_nxdl_files(paths):
        section = __add_section_from_nxdl(nxdl_file)
        if section is not None:
            sections.append(section)

    sections.sort(key=lambda x: x.name)

    for section in sections:
        package.section_definitions.append(section)
        if section.nx_category == 'application':
            nexus_section.sub_sections.append(
                SubSection(section_def=section, name=section.name))

    return package


nexus_metainfo_package: Optional[Package] = None  # pylint: disable=C0103


def init_nexus_metainfo():
    '''
    Initializes the metainfo package for the nexus definitions.
    '''
    global nexus_metainfo_package  # pylint: disable=global-statement

    if nexus_metainfo_package is not None:
        return

    # We take the application definitions and create a common parent section that allows
    # to include nexus in an EntryArchive.
    nexus_section = Section(validate=VALIDATE, name='NeXus')

    nexus_metainfo_package = __create_package_from_nxdl_directories(nexus_section)

    EntryArchive.nexus = SubSection(name='nexus', section_def=nexus_section)
    EntryArchive.nexus.init_metainfo()
    EntryArchive.m_def.sub_sections.append(EntryArchive.nexus)

    nexus_metainfo_package.section_definitions.append(nexus_section)

    # We need to initialize the metainfo definitions. This is usually done automatically,
    # when the metainfo schema is defined though MSection Python classes.
    nexus_metainfo_package.init_metainfo()

    # We skip the Python code generation for now and offer Python classes as variables
    # TO DO not necessary right now, could also be done case-by-case by the nexus parser
    python_module = sys.modules[__name__]
    for section in nexus_metainfo_package.section_definitions:  # pylint: disable=E1133
        setattr(python_module, section.name, section.section_cls)


init_nexus_metainfo()
