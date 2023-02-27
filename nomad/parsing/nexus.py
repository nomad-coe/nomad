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
import os.path
import xml.etree.ElementTree as ET
from typing import Optional

import numpy as np

from nexusutils.nexus import nexus as read_nexus  # pylint: disable=import-error
from nomad.datamodel import EntryArchive
from nomad.metainfo import MSection, nexus
from nomad.metainfo.util import MQuantity, resolve_variadic_name
from nomad.parsing import MatchingParser
from nomad.units import ureg
from nomad.utils import get_logger


def _to_group_name(nx_node: ET.Element):
    '''
    Normalise the given group name
    '''
    return nx_node.attrib.get('name', nx_node.attrib['type'][2:].upper())


# noinspection SpellCheckingInspection
def _to_section(
        hdf_name: Optional[str], nx_def: str, nx_node: Optional[ET.Element],
        current: MSection) -> MSection:
    '''
    Args:
        hdf_name : name of the hdf group/field/attribute (None for definition)
        nx_def : application definition
        nx_node : node in the nxdl.xml
        current : current section in which the new entry needs to be picked up from

    Note that if the new element did not exist, it will be created

    Returns:
        tuple: the new subsection

    The strict mapping is available between metainfo and nexus:
        Group <-> SubSection
        Field <-> Quantity
        Attribute <-> SubSection.Attribute or Quantity.Attribute

    If the given nxdl_node is a Group, return the corresponding Section.
    If the given nxdl_node is a Field, return the Section contains it.
    If the given nxdl_node is an Attribute, return the associated Section or the
    Section contains the associated Quantity.
    '''

    if hdf_name is None:
        nomad_def_name = nx_def
    elif nx_node.tag.endswith('group'):
        # it is a new group
        nomad_def_name = _to_group_name(nx_node)
    else:
        # no need to change section for quantities and attributes
        return current

    # for groups, get the definition from the package
    new_def = current.m_def.all_sub_sections[nomad_def_name]

    new_section: MSection = None  # type:ignore

    for section in current.m_get_sub_sections(new_def):
        if hdf_name is None or getattr(section, 'nx_name', None) == hdf_name:
            new_section = section
            break

    if new_section is None:
        current.m_create(new_def.section_def.section_cls)
        new_section = current.m_get_sub_section(new_def, -1)
        new_section.__dict__['nx_name'] = hdf_name

    return new_section


def _get_value(hdf_node):
    '''
    Get value from hdl5 node
    '''

    hdf_value = hdf_node[...]
    if str(hdf_value.dtype) == 'bool':
        val = bool(hdf_value)
    elif hdf_value.dtype.kind in 'iufc':
        val = hdf_value
    else:
        try:
            val = str(hdf_value.astype(str))
        except UnicodeDecodeError:
            val = str(hdf_node[()].decode())
    return val


class NexusParser(MatchingParser):
    '''
    NexusParser doc
    '''

    def __init__(self):
        super().__init__(
            metadata_path=os.path.join(os.path.dirname(__file__), 'metadata.yaml'),
            mainfile_mime_re=r'(application/.*)|(text/.*)',
            mainfile_name_re=r'.*\.nxs',
            supported_compressions=['gz', 'bz2', 'xz']
        )
        self.archive: Optional[EntryArchive] = None
        self.nx_root = None
        self._logger = None

    def _populate_data(self, depth: int, nx_path: list, nx_def: str, hdf_node, current: MSection):
        '''
        Populate attributes and fields
        '''

        if depth < len(nx_path):
            # it is an attribute of either field or group
            nx_attr = nx_path[depth]
            nx_parent: ET.Element = nx_path[depth - 1]

            if isinstance(nx_attr, str):
                if nx_attr != 'units':
                    # no need to handle units here
                    # as all quantities have flexible units
                    print(nx_attr)
            else:
                # get the name of parent (either field or group)
                # which will be used to set attribute
                # this is required by the syntax of metainfo mechanism
                # due to variadic/template quantity names
                parent_type = nx_parent.get('type').replace('NX', '').upper()
                parent_name = nx_parent.get('name', parent_type)  # type: ignore

                attr_name = nx_attr.get('name')
                # It could be 1D array, float or int
                attr_value = hdf_node.attrs[attr_name]
                if not isinstance(attr_value, str):
                    if isinstance(attr_value, np.ndarray):
                        attr_value = [value for value in attr_value]
                        if len(attr_value) == 1:
                            attr_value = attr_value[0]

                attr_name = attr_name + "__attribute"
                current = _to_section(attr_name, nx_def, nx_attr, current)

                try:
                    if nx_parent.tag.endswith('group'):
                        current.m_set_section_attribute(attr_name, attr_value)
                    else:
                        current.m_set_quantity_attribute(parent_name + "__field", attr_name, attr_value)
                except Exception as e:
                    self._logger.warning('Error while setting attribute.', target_name=attr_name, exe_info=str(e))
        else:
            # it is a field
            field = _get_value(hdf_node)

            # need to remove
            if hdf_node[...].dtype.kind in 'iufc' and isinstance(field, np.ndarray) and field.size > 1:
                field = np.array([np.mean(field), np.var(field), np.min(field), np.max(field)])

            # get the corresponding field name
            field_name = nx_path[-1].get('name') + "__field"
            metainfo_def = resolve_variadic_name(current.m_def.all_properties, field_name)

            # check if unit is given
            unit = hdf_node.attrs.get('units', None)

            pint_unit: Optional[ureg.Unit] = None

            if unit:
                try:
                    if unit != 'counts':
                        pint_unit = ureg.parse_units(unit)
                    else:
                        pint_unit = ureg.parse_units('1')
                    field = ureg.Quantity(field, pint_unit)
                except ValueError:
                    pass

            if metainfo_def.use_full_storage:
                field = MQuantity.wrap(field, hdf_node.name.split('/')[-1] + '__field')
            elif metainfo_def.unit is None and pint_unit is not None:
                metainfo_def.unit = pint_unit

            # may need to check if the given unit is in the allowable list

            try:
                current.m_set(metainfo_def, field)
            except Exception as e:
                self._logger.warning('Error while setting field.', target_name=field_name, exe_info=str(e))

    def __nexus_populate(self, params: dict, attr=None):  # pylint: disable=W0613
        '''
        Walks through name_list and generate nxdl nodes
        (hdf_info, nx_def, nx_path, val, logger) = params
        '''

        hdf_info: dict = params['hdf_info']
        nx_def: str = params['nxdef']
        nx_path: list = params['nxdl_path']

        hdf_path: str = hdf_info['hdf_path']
        hdf_node = hdf_info['hdf_node']

        if nx_path is None:
            return

        current: MSection = _to_section(None, nx_def, None, self.nx_root)
        depth: int = 1
        for name in hdf_path.split('/')[1:]:
            nx_node = nx_path[depth] if depth < len(nx_path) else name
            current = _to_section(name, nx_def, nx_node, current)
            depth += 1

        self._populate_data(depth, nx_path, nx_def, hdf_node, current)

    def parse(self, mainfile: str, archive: EntryArchive, logger=None, child_archives=None):
        self.archive = archive
        self.archive.m_create(nexus.NeXus)  # type: ignore # pylint: disable=no-member
        self.nx_root = self.archive.nexus
        self._logger = logger if logger else get_logger(__name__)

        nexus_helper = read_nexus.HandleNexus(logger, [mainfile])
        nexus_helper.process_nexus_master_file(self.__nexus_populate)

        if archive.metadata is None:
            return

        app_def: str = ''
        for var in dir(archive.nexus):
            if getattr(archive.nexus, var, None) is not None:
                app_def = var

        archive.metadata.entry_type = app_def
