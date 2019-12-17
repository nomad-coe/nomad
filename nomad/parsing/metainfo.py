# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from typing import List, Dict, Any
import json

from nomadcore.local_backend import ParserEvent

from nomad import utils
from nomad.metainfo import SubSection, Section, MSection, Quantity, MResource, Reference
from nomad.metainfo.legacy import LegacyMetainfoEnvironment

from .backend import LegacyParserBackend

# TODO nonOverlappingSection should reopen (?)
# TODO non repeat section should reopen existing (?)
# TODO references are set by index in old metainfo


class MetainfoBackend(LegacyParserBackend):
    """ A backend that uses the new metainfo to store all data. """

    def __init__(self, env: LegacyMetainfoEnvironment, logger=None):
        super().__init__(logger=logger)

        self.env = env.env
        self.legacy_env = env.legacy_info_env()
        self.resource = MResource()

        self.open_sections_by_def: Dict[Section, List[MSection]] = {}
        if logger is None:
            logger = utils.get_logger(__name__)
        self.logger = logger

    def metaInfoEnv(self):
        return self.legacy_env

    def resolve_context(self, context_uri: str):
        path = context_uri.strip('/').split('/')
        path.reverse()
        current = None
        while len(path) > 0:
            section = path.pop()
            index = int(path.pop())
            if current is None:
                section_def = self.env.resolve_definition(section, Section)
                if section_def is None:
                    raise KeyError
                current = self.resource.contents[index]
                if current.m_def != section_def:
                    raise KeyError
            else:
                sub_section_def = self.env.resolve_definition(section, SubSection)
                if sub_section_def is None:
                    raise KeyError
                current = current.m_get_sub_section(sub_section_def, index)

        return current

    def openContext(self, context_uri: str):
        """ Open existing archive data to introduce new data into an existing section. """
        resolved = self.resolve_context(context_uri)
        self.open_sections_by_def.setdefault(resolved.m_def, []).append(resolved)

    def closeContext(self, context_uri: str):
        """ Close priorly opened existing archive data again. """
        resolved = self.resolve_context(context_uri)
        self.open_sections_by_def.setdefault(resolved.m_def, []).remove(resolved)

    def openSection(self, name):
        """
        It will assume that there is a sub-section def with the given name.
        It will use the latest opened section of the sub-sections parent as the parent
        for the new section.
        An Exception will be known root sections, e.g. 'section_run'.
        """
        if name in ['section_run', 'section_entry_info']:
            section_def = self.env.resolve_definition(name, Section)
            sub_section = self.resource.create(section_def.section_cls)

        else:
            sub_section_def = self.env.resolve_definition(name, SubSection)
            parent_section_def = sub_section_def.m_parent_as(Section)
            section_def = sub_section_def.sub_section
            open_parent_sections = self.open_sections_by_def[parent_section_def]

            if len(open_parent_sections) == 0:
                raise KeyError('There is no open parent section %s' % parent_section_def)

            if len(open_parent_sections) > 1:
                raise KeyError('There is more then one open parent section %s' % parent_section_def)

            sub_section = open_parent_sections[-1].m_create(
                section_def.section_cls)

        self.open_sections_by_def.setdefault(section_def, []).append(sub_section)
        return sub_section.m_parent_index

    def get_open_section_for_quantity(self, name, g_index):
        """ Returns the open section that contains the quantity of the given name. """
        quantity_def = self.env.resolve_definition(name, Quantity)
        section_def = quantity_def.m_parent_as(Section)
        sections = self.open_sections_by_def.get(section_def, [])
        if g_index == -1:
            if len(sections) != 1:
                raise KeyError('Section %s with g_index %d not open' % (name, g_index))
            return sections[0], quantity_def
        else:
            for section in sections:
                if section.m_parent_index == g_index:
                    return section, quantity_def

        raise KeyError('Section %s with g_index %d not open' % (name, g_index))

    def closeSection(self, name, g_index):
        if name in ['section_run', 'section_entry_info']:
            section_def = self.env.resolve_definition(name, Section)

        else:
            sub_section_def = self.env.resolve_definition(name, SubSection)
            section_def = sub_section_def.sub_section

        sub_sections = self.open_sections_by_def.get(section_def, [])
        if g_index == -1:
            if len(sub_sections) != 1:
                raise KeyError('Section %s with g_index %d not open' % (name, g_index))
            del(sub_sections[0])
        else:
            to_delete = None
            for sub_section in sub_sections:
                if sub_section.m_parent_index == g_index:
                    to_delete = sub_section
            if to_delete is None:
                raise KeyError('Section %s with g_index %d not open' % (name, g_index))
            sub_sections.remove(to_delete)

    def openNonOverlappingSection(self, metaName):
        return self.openSection(metaName)

    def setSectionInfo(self, metaName, gIndex, references):
        """
        Sets info values of an open section references should be a dictionary with the
        gIndexes of the root sections this section refers to.
        """
        # TODO might be necessary to make references work?
        pass

    def closeNonOverlappingSection(self, name):
        return self.closeSection(name, -1)

    def openSections(self):
        """ Returns the sections that are still open as metaName, gIndex tuples. """
        for section_def, sub_sections in self.open_sections_by_def:
            for sub_section in sub_sections:
                yield section_def.name, sub_section.m_parent_index

    def addValue(self, name, value, g_index=-1):
        section, quantity_def = self.get_open_section_for_quantity(name, g_index)
        if isinstance(quantity_def.type, Reference):
            # quantity is a reference
            possible_targets = self.resource.all(
                quantity_def.type.target_section_def.section_cls)
            referenced_target = None
            for target in possible_targets:
                if target.m_parent_index == value:
                    referenced_target = target

            if referenced_target is None:
                raise ValueError('There is not section for the given reference index')

            value = referenced_target

        setattr(section, name, value)

    def addRealValue(self, name, value, g_index=-1):
        self.addValue(name, value, g_index)

    def addArray(self, name, shape, g_index=-1):
        """
        Adds an unannitialized array of the given shape for the given metaName.
        The gIndex is used to identify the right parent section.
        This is neccessary before array values can be set with :func:`setArrayValues`.
        """
        raise NotImplementedError()

    def setArrayValues(self, metaName, values, offset=None, gIndex=-1):
        """
        Adds values of the given numpy array to the last array added for the given
        metaName and parent gIndex.
        """
        raise NotImplementedError()

    def addArrayValues(self, name, values, g_index=-1, override: bool = False):
        """
        Adds an array with the given numpy array values for the given metaName and
        parent section gIndex. Override determines whether to rewrite exisiting values
        in the backend.
        """
        section, quantity_def = self.get_open_section_for_quantity(name, g_index)
        if isinstance(quantity_def.type, Reference):
            # quantity is a reference
            possible_targets = self.resource.all(
                quantity_def.type.target_section_def.section_cls)
            resolved_values = []
            for value in values:
                referenced_target = None
                for target in possible_targets:
                    if target.m_parent_index == value:
                        referenced_target = target

                if referenced_target is None:
                    raise ValueError('There is not section for the given reference index')

                resolved_values.append(referenced_target)
            values = resolved_values

        if not override:
            quantity_def = section.m_def.all_quantities[name]
            assert not section.m_is_set(quantity_def)
        setattr(section, name, values)

    # The following are extensions to the origin NOMAD-coe parser backend. And allow
    # access to existing data

    @property
    def data(self):
        raise NotImplementedError(
            'This method does not make sense in the context of the new metainfo.')

    def get_sections(self, meta_name: str, g_index: int = -1) -> List[int]:
        """ Return all gIndices for existing sections of the given meta_name and parent index. """
        section_def = self.env.resolve_definition(meta_name, Section)
        return [
            section.m_parent_index for section in self.resource.all(section_def.section_cls)
            if g_index == -1 or section.m_parent.m_parent_index == g_index]

    def get_value(self, meta_name: str, g_index=-1) -> Any:
        """
        Return the value set to the given meta_name in its parent section of the given index.
        An index of -1 (default) is only allowed if there is exactly one parent section.
        """
        try:
            quantity = self.env.resolve_definition(meta_name, Quantity)
        except KeyError:
            return None

        section_def = quantity.m_parent_as(Section)
        value = None

        for section in self.resource.all(section_def.section_cls):
            if section.m_parent_index == g_index:
                value = section.m_get(quantity)

            if section_def.name in ['section_run', 'section_entry_info']:
                value = section.m_get(quantity)

        if isinstance(quantity.type, Reference):
            return value.m_parent_index

        return value

    def add_mi2_section(self, section: MSection):
        self.resource.add(section)

    def get_mi2_section(self, section_def: Section):
        for content in self.resource.contents:
            if content.m_def == section_def:
                return content

    def write_json(self, out, pretty: bool = True, *args, **kwargs):
        if pretty:
            kwargs = dict(indent=2)
        else:
            kwargs = dict()
        json.dump([section.m_to_dict() for section in self.resource.contents], out, **kwargs)

    def traverse(self, *args, **kwargs):
        for content in self.resource.contents:
            for section in content.m_all_contents():
                yield section.m_def.name, ParserEvent.open_section, None
                for quantity_def in section.m_def.all_quantities.values():
                    if section.m_is_set(quantity_def):
                        if len(quantity_def.shape) == 0:
                            event = ParserEvent.add_value
                        else:
                            event = ParserEvent.add_array_value
                        yield quantity_def.name, event, section.m_get(quantity_def)
                yield section.m_def.name, ParserEvent.close_section, None
