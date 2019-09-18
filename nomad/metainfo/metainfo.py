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

from typing import Type, TypeVar, Union, Tuple, Iterable, List, Any, Dict
import sys


__module__ = sys.modules[__name__]
MObjectBound = TypeVar('MObjectBound', bound='MObject')


# Reflection

class MObjectMeta(type):

    def __new__(self, cls_name, bases, dct):
        cls = super().__new__(self, cls_name, bases, dct)
        init = getattr(cls, '__init_section_cls__')
        if init is not None:
            init()
        return cls


Content = Tuple[MObjectBound, Union[List[MObjectBound], MObjectBound], str, MObjectBound]


class MObject(metaclass=MObjectMeta):

    def __init__(self, **kwargs):
        self.m_section: 'Section' = kwargs.pop('m_section', getattr(self.__class__, 'm_section', None))  # TODO test with failure instead of None
        self.m_data = dict(**kwargs)

    @classmethod
    def __init_section_cls__(cls):
        if not hasattr(__module__, 'Quantity') or not hasattr(__module__, 'Section'):
            # no initialization during bootstrapping, will be done maunally
            return

        m_section = getattr(cls, 'm_section', None)
        if m_section is None:
            m_section = Section()
            setattr(cls, 'm_section', m_section)
        m_section.name = cls.__name__

        for name, value in cls.__dict__.items():
            if isinstance(value, Quantity):
                value.name = name
                m_section.m_add('Quantity', value)

    def m_create(self, section: Union['Section', Type[MObjectBound]], **kwargs) -> MObjectBound:
        """Creates a subsection and adds it this this section

        Args:
            section: The section definition of the subsection. It is either the
                definition itself, or the python class representing the section definition.
            **kwargs: Are used to initialize the subsection.

        Returns:
            The created subsection

        Raises:
            ValueError: If the given section is not a subsection of this section, or
                this given definition is not a section at all.
        """
        section_def: 'Section' = None
        if isinstance(section, type) and hasattr(section, 'm_section'):
            section_def = section.m_section
        elif isinstance(section, Section):
            section_def = section
        else:
            raise ValueError('Not a section definition')

        if section_def.parent != self.m_section:
            raise ValueError('Not a subsection')

        section_cls = section_def.section_cls
        section_instance = section_cls(**kwargs)

        if section_def.repeats:
            self.m_data.setdefault(section_def.name, []).append(section_instance)
        else:
            self.m_data[section_def.name] = section_instance

        return section_instance

    def m_set(self, name: str, value: Any):
        self.m_data[name] = value

    def m_add(self, name: str, value: Any):
        values = self.m_data.setdefault(name, [])
        values.append(value)

    def m_get(self, name: str) -> Any:
        try:
            return self.m_data[name]

        except KeyError:
            return self.m_section.attributes[name].default

    def m_delete(self, name: str):
        del self.m_data[name]

    def m_to_dict(self):
        pass

    def m_to_json(self):
        pass

    def m_all_contents(self) -> Iterable[Content]:
        for content in self.m_contents():
            for sub_content in content[0].m_all_contents():
                yield sub_content

            yield content

    def m_contents(self) -> Iterable[Content]:
        for name, attr in self.m_data.items():
            if isinstance(attr, list):
                for value in attr:
                    if isinstance(value, MObject):
                        yield value, attr, name, self

            elif isinstance(attr, MObject):
                yield value, value, name, self

    def __repr__(self):
        m_section_name = self.m_section.name
        name = ''
        if 'name' in self.m_data:
            name = self.m_data['name']

        return '%s:%s' % (name, m_section_name)


class Enum(list):
    pass


# M3

class Quantity(MObject):
    m_section: 'Section' = None
    name: 'Quantity' = None
    type: 'Quantity' = None
    shape: 'Quantity' = None

    __name = property(lambda self: self.m_data['name'])

    default = property(lambda self: None)

    def __get__(self, obj, type=None):
        return obj.m_get(self.__name)

    def __set__(self, obj, value):
        obj.m_set(self.__name, value)

    def __delete__(self, obj):
        obj.m_delete(self.__name)


class Section(MObject):
    m_section: 'Section' = None
    name: 'Quantity' = None
    repeats: 'Quantity' = None
    parent: 'Quantity' = None
    extends: 'Quantity' = None

    __all_instances: List['Section'] = []

    default = property(lambda self: [] if self.repeats else None)
    section_cls = property(lambda self: self.__class__)

    def __init__(self, **kwargs):
        # The mechanism that produces default values, depends on parent. Without setting
        # the parent default manually, an endless recursion will occur.
        kwargs.setdefault('parent', None)

        super().__init__(**kwargs)
        Section.__all_instances.append(self)

    # TODO cache
    @property
    def attributes(self) -> Dict[str, Union['Section', Quantity]]:
        """ All attribute (sub section and quantity) definitions. """

        attributes: Dict[str, Union[Section, Quantity]] = dict(**self.quantities)
        attributes.update(**self.sub_sections)
        return attributes

    # TODO cache
    @property
    def quantities(self) -> Dict[str, Quantity]:
        """ All quantity definition in the given section definition. """

        return {
            quantity.name: quantity
            for quantity in self.m_data.get('Quantity', [])}

    # TODO cache
    @property
    def sub_sections(self) -> Dict[str, 'Section']:
        """ All sub section definitions for this section definition. """

        return {
            sub_section.name: sub_section
            for sub_section in Section.__all_instances
            if sub_section.parent == self}


Section.m_section = Section(repeats=True, name='Section')
Section.m_section.m_section = Section.m_section

Section.name = Quantity(type=str, name='name')
Section.repeats = Quantity(type=bool, name='repeats')
Section.parent = Quantity(type=Section.m_section, name='parent')
Section.extends = Quantity(type=Section.m_section, shape=['0..*'], name='extends')

Quantity.m_section = Section(repeats=True, parent=Section.m_section, name='Quantity')
Quantity.name = Quantity(type=str, name='name')
Quantity.type = Quantity(type=Union[type, Enum, Section], name='type')
Quantity.shape = Quantity(type=Union[str, int], shape=['0..*'], name='shape')


class Package(MObject):
    m_section = Section()
    name = Quantity(type=str)


Section.m_section.parent = Package.m_section


class Definition(MObject):
    m_section = Section(extends=[Section.m_section, Quantity.m_section, Package.m_section])

    description = Quantity(type=str)
