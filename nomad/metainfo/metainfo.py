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

from typing import Type, TypeVar, Union, Tuple, Iterable, List, Any, Dict, cast
import sys


__module__ = sys.modules[__name__]
MObjectBound = TypeVar('MObjectBound', bound='MObject')

"""

Discussion:
-----------

"""


# Reflection

class Enum(list):
    pass


class MObjectMeta(type):

    def __new__(self, cls_name, bases, dct):
        cls = super().__new__(self, cls_name, bases, dct)
        init = getattr(cls, '__init_section_cls__')
        if init is not None:
            init()
        return cls


Content = Tuple[MObjectBound, Union[List[MObjectBound], MObjectBound], str, MObjectBound]
SectionDef = Union[str, 'Section', Type[MObjectBound]]


class MObject(metaclass=MObjectMeta):
    """Base class for all section objects on all meta-info levels.

    All metainfo objects instantiate classes that inherit from ``MObject``. Each
    section or quantity definition is an ``MObject``, each actual (meta-)data carrying
    section is an ``MObject``. This class consitutes the reflection interface of the
    meta-info, since it allows to manipulate sections (and therefore all meta-info data)
    without having to know the specific sub-class.

    It also carries all the data for each section. All sub-classes only define specific
    sections in terms of possible sub-sections and quantities. The data is managed here.

    The reflection insterface for reading and manipulating quantity values consists of
    Pythons build in ``getattr``, ``setattr``, and ``del``, as well as member functions
    :func:`m_add_value`, and :func:`m_add_values`.

    Sub-sections and parent sections can be read and manipulated with :data:`m_parent`,
    :func:`m_sub_section`, :func:`m_create`.

    ```
    system = run.m_create(System)
    assert system.m_parent == run
    assert run.m_sub_section(System, system.m_parent_index) == system
    ```

    Attributes:
        m_section: The section definition that defines this sections, its possible
            sub-sections and quantities.
        m_parent: The parent section instance that this section is a sub-section of.
        m_parent_index: For repeatable sections, parent keep a list of sub-sections for
            each section definition. This is the index of this section in the respective
            parent sub-section list.
        m_data: The dictionary that holds all data of this section. It keeps the quantity
            values and sub-section. It should only be read directly (and never manipulated)
            if you are know what you are doing. You should always use the reflection interface
            if possible.
    """

    def __init__(self, m_section: 'Section' = None, m_parent: 'MObject' = None, **kwargs):
        self.m_section: 'Section' = m_section
        self.m_parent: 'MObject' = m_parent
        self.m_parent_index = -1
        self.m_data = dict(**kwargs)

        if self.m_section is None:
            self.m_section = getattr(self.__class__, 'm_section', None)
        else:
            assert self.m_section == getattr(self.__class__, 'm_section', self.m_section), \
                'Section class and section definition must match'

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
        m_section.section_cls = cls

        for name, value in cls.__dict__.items():
            if isinstance(value, Quantity):
                value.name = name
                # manual manipulation of m_data due to bootstrapping
                m_section.m_data.setdefault('Quantity', []).append(value)

    @staticmethod
    def __type_check(definition: 'Quantity', value: Any, check_item: bool = False):
        """Checks if the value fits the given quantity in type and shape; raises
        ValueError if not."""

        def check_value(value):
            if isinstance(definition.type, Enum):
                if value not in definition.type:
                    raise ValueError('Not one of the enum values.')

            elif isinstance(definition.type, type):
                if not isinstance(value, definition.type):
                    raise ValueError('Value has wrong type.')

            elif isinstance(definition.type, Section):
                if not isinstance(value, MObject) or value.m_section != definition.type:
                    raise ValueError('The value is not a section of wrong section definition')

            else:
                raise Exception('Invalid quantity type: %s' % str(definition.type))

        shape = None
        try:
            shape = definition.shape
        except KeyError:
            pass

        if shape is None or len(shape) == 0 or check_item:
            check_value(value)

        elif len(shape) == 1:
            if not isinstance(value, list):
                raise ValueError('Wrong shape')

            for item in value:
                check_value(item)

        else:
            # TODO
            raise Exception('Higher shapes not implemented')

        # TODO check dimension

    def __resolve_section(self, definition: SectionDef) -> 'Section':
        """Resolves and checks the given section definition. """
        if isinstance(definition, str):
            section = self.m_section.sub_sections[definition]

        else:
            if isinstance(definition, type):
                section = getattr(definition, 'm_section')
            else:
                section = definition
            if section.name not in self.m_section.sub_sections:
                raise KeyError('Not a sub section.')

        return section

    def m_sub_section(self, definition: SectionDef, parent_index: int = -1) -> MObjectBound:
        """Returns the sub section for the given section definition and possible
           parent_index (for repeatable sections).

        Args:
            definition: The definition of the section.
            parent_index: The index of the desired section. This can be omitted for non
                repeatable sections. If omitted for repeatable sections a exception
                will be raised, if more then one sub-section exists. Likewise, if the given
                index is out of range.
        Raises:
            KeyError: If the definition is not for a sub section
            IndexError: If the given index is wrong, or if an index is given for a non
                repeatable section
        """
        section_def = self.__resolve_section(definition)

        m_data_value = self.m_data[section_def.name]

        if isinstance(m_data_value, list):
            m_data_values = m_data_value
            if parent_index == -1:
                if len(m_data_values) == 1:
                    return m_data_values[0]
                else:
                    raise IndexError()
            else:
                return m_data_values[parent_index]
        else:
            if parent_index != -1:
                raise IndexError('Not a repeatable sub section.')
            else:
                return m_data_value

    def m_create(self, definition: SectionDef, **kwargs) -> MObjectBound:
        """Creates a subsection and adds it this this section

        Args:
            section: The section definition of the subsection. It is either the
                definition itself, or the python class representing the section definition.
            **kwargs: Are used to initialize the subsection.

        Returns:
            The created subsection

        Raises:
            KeyError: If the given section is not a subsection of this section.
        """
        section_def: 'Section' = self.__resolve_section(definition)

        section_cls = section_def.section_cls
        section_instance = section_cls(m_section=section_def, m_parent=self, **kwargs)

        if section_def.repeats:
            m_data_sections = self.m_data.setdefault(section_def.name, [])
            section_index = len(m_data_sections)
            m_data_sections.append(section_instance)
            section_instance.m_parent_index = section_index
        else:
            self.m_data[section_def.name] = section_instance

        return cast(MObjectBound, section_instance)

    def __resolve_quantity(self, definition: Union[str, 'Quantity']) -> 'Quantity':
        """Resolves and checks the given quantity definition. """
        if isinstance(definition, str):
            quantity = self.m_section.quantities[definition]

        else:
            if definition.m_parent != self.m_section:
                raise KeyError('Quantity is not a quantity of this section.')
            quantity = definition

        return quantity

    def m_add(self, definition: Union[str, 'Quantity'], value: Any):
        """Adds the given value to the given quantity."""

        quantity = self.__resolve_quantity(definition)

        MObject.__type_check(quantity, value, check_item=True)

        m_data_values = self.m_data.setdefault(quantity.name, [])
        m_data_values.append(value)

    def m_add_values(self, definition: Union[str, 'Quantity'], values: Iterable[Any]):
        """Adds the given values to the given quantity."""

        quantity = self.__resolve_quantity(definition)

        for value in values:
            MObject.__type_check(quantity, value, check_item=True)

        m_data_values = self.m_data.setdefault(quantity.name, [])
        for value in values:
            m_data_values.append(value)

    def m_to_dict(self) -> Dict[str, Any]:
        """Returns the data of this section as a json serializeable dictionary. """
        pass

    def m_to_json(self):
        """Returns the data of this section as a json string. """
        pass

    def m_all_contents(self) -> Iterable[Content]:
        """Returns an iterable over all sub and sub subs sections. """
        for content in self.m_contents():
            for sub_content in content[0].m_all_contents():
                yield sub_content

            yield content

    def m_contents(self) -> Iterable[Content]:
        """Returns an iterable over all direct subs sections. """
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


# M3

class Quantity(MObject):
    m_section: 'Section' = None
    name: 'Quantity' = None
    type: 'Quantity' = None
    shape: 'Quantity' = None

    __name = property(lambda self: self.m_data['name'])

    default = property(lambda self: None)

    def __get__(self, obj, type=None):
        return obj.m_data[self.__name]

    def __set__(self, obj, value):
        MObject.__dict__['_MObject__type_check'].__get__(MObject)(self, value)
        obj.m_data[self.__name] = value

    def __delete__(self, obj):
        del obj.m_data[self.__name]


class Section(MObject):
    m_section: 'Section' = None
    section_cls: Type[MObject] = None
    name: 'Quantity' = None
    repeats: 'Quantity' = None
    parent: 'Quantity' = None
    extends: 'Quantity' = None

    __all_instances: List['Section'] = []

    default = property(lambda self: [] if self.repeats else None)

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
Section.m_section.section_cls = Section

Section.name = Quantity(type=str, name='name')
Section.repeats = Quantity(type=bool, name='repeats')
Section.parent = Quantity(type=Section.m_section, name='parent')
Section.extends = Quantity(type=Section.m_section, shape=['0..*'], name='extends')

Quantity.m_section = Section(repeats=True, parent=Section.m_section, name='Quantity')
Quantity.m_section.section_cls = Quantity
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
