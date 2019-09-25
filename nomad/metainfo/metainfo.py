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

"""
The NOMAD meta-info allows to define physics data quantities. These definitions are
necessary for all computer representations of respective data (e.g. in Python,
search engines, data-bases, and files).

This modules provides various Python interfaces for

- defining meta-info data
- to create and manipulate data that follows these definitions
- to (de-)serialize meta-info data in JSON (i.e. represent data in JSON formatted files)

Here is a simple example that demonstrates the definition of System related quantities:

.. code-block:: python

    class Run(MObject):
        pass

    class System(MObject):
        \"\"\"
        A system section includes all quantities that describe a single a simulated
        system (a.k.a. geometry).
        \"\"\"

        m_section = Section(repeats=True, parent=Run)

        n_atoms = Quantity(
            type=int, description='''
            A Defines the number of atoms in the system.
            ''')

        atom_labels = Quantity(type=Enum(ase.data.chemical_symbols), shape['n_atoms'])
        atom_positions = Quantity(type=float, shape=['n_atoms', 3], unit=Units.m)
        simulation_cell = Quantity(type=float, shape=[3, 3], unit=Units.m)
        pbc = Quantity(type=bool, shape=[3])

Here, we define a `section` called ``System``. The section mechanism allows to organize
related data into, well, sections. Sections form containment hierarchies. Here
containment is a parent-child (whole-part) relationship. In this example many ``Systems``,
are part of one ``Run``. Each ``System`` can contain values for the defined quantities:
``n_atoms``, ``atom_labels``, ``atom_positions``, ``simulation_cell``, and ``pbc``.
Quantities allow to state type, shape, and physics unit to specify possible quantity
values.

Here is an example, were we use the above definition to create, read, and manipulate
data that follows these definitions:

.. code-bock:: python

    run = Run()
    system = run.m_create(System)
    system.n_atoms = 3
    system.atom_labels = ['H', 'H', 'O']

    print(system.atom_labels)
    print(run.m_to_json(ident=2))

This last statement, will produce the following JSON:

.. code-block:: JSON

    {
        "m_section" = "Run",
        "System": [
            {
                "m_section" = "System",
                "m_parent_index" = 0,
                "n_atoms" = 3,
                "atom_labels" = [
                    "H",
                    "H",
                    "O"
                ]
            }
        ]
    }

This is the JSON representation, a serialized version of the Python representation in
the example above.

Sections can be extended with new quantities outside the original section definition.
This provides the key mechanism to extend commonly defined parts with (code) specific
quantities:

.. code-block:: Python

    class Method(nomad.metainfo.common.Method):
        x_vasp_incar_ALGO=Quantity(
            type=Enum(['Normal', 'VeryFast', ...]),
            links=['https://cms.mpi.univie.ac.at/wiki/index.php/ALGO'])
        \"\"\"
        A convenient option to specify the electronic minimisation algorithm (as of VASP.4.5)
        and/or to select the type of GW calculations.
        \"\"\"


All meta-info definitions and classes for meta-info data objects (i.e. section instances)
inherit from :class:` MObject`. This base-class provides common functions and attributes
for all meta-info data objects. Names of these common parts are prefixed with ``m_``
to distinguish them from user defined quantities. This also constitute's the `reflection`
interface (in addition to Python's build in ``getattr``, ``setattr``) that allows to
create and manipulate meta-info data, without prior program time knowledge of the underlying
definitions.

.. autoclass:: MObject

The following classes can be used to define and structure meta-info data:

- sections are defined by sub-classes :class:`MObject` and using :class:`Section` to
  populate the classattribute `m_section`
- quantities are defined by assigning classattributes of a section with :class:`Quantity`
  instances
- references (from one section to another) can be defined with quantities that use
  section definitions as type
- dimensions can use defined by simply using quantity names in shapes
- categories (former `abstract type definitions`) can be given in quantity definitions
  to assign quantities to additional specialization-generalization hierarchies

See the reference of classes :class:`Section` and :class:`Quantities` for details.

.. autoclass:: Section
.. autoclass:: Quantity
"""

from typing import Type, TypeVar, Union, Tuple, Iterable, List, Any, Dict, cast
import sys
import inspect
import re

from pint.unit import _Unit
from pint import UnitRegistry
import inflection

__module__ = sys.modules[__name__]
MObjectBound = TypeVar('MObjectBound', bound='MObject')


# Reflection

class Enum(list):
    """ Allows to define str types with values limited to a pre-set list of possible values. """
    pass


class DataType:
    """
    Allows to define custom data types that can be used in the meta-info.

    The metainfo supports most types out of the box. These includes the python build-in
    primitive types (int, bool, str, float, ...), references to sections, and enums.
    However, in some occasions you need to add custom data types.
    """
    def check_type(self, value):
        pass

    def normalize(self, value):
        return value

    def to_json_serializable(self, value):
        return value

    def from_json_serializable(self, value):
        return value


class Dimension(DataType):
    def check_type(self, value):
        if isinstance(value, int):
            return

        if isinstance(value, str):
            if value.isidentifier():
                return
            if re.match(r'(\d)\.\.(\d|\*)', value):
                return

        if isinstance(value, Section):
            return

        if isinstance(value, type) and hasattr(value, 'm_section'):
            return

        raise TypeError('%s is not a valid dimension' % str(value))
    # TODO


# TODO class Unit(DataType)
# TODO class MetainfoType(DataType)
# TODO class Datetime(DataType)


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

    .. code-block:: python

        system = run.m_create(System)
        assert system.m_parent == run
        assert run.m_sub_section(System, system.m_parent_index) == system

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

    m_section: 'Section' = None

    def __init__(self, m_section: 'Section' = None, m_parent: 'MObject' = None, _bs: bool = False, **kwargs):
        self.m_section: 'Section' = m_section
        self.m_parent: 'MObject' = m_parent
        self.m_parent_index = -1

        cls = self.__class__
        if self.m_section is None:
            self.m_section = cls.m_section

        if cls.m_section is not None:
            assert self.m_section == cls.m_section, \
                'Section class and section definition must match'

        self.m_data = dict(**kwargs)
        # TODO
        # self.m_data = {}
        # if _bs:
        #     self.m_data.update(**kwargs)
        # else:
        #     self.m_update(**kwargs)

    @classmethod
    def __init_section_cls__(cls):
        # only works after bootstrapping, since functionality is still missing
        if not all([hasattr(__module__, cls) for cls in ['Quantity', 'Section', 'sub_section']]):
            return

        # ensure that the m_section is defined
        m_section = cls.m_section
        if m_section is None and cls != MObject:
            m_section = Section()
            setattr(cls, 'm_section', m_section)

        # transfer name and description to m_section
        m_section.name = cls.__name__
        if cls.__doc__ is not None:
            m_section.description = inspect.cleandoc(cls.__doc__)
        m_section.section_cls = cls

        # add sub_section to parent section
        if m_section.parent is not None:
            sub_section_name = inflection.underscore(m_section.name)
            setattr(m_section.parent.section_cls, sub_section_name, sub_section(m_section))

        for name, attr in cls.__dict__.items():
            # transfer names and descriptions for quantities
            if isinstance(attr, Quantity):
                attr.name = name
                if attr.description is not None:
                    attr.description = inspect.cleandoc(attr.description)
                    attr.__doc__ = attr.description
                # manual manipulation of m_data due to bootstrapping
                m_section.m_data.setdefault('Quantity', []).append(attr)

            # set names and parent on sub-sections
            elif isinstance(attr, sub_section):
                attr.section_def.parent = m_section
                if attr.section_def.name is None:
                    attr.section_def.name = inflection.camelize(name)

    @staticmethod
    def m_type_check(definition: 'Quantity', value: Any, check_item: bool = False):
        """Checks if the value fits the given quantity in type and shape; raises
        TypeError if not."""

        if value is None and not check_item and definition.default is None:
            # Allow the default None value even if it would violate the type
            return

        def check_value(value):
            if isinstance(definition.type, Enum):
                if value not in definition.type:
                    raise TypeError('Not one of the enum values.')

            elif isinstance(definition.type, type):
                if not isinstance(value, definition.type):
                    raise TypeError('Value has wrong type.')

            elif isinstance(definition.type, Section):
                if not isinstance(value, MObject) or value.m_section != definition.type:
                    raise TypeError('The value is not a section of wrong section definition')

            else:
                # TODO
                # raise Exception('Invalid quantity type: %s' % str(definition.type))
                pass

        shape = None
        try:
            shape = definition.shape
        except KeyError:
            pass

        if shape is None or len(shape) == 0 or check_item:
            check_value(value)

        elif len(shape) == 1:
            if not isinstance(value, list):
                raise TypeError('Wrong shape')

            for item in value:
                check_value(item)

        else:
            # TODO
            raise Exception('Higher shapes not implemented')

        # TODO check dimension

    def _resolve_section(self, definition: SectionDef) -> 'Section':
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
        section_def = self._resolve_section(definition)

        m_data_value = self.m_data.get(section_def.name, None)

        if m_data_value is None:
            if section_def.repeats:
                m_data_value = []
            else:
                m_data_value = None

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

    def m_add_sub_section(self, sub_section: MObjectBound) -> MObjectBound:
        """Adds the given section instance as a sub section to this section."""

        section_def = sub_section.m_section

        if section_def.repeats:
            m_data_sections = self.m_data.setdefault(section_def.name, [])
            section_index = len(m_data_sections)
            m_data_sections.append(sub_section)
            sub_section.m_parent_index = section_index
        else:
            self.m_data[section_def.name] = sub_section

        return sub_section

    def m_create(self, definition: SectionDef, **kwargs) -> 'MObject':
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
        section_def: 'Section' = self._resolve_section(definition)

        section_cls = section_def.section_cls
        section_instance = section_cls(m_section=section_def, m_parent=self, **kwargs)

        return self.m_add_sub_section(section_instance)

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

        MObject.m_type_check(quantity, value, check_item=True)

        m_data_values = self.m_data.setdefault(quantity.name, [])
        m_data_values.append(value)

    def m_add_values(self, definition: Union[str, 'Quantity'], values: Iterable[Any]):
        """Adds the given values to the given quantity."""

        quantity = self.__resolve_quantity(definition)

        for value in values:
            MObject.m_type_check(quantity, value, check_item=True)

        m_data_values = self.m_data.setdefault(quantity.name, [])
        for value in values:
            m_data_values.append(value)

    def m_update(self, **kwargs):
        """ Updates all quantities and sub-sections with the given arguments. """
        for name, value in kwargs.items():
            attribute = self.m_section.attributes.get(name, None)
            if attribute is None:
                raise KeyError('%s is not an attribute of this section' % name)

            if isinstance(attribute, Section):
                if attribute.repeats:
                    if isinstance(value, List):
                        for item in value:
                            self.m_add_sub_section(item)
                    else:
                        raise TypeError('Sub section %s repeats, but no list was given' % attribute.name)
                else:
                    self.m_add_sub_section(item)

            else:
                setattr(self, name, value)

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


# M3, the definitions that are used to write definitions. These are the section definitions
# for sections Section and Quantity.They define themselves; i.e. the section definition
# for Section is the same section definition.
# Due to this circular nature (hen-egg-problem), the classes for sections Section and
# Quantity do only contain placeholder for their own section and quantity definitions.
# These placeholder are replaced, once the necessary classes are defined. This process
# is referred to as 'bootstrapping'.

class Definition(MObject):
    name: 'Quantity' = None
    description: 'Quantity' = None
    links: 'Quantity' = None


class Quantity(Definition):
    """Used to define quantities that store a certain piece of (meta-)data.

    Quantities are the basic building block with meta-info data. The Quantity class is
    used to define quantities within sections. A quantity definition
    is a (physics) quantity with name, type, shape, and potentially a unit.

    In Python terms, quantities are descriptors. Descriptors define how to get, set, and
    delete values for a object attribute. Meta-info descriptors ensure that
    type and shape fit the set values.
    """

    type: 'Quantity' = None
    shape: 'Quantity' = None
    unit: 'Quantity' = None
    default: 'Quantity' = None

    # TODO section = Quantity(type=Section), the section it belongs to
    # TODO synonym_for = Quantity(type=Quantity)
    # TODO derived_from = Quantity(type=Quantity, shape=['0..*'])
    # TODO categories = Quantity(type=Category, shape=['0..*'])
    # TODO converter = Quantity(type=Converter), a class with set of functions for
    #      normalizing, (de-)serializing values.

    # Some quantities of Quantity cannot be read as normal quantities due to bootstraping.
    # Those can be accessed internally through the following replacement properties that
    # read directly from m_data.
    __name = property(lambda self: self.m_data['name'])
    __default = property(lambda self: self.m_data.get('default', None))

    def __get__(self, obj, type=None):
        if obj is None:
            # class (def) attribute case
            return self

        # object (instance) attribute case
        try:
            return obj.m_data[self.__name]
        except KeyError:
            return self.__default

    def __set__(self, obj, value):
        if obj is None:
            # class (def) case
            raise KeyError('Cannot overwrite quantity definition. Only values can be set.')

        # object (instance) case
        MObject.m_type_check(self, value)
        obj.m_data[self.__name] = value

    def __delete__(self, obj):
        if obj is None:
            # class (def) case
            raise KeyError('Cannot delete quantity definition. Only values can be deleted.')

        # object (instance) case
        del obj.m_data[self.__name]


class Section(Definition):
    """Used to define section that organize meta-info data into containment hierarchies.

    Section definitions determine what quantities and sub-sections can appear in a section
    instance. A section instance itself can appear potentially many times in its parent
    section. See :data:`repeats` and :data:`parent`.

    In Python terms, sections are classes. Sub-sections and quantities are attribute of
    respective instantiating objects. For each section class there is a corresponding
    :class:`Section` instance that describes this class as a section. This instance
    is referred to as 'section definition' in contrast to the Python class that we call
    'section class'.
    """

    section_cls: Type[MObject] = None
    """ The section class that corresponse to this section definition. """

    repeats: 'Quantity' = None
    parent: 'Quantity' = None

    # TODO super = Quantity(type=Section, shape=['0..*']), inherit all quantity definition
    #      from the given sections, derived from Python base classes
    # TODO extends = Quantity(type=bool), denotes this section as a container for
    #      new quantities that belong to the base-class section definitions

    __all_instances: List['Section'] = []

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

    def add_quantity(self, quantity: Quantity):
        """
        Adds the given quantity to this section.

        Allows to add a quantity to a section definition outside the corresponding
        section class.

        .. code-block:: Python

        class System(MObject):
            pass

        System.m_section.add_quantity(Quantity(name='n_atoms', type=int))

        This will add the quantity definition to this section definition,
        and add the respective Python descriptor as an attribute to this class.
        """
        quantities = self.m_data.setdefault('Quantity', [])
        quantities.append(quantity)

        setattr(self.section_cls, quantity.name, quantity)


class Package(Definition):
    name: 'Quantity' = None
    """ The name of the package. """

    description: 'Quantity' = None
    """ A human readable description of this package. """


class sub_section:
    """ Allows to assign a section class as a sub-section to another section class. """

    def __init__(self, section: SectionDef, **kwargs):
        if isinstance(section, type):
            self.section_def = cast(MObject, section).m_section
        else:
            self.section_def = cast(Section, section)

    def __get__(self, obj: MObject, type=None) -> Union[MObject, Section]:
        if obj is None:
            # the class attribute case
            return self.section_def

        else:
            # the object attribute case
            m_data_value = obj.m_data.get(self.section_def.name, None)
            if m_data_value is None:
                if self.section_def.repeats:
                    m_data_value = []

            return m_data_value

    def __set__(self, obj: MObject, value: Union[MObject, List[MObject]]):
        raise NotImplementedError('Sub sections cannot be set directly. Use m_create.')

    def __delete__(self, obj):
        raise NotImplementedError('Sub sections cannot be deleted directly.')


# TODO Category(MObject), a definition kind for "abstract type declarations" and a
#      separate generalization/specialization/is-a hierarchy.


Section.m_section = Section(repeats=True, name='Section', _bs=True)
Section.m_section.m_section = Section.m_section
Section.m_section.section_cls = Section

Quantity.m_section = Section(repeats=True, parent=Section.m_section, name='Quantity', _bs=True)

Definition.name = Quantity(
    type=str, name='name', _bs=True, description='''
    The name of the quantity. Must be unique within a section.
    ''')
Definition.description = Quantity(
    type=str, name='description', _bs=True, description='''
    An optional human readable description.
    ''')
Definition.links = Quantity(
    type=str, shape=['0..*'], name='links', _bs=True, description='''
    A list of URLs to external resource that describe this definition.
    ''')

Section.repeats = Quantity(
    type=bool, name='repeats', default=False, _bs=True,
    description='''
    Wether instances of this section can occur repeatedly in the parent section.
    ''')
Section.parent = Quantity(
    type=Section.m_section, name='parent', _bs=True, description='''
    The section definition of parent sections. Only section instances of this definition
    can contain section instances of this definition.
    ''')

Quantity.m_section.section_cls = Quantity
Quantity.type = Quantity(
    type=Union[type, Enum, Section], name='type', _bs=True, description='''
    The type of the quantity.

    Can be one of the following:

    - none to support any value
    - a build-in primitive Python type, e.g. ``int``, ``str``
    - an instance of :class:`Enum`, e.g. ``Enum(['one', 'two', 'three'])
    - a instance of Section, i.e. a section definition. This will define a reference
    - a custom meta-info DataType

    In the NOMAD CoE meta-info this was basically the ``dTypeStr``.
    ''')
Quantity.shape = Quantity(
    type=Dimension, shape=['0..*'], name='shape', _bs=True, description='''
    The shape of the quantity that defines its dimensionality.

    A shape is a list, where each item defines a dimension. Each dimension can be:

    - an integer that defines the exact size of the dimension, e.g. ``[3]`` is the
      shape of a spacial vector
    - the name of an int typed quantity in the same section
    - a range specification as string build from a lower bound (i.e. int number),
      and an upper bound (int or ``*`` denoting arbitrary large), e.g. ``'0..*'``, ``'1..3'``
    ''')
Quantity.unit = Quantity(
    type=_Unit, _bs=True, description='''
    The optional physics unit for this quantity.

    Units are given in `pint` units. Pint is a Python package that defines units and
    their algebra. There is a default registry :data:`units` that you can use.
    Example units are: ``units.m``, ``units.m / units.s ** 2``.
    ''')
Quantity.default = Quantity(
    type=None, _bs=True, default=None, description='''
    The default value for this quantity.
    ''')

Package.m_section = Section(repeats=True, name='Package', _bs=True)
Package.m_section.parent = Package.m_section

Section.m_section.parent = Package.m_section

Package.__init_section_cls__()
Section.__init_section_cls__()
Quantity.__init_section_cls__()


units = UnitRegistry()
""" The default pint unit registry that should be used to give units to quantity definitions. """
