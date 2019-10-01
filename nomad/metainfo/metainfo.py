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

    class System(MSection):
        \"\"\"
        A system section includes all quantities that describe a single a simulated
        system (a.k.a. geometry).
        \"\"\"

        n_atoms = Quantity(
            type=int, description='''
            A Defines the number of atoms in the system.
            ''')

        atom_labels = Quantity(type=Enum(ase.data.chemical_symbols), shape['n_atoms'])
        atom_positions = Quantity(type=float, shape=['n_atoms', 3], unit=Units.m)
        simulation_cell = Quantity(type=float, shape=[3, 3], unit=Units.m)
        pbc = Quantity(type=bool, shape=[3])

    class Run(MSection):
        systems = SubSection(sub_section=System, repeats=True)

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
        "m_def" = "Run",
        "System": [
            {
                "m_def" = "System",
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
inherit from :class:` MSection`. This base-class provides common functions and properties
for all meta-info data objects. Names of these common parts are prefixed with ``m_``
to distinguish them from user defined quantities. This also constitute's the `reflection`
interface (in addition to Python's build in ``getattr``, ``setattr``) that allows to
create and manipulate meta-info data, without prior program time knowledge of the underlying
definitions.

.. autoclass:: MSection

The following classes can be used to define and structure meta-info data:

- sections are defined by sub-classes :class:`MSection` and using :class:`Section` to
  populate the classattribute `m_def`
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

# TODO validation

from typing import Type, TypeVar, Union, Tuple, Iterable, List, Any, Dict, cast
import sys
import inspect
import re
import json

import numpy as np
from pint.unit import _Unit
from pint import UnitRegistry

is_bootstrapping = True
MSectionBound = TypeVar('MSectionBound', bound='MSection')


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

        if isinstance(value, type) and hasattr(value, 'm_def'):
            return

        raise TypeError('%s is not a valid dimension' % str(value))
    # TODO


# TODO class Unit(DataType)
# TODO class MetainfoType(DataType)
# TODO class Datetime(DataType)


class MObjectMeta(type):

    def __new__(self, cls_name, bases, dct):
        cls = super().__new__(self, cls_name, bases, dct)
        init = getattr(cls, '__init_cls__')
        if init is not None and not is_bootstrapping:
            init()
        return cls


Content = Tuple[MSectionBound, Union[List[MSectionBound], MSectionBound], str, MSectionBound]

SectionDef = Union[str, 'Section', 'SubSection', Type[MSectionBound]]
""" Type for section definition references.

This can either be :

- the name of the section
- the section definition itself
- the definition of a sub section
- or the section definition Python class
"""


class MSection(metaclass=MObjectMeta):
    """Base class for all section instances on all meta-info levels.

    All metainfo objects instantiate classes that inherit from ``MSection``. Each
    section or quantity definition is an ``MSection``, each actual (meta-)data carrying
    section is an ``MSection``. This class consitutes the reflection interface of the
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
        m_def: The section definition that defines this sections, its possible
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

    m_def: 'Section' = None

    def __init__(self, m_def: 'Section' = None, m_parent: 'MSection' = None, _bs: bool = False, **kwargs):
        self.m_def: 'Section' = m_def
        self.m_parent: 'MSection' = m_parent
        self.m_parent_index = -1

        cls = self.__class__
        if self.m_def is None:
            self.m_def = cls.m_def

        if cls.m_def is not None:
            assert self.m_def == cls.m_def, \
                'Section class and section definition must match'

        self.m_annotations: Dict[str, Any] = {}
        self.m_data: Dict[str, Any] = {}
        for key, value in kwargs.items():
            if key.startswith('a_'):
                self.m_annotations[key[2:]] = value
            else:
                self.m_data[key] = value

        # TODO
        # self.m_data = {}
        # if _bs:
        #     self.m_data.update(**kwargs)
        # else:
        #     self.m_update(**kwargs)

    @classmethod
    def __init_cls__(cls):
        # ensure that the m_def is defined
        m_def = cls.m_def
        if m_def is None:
            m_def = Section()
            setattr(cls, 'm_def', m_def)

        # transfer name and description to m_def
        m_def.name = cls.__name__
        if cls.__doc__ is not None:
            m_def.description = inspect.cleandoc(cls.__doc__).strip()
        m_def.section_cls = cls

        for name, attr in cls.__dict__.items():
            # transfer names and descriptions for properties
            if isinstance(attr, Property):
                attr.name = name
                if attr.description is not None:
                    attr.description = inspect.cleandoc(attr.description).strip()
                    attr.__doc__ = attr.description

                # TODO manual manipulation of m_data due to bootstrapping
                if isinstance(attr, Quantity):
                    properties = m_def.m_data.setdefault('quantities', [])
                elif isinstance(attr, SubSection):
                    properties = m_def.m_data.setdefault('sub_sections', [])
                else:
                    raise NotImplementedError('Unknown property kind.')
                properties.append(attr)
                attr.m_parent = m_def
                attr.m_parent_index = len(properties) - 1

        # add section cls' section to the module's package
        module_name = cls.__module__
        pkg = Package.from_module(module_name)
        pkg.m_add_sub_section(cls.m_def)

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
                if not isinstance(value, MSection) or value.m_def != definition.type:
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

        else:
            if type(definition.type) == np.dtype:
                if len(shape) != len(value.shape):
                    raise TypeError('Wrong shape')
            else:
                if len(shape) == 1:
                    if not isinstance(value, list):
                        raise TypeError('Wrong shape')

                    for item in value:
                        check_value(item)

                else:
                    # TODO
                    # raise Exception('Higher shapes not implemented')
                    pass

        # TODO check dimension

    def _resolve_sub_section(self, definition: SectionDef) -> 'SubSection':
        """ Resolves and checks the given section definition. """

        if isinstance(definition, type):
            definition = getattr(definition, 'm_def', None)
            if definition is None:
                raise TypeError(
                    'The type/class %s is not definining a section, i.e. not derived from '
                    'MSection.' % str(definition))

        if isinstance(definition, Section):
            sub_section = self.m_def.all_sub_sections_by_section.get(definition, None)
            if sub_section is None:
                raise KeyError(
                    'The section %s is not a sub section of %s.' %
                    (definition.name, self.m_def.name))

        elif isinstance(definition, str):
            sub_section = self.m_def.all_sub_sections[definition]

        elif isinstance(definition, SubSection):
            sub_section = definition

        else:
            raise TypeError(
                '%s does not refer to a section definition. Either use the section '
                'definition, sub section definition, section class, or name.' %
                str(definition))

        if sub_section is None:
            raise KeyError(
                'The section %s is not a sub section of %s.' %
                (cast(Definition, definition).name, self.m_def.name))

        if sub_section.m_parent is not self.m_def:
            raise KeyError(
                'The section %s is not a sub section of %s.' %
                (cast(Definition, definition).name, self.m_def.name))

        return sub_section

    def m_sub_sections(self, definition: SectionDef) -> List[MSectionBound]:
        """Returns all sub sections for the given section definition

        Args:
            definition: The definition of the section.

        Raises:
            KeyError: If the definition is not for a sub section
        """
        sub_section = self._resolve_sub_section(definition)
        return getattr(self, sub_section.name)

    def m_sub_section(self, definition: SectionDef, parent_index: int = -1) -> MSectionBound:
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
        sub_section = self._resolve_sub_section(definition)

        m_data_value = getattr(self, sub_section.name)

        if m_data_value is None:
            if sub_section.repeats:
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

            return m_data_value

    def m_add_sub_section(self, sub_section: MSectionBound) -> MSectionBound:
        """Adds the given section instance as a sub section to this section."""

        sub_section_def = self._resolve_sub_section(sub_section.m_def.section_cls)
        sub_section.m_parent = self
        if sub_section_def.repeats:
            values = getattr(self, sub_section_def.name)
            sub_section.m_parent_index = len(values)
            values.append(sub_section)

        else:
            self.m_data[sub_section_def.name] = sub_section
            sub_section.m_parent_index = -1

        return sub_section

    # TODO this should work with the section constructor
    def m_create(self, definition: Type[MSectionBound], **kwargs) -> MSectionBound:
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
        sub_section: 'SubSection' = self._resolve_sub_section(definition)

        section_cls = sub_section.sub_section.section_cls
        section_instance = section_cls(m_def=section_cls.m_def, m_parent=self, **kwargs)

        return cast(MSectionBound, self.m_add_sub_section(section_instance))

    def __resolve_quantity(self, definition: Union[str, 'Quantity']) -> 'Quantity':
        """Resolves and checks the given quantity definition. """
        if isinstance(definition, str):
            quantity = self.m_def.all_quantities[definition]

        else:
            if definition.m_parent != self.m_def:
                raise KeyError('Quantity is not a quantity of this section.')
            quantity = definition

        return quantity

    def m_add(self, definition: Union[str, 'Quantity'], value: Any):
        """Adds the given value to the given quantity."""

        quantity = self.__resolve_quantity(definition)

        MSection.m_type_check(quantity, value, check_item=True)

        m_data_values = self.m_data.setdefault(quantity.name, [])
        m_data_values.append(value)

    def m_add_values(self, definition: Union[str, 'Quantity'], values: Iterable[Any]):
        """Adds the given values to the given quantity."""

        quantity = self.__resolve_quantity(definition)

        for value in values:
            MSection.m_type_check(quantity, value, check_item=True)

        m_data_values = self.m_data.setdefault(quantity.name, [])
        for value in values:
            m_data_values.append(value)

    def m_update(self, **kwargs):
        """ Updates all quantities and sub-sections with the given arguments. """
        for name, value in kwargs.items():
            prop = self.m_def.all_properties.get(name, None)
            if prop is None:
                raise KeyError('%s is not an attribute of this section' % name)

            if isinstance(prop, SubSection):
                if prop.repeats:
                    if isinstance(value, List):
                        for item in value:
                            self.m_add_sub_section(item)
                    else:
                        raise TypeError('Sub section %s repeats, but no list was given' % prop.name)
                else:
                    self.m_add_sub_section(item)

            else:
                setattr(self, name, value)

    def m_to_dict(self) -> Dict[str, Any]:
        """Returns the data of this section as a json serializeable dictionary. """

        def items() -> Iterable[Tuple[str, Any]]:
            yield 'm_def', self.m_def.name
            if self.m_parent_index != -1:
                yield 'm_parent_index', self.m_parent_index

            for name, sub_section in self.m_def.all_sub_sections.items():
                if name not in self.m_data:
                    continue

                if sub_section.repeats:
                    yield name, [item.m_to_dict() for item in self.m_data[name]]
                else:
                    yield name, self.m_data[name].m_to_dict()

            for name, quantity in self.m_def.all_quantities.items():
                if name in self.m_data:
                    value = getattr(self, name)
                    if hasattr(value, 'tolist'):
                        value = value.tolist()

                    # TODO
                    if isinstance(quantity.type, Section):
                        value = str(value)
                    # TODO
                    if isinstance(value, type):
                        value = str(value)
                    # TODO
                    if isinstance(value, np.dtype):
                        value = str(value)
                    # TODO
                    if isinstance(value, _Unit):
                        value = str(value)

                    yield name, value

        return {key: value for key, value in items()}

    @classmethod
    def m_from_dict(cls: Type[MSectionBound], dct: Dict[str, Any]) -> MSectionBound:
        section_def = cls.m_def

        # remove m_def and m_parent_index, they set themselves automatically
        assert section_def.name == dct.pop('m_def', None)
        dct.pop('m_parent_index', -1)

        def items():
            for name, sub_section_def in section_def.all_sub_sections.items():
                if name in dct:
                    sub_section_value = dct.pop(name)
                    if sub_section_def.repeats:
                        yield name, [
                            sub_section_def.sub_section.section_cls.m_from_dict(sub_section_dct)
                            for sub_section_dct in sub_section_value]
                    else:
                        yield name, sub_section_def.sub_section.section_cls.m_from_dict(sub_section_value)

            for key, value in dct.items():
                yield key, value

        dct = {key: value for key, value in items()}
        section_instance = cast(MSectionBound, section_def.section_cls())
        section_instance.m_update(**dct)
        return section_instance

    def m_to_json(self, **kwargs):
        """Returns the data of this section as a json string. """
        return json.dumps(self.m_to_dict(), **kwargs)

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
                    if isinstance(value, MSection):
                        yield value, attr, name, self

            elif isinstance(attr, MSection):
                yield value, value, name, self

    def __repr__(self):
        m_section_name = self.m_def.name
        name = ''
        if 'name' in self.m_data:
            name = self.m_data['name']

        return '%s:%s' % (name, m_section_name)


class MCategory(metaclass=MObjectMeta):

    m_def: 'Category' = None

    @classmethod
    def __init_cls__(cls):
        # ensure that the m_def is defined
        m_def = cls.m_def
        if m_def is None:
            m_def = Category()
            setattr(cls, 'm_def', m_def)

        # transfer name and description to m_def
        m_def.name = cls.__name__
        if cls.__doc__ is not None:
            m_def.description = inspect.cleandoc(cls.__doc__).strip()

        # add section cls' section to the module's package
        module_name = cls.__module__
        pkg = Package.from_module(module_name)
        pkg.m_add_sub_section(cls.m_def)


# M3, the definitions that are used to write definitions. These are the section definitions
# for sections Section and Quantity.They define themselves; i.e. the section definition
# for Section is the same section definition.
# Due to this circular nature (hen-egg-problem), the classes for sections Section and
# Quantity do only contain placeholder for their own section and quantity definitions.
# These placeholder are replaced, once the necessary classes are defined. This process
# is referred to as 'bootstrapping'.

_definition_change_counter = 0


class cached_property:
    """ A property that allows to cache the property value.

    The cache will be invalidated whenever a new definition is added. Once all definitions
    are loaded, the cache becomes stable and complex derived results become available
    instantaneous.
    """
    def __init__(self, f):
        self.__doc__ = getattr(f, "__doc__")
        self.f = f
        self.change = -1
        self.values: Dict[type(self), Any] = {}

    def __get__(self, obj, cls):
        if obj is None:
            return self

        global _definition_change_counter
        if self.change != _definition_change_counter:
            self.values = {}

        value = self.values.get(obj, None)
        if value is None:
            value = self.f(obj)
            self.values[obj] = value

        return value


class Definition(MSection):

    __all_definitions: Dict[Type[MSection], List[MSection]] = {}

    name: 'Quantity' = None
    description: 'Quantity' = None
    links: 'Quantity' = None
    categories: 'Quantity' = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        global _definition_change_counter
        _definition_change_counter += 1

        for cls in self.__class__.mro() + [self.__class__]:
            definitions = Definition.__all_definitions.setdefault(cls, [])
            definitions.append(self)

    @classmethod
    def all_definitions(cls: Type[MSectionBound]) -> Iterable[MSectionBound]:
        """ Returns all definitions of this definition class. """
        return cast(Iterable[MSectionBound], Definition.__all_definitions.get(cls, []))

    @cached_property
    def all_categories(self):
        """ All categories of this definition and its categories. """
        all_categories = list(self.categories)
        for category in self.categories:  # pylint: disable=not-an-iterable
            for super_category in category.all_categories:
                all_categories.append(super_category)

        return all_categories


class Property(Definition):
    pass


class Quantity(Property):
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
        if type(self.type) == np.dtype:
            if type(value) != np.ndarray:
                value = np.array(value, dtype=self.type)
            elif self.type != value.dtype:
                value = np.array(value, dtype=self.type)

        elif type(value) == np.ndarray:
            value = value.tolist()

        MSection.m_type_check(self, value)
        obj.m_data[self.__name] = value

    def __delete__(self, obj):
        if obj is None:
            # class (def) case
            raise KeyError('Cannot delete quantity definition. Only values can be deleted.')

        # object (instance) case
        del obj.m_data[self.__name]


class SubSection(Property):
    """ Allows to assign a section class as a sub-section to another section class. """

    sub_section: 'Quantity' = None
    repeats: 'Quantity' = None

    def __get__(self, obj: MSection, type=None) -> Union[MSection, 'Section']:
        if obj is None:
            # the class attribute case
            return self

        else:
            # the object attribute case
            m_data_value = obj.m_data.get(self.name, None)
            if m_data_value is None:
                if self.repeats:
                    m_data_value = []
                    obj.m_data[self.name] = m_data_value

            return m_data_value

    def __set__(self, obj: MSection, value: Union[MSection, List[MSection]]):
        raise NotImplementedError('Sub sections cannot be set directly. Use m_create.')

    def __delete__(self, obj):
        raise NotImplementedError('Sub sections cannot be deleted directly.')


class Section(Definition):
    """Used to define section that organize meta-info data into containment hierarchies.

    Section definitions determine what quantities and sub-sections can appear in a section
    instance.

    In Python terms, sections are classes. Sub-sections and quantities are attributes of
    respective instantiating objects. For each section class there is a corresponding
    :class:`Section` instance that describes this class as a section. This instance
    is referred to as 'section definition' in contrast to the Python class that we call
    'section class'.
    """

    section_cls: Type[MSection] = None
    """ The section class that corresponse to this section definition. """

    quantities: 'SubSection' = None
    sub_sections: 'SubSection' = None

    # TODO super = Quantity(type=Section, shape=['0..*']), inherit all quantity definition
    #      from the given sections, derived from Python base classes
    # TODO extends = Quantity(type=bool), denotes this section as a container for
    #      new quantities that belong to the base-class section definitions

    @cached_property
    def all_properties(self) -> Dict[str, Union['SubSection', Quantity]]:
        """ All attribute (sub section and quantity) definitions. """

        properties: Dict[str, Union[SubSection, Quantity]] = dict(**self.all_quantities)
        properties.update(**self.all_sub_sections)
        return properties

    @cached_property
    def all_quantities(self) -> Dict[str, Quantity]:
        """ All quantity definition in the given section definition. """

        return {
            quantity.name: quantity
            for quantity in self.m_data.get('quantities', [])}

    @cached_property
    def all_sub_sections(self) -> Dict[str, 'SubSection']:
        """ All sub section definitions for this section definition by name. """

        return {
            sub_section.name: sub_section
            for sub_section in self.m_data.get('sub_sections', [])}

    @cached_property
    def all_sub_sections_by_section(self) -> Dict['Section', 'SubSection']:
        """ All sub section definitions for this section definition by their section definition. """
        return {
            sub_section.sub_section: sub_section
            for sub_section in self.m_data.get('sub_sections', [])}


class Package(Definition):

    section_definitions: 'SubSection'
    category_definitions: 'SubSection'

    @staticmethod
    def from_module(module_name: str):
        module = sys.modules[module_name]

        pkg: 'Package' = getattr(module, 'm_package', None)
        if pkg is None:
            pkg = Package()
            setattr(module, 'm_package', pkg)

        pkg.name = module_name
        if pkg.description is None and module.__doc__ is not None:
            pkg.description = inspect.cleandoc(module.__doc__).strip()

        return pkg


class Category(Definition):
    """Can be used to define categories for definitions.

    Each definition, including categories themselves, can belong to a set of categories.
    Categories therefore form a hierarchy of concepts that definitions can belong to, i.e.
    they form a `is a` relationship.

    In the old meta-info this was known as `abstract types`.
    """

    @cached_property
    def definitions(self) -> Iterable[Definition]:
        """ All definitions that are directly or indirectly in this category. """
        return list([
            definition for definition in Definition.all_definitions()
            if self in definition.all_categories])


Section.m_def = Section(name='Section', _bs=True)
Section.m_def.m_def = Section.m_def
Section.m_def.section_cls = Section

Quantity.m_def = Section(name='Quantity', _bs=True)
SubSection.m_def = Section(name='SubSection', _bs=True)

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
Definition.categories = Quantity(
    type=Category.m_def, shape=['0..*'], default=[], name='categories', _bs=True,
    description='''
    The categories that this definition belongs to. See :class:`Category`.
    ''')

Section.quantities = SubSection(
    sub_section=Quantity.m_def, repeats=True,
    description='''The quantities of this section.''')

Section.sub_sections = SubSection(
    sub_section=SubSection.m_def, repeats=True,
    description='''The sub sections of this section.''')

SubSection.repeats = Quantity(
    type=bool, name='repeats', default=False, _bs=True,
    description='''Wether this sub section can appear only once or multiple times. ''')

SubSection.sub_section = Quantity(
    type=Section.m_def, name='sub_section', _bs=True, description='''
    The section definition for the sub section. Only section instances of this definition
    can be contained as sub sections.
    ''')

Quantity.m_def.section_cls = Quantity
Quantity.type = Quantity(
    type=Union[type, Enum, Section, np.dtype], name='type', _bs=True, description='''
    The type of the quantity.

    Can be one of the following:

    - none to support any value
    - a build-in primitive Python type, e.g. ``int``, ``str``
    - an instance of :class:`Enum`, e.g. ``Enum(['one', 'two', 'three'])
    - a instance of Section, i.e. a section definition. This will define a reference
    - a custom meta-info DataType
    - a numpy dtype,

    If set to a dtype, this quantity will use a numpy array to store values. It will use
    the given dtype. If not set, this quantity will use (nested) Python lists to store values.
    If values are set to the property, they will be converted to the respective
    representation.

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

Package.m_def = Section(name='Package', _bs=True)

Category.m_def = Section(name='Category', _bs=True)

Package.section_definitions = SubSection(
    sub_section=Section.m_def, name='section_definitions', repeats=True,
    description=''' The sections defined in this package. ''')

Package.category_definitions = SubSection(
    sub_section=Category.m_def, name='category_definitions', repeats=True,
    description=''' The categories defined in this package. ''')

is_bootstrapping = False

Package.__init_cls__()
Category.__init_cls__()
Section.__init_cls__()
SubSection.__init_cls__()
Quantity.__init_cls__()

units = UnitRegistry()
""" The default pint unit registry that should be used to give units to quantity definitions. """
