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

from typing import Type, TypeVar, Union, Tuple, Iterable, List, Any, Dict, Set, \
    Callable as TypingCallable, cast
from collections.abc import Iterable as IterableABC
import sys
import inspect
import re
import json
import itertools

import numpy as np
import pint
import pint.unit
import pint.quantity
import aniso8601
from datetime import datetime
import pytz
import docstring_parser
import flask_restplus.inputs


m_package: 'Package' = None

is_bootstrapping = True
MSectionBound = TypeVar('MSectionBound', bound='MSection')
T = TypeVar('T')


# Metainfo errors

class MetainfoError(Exception):
    ''' Metainfo related errors. '''
    pass


class DeriveError(MetainfoError):
    ''' An error occurred while computing a derived value. '''
    pass


class MetainfoReferenceError(MetainfoError):
    ''' An error indicating that a reference could not be resolved. '''
    pass


# Metainfo quantity data types

class MEnum():
    '''Allows to define str types with values limited to a pre-set list of possible values.'''
    def __init__(self, *args, **kwargs):
        # Supports one big list in place of args
        if len(args) == 1 and isinstance(args[0], list):
            args = args[0]

        # If non-named arguments are given, the default is to have them placed
        # into a dictionary with their string value as both the enum name and
        # the value.
        for arg in args:
            if arg in kwargs:
                raise ValueError("Duplicate value '{}' provided for enum".format(arg))
            kwargs[arg] = arg

        self._values = set(kwargs.values())  # For allowing constant time member check
        self._map = kwargs

    def __getattr__(self, attr):
        return self._map[attr]


class MProxy():
    ''' A placeholder object that acts as reference to a value that is not yet resolved.

    Attributes:
        url: The reference represented as an URL string.
    '''

    def __init__(self, url: str):
        self.url = url


class DataType:
    '''
    Allows to define custom data types that can be used in the meta-info.

    The metainfo supports the most types out of the box. These includes the python build-in
    primitive types (int, bool, str, float, ...), references to sections, and enums.
    However, in some occasions you need to add custom data types.

    This base class lets you customize various aspects of value treatment. This includes
    type checks and various value transformations. This allows to store values in the
    section differently from how the usermight set/get them, and it allows to have non
    serializeable values that are transformed on de-/serialization.
    '''
    def set_normalize(self, section: 'MSection', quantity_def: 'Quantity', value: Any) -> Any:
        ''' Transforms the given value before it is set and checks its type. '''
        return value

    def get_normalize(self, section: 'MSection', quantity_def: 'Quantity', value: Any) -> Any:
        ''' Transforms the given value when it is get. '''
        return value

    def serialize(self, section: 'MSection', quantity_def: 'Quantity', value: Any) -> Any:
        ''' Transforms the given value when making the section serializeable. '''
        return value

    def deserialize(self, section: 'MSection', quantity_def: 'Quantity', value: Any) -> Any:
        ''' Transforms the given value from its serializeable form. '''
        return value


range_re = re.compile(r'(\d)\.\.(\d|\*)')


class _Dimension(DataType):

    def set_normalize(self, section, quantity_def: 'Quantity', value):
        if isinstance(value, int):
            return value

        if isinstance(value, str):
            if value.isidentifier():
                return value
            if re.match(range_re, value):
                return value

        if isinstance(value, Section):
            return value

        if isinstance(value, type) and hasattr(value, 'm_def'):
            return value

        raise TypeError('%s is not a valid dimension' % str(value))

    @staticmethod
    def check_dimension(section, dimension, length):
        if isinstance(dimension, int):
            return dimension == length
        if isinstance(dimension, str):
            if dimension.isidentifier():
                return dimension == getattr(section, dimension)

            m = re.match(range_re, dimension)
            start = int(m.group(1))
            end = -1 if m.group(2) == '*' else int(m.group(2))
            return start <= length and (end == -1 or length <= end)


class _Unit(DataType):
    def set_normalize(self, section, quantity_def: 'Quantity', value):
        if isinstance(value, str):
            value = units.parse_units(value)

        elif not isinstance(value, pint.unit._Unit):
            raise TypeError('Units must be given as str or pint Unit instances.')

        return value

    def serialize(self, section, quantity_def: 'Quantity', value):
        return value.__str__()

    def deserialize(self, section, quantity_def: 'Quantity', value):
        return units.parse_units(value)


units = pint.UnitRegistry()
''' The default pint unit registry that should be used to give units to quantity definitions. '''


class _Callable(DataType):
    def serialize(self, section, quantity_def: 'Quantity', value):
        raise MetainfoError('Callables cannot be serialized')

    def deserialize(self, section, quantity_def: 'Quantity', value):
        raise MetainfoError('Callables cannot be serialized')


class _QuantityType(DataType):
    ''' Data type for defining the type of a metainfo quantity.

    A metainfo quantity type can be one of

    - python build-in primitives: int, float, bool, str
    - numpy dtypes, e.g. f, int32
    - a section definition to define references
    - an MEnum instance to use it's values as possible str values
    - a custom datatype, i.e. instance of :class:`DataType`
    - Any
    '''

    def set_normalize(self, section, quantity_def, value):
        if value in [str, int, float, bool]:
            return value

        if isinstance(value, MEnum):
            for enum_value in value._values:
                if not isinstance(enum_value, str):
                    raise TypeError('MEnum value %s is not a string.' % enum_value)
            return value

        if type(value) == np.dtype:
            return value

        if isinstance(value, Section):
            return value

        if isinstance(value, DataType):
            return value

        if value == Any:
            return value

        if isinstance(value, type):
            section = getattr(value, 'm_def', None)
            if section is not None:
                return Reference(section)

        raise MetainfoError(
            'Type %s of %s is not a valid metainfo quantity type' %
            (value, quantity_def))

    def serialize(self, section, quantity_def, value):
        if value is str or value is int or value is float or value is bool:
            return dict(type_kind='python', type_data=value.__name__)

        if isinstance(value, MEnum):
            return dict(type_kind='Enum', type_data=list(value))

        if type(value) == np.dtype:
            return dict(type_kind='numpy', type_data=str(value))

        if isinstance(value, Reference):
            return dict(type_kind='reference', type_data=value.target_section_def.m_path())

        if isinstance(value, DataType):
            module = value.__class__.__module__
            if module is None or module == str.__class__.__module__:
                type_data = value.__class__.__name__
            else:
                type_data = '%s.%s' % (module, value.__class__.__name__)

            return dict(type_kind='custom', type_data=type_data)

        if value == Any:
            return dict(type_kind='Any')

        raise MetainfoError(
            'Type %s of %s is not a valid metainfo quantity type' %
            (value, quantity_def))


class Reference(DataType):
    ''' Datatype used for reference quantities. '''

    def __init__(self, section_def: 'Section'):
        if not isinstance(section_def, Section):
            raise MetainfoError('%s is not a section definition.' % section_def)
        self.target_section_def = section_def

    def set_normalize(self, section: 'MSection', quantity_def: 'Quantity', value: Any) -> Any:
        if self.target_section_def.m_follows(Definition.m_def):
            # special case used in metainfo definitions, where we reference metainfo definitions
            # using their Python class. E.g. referencing a section definition using its
            # class instead of the object: Run vs. Run.m_def
            if isinstance(value, type):
                definition = getattr(value, 'm_def', None)
                if definition is not None and definition.m_follows(self.target_section_def):
                    return definition

        if isinstance(value, MProxy):
            return value

        if not isinstance(value, MSection):
            raise TypeError(
                'The value %s is not a section and can not be used as a reference.' % value)

        if not value.m_follows(self.target_section_def):
            raise TypeError(
                '%s is not a %s and therefore an invalid value of %s.' %
                (value, self.target_section_def, quantity_def))

        return value

    def resolve(self, section: 'MSection', quantity_def: 'Quantity', value: Any) -> 'MSection':
        return section.m_resolve(value.url)

    def get_normalize(self, section: 'MSection', quantity_def: 'Quantity', value: Any) -> Any:
        if isinstance(value, MProxy):
            resolved: 'MSection' = self.resolve(section, quantity_def, value)
            if resolved is None:
                raise ReferenceError('Could not resolve %s from %s.' % (value, section))

            return resolved

        return value

    def serialize(self, section: 'MSection', quantity_def: 'Quantity', value: Any) -> Any:
        return value.m_path()

    def deserialize(self, section: 'MSection', quantity_def: 'Quantity', value: Any) -> Any:
        return MProxy(value)


class _Datetime(DataType):

    def _parse(self, datetime_str: str) -> datetime:
        try:
            return aniso8601.parse_datetime(datetime_str)
        except ValueError:
            pass

        try:
            return aniso8601.parse_date(datetime_str)
        except ValueError:
            pass

        try:
            # TODO necessary?
            return flask_restplus.inputs.datetime_from_rfc822(datetime_str)
        except ValueError:
            pass

        raise TypeError('Invalid date literal "{0}"'.format(datetime_str))

    def _convert(self, value):
        if value is None:
            return None

        if isinstance(value, str):
            value = self._parse(value)

        elif isinstance(value, (int, float)):
            value = datetime.fromtimestamp(value)

        if not isinstance(value, datetime):
            raise TypeError('%s is not a datetime.' % value)

        return value

    def set_normalize(self, section: 'MSection', quantity_def: 'Quantity', value: Any) -> Any:
        return self._convert(value)

    def serialize(self, section: 'MSection', quantity_def: 'Quantity', value: Any) -> Any:
        if value is None:
            return None

        value.replace(tzinfo=pytz.utc)
        return value.isoformat()

    def deserialize(self, section: 'MSection', quantity_def: 'Quantity', value: Any) -> Any:
        return self._convert(value)


Dimension = _Dimension()
Unit = _Unit()
QuantityType = _QuantityType()
Callable = _Callable()
Datetime = _Datetime()


# Metainfo data storage and reflection interface

class MObjectMeta(type):

    def __new__(self, cls_name, bases, dct):
        do_init = dct.get('do_init', None)
        if do_init is not None:
            del(dct['do_init'])
        else:
            do_init = True

        cls = super().__new__(self, cls_name, bases, dct)

        init = getattr(cls, '__init_cls__')
        if init is not None and do_init and not is_bootstrapping:
            init()
        return cls


SectionDef = Union[str, 'Section', 'SubSection', Type[MSectionBound]]
''' Type for section definition references.

This can either be :

- the name of the section
- the section definition itself
- the definition of a sub section
- or the section definition Python class
'''


class MResource():
    '''Represents a collection of related metainfo data, i.e. a set of :class:`MSection` instances.

    MResource allows to keep related objects together and resolve sections of certain
    section definitions.
    '''
    def __init__(self):
        self.__data: Dict['Section', List['MSection']] = dict()
        self.contents: List['MSection'] = []

    def create(self, section_cls: Type[MSectionBound], *args, **kwargs) -> MSectionBound:
        ''' Create an instance of the given section class and adds it to this resource. '''
        result = section_cls(*args, **kwargs)
        self.add(result)
        return cast(MSectionBound, result)

    def add(self, section):
        section.m_resource = self
        self.__data.setdefault(section.m_def, []).append(section)
        if section.m_parent is None:
            self.contents.append(section)

    def remove(self, section):
        assert section.m_resource == self, 'Can only remove section from the resource that contains it.'
        section.m_resource = None
        self.__data.get(section.m_def).remove(section)
        if section.m_parent is not None:
            self.contents.remove(section)

    def all(self, section_cls: Type[MSectionBound]) -> List[MSectionBound]:
        ''' Returns all instances of the given section class in this resource. '''
        return cast(List[MSectionBound], self.__data.get(section_cls.m_def, []))

    def unload(self):
        ''' Breaks all references among the contain metainfo sections to allow GC. '''
        for collections in self.__data.values():
            for section in collections:
                section.m_parent = None
            collections.clear()

        # TODO break actual references via quantities


class MSection(metaclass=MObjectMeta):  # TODO find a way to make this a subclass of collections.abs.Mapping
    '''Base class for all section instances on all meta-info levels.

    All `section instances` indirectly instantiate the :class:`MSection` and therefore all
    members of :class:`MSection` are available on all `section instances`. :class:`MSection`
    provides many special attributes and functions (they all start with ``m_``) that allow
    to reflect on a `section's definition` and allow to manipulate the `section instance`
    without a priori knowledge of the `section defintion`.

    It also carries all the data for each section. All sub-classes only define specific
    sections in terms of possible sub-sections and quantities. The data is managed here.

    Attributes:
        m_def: The `section definition` that this `section instance` follows as a
            :class:`Section` object.

        m_parent:
            If this section is a sub-section, this references the parent section instance.

        m_parent_sub_section:
            If this section is a sub-section, this is the :class:`SubSection` that defines
            this relationship.

        m_parent_index:
            For repeatable sections, parent keep a list of sub-sections. This is the index
            of this section in the respective parent sub-section list.

        m_resource: The :class:`MResource` that contains and manages this section.

    '''

    m_def: 'Section' = None

    def __init__(
            self, m_def: 'Section' = None, m_resource: MResource = None, **kwargs):

        self.m_def: 'Section' = m_def
        self.m_parent: 'MSection' = None
        self.m_parent_sub_section: 'SubSection' = None
        self.m_parent_index = -1
        self.m_resource = m_resource

        # get missing m_def from class
        cls = self.__class__
        if self.m_def is None:
            self.m_def = cls.m_def

        # check m_def
        if cls.m_def is not None:
            if self.m_def != cls.m_def:
                MetainfoError('Section class and section definition must match.')

            if self.m_def.extends_base_section:
                MetainfoError('Section extends another section and cannot be instantiated.')

        else:
            if not is_bootstrapping:
                MetainfoError('Section has not m_def.')

        # get annotations from kwargs
        self.m_annotations: Dict[str, Any] = {}
        other_kwargs = {}
        for key, value in kwargs.items():
            if key.startswith('a_'):
                self.m_annotations[key[2:]] = value
            else:
                other_kwargs[key] = value

        # get additional annotations from the section definition
        if not is_bootstrapping:
            for section_annotation in self.m_def.m_x(SectionAnnotation, as_list=True):
                for name, annotation in section_annotation.new(self).items():
                    self.m_annotations[name] = annotation

        # add annotation attributes for names annotations
        for annotation_name, annotation in self.m_annotations.items():
            setattr(self, 'a_%s' % annotation_name, annotation)

        # set remaining kwargs
        if is_bootstrapping:
            self.__dict__.update(**other_kwargs)  # type: ignore
        else:
            self.m_update(**other_kwargs)

    @classmethod
    def __init_cls__(cls):
        # ensure that the m_def is defined
        m_def = cls.m_def
        if m_def is None:
            m_def = Section()
            setattr(cls, 'm_def', m_def)

        # transfer name m_def
        m_def.name = cls.__name__
        m_def.section_cls = cls

        # add base sections
        if m_def.extends_base_section:
            base_sections_count = len(cls.__bases__)
            if base_sections_count == 0:
                raise MetainfoError(
                    'Section %s extend the base section, but has no base section.' % m_def)

            if base_sections_count > 1:
                raise MetainfoError(
                    'Section %s extend the base section, but has more than one base section' % m_def)

            base_section_cls = cls.__bases__[0]
            base_section = getattr(base_section_cls, 'm_def', None)
            if base_section is None:
                raise MetainfoError(
                    'The base section of %s is not a section class.' % m_def)

            for name, attr in cls.__dict__.items():
                if isinstance(attr, Property):
                    setattr(base_section_cls, name, attr)

            section_to_add_properties_to = base_section
            base_sections = [base_section]
        else:
            base_sections: List[Section] = []
            for base_cls in cls.__bases__:
                if base_cls != MSection:
                    base_section = getattr(base_cls, 'm_def')
                    if base_section is None:
                        raise TypeError(
                            'Section defining classes must have MSection or a decendant as '
                            'base classes.')

                    base_sections.append(base_section)

            section_to_add_properties_to = m_def

        m_def.m_set(Section.base_sections, base_sections)

        # transfer names, descriptions, constraints, event_handlers
        constraints: Set[str] = set()
        event_handlers: Set[Callable] = set(m_def.event_handlers)
        for name, attr in cls.__dict__.items():
            # transfer names and descriptions for properties, init properties
            if isinstance(attr, Property):
                attr.name = name
                if attr.description is not None:
                    attr.description = inspect.cleandoc(attr.description).strip()
                    attr.__doc__ = attr.description

                if isinstance(attr, Quantity):
                    section_to_add_properties_to.m_add_sub_section(Section.quantities, attr)
                elif isinstance(attr, SubSection):
                    section_to_add_properties_to.m_add_sub_section(Section.sub_sections, attr)
                else:
                    raise NotImplementedError('Unknown property kind.')

            if inspect.isfunction(attr):
                method_name = attr.__name__

                # transfer constraints
                if method_name.startswith('c_'):
                    constraint = method_name[2:]
                    constraints.add(constraint)

                # register event_handlers from event_handler methods
                if method_name.startswith('on_set') or method_name.startswith('on_add_sub_section'):
                    if attr not in event_handlers:
                        event_handlers.add(attr)

        # add handler and constraints from base sections
        for base_section in m_def.all_base_sections:
            for constraint in base_section.constraints:
                constraints.add(constraint)
            for event_handler in base_section.event_handlers:
                event_handlers.add(event_handler)

        m_def.constraints = list(constraints)
        m_def.event_handlers = list(event_handlers)

        # add section cls' section to the module's package
        module_name = cls.__module__
        pkg = Package.from_module(module_name)
        pkg.m_add_sub_section(Package.section_definitions, cls.m_def)

        # apply_google_docstrings
        # Parses the google doc string of the given class and properly updates the
        # definition descriptions.

        # This allows to document quantities and sub-sections with 'Args:' in the section
        # class. It will remove the 'Args' section from the section definition and will
        # set the respective pieces to the quantity and sub-section descriptions.
        docstring = cls.__doc__
        if docstring is not None:
            parsed_docstring = docstring_parser.parse(docstring)
            short = parsed_docstring.short_description
            dsc = parsed_docstring.long_description

            if short and dsc:
                description = '%s %s' % (short.strip(), dsc.strip())
            elif short:
                description = short.strip()
            elif dsc:
                description = dsc.strip()
            else:
                description = None

            if m_def.description is None:
                m_def.description = description

            for param in parsed_docstring.params:
                prop = m_def.all_properties.get(param.arg_name)
                if prop is not None:
                    if prop.description is None:
                        prop.description = param.description

        m_def.__init_metainfo__()

    def __getattr__(self, name):
        # This will make mypy and pylint ignore 'missing' dynamic attributes and functions
        # and wrong types of those.
        # Ideally we have a plugin for both that add the corrent type info
        raise AttributeError(name)

    def __check_np(self, quantity_def: 'Quantity', value: np.ndarray) -> np.ndarray:
        # TODO this feels expensive, first check, then possible convert very often?
        # if quantity_ref.type != value.dtype:
        #     raise MetainfoError(
        #         'Quantity dtype %s and value dtype %s do not match.' %
        #         (quantity_ref.type, value.dtype))

        return value

    def __set_normalize(self, quantity_def: 'Quantity', value: Any) -> Any:

        if isinstance(quantity_def.type, DataType):
            return quantity_def.type.set_normalize(self, quantity_def, value)

        elif isinstance(quantity_def.type, Section):
            if isinstance(value, MProxy):
                return value

            if not isinstance(value, MSection):
                raise TypeError(
                    'The value %s for reference quantity %s is not a section instance.' %
                    (value, quantity_def))

            if not value.m_follows(quantity_def.type):
                raise TypeError(
                    'The value %s for quantity %s does not follow %s' %
                    (value, quantity_def, quantity_def.type))

        elif isinstance(quantity_def.type, MEnum):
            if value not in quantity_def.type._values:
                raise TypeError(
                    'The value %s is not an enum value for quantity %s.' %
                    (value, quantity_def))

        elif quantity_def.type == Any:
            pass

        elif quantity_def.type == str and type(value) == np.str_:
            return str(value)

        elif quantity_def.type == bool and type(value) == np.bool_:
            return bool(value)

        else:
            if value is not None and type(value) != quantity_def.type:
                raise TypeError(
                    'The value %s with type %s for quantity %s is not of type %s' %
                    (value, type(value), quantity_def, quantity_def.type))

        return value

    def __resolve_synonym(self, quantity_def: 'Quantity') -> 'Quantity':
        if quantity_def.synonym_for is not None:
            return self.m_def.all_quantities[quantity_def.synonym_for]
        return quantity_def

    def __to_np(self, quantity_def: 'Quantity', value):
        if isinstance(value, pint.quantity._Quantity):
            if quantity_def.unit is None:
                raise MetainfoError(
                    'The quantity %s has not a unit, but value %s has.' %
                    (quantity_def, value))
            value = value.to(quantity_def.unit).magnitude

        if len(quantity_def.shape) > 0 and type(value) != np.ndarray:
            try:
                value = np.asarray(value)
            except TypeError:
                raise TypeError(
                    'Could not convert value %s of %s to a numpy array' %
                    (value, quantity_def))
        elif type(value) != quantity_def.type.type:
            try:
                value = quantity_def.type.type(value)
            except TypeError:
                raise TypeError(
                    'Could not convert value %s of %s to a numpy scalar' %
                    (value, quantity_def))

        return self.__check_np(quantity_def, value)

    def m_set(self, quantity_def: 'Quantity', value: Any) -> None:
        ''' Set the given value for the given quantity. '''
        quantity_def = self.__resolve_synonym(quantity_def)

        if quantity_def.derived is not None:
            raise MetainfoError('The quantity %s is derived and cannot be set.' % quantity_def)

        if type(quantity_def.type) == np.dtype:
            value = self.__to_np(quantity_def, value)

        else:
            dimensions = len(quantity_def.shape)
            if dimensions == 0:
                value = self.__set_normalize(quantity_def, value)

            elif dimensions == 1:
                if type(value) == str or not isinstance(value, IterableABC):
                    raise TypeError(
                        'The shape of %s requires an iterable value, but %s is not iterable.' %
                        (quantity_def, value))

                value = list(self.__set_normalize(quantity_def, item) for item in value)

            else:
                raise MetainfoError(
                    'Only numpy arrays and dtypes can be used for higher dimensional '
                    'quantities.')

        self.__dict__[quantity_def.name] = value

        for handler in self.m_def.event_handlers:
            if handler.__name__.startswith('on_set'):
                handler(self, quantity_def, value)

    def m_get(self, quantity_def: 'Quantity') -> Any:
        ''' Retrieve the given value for the given quantity. '''
        return quantity_def.__get__(self, Quantity)

    def m_is_set(self, quantity_def: 'Quantity') -> bool:
        ''' True if the given quantity is set. '''
        quantity_def = self.__resolve_synonym(quantity_def)
        if quantity_def.derived is not None:
            return True

        return quantity_def.name in self.__dict__

    def m_add_values(self, quantity_def: 'Quantity', values: Any, offset: int) -> None:
        ''' Add (partial) values for the given quantity of higher dimensionality. '''
        # TODO
        raise NotImplementedError()

    def m_add_sub_section(self, sub_section_def: 'SubSection', sub_section: 'MSection') -> None:
        ''' Adds the given section instance as a sub section of the given sub section definition. '''

        parent_index = -1
        if sub_section_def.repeats:
            parent_index = self.m_sub_section_count(sub_section_def)

        else:
            old_sub_section = self.__dict__.get(sub_section_def.name)
            if old_sub_section is not None:
                old_sub_section.m_parent = None
                old_sub_section.m_parent_sub_section = None
                old_sub_section.m_parent_index = -1
                if self.m_resource is not None:
                    self.m_resource.remove(sub_section)

        sub_section.m_parent = self
        sub_section.m_parent_sub_section = sub_section_def
        sub_section.m_parent_index = parent_index
        if sub_section.m_resource is not None:
            sub_section.m_resource.remove(sub_section)
        if self.m_resource is not None:
            self.m_resource.add(sub_section)

        sub_section_name = sub_section_def.name
        if sub_section_def.repeats:
            sub_section_lst = self.__dict__.get(sub_section_name)
            if sub_section_lst is None:
                sub_section_lst = self.__dict__.setdefault(sub_section_name, [])

            sub_section_lst.append(sub_section)

        else:
            self.__dict__[sub_section_name] = sub_section

        for handler in self.m_def.event_handlers:
            if handler.__name__.startswith('on_add_sub_section'):
                handler(self, sub_section_def, sub_section)

    def m_remove_sub_section(self, sub_section_def: 'SubSection', index: int) -> None:
        ''' Removes the exiting section for a non repeatable sub section '''
        if sub_section_def.repeats:
            del(self.__dict__[sub_section_def.name][index])

        elif sub_section_def.name in self.__dict__:
            del(self.__dict__[sub_section_def.name])

    def m_get_sub_section(self, sub_section_def: 'SubSection', index: int) -> 'MSection':
        ''' Retrieves a single sub section of the given sub section definition. '''
        if sub_section_def.repeats:
            return self.__dict__[sub_section_def.name][index]

        else:
            return self.__dict__.get(sub_section_def.name, None)

    def m_get_sub_sections(self, sub_section_def: 'SubSection') -> Iterable['MSection']:
        ''' Retrieves  all sub sections of the given sub section definition. '''
        try:
            return self.__dict__[sub_section_def.name]
        except KeyError:
            return []

    def m_sub_section_count(self, sub_section_def: 'SubSection') -> int:
        ''' Returns the number of sub sections for the given sub section definition. '''
        try:
            value = self.__dict__[sub_section_def.name]
            if sub_section_def.repeats:
                return len(value)
            else:
                return 1
        except KeyError:
            return 0

    def m_create(
            self, section_cls: Type[MSectionBound], sub_section_def: 'SubSection' = None,
            **kwargs) -> MSectionBound:
        ''' Creates a section instance and adds it to this section provided there is a
        corresponding sub section.

        Args:
            section_cls: The section class for the sub-secton to create
            sub_section_def: If there are multiple sub-sections for the given class,
                this must be used to explicitely state the sub-section definition.
        '''

        section_def = section_cls.m_def
        sub_section_defs = self.m_def.all_sub_sections_by_section.get(section_def, [])
        n_sub_section_defs = len(sub_section_defs)
        if n_sub_section_defs == 0:
            raise TypeError('There is no sub section to hold a %s in %s.' % (section_def, self.m_def))

        if n_sub_section_defs > 1 and sub_section_def is None:
            raise MetainfoError(
                'There are multiple sub section to hold a %s in %s, '
                'but no sub-section was explicitely given.' % (section_def, self.m_def))

        if sub_section_def is not None and sub_section_def not in sub_section_defs:
            raise MetainfoError(
                'The given sub-section class %s does not match the given sub-section '
                'definition %s.' % (section_cls, sub_section_def))

        if sub_section_def is None:
            sub_section_def = sub_section_defs[0]

        sub_section = section_cls(**kwargs)
        self.m_add_sub_section(sub_section_def, sub_section)

        return cast(MSectionBound, sub_section)

    def m_update(self, safe: bool = True, **kwargs):
        ''' Updates all quantities and sub-sections with the given arguments. '''
        if safe:
            for name, value in kwargs.items():
                prop = self.m_def.all_properties.get(name, None)
                if prop is None:
                    raise KeyError('%s is not an attribute of this section %s' % (name, self))

                if isinstance(prop, SubSection):
                    if prop.repeats:
                        if isinstance(value, List):
                            for item in value:
                                self.m_add_sub_section(prop, item)
                        else:
                            raise TypeError('Sub section %s repeats, but no list was given' % prop.name)
                    else:
                        self.m_add_sub_section(prop, item)

                else:
                    self.m_set(prop, value)

        else:
            self.__dict__.update(**kwargs)

    def m_as(self, section_cls: Type[MSectionBound]) -> MSectionBound:
        ''' 'Casts' this section to the given extending sections. '''
        return cast(MSectionBound, self)

    def m_follows(self, definition: 'Section') -> bool:
        ''' Determines if this section's definition is or is derived from the given definition. '''
        return self.m_def == definition or definition in self.m_def.all_base_sections

    def m_to_dict(self, with_meta: bool = False, include_defaults: bool = False) -> Dict[str, Any]:
        '''
        Returns the data of this section as a json serializeable dictionary.

        Arguments:
            with_meta: Include information about the section definition and the sections
                position in its parent.
            include_defaults: Include default values of unset quantities.
        '''

        def items() -> Iterable[Tuple[str, Any]]:
            # metadata
            if with_meta:
                yield 'm_def', self.m_def.name
                if self.m_parent_index != -1:
                    yield 'm_parent_index', self.m_parent_index
                if self.m_parent_sub_section is not None:
                    yield 'm_parent_sub_section', self.m_parent_sub_section.name

            # quantities
            for name, quantity in self.m_def.all_quantities.items():
                if quantity.virtual:
                    continue

                is_set = self.m_is_set(quantity)
                if not is_set:
                    if not include_defaults or not quantity.m_is_set(Quantity.default):
                        continue

                quantity_type = quantity.type

                serialize: TypingCallable[[Any], Any] = str
                if isinstance(quantity_type, Reference):

                    def reference_serialize(value):
                        if isinstance(value, MProxy):
                            return value.url
                        else:
                            return quantity_type.serialize(self, quantity, value)

                    serialize = reference_serialize

                elif isinstance(quantity_type, DataType):

                    def data_type_serialize(value):
                        return quantity_type.serialize(self, quantity, value)

                    serialize = data_type_serialize

                elif quantity_type in [str, int, float, bool]:
                    serialize = quantity_type

                elif type(quantity_type) == np.dtype:
                    pass

                elif isinstance(quantity_type, MEnum):
                    pass

                elif quantity_type == Any:
                    def _serialize(value: Any):
                        if type(value) not in [str, int, float, bool, list, type(None)]:
                            raise MetainfoError(
                                'Only python primitives are allowed for Any typed non '
                                'virtual quantities: %s of quantity %s in section %s' %
                                (value, quantity, self))

                        return value

                    serialize = _serialize

                else:
                    raise MetainfoError(
                        'Do not know how to serialize data with type %s for quantity %s' %
                        (quantity_type, quantity))

                if is_set:
                    value = self.__dict__[name]
                else:
                    value = quantity.default

                if type(quantity_type) == np.dtype:
                    serializable_value = value.tolist()

                else:
                    if len(quantity.shape) == 0:
                        serializable_value = serialize(value)
                    elif len(quantity.shape) == 1:
                        serializable_value = [serialize(i) for i in value]
                    else:
                        raise NotImplementedError('Higher shapes (%s) not supported: %s' % (quantity.shape, quantity))

                yield name, serializable_value

            # sub sections
            for name, sub_section_def in self.m_def.all_sub_sections.items():
                if sub_section_def.repeats:
                    if self.m_sub_section_count(sub_section_def) > 0:
                        yield name, [
                            item.m_to_dict(with_meta=with_meta, include_defaults=include_defaults)
                            for item in self.m_get_sub_sections(sub_section_def)]
                else:
                    sub_section = self.m_get_sub_section(sub_section_def, -1)
                    if sub_section is not None:
                        yield name, sub_section.m_to_dict(with_meta=with_meta, include_defaults=include_defaults)

        return {key: value for key, value in items()}

    @classmethod
    def m_from_dict(cls: Type[MSectionBound], dct: Dict[str, Any]) -> MSectionBound:
        ''' Creates a section from the given serializable data dictionary.

        This is the 'opposite' of :func:`m_to_dict`. It takes a deserialised dict, e.g
        loaded from JSON, and turns it into a proper section, i.e. instance of the given
        section class.
        '''

        section_def = cls.m_def

        # remove m_def, m_parent_index, m_parent_sub_section metadata,
        # they set themselves automatically
        dct.pop('m_def', None)
        dct.pop('m_parent_index', None)
        dct.pop('m_parent_sub_section', None)

        section = cls()

        for name, sub_section_def in section_def.all_sub_sections.items():
            if name in dct:
                sub_section_value = dct.pop(name)
                if sub_section_def.repeats:
                    for sub_section_dct in sub_section_value:
                        sub_section = sub_section_def.sub_section.section_cls.m_from_dict(sub_section_dct)
                        section.m_add_sub_section(sub_section_def, sub_section)

                else:
                    sub_section = sub_section_def.sub_section.section_cls.m_from_dict(sub_section_value)
                    section.m_add_sub_section(sub_section_def, sub_section)

        for name, quantity_def in section_def.all_quantities.items():
            if name in dct:
                quantity_value = dct[name]

                if type(quantity_def.type) == np.dtype:
                    quantity_value = np.asarray(quantity_value)

                if isinstance(quantity_def.type, DataType):
                    dimensions = len(quantity_def.shape)
                    if dimensions == 0:
                        quantity_value = quantity_def.type.deserialize(
                            section, quantity_def, quantity_value)
                    elif dimensions == 1:
                        quantity_value = list(
                            quantity_def.type.deserialize(section, quantity_def, item)
                            for item in quantity_value)
                    else:
                        raise MetainfoError(
                            'Only numpy quantities can have more than 1 dimension.')

                section.__dict__[name] = quantity_value  # type: ignore

        return section

    def m_to_json(self, **kwargs):
        ''' Returns the data of this section as a json string. '''
        return json.dumps(self.m_to_dict(), **kwargs)

    def m_all_contents(self) -> Iterable['MSection']:
        ''' Returns an iterable over all sub and sub subs sections. '''
        for content in self.m_contents():
            for sub_content in content.m_all_contents():
                yield sub_content

            yield content

    def m_contents(self) -> Iterable['MSection']:
        ''' Returns an iterable over all direct subs sections. '''
        for sub_section_def in self.m_def.all_sub_sections.values():
            if sub_section_def.repeats:
                index = 0
                for sub_section in self.m_get_sub_sections(sub_section_def):
                    yield sub_section
                    index += 1

            else:
                sub_section = self.m_get_sub_section(sub_section_def, -1)
                yield sub_section

    def m_path(self, quantity_def: 'Quantity' = None) -> str:
        ''' Returns the path of this section or the given quantity within the section hierarchy. '''
        if self.m_parent is None:
            return '/'

        if self.m_parent_index == -1:
            segment = self.m_parent_sub_section.name
        else:
            segment = '%s/%d' % (self.m_parent_sub_section.name, self.m_parent_index)

        if quantity_def is not None:
            segment = '%s/%s' % (segment, quantity_def.name)

        return '%s/%s' % (self.m_parent.m_path().rstrip('/'), segment)

    def m_root(self, cls: Type[MSectionBound] = None) -> MSectionBound:
        ''' Returns the first parent of the parent section that has no parent; the root. '''
        if self.m_parent is None:
            return cast(MSectionBound, self)
        else:
            return self.m_parent.m_root(cls)

    def m_parent_as(self, cls: Type[MSectionBound] = None) -> MSectionBound:
        ''' Returns the parent section with the given section class type. '''
        return cast(MSectionBound, self.m_parent)

    def m_resolve(self, path: str, cls: Type[MSectionBound] = None) -> MSectionBound:
        '''
        Resolves the given path or dotted quantity name using this section as context and
        returns the sub_section or value.
        '''
        if path.startswith('/'):
            context: 'MSection' = self.m_root()
        else:
            context = self

        path_stack = path.strip('/').split('/')
        path_stack.reverse()
        while len(path_stack) > 0:
            prop_name = path_stack.pop()
            prop_def = context.m_def.all_properties.get(prop_name, None)

            if prop_def is None:
                raise ReferenceError(
                    'Could not resolve %s, property %s does not exist in %s' %
                    (path, prop_name, context.m_def))

            if isinstance(prop_def, SubSection):
                if prop_def.repeats:
                    try:
                        index = int(path_stack.pop())
                    except ValueError:
                        raise ReferenceError(
                            'Could not resolve %s, %s repeats but there is no index in the path' %
                            (path, prop_name))

                    try:
                        context = context.m_get_sub_section(prop_def, index)
                    except Exception:
                        raise ReferenceError(
                            'Could not resolve %s, there is no sub section for %s at %d' %
                            (path, prop_name, index))

                else:
                    context = context.m_get_sub_section(prop_def, -1)
                    if context is None:
                        raise ReferenceError(
                            'Could not resolve %s, there is no sub section for %s' %
                            (path, prop_name))

            elif isinstance(prop_def, Quantity):
                if len(path_stack) > 0:
                    raise ReferenceError(
                        'Could not resolve %s, %s is not a sub section' % (path, prop_name))

                return context.m_get(prop_def)

        return cast(MSectionBound, context)

    def m_x(self, *args, **kwargs):
        # TODO remove
        return self.m_get_annotations(*args, **kwargs)

    def m_get_annotations(self, key: Union[str, type], default=None, as_list: bool = False):
        '''
        Convinience method to get annotations

        Arguments:
            key: Either the optional annoation name or an annotation class. In the first
                case the annotation is returned, regardless of its type. In the second
                case, all names and list for names are iterated and all annotations of the
                given class are returned.
            default: The default, if no annotation is found. None is  the default default.
            as_list: Returns a list, no matter how many annoations have been found.
        '''
        if isinstance(key, str):
            value = self.m_annotations.get(key, default)
            if as_list and not isinstance(value, (list, tuple)):
                return [value]
            else:
                return value

        elif isinstance(key, type):
            result_list = []
            for values in self.m_annotations.values():
                if isinstance(values, (tuple, list)):
                    for value in values:
                        if isinstance(value, key):
                            result_list.append(value)
                elif isinstance(values, key):
                    result_list.append(values)

            result_list_len = len(result_list)
            if not as_list:
                if result_list_len == 1:
                    return result_list[0]
                elif result_list_len == 0:
                    return default

            return result_list

        raise TypeError('Key must be str or annotation class.')

    def __validate_shape(self, quantity_def: 'Quantity', value):
        if quantity_def == Quantity.default:
            return True

        quantity_shape = quantity_def.shape

        if type(value) == np.ndarray:
            value_shape = value.shape
        if isinstance(value, list) and not isinstance(value, MEnum):
            value_shape = [len(value)]
        else:
            value_shape = []

        if len(value_shape) != len(quantity_shape):
            return False

        for i in range(0, len(value_shape)):
            if not _Dimension.check_dimension(self, quantity_shape[i], value_shape[i]):
                return False

        return True

    def m_validate(self):
        ''' Evaluates all constraints and shapes of this section and returns a list of errors. '''
        errors: List[str] = []
        for constraint_name in self.m_def.constraints:
            constraint = getattr(self, 'c_%s' % constraint_name, None)
            if constraint is None:
                raise MetainfoError(
                    'Could not find implementation for contraint %s of section %s.' %
                    (constraint_name, self.m_def))

            try:
                constraint()
            except AssertionError as e:
                error_str = str(e).strip()
                if error_str == '':
                    error_str = 'Constraint %s violated.' % constraint_name
                errors.append(error_str)

        for quantity in self.m_def.all_quantities.values():
            if self.m_is_set(quantity):
                if not self.__validate_shape(quantity, self.m_get(quantity)):
                    errors.append(
                        'The shape of quantity %s does not match its value.' % quantity)

        return errors

    def m_all_validate(self):
        ''' Evaluates all constraints in the whole section hierarchy, incl. this section. '''
        errors: List[str] = []
        for section in itertools.chain([self], self.m_all_contents()):
            for error in section.m_validate():
                errors.append(error)

        return errors

    def __repr__(self):
        m_section_name = self.m_def.name
        # name_quantity_def = self.m_def.all_quantities.get('name', None)
        # if name_quantity_def is not None:
        #     name = self.m_get(name_quantity_def)
        try:
            name = self.__dict__['name']
        except KeyError:
            name = '<noname>'

        return '%s:%s' % (name, m_section_name)

    def __getitem__(self, key):
        key = key.replace('.', '/')
        return self.m_resolve(key)

    def __iter__(self):
        return self.m_def.all_properties.__iter__()

    def __len__(self):
        return len(self.m_def.all_properties)


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
        pkg.m_add_sub_section(Package.category_definitions, cls.m_def)


# Metainfo M3 (i.e. definitions of definitions)

class Definition(MSection):
    ''' A common base for all metainfo definitions.

    All metainfo `definitions` (sections, quantities, sub-sections, packages, ...) share
    some common attributes. These are defined in a common base: all
    metainfo items extend this common base and inherit from ``Definition``.

    Args:
        name: Each `definition` has a name. Names have to be valid Python identifier.
            They can contain letters, numbers and _, but must not start with a number.
            This also qualifies them as identifier in most storage formats, databases,
            makes them URL safe, etc.

            Names must be unique within the :class:`Package` or :class:`Section` that
            this definition is part of.

        description: The description can be an arbitrary human readable text that explains
            what this definition is about.

        links: Each definition can be accompanied by a list of URLs. These should point
            to resources that further explain the definition.

        categories: All metainfo definitions can be put into one or more `categories`.
            Categories allow to organize the definitions themselves. It is different from
            sections, which organize the data (e.g. quantity values) and not the definitions
            of data (e.g. quantities definitions). See :ref:`metainfo-categories` for more details.

    Additional helper functions for `definitions`:

    .. automethod:: all_definitions
    '''

    __all_definitions: Dict[Type[MSection], List[MSection]] = {}

    name: 'Quantity' = None
    description: 'Quantity' = None
    links: 'Quantity' = None
    categories: 'Quantity' = None

    def __init__(self, *args, **kwargs):
        self.all_categories: Set[Category] = set()

        super().__init__(*args, **kwargs)

        for cls in self.__class__.mro() + [self.__class__]:
            definitions = Definition.__all_definitions.setdefault(cls, [])
            definitions.append(self)

    def __init_metainfo__(self):
        '''
        An initialization method that is called after the class context of the definition
        has been initialized. For example it is called on all quantities of a section
        class after the class was created. If metainfo definitions are created without
        a class context, this method must be called manually on all definitions.
        '''

        # initialize definition annotations
        for annotation in self.m_x(DefinitionAnnotation, as_list=True):
            annotation.init_annotation(self)

    @classmethod
    def all_definitions(cls: Type[MSectionBound]) -> Iterable[MSectionBound]:
        ''' Class method that returns all definitions of this class.

        This can be used to get a list of all globally available `defintions` or a certain
        kind. E.g. to get all `quantities`: ``Quantity.all_definitions()``.
        '''
        return cast(Iterable[MSectionBound], Definition.__all_definitions.get(cls, []))

    def qualified_name(self):
        names = []
        current = self
        while current is not None and current.m_follows(Definition.m_def):
            names.append(current.name)
            current = current.m_parent

        return '.'.join(reversed(names))

    def on_set(self, quantity_def, value):
        if quantity_def == Definition.categories:
            for category in value:
                category.definitions.add(self)
                self.all_categories.add(category)
                for base_category in category.all_categories:
                    self.all_categories.add(base_category)


class Property(Definition):
    ''' A common base-class for section properties: sub sections and quantities. '''
    pass


class Quantity(Property):
    ''' Definition of an atomic piece of data.

    Quantity definitions are the main building block of meta-info schemas. Each quantity
    represents a single piece of data.

    To define quantities, use objects of this class as classattribute values in
    `section classes`. The name of a quantity is automatically taken from its `section class`
    attribute. You can provide all other attributes to the constructor with keyword arguments

    See :ref:`metainfo-sections` to learn about `section classes`.
    In Python terms, ``Quantity`` is a descriptor. Descriptors define how to get and
    set attributes in a Python object. This allows us to use sections like regular
    Python objects and quantity like regular Python attributes.

    Beyond basic :class:`Definition` attributes, Quantities are defined with the following
    attributes.

    Args:
        type:
            Defines the datatype of quantity values. This is the type of individual elements
            in a potentially complex shape. If you define a list of integers for example,
            the `shape` would be list and the `type` integer:
            ``Quantity(type=int, shape=['0..*'])``.

            The `type` can be one of:

            - a build-in primitive Python type: ``int``, ``str``, ``bool``, ``float``
            - an instance of :class:`MEnum`, e.g. ``MEnum('one', 'two', 'three')``
            - a section to define references to other sections as quantity values
            - a custom meta-info :class:`DataType`, see :ref:`metainfo-custom-types`
            - a numpy `dtype`, e.g. ``np.dtype('float32')``
            - ``typing.Any`` to support any value

            If set to `dtype`, this quantity will use a numpy array or scalar to store values
            internally. If a regular (nested) Python list or Python scalar is given, it will
            be automatically converted. The given `dtype` will be used in the numpy value.

            To define a reference, either a `section class` or instance of :class:`Section`
            can be given. See :ref:`metainfo-sections` for details. Instances of the given section
            constitute valid values for this type. Upon serialization, references section
            instance will represented with metainfo URLs. See :ref:`metainfo-urls`.

            For quantities with more than one dimension, only numpy arrays and `dtypes`
            are allowed.

        shape:
            The shape of the quantity. It defines its dimensionality.

            A shape is a list, where each item defines one dimension.
            Each dimension can be:

            - an integer that defines the exact size of the dimension, e.g. ``[3]`` is the
              shape of a 3D spacial vector
            - a string that specifies a possible range, e.g. ``0..*``, ``1..*``, ``3..6``
            - the name of an int typed and shapeless quantity in the same section which
              values define the length of this dimension, e.g. ``number_of_atoms`` defines
              the length of ``atom_positions``

            Range specifications define lower and upper bounds for the possible dimension
            length. The ``*`` can be used to denote an arbitrarily high upper bound.

            Quantities with dimensionality (length of the shape) higher than 1, must be
            numpy arrays. Theire type must be a `dtype`.

        unit:
            The physics unit for this quantity. It is optional.

            Units are represented with the pint_ Python package. Pint defines units and
            their algebra. You can either use `pint` units directly, e.g. ``units.m / units.s``.
            The metainfo provides a preconfigured `pint` unit registry :py:data:`units`.
            You can also provide the unit as `pint` parsable string, e.g. ``'meter / seconds'`` or
            ``'m/s'``.

        default:
            The default value for this quantity. The value must match type and shape.

            Be careful with a default value like ``[]`` as it will be the default value for
            all occurences of this quantity.

        synonym_for:
            The name of a quantity defined in the same section as string. This will make
            this quantity a synonym for the other quantity. All other properties (type,
            shape, unit, etc.) are ignored. Getting or setting from/to this quantity will
            be delegated to the other quantity. Synonyms are always virtual.

        derived:
            A Python callable that takes the containing section as input and outputs the
            value for this quantity. This quantity cannot be set directly, its value
            is only derived by the given callable. The callable is executed when this
            quantity is get. Derived quantities are always virtual.

        virtual:
            A boolean that determines if this quantity is virtual. Virtual quantities can
            be get/set like regular quantities, but their values are not (de-)serialized,
            hence never permanently stored.

        is_scalar:
            Derived quantity that is True, iff this quantity has shape of length 0
        '''

    type: 'Quantity' = None
    shape: 'Quantity' = None
    unit: 'Quantity' = None
    default: 'Quantity' = None
    synonym_for: 'Quantity' = None
    derived: 'Quantity' = None
    virtual: 'Quantity' = None
    is_scalar: 'Quantity' = None

    # TODO derived_from = Quantity(type=Quantity, shape=['0..*'])
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._synonym = None

    def __init_metainfo__(self):
        super().__init_metainfo__()

        if self.derived is not None:
            self.virtual = True

        if self.synonym_for is not None:
            self._synonym = self.m_parent.all_quantities[self.synonym_for]

        # replace the quantity implementation with an optimized version for the most
        # primitive quantities if applicable
        is_primitive = not self.synonym_for and not self.derived
        is_primitive = is_primitive and len(self.shape) <= 1
        is_primitive = is_primitive and self.type in [str, bool, float, int]
        if is_primitive:
            self._default = self.default
            self._name = self.name
            self._type = self.type
            self._list = len(self.shape) == 1
            self.__class__ = PrimitiveQuantity

    def __get__(self, obj, cls):
        try:
            value = obj.__dict__[self.name]

        except KeyError:
            if self._synonym:
                return obj.m_get(self._synonym)

            if self.derived is not None:
                try:
                    return self.derived(obj)  # pylint: disable=not-callable
                except Exception as e:
                    raise DeriveError('Could not derive value for %s: %s' % (self, str(e)))
            value = self.default

        except AttributeError:
            # class (def) attribute case, because obj is None
            return self

        if value is None:
            return value

        if isinstance(self.type, DataType) and self.type.get_normalize != DataType.get_normalize:
            dimensions = len(self.shape)
            if dimensions == 0:
                value = self.type.get_normalize(obj, self, value)

            elif dimensions == 1:
                value = list(
                    self.type.get_normalize(obj, self, item)
                    for item in value)

            else:
                raise MetainfoError(
                    'Only numpy arrays and dtypes can be used for higher dimensional '
                    'quantities.')

            if isinstance(self.type, Reference):
                # save the resolved values for the next access to avoid re-resolve
                self.__dict__[self.name] = value

        elif type(self.type) == np.dtype:
            if self.unit is not None:
                value = value * self.unit

        return value

    def __set__(self, obj, value):
        if obj is None:
            # class (def) case
            raise KeyError('Cannot overwrite quantity definition. Only values can be set.')

        # object (instance) case
        obj.m_set(self, value)

    def __delete__(self, obj):
        if obj is None:
            # class (def) case
            raise KeyError('Cannot delete quantity definition. Only values can be deleted.')

        # object (instance) case
        raise NotImplementedError('Deleting quantity values is not supported.')

    def c_dimensions(self):
        for dimension in self.shape:  # pylint: disable=not-an-iterable
            if isinstance(dimension, str):
                if dimension.isidentifier():
                    dim_quantity = self.m_parent.all_quantities.get(dimension, None)

                    assert dim_quantity is not None, \
                        'Dimensions (%s) must be quantities of the same section (%s).' % (
                            dimension, self.m_parent)

                    assert len(dim_quantity.shape) == 0 and \
                        dim_quantity.type in [int, np.int16, np.int32, np.int8, np.uint8], \
                        'Dimensions (%s) must be shapeless (%s) and int (%s) typed.' % (
                            dimension, dim_quantity.shape, dim_quantity.type)

    def c_higher_shapes_require_dtype(self):
        if len(self.shape) > 1:
            assert type(self.type) == np.dtype, \
                'Higher dimensional quantities (%s) need a dtype and will be treated as ' \
                'numpy arrays.' % self


class DirectQuantity(Quantity):
    ''' Used for quantities that would cause indefinite loops due to bootstrapping. '''

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._name = kwargs.get('name')
        self._default = kwargs.get('default')

    def __get__(self, obj, cls):
        try:
            return obj.__dict__[self._name]
        except KeyError:
            return self._default

        except AttributeError:
            # class (def) attribute case
            return self

    def __set__(self, obj, value):
        if obj is None:
            # class (def) case
            raise KeyError('Cannot overwrite quantity definition. Only values can be set.')

        # object (instance) case
        obj.__dict__[self._name] = value


class PrimitiveQuantity(Quantity):
    ''' An optimized replacement for Quantity suitable for primitive properties. '''
    def __get__(self, obj, cls):
        try:
            return obj.__dict__[self._name]
        except KeyError:
            return self._default
        except AttributeError:
            return self

    def __set__(self, obj, value):
        if self._list:
            if not isinstance(value, list):
                raise TypeError(
                    'The value %s for quantity %s has not shape %s' %
                    (value, self, self.shape))

            if any(v is not None and type(v) != self._type for v in value):
                raise TypeError(
                    'The value %s with type %s for quantity %s is not of type %s' %
                    (value, type(value), self, self.type))

        elif value is not None and type(value) != self._type:
            raise TypeError(
                'The value %s with type %s for quantity %s is not of type %s' %
                (value, type(value), self, self.type))

        try:
            obj.__dict__[self._name] = value
        except AttributeError:
            raise KeyError('Cannot overwrite quantity definition. Only values can be set.')


class SubSection(Property):
    ''' Defines what sections can appear as sub-sections of another section.

    Like quantities, sub-sections are defined in a `section class` as attributes
    of this class. An like quantities, each sub-section definition becomes a property of
    the corresponding `section definition` (parent). A sub-section definition references
    another `section definition` as the sub-section (child). As a consequence, parent
    `section instances` can contain child `section instances` as sub-sections.

    Contrary to the old NOMAD metainfo, we distinguish between sub-section the section
    and sub-section the property. This allows to use on child `section definition` as
    sub-section of many different parent `section definitions`.

    Args:
        sub_section: A :class:`Section` or Python class object for a `section class`. This
            will be the child `section definition`. The defining section the child
            `section definition`.

        repeats: A boolean that determines wether this sub-section can appear multiple
            times in the parent section.
    '''

    sub_section: 'Quantity' = None
    repeats: 'Quantity' = None

    def __get__(self, obj, type=None):
        if obj is None:
            # the class attribute case
            return self

        else:
            # the object attribute case
            if self.repeats:
                return obj.m_get_sub_sections(self)
            else:
                return obj.m_get_sub_section(self, -1)

    def __set__(self, obj, value):
        if obj is None:
            raise NotImplementedError()

        if self.repeats:
            raise NotImplementedError('Cannot set a repeating sub section use m_create or m_add_sub_section.')

        if value is None:
            obj.m_remove_sub_section(self, -1)
        else:
            obj.m_add_sub_section(self, value)

    def __delete__(self, obj):
        raise NotImplementedError('Deleting sub sections is not supported.')


class Section(Definition):
    ''' Sections define blocks of related quantities and allows hierarchical data.

    Section definitions determine what quantities and sub-sections can appear in a
    following section instance.

    Args:
        quantities:
            The quantities definitions of this section definition as list of :class:`Quantity`.
            Will be automatically set from the `section class`.

        sub_sections:
            The sub-section definitions of this section definition as list of :class:`SubSection`.
            Will be automatically set from the `section class`.

        base_sections:
            A list of `section definitions` (:class:`Section`). By default this definition will
            inherit all quantity and sub section definitions from the given section definitions.
            This behavior might be altered with ``extends_base_section``.

        extends_base_section:
            If True, this definition must have exactly one ``base_sections``.
            Instead of inheriting properties, he quantity and sub-section definitions
            of this section will be added to the base section.

            This allows to add further properties to an existing section definition.
            To use such extension on section instances in a type-safe manner
            :py:func:`MSection.m_as` can be used to cast the base section to the extending
            section.

        constraints:
            Constraints are rules that a section must fulfil to be valid. This allows to implement
            semantic checks that go behind mere type or shape checks. This quantity takes
            the names of constraints as string. Constraints have to be implemented as methods of the
            section class. These constraints methods must be named ``c_<constraint name>``
            and have no additional parameters. They can raise :class:`ConstraintVialated` or
            an AssertionError to indicate that the constraint is not fulfilled for the ``self``
            section. This quantity will be set automatically from all ``c_`` methods in the
            respective section class. To run validation of a section use :py:meth:`MSection.m_validate`.

        event_handlers:
            Event handler are functions that get called when the section data is changed.
            There are two types of events: ``set`` and ``add_sub_section``. The handler type
            is determined by the handler (i.e. function) name: ``on_set`` and ``on_add_sub_section``.
            The handler arguments correspond to :py:meth:`MSection.m_set` (section, quantity_def, value) and
            :py:meth:`MSection.m_add_sub_section` (section, sub_section_def, sub_section).
            Handler are called after the respective action was performed. This quantity is
            automatically populated with handler from the section classes methods. If there
            is a method ``on_set`` or ``on_add_sub_section``, it will be added as handler.

        section_cls:
            A helper attribute that gives the `section class` as a Python class object.

        all_base_sections:
            A helper attribute that gives direct and indirect base sections.

        all_properties:
            A helper attribute that gives all properties (sub section and quantity) definitions
            including inherited properties as a dictionary with names and definitions.

        all_quantities:
            A helper attribute that gives all quantity definition including inherited ones
            as a dictionary that maps names (strings) to :class:`Quantity`.

        all_sub_sections:
            A helper attribute that gives all sub-section definition including inherited ones
            as a dictionary that maps names (strings) to :class:`SubSection`.

        all_sub_sections_by_section:
            A helper attribute that gives all sub-section definition including inherited ones
            as a dictionary that maps section classes (i.e. Python class objects) to
            lists of :class:`SubSection`.

        parent_section_sub_section_defs:
            A helper attribute that gives all sub-section definitions that this section
            is used in.
    '''

    section_cls: Type[MSection] = None

    quantities: 'SubSection' = None
    sub_sections: 'SubSection' = None

    base_sections: 'Quantity' = None
    extends_base_section: 'Quantity' = None
    constraints: 'Quantity' = None
    event_handlers: 'Quantity' = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.all_base_sections: Set['Section'] = set()
        self.all_properties: Dict[str, Union['SubSection', Quantity]] = dict()
        self.all_quantities: Dict[str, Quantity] = dict()
        self.all_sub_sections: Dict[str, SubSection] = dict()
        self.all_sub_sections_by_section: Dict['Section', List['SubSection']] = dict()
        self.parent_section_sub_section_defs: List['SubSection'] = list()

    def __init_metainfo__(self):
        # initialize properties
        for property in self.all_properties.values():
            property.__init_metainfo__()

        super().__init_metainfo__()

        # validate
        def validate(definition):
            errors = definition.m_all_validate()
            if len(errors) > 0:
                raise MetainfoError(
                    '%s. The section definition %s violates %d more constraints' %
                    (str(errors[0]).strip('.'), definition, len(errors) - 1))

        if self.extends_base_section:
            validate(self.m_get(Section.base_sections)[0])
        else:
            validate(self)

    def on_add_sub_section(self, sub_section_def, sub_section):
        if sub_section_def == Section.quantities:
            self.all_properties[sub_section.name] = sub_section
            self.all_quantities[sub_section.name] = sub_section

        if sub_section_def == Section.sub_sections:
            self.all_properties[sub_section.name] = sub_section
            self.all_sub_sections[sub_section.name] = sub_section
            self.all_sub_sections_by_section.setdefault(
                sub_section.sub_section, []).append(sub_section)

            if isinstance(sub_section, SubSection):
                sub_section.sub_section.parent_section_sub_section_defs.append(sub_section)

    def on_set(self, quantity_def, value):
        if quantity_def == Section.base_sections:
            for base_section in value:
                self.all_base_sections.add(base_section)
                for base_base_section in base_section.all_base_sections:
                    self.all_base_sections.add(base_base_section)

                if not self.extends_base_section:
                    self.all_properties.update(**base_section.all_properties)
                    self.all_quantities.update(**base_section.all_quantities)
                    self.all_sub_sections.update(**base_section.all_sub_sections)
                    self.all_sub_sections_by_section.update(**base_section.all_sub_sections_by_section)

    def c_unique_names(self):
        # start with the names of all base_sections
        names: Set[str] = set(
            name
            for base in self.all_base_sections
            for name in base.all_properties.keys())

        for def_list in [self.quantities, self.sub_sections]:
            for definition in def_list:
                assert definition.name not in names, 'All names in a section must be unique. ' \
                    'Name %s of %s in %s already exists in %s.' % (definition.name, definition, definition.m_parent, self)
                names.add(definition.name)


class Package(Definition):
    ''' Packages organize metainfo defintions alongside Python modules

    Each Python module with metainfo Definition (explicitely or implicitely) has a member
    ``m_package`` with an instance of this class. Definitions (categories, sections) in
    Python modules are automatically added to the module's :class:`Package`.
    Packages are not nested and rather have the fully qualitied Python module name as
    name.

    This allows to inspect all definitions in a Python module and automatically puts
    module name and docstring as :class:`Package` name and description.

    Besides the regular :class:`Defintion` attributes, packages can have the following
    attributes:

    Args:
        section_definitions: All `section definitions` in this package as :class:`Section`
            objects.

        category_definitions: All `category definitions` in this package as :class:`Category`
            objects.

        all_definitions: A helper attribute that provides all section definitions
            by name.
    '''

    section_definitions: 'SubSection' = None
    category_definitions: 'SubSection' = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.all_definitions: Dict[str, 'Section'] = dict()

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

    def on_add_sub_section(self, sub_section_def, sub_section):
        if sub_section_def in [Package.section_definitions, Package.category_definitions]:
            self.all_definitions[sub_section.name] = sub_section


class Category(Definition):
    ''' Categories allow to organize metainfo definitions (not metainfo data like sections do)

    Each definition, including categories themselves, can belong to a set of categories.
    Categories therefore form a hierarchy of concepts that definitions can belong to, i.e.
    they form a `is a` relationship.

    Args:
        definitions: A helper attribute that gives all definitions that are directly or
            indirectly in this category.
    '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.definitions: Set[Definition] = set()


class Annotation:
    ''' Base class for annotations. '''
    pass


class DefinitionAnnotation(Annotation):
    ''' Base class for annotations for definitions. '''

    def __init__(self):
        self.definition: Definition = None

    def init_annotation(self, definition: Definition):
        self.definition = definition


class SectionAnnotation(DefinitionAnnotation):
    '''
    Special annotation class for section definition that allows to auto add annotations
    to section instances.
    '''

    def new(self, section) -> Dict[str, Any]:
        return {}


Section.m_def = Section(name='Section')
Section.m_def.m_def = Section.m_def
Section.m_def.section_cls = Section

Definition.m_def = Section(name='Definition')
Property.m_def = Section(name='Property')
Quantity.m_def = Section(name='Quantity')
SubSection.m_def = Section(name='SubSection')
Category.m_def = Section(name='Category')
Package.m_def = Section(name='Package')

Definition.name = DirectQuantity(type=str, name='name')
Definition.description = Quantity(type=str, name='description')
Definition.links = Quantity(type=str, shape=['0..*'], name='links')
Definition.categories = Quantity(
    type=Reference(Category.m_def), shape=['0..*'], default=[], name='categories')

Section.quantities = SubSection(
    sub_section=Quantity.m_def, name='quantities', repeats=True)

Section.sub_sections = SubSection(
    sub_section=SubSection.m_def, name='sub_sections', repeats=True)
Section.base_sections = Quantity(
    type=Reference(Section.m_def), shape=['0..*'], default=[], name='base_sections')
Section.extends_base_section = Quantity(type=bool, default=False, name='extends_base_section')
Section.constraints = Quantity(type=str, shape=['0..*'], default=[], name='constraints')
Section.event_handlers = Quantity(
    type=Callable, shape=['0..*'], name='event_handlers', virtual=True, default=[])

SubSection.repeats = Quantity(type=bool, name='repeats', default=False)

SubSection.sub_section = Quantity(type=Reference(Section.m_def), name='sub_section')

Quantity.m_def.section_cls = Quantity
Quantity.type = DirectQuantity(type=QuantityType, name='type')
Quantity.shape = DirectQuantity(type=Dimension, shape=['0..*'], name='shape', default=[])
Quantity.unit = Quantity(type=Unit, name='unit')
Quantity.default = DirectQuantity(type=Any, default=None, name='default')
Quantity.synonym_for = DirectQuantity(type=str, name='synonym_for')
Quantity.derived = DirectQuantity(type=Callable, default=None, name='derived', virtual=True)
Quantity.virtual = DirectQuantity(type=bool, default=False, name='virtual')
Quantity.is_scalar = Quantity(
    type=bool, name='is_scalar', derived=lambda quantity: len(quantity.shape) == 0)

Package.section_definitions = SubSection(
    sub_section=Section.m_def, name='section_definitions', repeats=True)

Package.category_definitions = SubSection(
    sub_section=Category.m_def, name='category_definitions', repeats=True)

is_bootstrapping = False

Section.m_def.event_handlers = [Section.on_add_sub_section, Section.on_set]

Definition.__init_cls__()
Property.__init_cls__()
Package.__init_cls__()
Section.__init_cls__()
Category.__init_cls__()
Quantity.__init_cls__()
SubSection.__init_cls__()


class Environment(MSection):
    ''' Environments allow to manage many metainfo packages and quickly access all definitions.

    Environments provide a name-table for large-sets of metainfo definitions that span
    multiple packages. It provides various functions to resolve metainfo definitions by
    their names, legacy names, and qualified names.

    Args:
        packages: Packages in this environment.
    '''

    packages = SubSection(sub_section=Package, repeats=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.all_definitions_by_name: Dict[str, List[Definition]] = dict()
        self.all_packages: Dict[str, Package] = dict()

    def resolve_definitions(  # type: ignore
            self, name: str, cls: Type[MSectionBound] = Definition) -> List[MSectionBound]:

        return [
            cast(MSectionBound, definition)
            for definition in self.all_definitions_by_name.get(name, [])
            if isinstance(definition, cls)]

    def resolve_definition(  # type: ignore
            self, name, cls: Type[MSectionBound] = Definition) -> MSectionBound:

        defs = self.resolve_definitions(name, cls)
        if len(defs) == 1:
            return defs[0]
        elif len(defs) > 1:
            raise KeyError('Could not uniquely identify %s' % name)
        else:
            raise KeyError('Could not resolve %s' % name)

    def on_add_sub_section(self, sub_section_def: SubSection, sub_section: MSection):
        if sub_section_def == Environment.packages:
            package = sub_section.m_as(Package)
            self.all_packages[package.name] = package
            for definition in package.m_all_contents():
                if isinstance(definition, Definition):
                    definitions = self.all_definitions_by_name.setdefault(definition.name, [])
                    definitions.append(definition)
