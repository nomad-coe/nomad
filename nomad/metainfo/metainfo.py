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

from typing import Type, TypeVar, Union, Tuple, Iterable, List, Any, Dict, Set, \
    Callable as TypingCallable, cast
from collections.abc import Iterable as IterableABC, Sequence
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
import jmespath

from nomad.units import ureg as units


m_package: 'Package' = None

is_bootstrapping = True
MSectionBound = TypeVar('MSectionBound', bound='MSection')
T = TypeVar('T')

# Make pylint believe all bootstrap quantities are actual properties even though
# we have to initialize them to None due to bootstrapping
_placeholder_quantity: 'Quantity' = property()  # type: ignore
if True:
    _placeholder_quantity: 'Quantity' = None  # type: ignore

_primitive_types = {str: str, int: int, float: float, bool: bool, np.bool_: bool}


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

class MEnum(Sequence):
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

        self._list = list(kwargs.values())
        self._values = set(kwargs.values())  # For allowing constant time member check
        self._map = kwargs

    def __getattr__(self, attr):
        return self._map[attr]

    def __getitem__(self, index):
        return self._list[index]

    def __len__(self):
        return len(self._list)


class MProxy():
    ''' A placeholder object that acts as reference to a value that is not yet resolved.

    Attributes:
        url: The reference represented as an URL string.
    '''

    def __init__(
            self, m_proxy_value: Union[str, int, dict], m_proxy_section: 'MSection' = None,
            m_proxy_quantity: 'Quantity' = None):
        self.m_proxy_value = m_proxy_value
        self.m_proxy_section = m_proxy_section
        self.m_proxy_resolved = None
        self.m_proxy_quantity = m_proxy_quantity

    def m_proxy_resolve(self):
        if self.m_proxy_section and self.m_proxy_quantity and not self.m_proxy_resolved:
            self.m_proxy_resolved = self.m_proxy_quantity.type.resolve(self)

        if self.m_proxy_resolved is not None and isinstance(self, MProxy):
            setattr(self, '__class__', self.m_proxy_resolved.__class__)
            self.__dict__.update(**self.m_proxy_resolved.__dict__)

        return self.m_proxy_resolved

    def __getattr__(self, key):
        if self.m_proxy_resolve() is not None:
            return getattr(self.m_proxy_resolved, key)

        raise ReferenceError('could not resolve %s' % self.m_proxy_value)


class SectionProxy(MProxy):
    def m_proxy_resolve(self):
        if self.m_proxy_section and not self.m_proxy_resolved:
            root = self.m_proxy_section
            while root.m_parent is not None and not isinstance(root, Package):
                root = root.m_parent

            if isinstance(root, Package):
                self.m_proxy_resolved = root.all_definitions.get(self.m_proxy_value)

            if self.m_proxy_resolved is None:
                raise ReferenceError('could not resolve %s' % self.m_proxy_value)

        return self.m_proxy_resolved


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

        if isinstance(value, str):
            # TODO raise a warning or allow this?
            # In the old metainfo there are cases where an expression is used
            # that is later evaluated in the parser
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
        if value in _primitive_types:
            return value

        if isinstance(value, MEnum):
            for enum_value in value._values:
                if not isinstance(enum_value, str):
                    raise TypeError('MEnum value %s is not a string.' % enum_value)
            return value

        if isinstance(value, np.dtype):
            return value

        if isinstance(value, Section):
            return value

        if isinstance(value, Reference) and isinstance(value.target_section_def, MProxy):
            proxy = value.target_section_def
            proxy.m_proxy_section = section
            proxy.m_proxy_quantity = quantity_def
            return value

        if isinstance(value, DataType):
            return value

        if value == Any:
            return value

        if isinstance(value, type):
            section = getattr(value, 'm_def', None)
            if section is not None:
                return Reference(section)

        if isinstance(value, Quantity):
            return QuantityReference(value)

        if isinstance(value, MProxy):
            value.m_proxy_section = section
            value.m_proxy_quantity = quantity_def
            return value

        raise MetainfoError(
            'Type %s of %s is not a valid metainfo quantity type' %
            (value, quantity_def))

    def serialize(self, section, quantity_def, value):
        if value in _primitive_types:
            return dict(type_kind='python', type_data=value.__name__)

        if isinstance(value, MEnum):
            return dict(type_kind='Enum', type_data=list(value))

        if isinstance(value, np.dtype):
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

    def __init__(self, section_def: Union['Section', 'SectionProxy']):
        self.target_section_def = section_def

    def resolve(self, proxy) -> 'MSection':
        '''
        Resolve the given proxy. The proxy is guaranteed to have a context and
        will to be not yet resolved.
        '''
        return proxy.m_proxy_section.m_resolve(proxy.m_proxy_value)

    def set_normalize(self, section: 'MSection', quantity_def: 'Quantity', value: Any) -> Any:
        if isinstance(self.target_section_def, MProxy):
            proxy = self.target_section_def
            proxy.m_proxy_section = section.m_def
            proxy.m_proxy_quantity = quantity_def
            self.target_section_def = proxy.m_proxy_resolve()

        if self.target_section_def.m_follows(Definition.m_def):
            # special case used in metainfo definitions, where we reference metainfo definitions
            # using their Python class. E.g. referencing a section definition using its
            # class instead of the object: Run vs. Run.m_def
            if isinstance(value, type):
                definition = getattr(value, 'm_def', None)
                if definition is not None and definition.m_follows(self.target_section_def):
                    return definition

        if isinstance(value, (str, int, dict)):
            return MProxy(value, m_proxy_section=section, m_proxy_quantity=quantity_def)

        if isinstance(value, MProxy):
            value.m_proxy_section = section
            value.m_proxy_quantity = quantity_def
            return value

        if not isinstance(value, MSection):
            raise TypeError(
                'The value %s is not a section and can not be used as a reference.' % value)

        if not value.m_follows(self.target_section_def):  # type: ignore
            raise TypeError(
                '%s is not a %s and therefore an invalid value of %s.' %
                (value, self.target_section_def, quantity_def))

        return value

    def serialize(self, section: 'MSection', quantity_def: 'Quantity', value: Any) -> Any:
        return value.m_path()

    def deserialize(self, section: 'MSection', quantity_def: 'Quantity', value: Any) -> Any:
        return MProxy(value, m_proxy_section=section, m_proxy_quantity=quantity_def)


class QuantityReference(Reference):
    ''' Datatype used for reference quantities that reference other quantities. '''

    def __init__(self, quantity_def: Union['Quantity']):
        super().__init__(cast(Section, quantity_def.m_parent))
        self.target_quantity_def = quantity_def

    def get_normalize(self, section: 'MSection', quantity_def: 'Quantity', value: Any) -> Any:
        section = super().get_normalize(section, quantity_def, value)
        return getattr(section, self.target_quantity_def.name)

    def serialize(self, section: 'MSection', quantity_def: 'Quantity', value: Any) -> Any:
        section_path = super().serialize(section, quantity_def, value)
        return f'{section_path}/{self.target_quantity_def.name}'

    def deserialize(self, section: 'MSection', quantity_def: 'Quantity', value: Any) -> Any:
        section_path = value.rsplit('/', 1)[0]
        return MProxy(section_path, m_proxy_section=section, m_proxy_quantity=quantity_def)


class _Datetime(DataType):

    def _parse(self, datetime_str: str) -> datetime:
        try:
            return aniso8601.parse_datetime(datetime_str)
        except ValueError:
            pass

        try:
            date = aniso8601.parse_date(datetime_str)
            if isinstance(date, datetime):
                return date
        except ValueError as e:
            pass

        try:
            # TODO necessary?
            import flask_restplus.inputs
            return flask_restplus.inputs.datetime_from_rfc822(datetime_str)
        except ValueError:
            pass

        try:
            return datetime.strptime(datetime_str, '%Y-%m-%d %H:%M:%S.%f')
        except ValueError:
            pass

        try:
            return datetime.strptime(datetime_str, '%Y-%m-%d %H:%M:%S')
        except ValueError:
            pass

        try:
            return datetime.strptime(datetime_str, '%Y-%m-%d')
        except ValueError:
            pass

        raise TypeError('Invalid date literal %s' % datetime_str)

    def _convert(self, value):
        if value is None:
            return None

        if isinstance(value, str):
            value = self._parse(value)

        elif isinstance(value, (int, float)):
            value = datetime.fromtimestamp(value)

        elif isinstance(value, pint.Quantity):
            value = datetime.fromtimestamp(value.magnitude)

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


class _JSON(DataType):
    pass


class _Capitalized(DataType):
    def set_normalize(self, section: 'MSection', quantity_def: 'Quantity', value: Any) -> Any:
        if value is not None and len(value) >= 1:
            return value[0].capitalize() + value[1:]

        return value


Dimension = _Dimension()
Unit = _Unit()
QuantityType = _QuantityType()
Callable = _Callable()
Datetime = _Datetime()
JSON = _JSON()
Capitalized = _Capitalized()


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


def constraint(warning):
    ''' A decorator for methods implementing constraints. '''
    f = None
    if not isinstance(warning, bool):
        f = warning
        warning = False

    def decorator(f):
        setattr(f, 'm_constraint', True)
        setattr(f, 'm_warning', warning)
        return f

    if f is None:
        return decorator
    else:
        return decorator(f)


class MResource():
    '''
    Represents a collection of related metainfo data, i.e. a set of :class:`MSection` instances.
    '''

    def __init__(self, logger=None):
        self.__data: Dict['Section', List['MSection']] = dict()
        self.contents: List['MSection'] = []
        self.logger = logger

    def create(self, section_cls: Type[MSectionBound], *args, **kwargs) -> MSectionBound:
        '''
        Create an instance of the given section class and adds it to this resource as
        a root section. The m_parent_index will be set sequentially among root sections of
        the same section definition starting with 0.
        '''
        index = 0
        for content in self.contents:
            if content.m_follows(section_cls.m_def):
                index = max(index, content.m_parent_index + 1)

        result = section_cls(*args, **kwargs)
        result.m_parent_index = index

        self.add(result)
        return cast(MSectionBound, result)

    def add(self, section):
        '''
        Add the given section to this resource. Will also add all its contents to the
        resource and make all contest available for :func:`all`. Will also remove
        all contents from possible other resources. A section can only be contained in
        one resource at a time.

        This is potentially expensive. Do not add a section that already has a deep tree
        of sub-sections. Ideally, add the root section first. If you create sub sections
        afterwards, they will be automatically added to this resource.
        '''
        if section.m_resource is not None:
            section.m_resource.remove(section)

        for content in section.m_all_contents(include_self=True):
            content.m_resource = self
            self.__data.setdefault(content.m_def, []).append(content)

        if section.m_parent is None:
            self.contents.append(section)

    def remove(self, section):
        assert section.m_resource == self, 'Can only remove section from the resource that contains it.'
        section.m_resource = None
        self.__data.get(section.m_def).remove(section)
        if section.m_parent is None:
            self.contents.remove(section)

    def all(self, section_cls: Type[MSectionBound]) -> List[MSectionBound]:
        ''' Returns all instances of the given section class in this resource. '''
        return cast(List[MSectionBound], self.__data.get(section_cls.m_def, []))

    def unload(self):
        ''' Breaks all references among the contain metainfo sections to allow GC. '''
        for collections in self.__data.values():
            for section in collections:
                section.m_parent = None
                section.__dict__.clear()
            collections.clear()

    def m_to_dict(self, filter: TypingCallable[['MSection'], bool] = None):
        if filter is None:
            def filter(_):  # pylint: disable=function-redefined
                return True

        return {
            section.m_def.name: section.m_to_dict()
            for section in self.contents
            if filter(section)}

    def warning(self, *args, **kwargs):
        if self.logger is not None:
            self.logger.warn(*args, **kwargs)


class MSection(metaclass=MObjectMeta):  # TODO find a way to make this a subclass of collections.abs.Mapping
    '''
    The base-class for all *section defining classes* and respectively the base-class
    for all section objects.

    While we use *section classes* to *define sections*, it is important to note that
    the *section class* is a different Python object than the actual *section definition*.
    For each *section class* (a Python class), we automatically generate a *section definition*
    Python object that instantiates :class:`Section`. :class:`MSection` and :class:`Section`
    are completely different classes. :class:`MSection` is used as a base-class for all
    *section defining* classes and :class:`Section` is a *section class* that defines the
    section `Section`.

    Attributes:
        m_def: Each *section class* (and also *section instance*) has a build-in
            property ``m_def`` that refers to the actual *section definition*. While this defined
            automatically, you can do it manually to provide additional characteristics that cannot
            be covered in a Python class definition.

    All `section instances` indirectly instantiate the :class:`MSection` and therefore all
    members of :class:`MSection` are available on all `section instances`. :class:`MSection`
    provides many special attributes and functions (they all start with ``m_``) that allow
    to reflect on a `section's definition` and allow to manipulate the `section instance`
    without a priori knowledge of the `section defintion`.

    .. automethod:: m_set
    .. automethod:: m_get
    .. automethod:: m_add_values
    .. automethod:: m_get_sub_section
    .. automethod:: m_get_sub_sections
    .. automethod:: m_create
    .. automethod:: m_add_sub_section
    .. automethod:: m_remove_sub_section

    There are some specific attributes for section instances that are sub-sections of
    another section. While sub-sections are directly accessible from the containing
    section by using the Python property that represents the sub-section (e.g.
    `run.section_system`), there is also a way to navigate from the sub-section to
    the containing section (`parent section`) using these Python properties:

    Attributes:
        m_parent:
            If this section is a sub-section, this references the parent section instance.

        m_parent_sub_section:
            If this section is a sub-section, this is the :class:`SubSection` that defines
            this relationship.

        m_parent_index:
            For repeatable sections, parent keep a list of sub-sections. This is the index
            of this section in the respective parent sub-section list.

        m_resource: The :class:`MResource` that contains and manages this section.

    Often some general tasks have to be performed on a whole tree of sections without
    knowing about the definitions in advance. The following methods allow to access
    sub-sections reflectively.

    .. automethod:: m_traverse
    .. automethod:: m_all_contents
    .. automethod:: m_contents
    .. automethod:: m_xpath

    Each section and all its quantities and contents can be transformed into a general
    JSON-serializable Python dictionary. Similarely, a section can be instantiated from
    such a Python dictionary. This allows to save and load sections to JSON-files or
    by other compatible means (e.g. document databases, binary JSON flavours).

    .. automethod:: m_to_dict
    .. automethod:: m_from_dict
    .. automethod:: m_update_from_dict
    .. automethod:: m_to_json
    '''

    m_def: 'Section' = None

    def __init__(
            self, m_def: 'Section' = None, m_resource: MResource = None, **kwargs):

        self.m_def: 'Section' = m_def
        self.m_parent: 'MSection' = None
        self.m_parent_sub_section: 'SubSection' = None
        self.m_parent_index = -1
        self.m_resource = m_resource
        self.m_mod_count = 0
        self.m_cache: dict = {}  # Dictionary for caching temporary values that are not persisted to the Archive

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
            for section_annotation in self.m_def.m_get_annotations(SectionAnnotation, as_list=True):
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
        m_def = cls.__dict__.get('m_def')  # do not accedentally get the m_def from a potential base section
        if m_def is None:
            m_def = Section()
            setattr(cls, 'm_def', m_def)

        # Use class name if name is not explicitly defined
        if m_def.name is None:
            m_def.name = cls.__name__
        m_def._section_cls = cls

        # add base sections
        base_sections: List[Section] = []
        for base_cls in cls.__bases__:
            if base_cls != MSection:
                base_section = getattr(base_cls, 'm_def')
                if base_section is None:
                    raise TypeError(
                        'Section defining classes must have MSection or a decendant as '
                        'base classes.')
                base_sections.append(base_section)

        m_def.m_set(Section.base_sections, base_sections)

        # transfer names, descriptions, constraints, event_handlers
        constraints: Set[str] = set()
        event_handlers: Set[Callable] = set(m_def.event_handlers)
        for name, attr in cls.__dict__.items():
            # transfer names and descriptions for properties, init properties
            if isinstance(attr, Property):
                attr.name = name
                if attr.description is not None:
                    description = inspect.cleandoc(attr.description)
                    description = description.strip()
                    description = re.sub(
                        r'\(https?://[^\)]*\)',
                        lambda m: re.sub(r'\n', '', m.group(0)),
                        description)
                    attr.description = description
                    attr.__doc__ = attr.description

                if isinstance(attr, Quantity):
                    m_def.m_add_sub_section(Section.quantities, attr)
                elif isinstance(attr, SubSection):
                    m_def.m_add_sub_section(Section.sub_sections, attr)
                else:
                    raise NotImplementedError('Unknown property kind.')

            if inspect.isfunction(attr):
                method_name = attr.__name__

                # transfer constraints
                if getattr(attr, 'm_constraint', False):
                    constraint = method_name
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

        for content in m_def.m_all_contents(depth_first=True, include_self=True):
            cast(Definition, content).__init_metainfo__()

    def __setattr__(self, name, value):
        all_aliases = None
        if self.m_def is not None:
            all_aliases = self.m_def.all_aliases

        if all_aliases is not None and name in self.m_def.all_aliases:
            name = self.m_def.all_aliases[name].name

        return super().__setattr__(name, value)

    def __getattr__(self, name):
        # The existence of __getattr__ will make mypy and pylint ignore 'missing' dynamic
        # attributes and functions and wrong types of those.
        # Ideally we have a plugin for both that add the corrent type info

        if name in self.m_def.all_aliases:
            return getattr(self, self.m_def.all_aliases[name].name)

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

    def __to_np(self, quantity_def: 'Quantity', value):
        if isinstance(value, pint.quantity._Quantity):
            if quantity_def.unit is None:
                raise MetainfoError(
                    'The quantity %s has not a unit, but value %s has.' %
                    (quantity_def, value))

            if type(value.magnitude) == np.ndarray and quantity_def.type != value.dtype:
                value = value.astype(quantity_def.type)

            value = value.to(quantity_def.unit).magnitude

        if type(value) != np.ndarray:
            if len(quantity_def.shape) > 0:
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
        self.m_mod_count += 1

        if quantity_def.derived is not None:
            raise MetainfoError('The quantity %s is derived and cannot be set.' % quantity_def)

        if value is None:
            # This implements the implicit "unset" semantics of assigned None as a
            # value
            self.__dict__.pop(quantity_def.name, None)
            return

        if isinstance(quantity_def.type, np.dtype):
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
                    'quantities: %s' % quantity_def)

        self.__dict__[quantity_def.name] = value

        for handler in self.m_def.event_handlers:
            if handler.__name__.startswith('on_set'):
                handler(self, quantity_def, value)

    def m_get(self, quantity_def: 'Quantity') -> Any:
        ''' Retrieve the given value for the given quantity. '''
        return quantity_def.__get__(self, Quantity)

    def m_is_set(self, quantity_def: 'Quantity') -> bool:
        ''' True if the given quantity is set. '''
        if quantity_def.derived is not None:
            return True

        return quantity_def.name in self.__dict__

    def m_add_values(self, quantity_def: 'Quantity', values: Any, offset: int) -> None:
        ''' Add (partial) values for the given quantity of higher dimensionality. '''
        # TODO
        raise NotImplementedError()

    def m_add_sub_section(self, sub_section_def: 'SubSection', sub_section: 'MSection') -> None:
        ''' Adds the given section instance as a sub section of the given sub section definition. '''

        self.m_mod_count += 1
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

        if sub_section is not None:
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
        self.m_mod_count += 1

        if sub_section_def.repeats:
            sub_section = self.__dict__[sub_section_def.name][index]
            del(self.__dict__[sub_section_def.name][index])

        elif sub_section_def.name in self.__dict__:
            sub_section = self.__dict__[sub_section_def.name]
            del(self.__dict__[sub_section_def.name])

        if sub_section.m_resource is not None:
            sub_section.m_resource.remove(sub_section)

    def m_get_sub_section(self, sub_section_def: 'SubSection', index: int) -> 'MSection':
        ''' Retrieves a single sub section of the given sub section definition. '''
        if sub_section_def.repeats:
            return self.__dict__[sub_section_def.name][index]

        else:
            return self.__dict__.get(sub_section_def.name, None)

    def m_get_sub_sections(self, sub_section_def: 'SubSection') -> List['MSection']:
        ''' Retrieves  all sub sections of the given sub section definition. '''
        try:
            if sub_section_def.repeats:
                return self.__dict__[sub_section_def.name]
            else:
                return [self.__dict__[sub_section_def.name]]
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
        self.m_mod_count += 1
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
                        self.m_add_sub_section(prop, value)

                else:
                    self.m_set(prop, value)

        else:
            self.__dict__.update(**kwargs)

    def m_as(self, section_cls: Type[MSectionBound]) -> MSectionBound:
        ''' 'Casts' this section to the given extending sections. '''
        return cast(MSectionBound, self)

    def m_follows(self, definition: 'Section') -> bool:
        ''' Determines if this section's definition is or is derived from the given definition. '''
        if not isinstance(definition, Section):
            raise TypeError('%s is of class Section' % definition)
        return self.m_def == definition or definition in self.m_def.all_base_sections

    def m_to_dict(
            self, with_meta: bool = False,
            include_defaults: bool = False,
            include_derived: bool = False,
            categories: List[Union['Category', Type['MCategory']]] = None,
            partial: TypingCallable[['Definition', 'MSection'], bool] = None) -> Dict[str, Any]:
        '''
        Returns the data of this section as a json serializeable dictionary.

        Arguments:
            with_meta: Include information about the section definition and the sections
                position in its parent.
            include_defaults: Include default values of unset quantities.
            include_derived: Include values of derived quantities.
            categories: A list of category classes or category definitions that is used
                to filter the included quantities and sub sections. Only applied to
                properties of this section, not on sub-sections. Is overwritten
                by partial.
            partial: A function that determines if a definition should be included in
                the output dictionary. Takes a definition and the containing section
                as arguments. Two default functions can be used by providing a
                string instead:

                - 'mongo': Only include quantities that have an a_mongo
                  annotation.
                - 'es': Only include quantities that have an a_elastic or
                  an an a_search annotation.

                Partial is applied recursively on sub-sections. Overrides
                categories.
        '''
        # determine partial for sub-sections and partial based on categories
        if partial is not None:
            if partial == "es":
                partial = lambda d, s: hasattr(d, "a_search") or hasattr(d, "a_search")
            if partial == "mongo":
                partial = lambda d, s: hasattr(d, "a_mongo")
            child_partial = partial
        else:
            if categories is None:
                partial = lambda *args, **kwargs: True
                child_partial = lambda *args, **kwargs: True

            else:
                category_defs: List[Category] = None
                if categories is not None:
                    category_defs = []
                    for category in categories:
                        if issubclass(category, MCategory):  # type: ignore
                            category_defs.append(category.m_def)  # type: ignore
                        elif isinstance(category, Category):
                            category_defs.append(category)
                        else:
                            raise TypeError('%s is not a category' % category)

                partial = lambda definition, *args, **kwargs: any(
                    definition in category.get_all_definitions()
                    for category in category_defs)
                child_partial = lambda *args, **kwargs: True

        def serialize_quantity(quantity, is_set, is_derived):
            quantity_type = quantity.type

            serialize: TypingCallable[[Any], Any] = str
            if isinstance(quantity_type, Reference):

                def reference_serialize(value):
                    if isinstance(value, MProxy):
                        if value.m_proxy_resolved is not None:
                            return quantity_type.serialize(self, quantity, value)
                        else:
                            return value.m_proxy_value
                    else:
                        return quantity_type.serialize(self, quantity, value)
                serialize = reference_serialize

            elif isinstance(quantity_type, DataType):

                def data_type_serialize(value):
                    return quantity_type.serialize(self, quantity, value)

                serialize = data_type_serialize

            elif quantity_type in _primitive_types:
                serialize = _primitive_types[quantity_type]

            elif isinstance(quantity_type, np.dtype):
                is_scalar = quantity.is_scalar

                def serialize_dtype(value):
                    if isinstance(value, np.ndarray):
                        if is_scalar:
                            self.m_warning('numpy quantity has wrong shape', quantity=str(quantity))

                        return value.tolist()

                    else:
                        if not is_scalar:
                            self.m_warning('numpy quantity has wrong shape', quantity=str(quantity))

                        return value.item()

                serialize = serialize_dtype

            elif isinstance(quantity_type, MEnum):
                serialize = str

            elif quantity_type == Any:
                def _serialize(value: Any):
                    if type(value) not in [str, int, float, bool, np.bool_, list, dict, type(None)]:
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
                value = self.__dict__[quantity.name]
            elif is_derived:
                value = quantity.derived(self)
            else:
                value = quantity.default

            if isinstance(quantity_type, np.dtype):
                return serialize(value)
            elif len(quantity.shape) == 0:
                return serialize(value)
            elif len(quantity.shape) == 1:
                return [serialize(i) for i in value]
            else:
                raise NotImplementedError('Higher shapes (%s) not supported: %s' % (quantity.shape, quantity))

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
                if not partial(quantity, self):
                    continue

                try:
                    if quantity.virtual:
                        if include_derived and quantity.derived is not None:
                            yield name, serialize_quantity(quantity, False, True)
                        continue

                    is_set = self.m_is_set(quantity)
                    if not is_set:
                        if not include_defaults or not quantity.m_is_set(Quantity.default):
                            continue

                    yield name, serialize_quantity(quantity, is_set, False)

                except ValueError as e:
                    import traceback
                    traceback.print_exc()
                    raise ValueError('Value error (%s) for %s' % (str(e), quantity))

            # sub sections
            for name, sub_section_def in self.m_def.all_sub_sections.items():
                if not partial(sub_section_def, self):
                    continue

                if sub_section_def.repeats:
                    if self.m_sub_section_count(sub_section_def) > 0:
                        yield name, [
                            None if item is None else item.m_to_dict(
                                with_meta=with_meta,
                                include_defaults=include_defaults,
                                include_derived=include_derived,
                                partial=child_partial)
                            for item in self.m_get_sub_sections(sub_section_def)]
                else:
                    sub_section = self.m_get_sub_section(sub_section_def, -1)
                    if sub_section is not None:
                        yield name, sub_section.m_to_dict(
                            with_meta=with_meta,
                            include_defaults=include_defaults,
                            include_derived=include_derived,
                            partial=child_partial)

        return {key: value for key, value in items()}

    def m_update_from_dict(self, dct: Dict[str, Any]) -> None:
        '''
        Updates this section with the serialized data from the given dict, e.g. data
        produced by :func:`m_to_dict`.
        '''
        section_def = self.m_def
        section = self

        for name, sub_section_def in section_def.all_sub_sections.items():
            if name in dct:
                sub_section_value = dct.get(name)
                if sub_section_def.repeats:
                    for sub_section_dct in sub_section_value:
                        if sub_section_dct is None:
                            sub_section = None
                        else:
                            sub_section = sub_section_def.sub_section.section_cls.m_from_dict(sub_section_dct)
                        section.m_add_sub_section(sub_section_def, sub_section)

                else:
                    sub_section = sub_section_def.sub_section.section_cls.m_from_dict(sub_section_value)
                    section.m_add_sub_section(sub_section_def, sub_section)

        for name, quantity_def in section_def.all_quantities.items():
            if name in dct:
                quantity_value = dct[name]

                if isinstance(quantity_def.type, np.dtype):
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

    @classmethod
    def m_from_dict(cls: Type[MSectionBound], dct: Dict[str, Any]) -> MSectionBound:
        ''' Creates a section from the given serializable data dictionary.

        This is the 'opposite' of :func:`m_to_dict`. It takes a deserialized dict, e.g
        loaded from JSON, and turns it into a proper section, i.e. instance of the given
        section class.
        '''
        section = cls()
        section.m_update_from_dict(dct)
        return section

    def m_to_json(self, **kwargs):
        ''' Returns the data of this section as a json string. '''
        return json.dumps(self.m_to_dict(), **kwargs)

    def m_all_contents(
            self, depth_first: bool = False, include_self: bool = False,
            stop: TypingCallable[['MSection'], bool] = None) -> Iterable['MSection']:
        '''
        Returns an iterable over all sub and sub subs sections.

        Arguments:
            depth_first: A boolean indicating that children should be returned before
                parents.
            include_self: A boolean indicating that the results should contain this section.
            stop: A predicate that determines if the traversal should be stopped or if
                children should be returned. The sections for which this returns True
                are included in the results.
        '''
        if include_self and not depth_first:
            yield self

        if stop is None or not stop(self):
            for content in self.m_contents():
                if not depth_first:
                    yield content

                for sub_content in content.m_all_contents(depth_first=depth_first, stop=stop):
                    yield sub_content

                if depth_first:
                    yield content

        if include_self and depth_first:
            yield self

    def m_traverse(self):
        '''
        Performs a depth-first traversal and yield tuples of section, property def,
        parent index for all set properties.
        '''
        for key in self.__dict__:
            property_def = self.m_def.all_properties.get(key)
            if property_def is None:
                continue

            if isinstance(property_def, SubSection):
                for sub_section in self.m_get_sub_sections(property_def):
                    for i in sub_section.m_traverse():
                        yield i

                    yield self, property_def, sub_section.m_parent_index

            else:
                yield self, property_def, -1

    def m_pretty_print(self, indent=None):
        ''' Pretty prints the containment hierarchy '''
        if indent is None:
            indent = []
        if len(indent) == 0:
            indent_str = ''
        else:
            indent_str = ''.join(['   ' if last else ' │ ' for last in indent[:-1]])
            indent_str += ' └╴' if indent[-1] else ' ├╴'
        print(indent_str + str(self))

        contents = list(self.m_contents())
        length = len(contents)
        for i, content in enumerate(contents):
            content.m_pretty_print(indent + [i == length - 1])

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
                if sub_section is not None:
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

    def m_resolved(self):
        '''
        Returns the original resolved object, if this instance used to be a proxy.

        For most purposes a resolved proxy is equal to the section it was resolved to.
        The exception are hashes. So if you want to use a potential former proxy in
        a hash table and make it really equal to the section it was resolved to, use
        the result of this method instead of the section/proxy itself.
        '''
        return getattr(self, 'm_proxy_resolved', self)

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
                    if len(path_stack) == 0:
                        return context.m_get_sub_sections(prop_def)  # type: ignore

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
                        'Could not resolve %s, %s is not a property' % (path, prop_name))

                if not context.m_is_set(prop_def):
                    raise ReferenceError(
                        'Could not resolve %s, %s is not set on %s' % (path, prop_name, context))

                return context.m_get(prop_def)

        return cast(MSectionBound, context)

    def m_get_annotations(self, key: Union[str, type], default=None, as_list: bool = False):
        '''
        Convenience method to get annotations

        Arguments:
            key: Either the optional annotation name or an annotation class. In the first
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

    def m_validate(self) -> Tuple[List[str], List[str]]:
        ''' Evaluates all constraints and shapes of this section and returns a list of errors. '''
        errors: List[str] = []
        warnings: List[str] = []
        for constraint_name in self.m_def.constraints:
            constraint = getattr(self, constraint_name, None)
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
                if getattr(constraint, 'm_warning', False):
                    warnings.append(error_str)
                else:
                    errors.append(error_str)

        for quantity in self.m_def.all_quantities.values():
            if self.m_is_set(quantity) and not quantity.derived:
                if not self.__validate_shape(quantity, self.m_get(quantity)):
                    errors.append(
                        'The shape of quantity %s does not match its value.' % quantity)

        return errors, warnings

    def m_copy(self: MSectionBound, deep=False, parent=None) -> MSectionBound:
        # TODO this a shallow copy, but should be a deep copy
        copy = self.m_def.section_cls()
        copy.__dict__.update(**self.__dict__)
        copy.m_parent = parent
        copy.m_parent_index = -1 if parent is None else self.m_parent_index
        copy.m_parent_sub_section = None if parent is None else self.m_parent_sub_section
        copy.m_resource = None

        if deep:
            for sub_section_def in self.m_def.all_sub_sections.values():
                sub_sections_copy = [
                    sub_section.m_copy(deep=True, parent=copy)
                    for sub_section in self.m_get_sub_sections(sub_section_def)]

                if sub_section_def.repeats:
                    copy.__dict__[sub_section_def.name] = sub_sections_copy
                else:
                    if len(sub_sections_copy) == 1:
                        copy.__dict__[sub_section_def.name] = sub_sections_copy[0]
                    else:
                        copy.__dict__[sub_section_def.name] = None

        return cast(MSectionBound, copy)

    def m_all_validate(self):
        ''' Evaluates all constraints in the whole section hierarchy, incl. this section. '''
        errors: List[str] = []
        warnings: List[str] = []
        for section in itertools.chain([self], self.m_all_contents()):
            more_errors, more_warnings = section.m_validate()
            errors.extend(more_errors)
            warnings.extend(more_warnings)

        return errors, warnings

    def m_warning(self, *args, **kwargs):
        if self.m_resource is not None:
            self.m_resource.warning(*args, **kwargs)

    def __repr__(self):
        m_section_name = self.m_def.name
        # name_quantity_def = self.m_def.all_quantities.get('name', None)
        # if name_quantity_def is not None:
        #     name = self.m_get(name_quantity_def)
        try:
            name = self.__dict__['name']
            main = '%s:%s' % (name, m_section_name)
        except KeyError:
            main = m_section_name

        more = ''
        props = [
            prop
            for prop in self.m_def.all_properties
            if prop in self.__dict__]

        if len(props) > 10:
            more = ', +%d more properties' % (len(props) - 10)

        return '%s(%s%s)' % (main, ', '.join(props[0:10]), more)

    def __getitem__(self, key):
        try:
            key = key.replace('.', '/')
            return self.m_resolve(key)
        except ReferenceError:
            raise KeyError(key)

    def __iter__(self):
        return self.m_def.all_properties.__iter__()

    def __len__(self):
        return len(self.m_def.all_properties)

    def get(self, key):
        return self.__dict__.get(key, None)

    def values(self):
        return {key: val for key, val in self.__dict__.items() if not key.startswith('m_')}.values()

    def m_xpath(self, expression: str):
        '''
        Provides an interface to jmespath search functionality.

        Arguments:
            expression: A string compatible with the jmespath specs representing the
                search. See https://jmespath.org/ for complete description.

        .. code-block:: python

            metainfo_section.m_xpath('code_name')
            metainfo_section.m_xpath('systems[-1].system_type')
            metainfo_section.m_xpath('sccs[0].system.atom_labels')
            metainfo_section.m_xpath('systems[?system_type == `molecule`].atom_labels')
            metainfo_section.m_xpath('sccs[?energy_total < `1.0E-23`].system')

        '''
        def to_dict(entries):
            if not isinstance(entries, list):
                try:
                    entries = entries.m_to_dict()
                except Exception:
                    pass
                return entries
            else:
                return [to_dict(entry) for entry in entries]

        result = jmespath.search(expression, self)
        return to_dict(result)


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
    '''
    :class:`Definition` is the common base class for all metainfo definitions.

    All metainfo `definitions` (sections, quantities, sub-sections, packages, ...) share
    some common properties.

    Attributes:
        name: Each `definition` has a name. Names have to be valid Python identifier.
            They can contain letters, numbers and _, but must not start with a number.
            This also qualifies them as identifier in most storage formats, databases,
            makes them URL safe, etc.

            Names must be unique within the :class:`Package` or :class:`Section` that
            this definition is part of.

            By convention, we use capitalized `CamelCase` identifier to refer to *sections
            definitions* (i.e. section definitions are represented by Python classes),
            lower case `snake_case` identifier for variables that hold *sections*, and for
            *properties* (i.e. fields in a Python class) we typically use lower
            case `snake_case` identifier. Sub-sections are often prefixed with ``section_``
            to clearly separate sub-sections from quantities.

            Generally, you do not have to set this attribute manually, it will be derived
            from Python identifiers automatically.

        description: The description can be an arbitrary human readable text that explains
            what a definition is about. For section definitions you do not have to set
            this manually as it will be derived from the classes doc string. Quantity and
            sub-section descriptions can also be taken from the containing section class'
            doc-string ``Attributes:`` section.

        links: Each definition can be accompanied by a list of URLs. These should point
            to resources that further explain the definition.

        aliases: A list of alternative names. For quantities and sub-sections these
            can be used to access the respective property with a different name from
            its containing section.

        deprecated: If set this definition is marked deprecated. The value should be a
            string that describes how to replace the deprecated definition.

        categories: All metainfo definitions can be put into one or more `categories`.
            Categories allow to organize the definitions themselves. It is different from
            sections, which organize the data (e.g. quantity values) and not the definitions
            of data (e.g. quantities definitions). See :ref:`metainfo-categories` for more
            details.
    '''

    name: 'Quantity' = _placeholder_quantity
    description: 'Quantity' = _placeholder_quantity
    links: 'Quantity' = _placeholder_quantity
    categories: 'Quantity' = _placeholder_quantity
    deprecated: 'Quantity' = _placeholder_quantity
    aliases: 'Quantity' = _placeholder_quantity

    def __init_metainfo__(self):
        '''
        An initialization method that is called after the class context of the definition
        has been initialized. For example it is called on all quantities of a section
        class after the class was created. If metainfo definitions are created without
        a class context, this method must be called manually on all definitions.
        '''

        # initialize definition annotations
        for annotation in self.m_get_annotations(DefinitionAnnotation, as_list=True):
            annotation.init_annotation(self)

    def init_metainfo(self):
        '''
        Calls __init_metainfo__ on all its children. This is necessary if the
        package, section was created without corresponding python classes
        packages, etc.
        '''
        for content in self.m_all_contents(depth_first=True):
            content.__init_metainfo__()

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

    def __repr__(self):
        return '%s:%s' % (self.qualified_name(), self.m_def.name)


class Property(Definition):
    ''' A common base-class for section properties: sub sections and quantities. '''
    pass


class Quantity(Property):
    '''
    To define quantities, instantiate :class:`Quantity` as a classattribute values in
    a `section classes`. The name of a quantity is automatically taken from its `section class`
    attribute. You can provide all other attributes to the constructor with keyword arguments

    See :ref:`metainfo-sections` to learn about `section classes`.
    In Python terms, ``Quantity`` is a descriptor. Descriptors define how to get and
    set attributes in a Python object. This allows us to use sections like regular
    Python objects and quantity like regular Python attributes.

    Each quantity must define a basic data type and a shape. The values of a quantity must
    fulfil the given type. The default shape is a single value. Quantities can also have
    physical units. Units are applied to all values.

    Attributes:
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

        is_scalar:
            Derived quantity that is True, iff this quantity has shape of length 0

        unit:
            The physics unit for this quantity. It is optional.

            Units are represented with the Pint Python package. Pint defines units and
            their algebra. You can either use *pint* units directly, e.g. ``units.m / units.s``.
            The metainfo provides a preconfigured *pint* unit registry :py:data:`ureg`.
            You can also provide the unit as *pint* parsable string, e.g. ``'meter / seconds'`` or
            ``'m/s'``.

        default:
            The default value for this quantity. The value must match type and shape.

            Be careful with a default value like ``[]`` as it will be the default value for
            all occurrences of this quantity.

    Quantities are mapped to Python properties on all section objects that instantiate
    the Python class/section definition that has this quantity. This means quantity values
    can be read and set like normal Python attributes.

    In some cases it might be desirable to have virtual and read only quantities that are
    not real quantities used for storing values, but rather define an interface to other
    quantities. Examples for this are `synonyms` and `derived` quantities.

    Attributes:
        derived:
            A Python callable that takes the containing section as input and outputs the
            value for this quantity. This quantity cannot be set directly, its value
            is only derived by the given callable. The callable is executed when this
            quantity is get. Derived quantities are always virtual.

        cached:
            A bool indicating that derived values should be cached unless the underlying
            section has changed.

        virtual:
            A boolean that determines if this quantity is virtual. Virtual quantities can
            be get/set like regular quantities, but their values are not (de-)serialized,
            hence never permanently stored.
    '''

    type: 'Quantity' = _placeholder_quantity
    shape: 'Quantity' = _placeholder_quantity
    unit: 'Quantity' = _placeholder_quantity
    default: 'Quantity' = _placeholder_quantity
    derived: 'Quantity' = _placeholder_quantity
    cached: 'Quantity' = _placeholder_quantity
    virtual: 'Quantity' = _placeholder_quantity
    is_scalar: 'Quantity' = _placeholder_quantity

    # TODO derived_from = Quantity(type=Quantity, shape=['0..*'])
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __init_metainfo__(self):
        super().__init_metainfo__()

        if self.derived is not None:
            self.virtual = True

        # replace the quantity implementation with an optimized version for the most
        # primitive quantities if applicable
        is_primitive = not self.derived
        is_primitive = is_primitive and len(self.shape) <= 1
        is_primitive = is_primitive and self.type in [str, bool, float, int]
        is_primitive = is_primitive and not isinstance(self.type, np.dtype)
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
            if self.derived is not None:
                try:
                    if self.cached:
                        cached = obj.__dict__.setdefault(self.name + '_cached', [-1, None])
                        if cached[0] != obj.m_mod_count:
                            cached[0] = obj.m_mod_count
                            cached[1] = self.derived(obj)  # pylint: disable=not-callable
                        return cached[1]
                    else:
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

        elif isinstance(self.type, np.dtype):
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

    @constraint(warning=True)
    def dimensions(self):
        for dimension in self.shape:
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

    @constraint
    def higher_shapes_require_dtype(self):
        if len(self.shape) > 1:
            assert isinstance(self.type, np.dtype), \
                'Higher dimensional quantities (%s) need a dtype and will be treated as ' \
                'numpy arrays.' % self


def derived(**kwargs):
    def decorator(f) -> Quantity:
        if 'name' not in kwargs:
            kwargs['name'] = f.__name__
        if 'description' not in kwargs:
            kwargs['description'] = f.__doc__.strip() if f.__doc__ is not None else None
        if 'type' not in kwargs:
            kwargs['type'] = Any

        return Quantity(derived=f, **kwargs)
    return decorator


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
        obj.m_mod_count += 1
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
        obj.m_mod_count += 1

        if value is None:
            obj.__dict__.pop(self.name, None)
            return

        if self._list:
            if not isinstance(value, list):
                if hasattr(value, 'tolist'):
                    value = value.tolist()

                else:
                    raise TypeError(
                        'The value %s for quantity %s has not shape %s' %
                        (value, self, self.shape))

            if any(v is not None and type(v) != self._type for v in value):
                raise TypeError(
                    'The value %s with type %s for quantity %s is not of type %s' %
                    (value, type(value), self, self.type))

        elif type(value) != self._type:
            raise TypeError(
                'The value %s with type %s for quantity %s is not of type %s' %
                (value, type(value), self, self.type))

        try:
            obj.__dict__[self._name] = value
        except AttributeError:
            raise KeyError('Cannot overwrite quantity definition. Only values can be set.')


class SubSection(Property):
    '''
    Like quantities, sub-sections are defined in a `section class` as attributes
    of this class. An like quantities, each sub-section definition becomes a property of
    the corresponding `section definition` (parent). A sub-section definition references
    another `section definition` as the sub-section (child). As a consequence, parent
    `section instances` can contain child `section instances` as sub-sections.

    Contrary to the old NOMAD metainfo, we distinguish between sub-section the section
    and sub-section the property. This allows to use on child `section definition` as
    sub-section of many different parent `section definitions`.

    Attributes:
        sub_section: A :class:`Section` or Python class object for a `section class`. This
            will be the child `section definition`. The defining section the child
            `section definition`.

        repeats: A boolean that determines wether this sub-section can appear multiple
            times in the parent section.
    '''

    sub_section: 'Quantity' = _placeholder_quantity
    repeats: 'Quantity' = _placeholder_quantity

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
    '''
    Instances of the class :class:`Section` are created by writing Python classes
    that extend :class:`MSection` like this:

    .. code-block:: python

        class SectionName(BaseSection):
            \'\'\' Section description \'\'\'
            m_def = Section(**section_attributes)

            quantity_name = Quantity(**quantity_attributes)
            sub_section_name = SubSection(**sub_section_attributes)

    We call such classes *section classes*. They are not the *section definition*, but just
    representation of it in Python syntax. The *section definition* (in instance of :class:`Section`)
    will be created for each of these classes and stored in the ``m_def`` property. See
    :ref:`metainfo-reflection` for more details.

    Most of the attributes for a :class:`Section` instance will be set automatically from
    the section class:

    Attributes:
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

            If there are no base sections to define, you have to use :class:`MSection`.

    The Metainfo supports two inheritance mechanism. By default it behaves like regular
    Python inheritance and the class inherits all its base
    classes' properties. The other mode (enabled via ``extends_base_section=True``), will
    add all sub-class properties to the base-class. This is used throughout the NOMAD metainfo
    to add code-specific metadata to common section definitions. Here is an example:

    .. code-block:: python

        class Method(MSection):
            code_name = Quantity(str)

        class VASPMethod(Method):
            m_def = Section(extends_base_section=True)
            x_vasp_some_incar_parameter = Quantity(str)

        method = Method()
        methid.x_vasp_same_incar_parameter = 'value'

    In this example, the section class ``VASPMethod`` defines a section definition that inherits
    from section definition ``Method``. The quantity `x_vasp_some_incar_parameter` will
    be added to `Method` and can be used in regular `Method` instances.

    The following :class:`Section` attributes maniputlate the inheritance semantics:

    Attributes:
        extends_base_section:
            If True, this definition must have exactly one ``base_sections``.
            Instead of inheriting properties, the quantity and sub-section definitions
            of this section will be added to the base section.

            This allows to add further properties to an existing section definition.
            To use such extension on section instances in a type-safe manner
            :py:func:`MSection.m_as` can be used to cast the base section to the extending
            section.

        extending_sections:
            A list of `section definitions` (:class:`Section`). These are those sections
            that add their properties to this section via :attr:`extends_base_section`.
            This quantity will be set automatically.


    Besides defining quantities and sub-sections, a section definition can also provide
    constraints that are used to validate a section and its quantities and sub-sections.
    Constraints allow to define more specific data structures beyond types and shapes.
    But constraints are not enforced automatically, sections have to be explicitly
    validated in order to evaluate constraints.

    Constrains can be defined via methods with the ``constraint`` decorator:

    .. code-block:: python

        class System(MSection):
            lattice = Quantity(float, shape=[3, 3], unit='meter')

            @constraint
            def non_empty_lattice(self):
                assert np.abs(np.linalg.det(self.lattice.magnitude)) > 0

        system = System()
        system.m_validate()

    Attributes:
        constraints:
            Constraints are rules that a section must fulfil to be valid. This allows to implement
            semantic checks that go behind mere type or shape checks. This quantity takes
            the names of constraints as string. Constraints have to be implemented as methods
            with the :func:`constraint` decorator. They can raise :class:`ConstraintVialated`
            or an AssertionError to indicate that the constraint is not fulfilled for the ``self``
            section. This quantity will be set automatically from all constraint methods in the
            respective section class. To run validation of a section use :py:meth:`MSection.m_validate`.

    Other attributes and helper properties:

    Attributes:
        section_cls:
            A helper attribute that gives the `section class` as a Python class object.

        inherited_sections:
            A helper attribute that gives direct and indirect base sections and extending
            sections including this section. These are all sections that this sections
            gets its properties from.

        all_base_sections:
            A helper attribute that gives direct and indirect base sections.

        all_properties:
            A helper attribute that gives all properties (sub section and quantity) definitions
            including inherited properties and properties from extending sections as a
            dictionary with names and definitions.

        all_quantities:
            A helper attribute that gives all quantity definition including inherited ones
            and ones from extending sections as a dictionary that maps names (strings)
            to :class:`Quantity`.

        all_sub_sections:
            A helper attribute that gives all sub-section definition including inherited ones
            and ones from extending sections as a dictionary that maps names (strings)
            to :class:`SubSection`.

        all_sub_sections_by_section:
            A helper attribute that gives all sub-section definition including inherited ones
            and ones from extending sections as a dictionary that maps section classes
            (i.e. Python class objects) to lists of :class:`SubSection`.

        all_aliases:
            A helper attribute that gives all aliases for all properties including
            inherited properties and properties form extending sections as a
            dictionary with aliases and the definitions.

        event_handlers:
            Event handler are functions that get called when the section data is changed.
            There are two types of events: ``set`` and ``add_sub_section``. The handler type
            is determined by the handler (i.e. function) name: ``on_set`` and ``on_add_sub_section``.
            The handler arguments correspond to :py:meth:`MSection.m_set` (section, quantity_def, value) and
            :py:meth:`MSection.m_add_sub_section` (section, sub_section_def, sub_section).
            Handler are called after the respective action was performed. This quantity is
            automatically populated with handler from the section classes methods. If there
            is a method ``on_set`` or ``on_add_sub_section``, it will be added as handler.

        errors:
            A list of errors. These issues prevent the section definition from being usable.

        warnings:
            A list of warnings. These still allow to use the section definition.
    '''

    quantities: 'SubSection' = None
    sub_sections: 'SubSection' = None

    base_sections: 'Quantity' = _placeholder_quantity
    extending_sections: 'Quantity' = _placeholder_quantity
    extends_base_section: 'Quantity' = _placeholder_quantity
    constraints: 'Quantity' = _placeholder_quantity
    event_handlers: 'Quantity' = _placeholder_quantity

    inherited_sections: 'Quantity' = _placeholder_quantity
    all_base_sections: 'Quantity' = _placeholder_quantity
    all_properties: 'Quantity' = _placeholder_quantity
    all_quantities: 'Quantity' = _placeholder_quantity
    all_sub_sections: 'Quantity' = _placeholder_quantity
    all_sub_sections_by_section: 'Quantity' = _placeholder_quantity
    all_aliases: 'Quantity' = _placeholder_quantity

    def __init__(self, *args, validate: bool = True, **kwargs):
        self._section_cls: Type[MSection] = None

        super().__init__(*args, **kwargs)
        self.validate = validate
        self.errors: List[str] = []
        self.warnings: List[str] = []

    @property
    def section_cls(self) -> Type[MSection]:
        if self._section_cls is None:
            # Create a section class if this does not exist. This happens if the section
            # is not created through a class definition.
            attrs = {
                prop.name: prop
                for prop in self.quantities + self.sub_sections}
            attrs.update(m_def=self, do_init=False)
            self._section_cls = type(self.name, (MSection,), attrs)

        return self._section_cls

    def __init_metainfo__(self):
        super().__init_metainfo__()

        if self.extends_base_section:
            base_sections_count = len(self.base_sections)
            if base_sections_count == 0:
                raise MetainfoError(
                    'Section %s extend the base section, but has no base section.' % self)

            if base_sections_count > 1:
                raise MetainfoError(
                    'Section %s extend the base section, but has more than one base section' % self)

            base_section = self.base_sections[0]
            for name, attr in self.section_cls.__dict__.items():
                if isinstance(attr, Property):
                    setattr(base_section.section_cls, name, attr)

            base_section.extending_sections = base_section.extending_sections + [self]

        # validate
        def validate(definition):
            errors, warnings = definition.m_all_validate()
            self.errors.extend(errors)
            self.warnings.extend(warnings)

        if not self.validate:
            return

        if self.extends_base_section:
            validate(self.m_get(Section.base_sections)[0])
        else:
            validate(self)

        if len(self.errors) > 0:
            raise MetainfoError(
                '%s. The section definition %s violates %d more constraints' %
                (str(self.errors[0]).strip('.'), self, len(self.errors) - 1))

    @constraint
    def unique_names(self):
        # start with the names of all base_sections
        names: Set[str] = set()
        for base in list(self.all_base_sections) + self.extending_sections:
            for quantity in base.quantities + base.sub_sections:
                for alias in quantity.aliases:
                    names.add(alias)
                names.add(quantity.name)

        for def_list in [self.quantities, self.sub_sections]:
            for definition in def_list:
                assert definition.name not in names, 'All names in a section must be unique. ' \
                    'Name %s of %s in %s already exists in %s.' % (definition.name, definition, definition.m_parent, self)
                names.add(definition.name)
                for alias in definition.aliases:
                    assert alias not in names, 'All names (incl. aliases) in a section must be unique. ' \
                        'Alias %s of %s in %s already exists in %s.' % (alias, definition, definition.m_parent, self)
                    names.add(alias)


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

    Attributes:
        section_definitions: All `section definitions` in this package as :class:`Section`
            objects.

        category_definitions: All `category definitions` in this package as :class:`Category`
            objects.

        all_definitions: A helper attribute that provides all section definitions
            by name.
    '''

    section_definitions: 'SubSection' = None
    category_definitions: 'SubSection' = None

    all_definitions: 'Quantity' = _placeholder_quantity
    dependencies: 'Quantity' = _placeholder_quantity

    registry: Dict[str, 'Package'] = {}
    ''' A static member that holds all currently known packages. '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __init_metainfo__(self):
        super().__init_metainfo__()

        Package.registry[self.name] = self

        # access potential SectionProxies to resolve them
        for content in self.m_all_contents():
            if isinstance(content, Quantity):
                if isinstance(content.type, MProxy):
                    content.type.m_proxy_resolve()
            elif isinstance(content, SubSection):
                if isinstance(content.sub_section, MProxy):
                    content.sub_section.m_proxy_resolve()

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
    ''' Categories allow to organize metainfo definitions (not metainfo data like sections do)

    Each definition, including categories themselves, can belong to a set of categories.
    Categories therefore form a hierarchy of concepts that definitions can belong to, i.e.
    they form a `is a` relationship.

    Attributes:
        definitions: A helper attribute that gives all definitions that are directly or
            indirectly in this category.
    '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.definitions: Set[Definition] = set()

    def get_all_definitions(self, definitions: Set[Definition] = None) -> Set[Definition]:
        '''
        Helper method that collects all non category definitions, including those
        in categories of this category.
        '''
        if definitions is None:
            definitions = set()

        for definition in self.definitions:
            if isinstance(definition, MCategory):
                definition = definition.m_def

            if isinstance(definition, Category):
                definition.get_all_definitions(definitions)

            if definition not in definitions:
                definitions.add(definition)

        return definitions


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
Section.m_def._section_cls = Section

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
Definition.deprecated = Quantity(type=str, name='deprecated')
Definition.aliases = Quantity(type=str, shape=['0..*'], default=[], name='aliases')

Section.quantities = SubSection(
    sub_section=Quantity.m_def, name='quantities', repeats=True)

Section.sub_sections = SubSection(
    sub_section=SubSection.m_def, name='sub_sections', repeats=True)
Section.base_sections = Quantity(
    type=Reference(Section.m_def), shape=['0..*'], default=[], name='base_sections')
Section.extending_sections = Quantity(
    type=Reference(Section.m_def), shape=['0..*'], default=[], name='extending_sections')
Section.extends_base_section = Quantity(type=bool, default=False, name='extends_base_section')
Section.constraints = Quantity(type=str, shape=['0..*'], default=[], name='constraints')
Section.event_handlers = Quantity(
    type=Callable, shape=['0..*'], name='event_handlers', virtual=True, default=[])


@derived(cached=True)
def inherited_sections(self) -> Set[Section]:
    result: Set[Section] = set()
    result.add(self)
    for base_section in self.base_sections:
        result.add(base_section)
        for base_base_section in base_section.all_base_sections:
            result.add(base_base_section)
    for extending_section in self.extending_sections:
        result.add(extending_section)
    return result


@derived(cached=True)
def all_base_sections(self) -> Set[Section]:
    result: Set[Section] = set()
    for base_section in self.base_sections:
        result.add(base_section)
        for base_base_section in base_section.all_base_sections:
            result.add(base_base_section)
    return result


@derived(cached=True)
def all_properties(self) -> Dict[str, Union[SubSection, Quantity]]:
    result: Dict[str, Union[SubSection, Quantity]] = dict()
    for section in self.inherited_sections:
        for definition in section.quantities + section.sub_sections:
            result[definition.name] = definition
    return result


@derived(cached=True)
def all_quantities(self) -> Dict[str, Quantity]:
    result: Dict[str, Quantity] = dict()
    for section in self.inherited_sections:
        for definition in section.quantities:
            result[definition.name] = definition
    return result


@derived(cached=True)
def all_sub_sections(self) -> Dict[str, SubSection]:
    result: Dict[str, SubSection] = dict()
    for section in self.inherited_sections:
        for definition in section.sub_sections:
            result[definition.name] = definition
    return result


@derived(cached=True)
def all_sub_sections_by_section(self) -> Dict[Section, List[SubSection]]:
    result: Dict[Section, List[SubSection]] = dict()
    for section in self.inherited_sections:
        for definition in section.sub_sections:
            sub_sections = result.setdefault(definition.sub_section.m_resolved(), [])
            sub_sections.append(definition)
    return result


@derived(cached=True)
def all_aliases(self) -> Dict[str, Union[SubSection, Quantity]]:
    result: Dict[str, Union[SubSection, Quantity]] = dict()
    for section in self.inherited_sections:
        for definition in section.quantities + section.sub_sections:
            for alias in definition.aliases:
                result[alias] = definition
            result[definition.name] = definition

    return result


Section.inherited_sections = inherited_sections
Section.all_base_sections = all_base_sections
Section.all_properties = all_properties
Section.all_quantities = all_quantities
Section.all_sub_sections = all_sub_sections
Section.all_sub_sections_by_section = all_sub_sections_by_section
Section.all_aliases = all_aliases


SubSection.repeats = Quantity(type=bool, name='repeats', default=False)

SubSection.sub_section = Quantity(type=Reference(Section.m_def), name='sub_section')

Quantity.m_def._section_cls = Quantity
Quantity.type = DirectQuantity(type=QuantityType, name='type')
Quantity.shape = DirectQuantity(type=Dimension, shape=['0..*'], name='shape', default=[])
Quantity.unit = Quantity(type=Unit, name='unit')
Quantity.default = DirectQuantity(type=Any, default=None, name='default')
Quantity.derived = DirectQuantity(type=Callable, default=None, name='derived', virtual=True)
Quantity.virtual = DirectQuantity(type=bool, default=False, name='virtual')
Quantity.is_scalar = Quantity(
    type=bool, name='is_scalar', derived=lambda quantity: len(quantity.shape) == 0)
Quantity.cached = Quantity(type=bool, name='cached', default=False)

Package.section_definitions = SubSection(
    sub_section=Section.m_def, name='section_definitions', repeats=True)

Package.category_definitions = SubSection(
    sub_section=Category.m_def, name='category_definitions', repeats=True)


@derived(cached=True)
def all_definitions(self):
    all_definitions: Dict[str, Definition] = dict()
    for sub_section_def in [Package.section_definitions, Package.category_definitions]:
        for definition in self.m_get_sub_sections(sub_section_def):
            all_definitions[definition.name] = definition
            for alias in definition.aliases:
                all_definitions[alias] = definition
    return all_definitions


@derived(cached=True)
def dependencies(self):
    '''
    All packages which have definitions that definitions from this package need. Being
    'needed' includes categories, base sections, and referenced definitions.
    '''
    dependencies: Set[Package] = set()
    for content in self.m_all_contents():
        to_add = None
        if isinstance(content, Definition):
            for category in content.categories:
                if category.m_parent != self:
                    to_add = category.m_parent

        if isinstance(content, Quantity):
            if isinstance(content.type, Reference):
                if content.type.target_section_def.m_parent != self:
                    to_add = content.type.target_section_def.m_parent

        if isinstance(content, Section):
            for section in content.base_sections:
                if section.m_parent != self:
                    to_add = section.m_parent

        more_dependencies = []
        if to_add is not None:
            more_dependencies.append(to_add)
        while len(more_dependencies) > 0:
            dependency = more_dependencies.pop()
            if dependency not in dependencies:
                dependencies.add(dependency)
                more_dependencies.extend(dependency.dependencies)

    return dependencies


Package.all_definitions = all_definitions
Package.dependencies = dependencies

is_bootstrapping = False

Definition.__init_cls__()
Property.__init_cls__()
Section.__init_cls__()
Package.__init_cls__()
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

    @derived(cached=True)
    def all_definitions_by_name(self):
        all_definitions_by_name: Dict[str, List[Definition]] = dict()
        for definition in self.m_all_contents():
            if isinstance(definition, Definition):
                for name in [definition.name] + definition.aliases:
                    definitions = all_definitions_by_name.setdefault(name, [])
                    assert definition not in definitions, '%s must be unique' % definitions
                    definitions.append(definition)

        return all_definitions_by_name

    def resolve_definitions(
            self, name: str, section_cls: Type[MSectionBound],
            filter: TypingCallable[[MSection], bool] = None) -> List[MSectionBound]:

        return [
            definition
            for definition in self.all_definitions_by_name.get(name, [])
            if isinstance(definition, section_cls)
            if not (isinstance(definition, Section) and definition.extends_base_section)
            if filter is None or filter(definition)]

    def resolve_definition(
            self, name, section_cls: Type[MSectionBound],
            filter: TypingCallable[[MSection], bool] = None) -> MSectionBound:

        defs = self.resolve_definitions(name, section_cls, filter=filter)
        if len(defs) == 1:
            return defs[0]
        elif len(defs) > 1:
            raise KeyError('Could not uniquely identify %s, candidates are %s' % (name, defs))
        else:
            raise KeyError('Could not resolve %s' % name)
