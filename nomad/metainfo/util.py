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

import email.utils
import hashlib
import os
import re
from dataclasses import dataclass
from datetime import date, datetime
from difflib import SequenceMatcher
from functools import reduce
from typing import Any, Dict, Optional, Sequence, Tuple, Union
from urllib.parse import SplitResult, urlsplit, urlunsplit

import aniso8601
import numpy as np
import pandas as pd
import pint
import pytz

from nomad.units import ureg

__hash_method = 'sha1'  # choose from hashlib.algorithms_guaranteed
_delta_symbols = {'delta_', 'Δ'}


@dataclass(frozen=True)
class MRegEx:
    # matches the range of indices, e.g., 1..3, 0..*
    index_range = re.compile(r'(\d)\.\.(\d|\*)')
    # matches the reserved name
    reserved_name = re.compile(r'^(m_|a_|_+).*$')
    # matches for example
    # Python package/module name: nomad.metainfo.section
    # Python name + 40 digits id: nomad.metainfo.section@1a2b3c...
    python_definition = re.compile(r'^\w*(\.\w*)*(@\w{40})?$')
    # matches url
    url = re.compile(
        r'^(?:http|ftp)s?://'
        r'(?:(?:[A-Z0-9](?:[A-Z0-9-]{0,61}[A-Z0-9])?\.)+(?:[A-Z]{2,6}\.?|[A-Z0-9-]{2,}\.?)|'
        r'localhost|'
        r'\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})'
        r'(?::\d+)?'
        r'(?:/?|[/?]\S+)$', re.IGNORECASE)
    complex_str = re.compile(
        r'^(?=[iIjJ.\d+-])([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?(?![iIjJ.\d]))?'
        r'([+-]?(?:(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)?[iIjJ])?$')


def normalize_complex(value, complex_type, to_unit: Union[str, ureg.Unit, None]):
    '''
    Try to convert a given value to a complex number.
    '''

    def __check_precision(_type):
        if isinstance(_type, type(None)):
            return

        precision_error = ValueError(
            f'Cannot type {_type.__name__} to complex number of type {complex_type.__name__} '
            f'due to possibility of loss of precision.')

        def __check_unix():
            if os.name != 'nt' and _type in (np.float128, np.complex256):
                raise precision_error

        if complex_type in (np.complex128, complex):  # 64-bit complex
            if _type in (np.int64, np.uint64):
                raise precision_error
            __check_unix()
        elif complex_type == np.complex64:  # 32-bit complex
            if _type in (int, float, np.int32, np.int64, np.uint32, np.uint64, np.float64, np.complex128):
                raise precision_error
            __check_unix()

    if isinstance(value, pint.Quantity):
        scaled: np.ndarray = value.to(to_unit).magnitude if to_unit else value.magnitude
        return normalize_complex(scaled, complex_type, None)

    # a list of complex numbers represented by int, float or str
    if isinstance(value, list):
        normalized = [normalize_complex(v, complex_type, to_unit) for v in value]
        return normalized if complex_type == complex else np.array(normalized, dtype=complex_type)

    # complex or real part only
    if type(value) in MTypes.num:
        __check_precision(type(value))
        return complex_type(value)

    # np array
    if isinstance(value, np.ndarray):
        __check_precision(value.dtype.type)
        return value.astype(complex_type)

    # dict representation of complex number
    if isinstance(value, dict):
        real = value.get('re')
        imag = value.get('im')
        assert real is not None or imag is not None, 'Cannot convert an empty dict to complex number.'

        def __combine(_real, _imag):
            _real_list: bool = isinstance(_real, list)
            _imag_list: bool = isinstance(_imag, list)
            if _real_list or _real_list:
                if _real is None:
                    return [__combine(None, i) for i in _imag]
                if _imag is None:
                    return [__combine(r, None) for r in _real]
                # leverage short-circuit evaluation, do not change order
                if _real_list and _imag_list and len(_real) == len(_imag):
                    return [__combine(r, i) for r, i in zip(_real, _imag)]

                raise ValueError('Cannot combine real and imaginary parts of complex numbers.')

            __check_precision(type(_real))
            __check_precision(type(_imag))
            if _real is None:
                return complex_type(_imag) * 1j
            if _imag is None:
                return complex_type(_real)
            return complex_type(_real) + complex_type(_imag) * 1j

        combined = __combine(real, imag)
        return combined if complex_type == complex else np.array(combined, dtype=complex_type)

    # a string, '1+2j'
    # one of 'i', 'I', 'j', 'J' can be used to represent the imaginary unit
    if isinstance(value, str):
        match = MRegEx.complex_str.match(value)
        if match is not None:
            return complex_type(reduce(lambda a, b: a.replace(b, 'j'), 'iIJ', value))

    raise ValueError(f'Cannot convert {value} to complex number.')


def serialize_complex(value):
    '''
    Convert complex number to string.
    '''
    # scalar
    if type(value) in MTypes.complex:
        return {'re': value.real, 'im': value.imag}

    # 1D
    if isinstance(value, (list, tuple)):
        return {'re': [v.real for v in value], 'im': [v.imag for v in value]}

    # ND
    if isinstance(value, np.ndarray):
        return {'re': value.real.tolist(), 'im': value.imag.tolist()}

    raise ValueError(f'Cannot serialize {value}.')


@dataclass(frozen=True)
class MTypes:
    # todo: account for bytes which cannot be naturally serialized to JSON
    primitive = {
        str: lambda v: None if v is None else str(v),
        int: lambda v: None if v is None else int(v),
        float: lambda v: None if v is None else float(v),
        complex: lambda v: None if v is None else complex(v),
        bool: lambda v: None if v is None else bool(v),
        np.bool_: lambda v: None if v is None else bool(v)}

    primitive_name = {v.__name__: v for v in primitive} | {'string': str, 'boolean': bool}

    int_numpy = {np.int8, np.int16, np.int32, np.int64, np.uint8, np.uint16, np.uint32, np.uint64}
    int_python = {int}
    int = int_python | int_numpy
    float_numpy = {np.float16, np.float32, np.float64} | (set() if os.name == 'nt' else {np.float128})
    complex_numpy = {np.complex64, np.complex128} | (set() if os.name == 'nt' else {np.complex256})
    float_python = {float}
    complex_python = {complex}
    float = float_python | float_numpy
    complex = complex_python | complex_numpy
    num_numpy = int_numpy | float_numpy | complex_numpy
    num_python = int_python | float_python | complex_python
    num = num_python | num_numpy
    str_numpy = {np.str_}
    bool_numpy = {np.bool_}
    numpy = num_numpy | str_numpy | bool_numpy


class MEnum(Sequence):
    '''
    Allows to define string types with values limited to a pre-set list of possible values.

    The allowed values can be provided as a list of strings, the keys of which will be identical to values.
    Alternatively, they can be provided as key-value pairs.

    For example:
        some_variable = MEnum(['a', 'b', 'c'])
        some_variable = MEnum(a='a', b='b', c='c')

    The values are stored in __dict__ and can be accessed as attributes:
        some_variable.a # gives 'a'

    For description of each possible value, it can be organized into a dictionary.

    For example:
        some_variable = MEnum(['a', 'b', 'c'], m_descriptions={'a': 'first', 'b': 'second', 'c': 'third'})
    '''

    def __init__(self, *args, **kwargs):
        # Supports one big list in place of args
        if len(args) == 1 and isinstance(args[0], list):
            args = args[0]

        self._descriptions: Dict[str, str] = {}
        if 'm_descriptions' in kwargs:
            self._descriptions = kwargs.pop('m_descriptions')

        # If non-named arguments are given, the default is to have them placed
        # into a dictionary with their string value as both the enum name and
        # the value.
        for arg in args:
            if arg in kwargs:
                raise ValueError(f"Duplicate value '{arg}' provided for enum")
            kwargs[arg] = arg

        self._list = list(kwargs.values())
        self._values = set(kwargs.values())  # For allowing constant time member check

        for enum_value in self._values:
            if not isinstance(enum_value, str):
                raise TypeError(f'MEnum value {enum_value} is not a string.')

        self.__dict__.update(kwargs)

    def set_description(self, value: str, description: str):
        if value not in self._values:
            raise ValueError(f'{value} is not a value of this MEnum')
        self._descriptions[value] = description

    def get_description(self, value: str) -> str:
        if value not in self._values:
            raise ValueError(f'{value} is not a value of this MEnum')
        return self._descriptions.get(value, '')

    def get_all_descriptions(self) -> Dict[str, str]:
        return self._descriptions

    def get_all_values(self) -> set:
        return self._values

    # no need to implement __getattr__ as all attributes are stored in the __dict__
    # def __getattr__(self, attr):
    #     pass

    def __getitem__(self, index):
        return self._list[index]

    def __len__(self):
        return len(self._list)


class MQuantity:
    '''
    A simple wrapper to represent complex quantities that may have multiple values,
    additional attributes, and more.
    '''
    name: str = None
    value: Any = None
    unit: Optional[pint.Unit] = None
    original_unit: Optional[pint.Unit] = None
    attributes: dict = None

    def __init__(
            self,
            in_name: Optional[str],
            in_value: Any,
            in_unit: Optional[pint.Unit] = None,
            in_attributes: Optional[dict] = None):
        '''
        The validation of value/unit/attribute is performed at 'MSection' level.
        '''
        self.name = in_name
        if self.name:
            assert isinstance(self.name, str), 'Name must be a string'

        self.unit = None
        if isinstance(in_value, pint.Quantity):
            self.value = in_value.m  # magnitude
            self.unit = in_value.u  # unit
            assert in_unit is None, f'Unit is already defined in the value {in_value}'
        else:
            # the input argument is not a pint quantity
            # the unit is set to None
            self.value = in_value
            if isinstance(in_unit, pint.Unit):
                self.unit = in_unit
            elif isinstance(in_unit, str):
                self.unit = ureg.parse_units(in_unit)

        self.original_unit = self.unit

        self.attributes: dict = {}
        if in_attributes is not None:
            self.attributes.update(**in_attributes)
            self.__dict__.update(**in_attributes)

    @staticmethod
    def wrap(in_value: Any, in_name: Optional[str] = None):
        '''
        Syntax sugar to wrap a value into a MQuantity. The name is optional.

        This would be useful for non-variadic primitive quantities with additional attributes.
        '''
        return MQuantity(in_name, in_value)

    def __repr__(self):
        return self.name if self.name else 'Unnamed quantity'

    def m_set_attribute(self, name, value):
        '''
        Validation is done outside this container
        '''
        self.attributes[name] = value


class MSubSectionList(list):
    def __init__(self, section, sub_section_def):
        self.section = section
        self.sub_section_def = sub_section_def
        super().__init__()

    def __setitem__(self, key, value):
        raise NotImplementedError('You can only append subsections.')

    def __delitem__(self, key):
        old_value = self[key]
        list.__delitem__(self, key)
        for index in range(key, len(self)):
            self[index].m_parent_index = index

        # noinspection PyProtectedMember
        self.section._on_remove_sub_section(self.sub_section_def, old_value)

    def __setslice__(self, i, j, sequence):
        raise NotImplementedError('You can only append subsections.')

    def __delslice__(self, i, j):
        raise NotImplementedError('You can only append subsections.')

    def append(self, value):
        list.append(self, value)
        if value is not None:
            # noinspection PyProtectedMember
            self.section._on_add_sub_section(self.sub_section_def, value, len(self) - 1)

    def pop(self, index=...):
        raise NotImplementedError('You can only append subsections.')

    def extend(self, new_value):
        start_index = len(self)
        list.extend(self, new_value)
        for index, value in enumerate(new_value):
            # noinspection PyProtectedMember
            self.section._on_add_sub_section(self.sub_section_def, value, start_index + index)

    def insert(self, i, element):
        raise NotImplementedError('You can only append subsections.')

    def remove(self, element):
        raise NotImplementedError('You can only append subsections.')

    def reverse(self):
        raise NotImplementedError('You can only append subsections.')

    def sort(self, *, key=..., reverse=...):
        raise NotImplementedError('You can only append subsections.')

    def clear(self):
        old_values = list(self)
        list.clear(self)
        for old_value in old_values:
            # noinspection PyProtectedMember
            self.section._on_remove_sub_section(self.sub_section_def, old_value)


@dataclass
class ReferenceURL:
    fragment: str
    archive_url: str
    url_parts: SplitResult

    def __init__(self, url: str):
        if '#' not in url:
            url = f'#{url}'

        self.url_parts = urlsplit(url)
        archive_url = urlunsplit(self.url_parts[0:4] + ('',))
        self.archive_url = None if archive_url is None else archive_url
        self.fragment = self.url_parts.fragment


class Annotation:
    ''' Base class for annotations. '''

    def m_to_dict(self):
        '''
        Returns a JSON serializable representation that is used for exporting the
        annotation to JSON.
        '''
        return str(self.__class__.__name__)


class DefinitionAnnotation(Annotation):
    ''' Base class for annotations for definitions. '''

    def __init__(self):
        self.definition = None

    def init_annotation(self, definition):
        self.definition = definition


class SectionAnnotation(DefinitionAnnotation):
    '''
    Special annotation class for section definition that allows to auto add annotations
    to section instances.
    '''

    def new(self, section) -> Dict[str, Any]:
        return {}


def to_dict(entries):
    if isinstance(entries, list):
        return [to_dict(entry) for entry in entries]

    # noinspection PyBroadException
    try:
        entries = entries.m_to_dict()
    except Exception:
        pass

    return entries


def convert_to(from_magnitude, from_unit: Optional[ureg.Unit], to_unit: Optional[ureg.Unit]):
    '''
    Convert a magnitude from one unit to another.

    Arguments:
        from_magnitude: the magnitude to be converted
        from_unit: the unit of the magnitude
        to_unit: the unit to convert to

    Return:
        the converted magnitude
    '''

    if to_unit is None:
        return from_magnitude

    from_quantity: ureg.Quantity = from_magnitude * from_unit

    return from_quantity.to(to_unit).m


def __similarity_match(candidates: list, name: str):
    '''
    Use similarity to find the best match for a name.
    '''
    similarity: list = [SequenceMatcher(None, v.name.upper(), name.upper()).ratio() for v in candidates]

    return candidates[similarity.index(max(similarity))]


def resolve_variadic_name(definitions: dict, name: str, hint: Optional[str] = None):
    '''
    For properties with variadic names, it is necessary to check all possible definitions
    in the schema to find the unique and correct definition that matches the naming pattern.

    In the schema defines a property with the name 'FOO_bar', implying the prefix 'FOO' is
    merely a placeholder, the actual name in the data can be anything, such as 'a_bar' or 'b_bar'.

    This method checks each definition name by replacing the placeholder with '.*' and then check if
    the property name matches the pattern. If it does, it returns the corresponding definition.

    For example, the definition name 'FOO_bar' will be replaced by '.*_bar', which further matches
    'a_bar', 'aa_bar', etc.

    In case of multiple quantities with identical template/variadic patterns, the following strategy
    is used:
        1. Check all quantities and collect all qualified quantities that match the naming pattern
            in a candidate list.
        2. Use the optionally provided hint string, which shall be one of attribute names of the desired
            quantity. Check all candidates if this attribute exists. The existence of a hint attribute
            prioritize this quantity, and it will be put into a prioritized list.
        3. If the prioritized candidate list contains multiple matches, use name similarity determine
            which to be used.
        4. If no hint is provided, or no candidate has the hint attribute, check all quantities in the
            first candidate list and use name similarity to determine which to be used.

    '''

    # check the exact name match
    if name in definitions:
        return definitions[name]

    # check naming pattern match
    candidates: list = []
    for definition in set(definitions.values()):
        if not definition.variable:
            continue

        name_pattern = re.sub(r'^([a-z0-9_]*)[A-Z0-9]+([a-z0-9_]*)$', r'\1[a-z0-9]+\2', definition.name)
        if re.match(name_pattern, name):
            candidates.append(definition)

    if len(candidates) == 0:
        raise ValueError(f'Cannot find a proper definition for name {name}')

    if len(candidates) == 1:
        return candidates[0]

    hinted_candidates: list = []
    if hint is not None:
        for definition in candidates:
            try:
                if resolve_variadic_name(definition.all_attributes, hint):
                    hinted_candidates.append(definition)
            except ValueError:
                pass

    if len(hinted_candidates) == 1:
        return hinted_candidates[0]

    # multiple matches, check similarity
    if len(hinted_candidates) > 1:
        return __similarity_match(hinted_candidates, name)

    return __similarity_match(candidates, name)


def retrieve_attribute(section, definition, attr_name: str) -> tuple:
    '''
    Retrieve the attribute of a definition by its name.
    In the case of variadic/template name, the name is also resolved by checking naming pattern.
    '''
    from nomad.metainfo.metainfo import Definition

    # find the section or quantity where attribute is defined
    if definition is None:
        tgt_def = section
    elif isinstance(definition, Definition):
        tgt_def = definition
    else:
        tgt_def = resolve_variadic_name(section.all_quantities, definition)
    if tgt_def is None:
        raise ValueError(f'Cannot find the definition by the given {definition}')

    # find the corresponding attribute
    tgt_attr = resolve_variadic_name(tgt_def.all_attributes, attr_name)
    if tgt_attr is None:
        raise ValueError('The given attribute name is not found in the given property.')

    return tgt_def, tgt_attr


def validate_allowable_unit(
        dimensionality: Optional[str],
        allowable_list: Union[str, list, pint.Unit, pint.Quantity]) -> bool:
    '''
    For a given list of units, e.g., ['m', 'cm', 'mm'], and a target NX unit token such as 'NX_LENGTH',
    this function checks the compatibility of the target unit with the list of units.

    Returns:
        True if ALL units are compatible with the unit token (dimensionality).
        False if at least one unit cannot be represented by the unit token (dimensionality).
    '''
    if not dimensionality:
        return True

    if isinstance(allowable_list, str):
        if dimensionality in ('1', 'dimensionless'):
            return ureg.Quantity(1, allowable_list).dimensionless

        try:
            return ureg.Quantity(1, allowable_list).check(dimensionality)
        except KeyError:
            return False

    if isinstance(allowable_list, (pint.Unit, pint.Quantity)):
        if dimensionality in ('1', 'dimensionless'):
            return allowable_list.dimensionless

        return allowable_list.dimensionality == dimensionality

    for unit in allowable_list:
        if not validate_allowable_unit(dimensionality, unit):
            return False

    return True


def default_hash():
    '''
    Returns a hash object using the designated hash algorithm.
    '''
    return hashlib.new(__hash_method)


def split_python_definition(definition_with_id: str) -> Tuple[list, Optional[str]]:
    '''
    Split a Python type name into names and an optional id.

    Example:
        my_package.my_section@my_id  ==> (['my_package', 'my_section'], 'my_id')

        my_package.my_section       ==> (['my_package', 'my_section'], None)
    '''
    if '@' not in definition_with_id:
        return definition_with_id.split('.'), None

    definition_names, definition_id = definition_with_id.split('@')
    return definition_names.split('.'), definition_id


def check_dimensionality(quantity_def, unit: Optional[pint.Unit]) -> None:
    if quantity_def is None or unit is None:
        return

    dimensionality = getattr(quantity_def, 'dimensionality', None)

    if dimensionality is None:  # not set, do not validate
        return

    if dimensionality in ('dimensionless', '1') and unit.dimensionless:  # dimensionless
        return

    if dimensionality == 'transformation':
        # todo: check transformation dimensionality
        return

    if ureg.Quantity(1 * unit).check(dimensionality):  # dimensional
        return

    raise TypeError(f'Dimensionality {dimensionality} is not met by unit {unit}')


def check_unit(unit: Union[str, pint.Unit]) -> None:
    '''Check that the unit is valid.
    '''
    if isinstance(unit, str):
        unit_str = unit
    elif isinstance(unit, pint.Unit):
        unit_str = str(unit)
    else:
        raise TypeError('Units must be given as str or pint Unit instances.')

    # Explicitly providing a Pint delta-unit is not currently allowed.
    # Implicit conversions are fine as MathJS on the frontend supports them.
    if any(x in unit_str for x in _delta_symbols):
        raise TypeError('Explicit Pint "delta"-units are not yet supported.')


def to_section_def(section_def):
    '''
    Resolves duck-typing for values that are section definitions or section classes to
    section definition.
    '''
    return section_def.m_def if isinstance(section_def, type) else section_def  # type: ignore


def to_numpy(np_type, shape: list, unit: Optional[pint.Unit], definition, value: Any):
    check_dimensionality(definition, unit)

    if isinstance(value, pint.Quantity):
        # if flexible unit is set, do not check unit in the definition
        # it will be handled specially
        # the stored unit would not be serialized
        flexible_unit = getattr(definition, 'flexible_unit', False)

        if not flexible_unit and unit is None:
            raise TypeError(f'The quantity {definition} does not have a unit, but value {value} does.')

        if type(value.magnitude) == np.ndarray and np_type != value.dtype:
            value = value.astype(np_type)

        if not flexible_unit:
            value = value.to(unit).magnitude
        else:
            value = value.magnitude

    if isinstance(value, pd.DataFrame):
        try:
            value = value.to_numpy()
        except AttributeError:
            raise AttributeError(
                f'Could not convert value {value} of type pandas.Dataframe to a numpy array')

    if np_type in MTypes.complex:
        value = normalize_complex(value, np_type, unit)

    if type(value) != np.ndarray:
        if len(shape) > 0:
            try:
                value = np.asarray(value, dtype=np_type)
            except TypeError:
                raise TypeError(f'Could not convert value {value} of {definition} to a numpy array')
        elif type(value) != np_type:
            try:
                value = np_type(value)
            except TypeError:
                raise TypeError(f'Could not convert value {value} of {definition} to a numpy scalar')
    elif value.dtype != np_type and np_type in MTypes.complex:
        try:
            value = value.astype(np_type)
        except TypeError:
            raise TypeError(f'Could not convert value {value} of {definition} to a numpy array')

    return value


def __validate_shape(section, dimension: Union[str, int], length: int) -> bool:
    if isinstance(dimension, int):
        return dimension == length

    if not isinstance(dimension, str):
        raise TypeError(f'Invalid dimension type {type(dimension)}')

    if dimension.isidentifier():
        return dimension == getattr(section, dimension)

    m = re.match(MRegEx.index_range, dimension)
    start = int(m.group(1))
    end = -1 if m.group(2) == '*' else int(m.group(2))
    return start <= length and (end == -1 or length <= end)


def validate_shape(section, quantity_def, value: Any) -> bool:
    quantity_shape: list = quantity_def.shape

    if type(value) == np.ndarray:
        value_shape = value.shape
    elif isinstance(value, list) and not isinstance(value, MEnum):
        value_shape = (len(value),)
    else:
        value_shape = ()

    if len(value_shape) != len(quantity_shape):
        return False

    return all(__validate_shape(section, x, y) for x, y in zip(quantity_shape, value_shape))


def dict_to_named_list(data) -> list:
    if not isinstance(data, dict):
        return data

    results: list = []
    for key, value in data.items():
        if value is None:
            value = {}
        value.update(dict(name=key))
        results.append(value)
    return results


def validate_url(url_str: str) -> Optional[str]:
    if url_str is None:
        return None

    if not isinstance(url_str, str):
        raise TypeError('Links need to be given as URL strings')
    if re.match(MRegEx.url, url_str) is None:
        raise ValueError('The given URL is not valid')

    return url_str


def __parse_datetime(datetime_str: str) -> datetime:
    # removing trailing spaces and replacing the potential white space between date and time with char "T"
    if datetime_str[0].isdigit():
        datetime_str = datetime_str.strip().replace(' ', 'T')

    try:
        return aniso8601.parse_datetime(datetime_str)
    except ValueError:
        pass

    try:
        date_value = aniso8601.parse_date(datetime_str)
        if isinstance(date_value, datetime):
            return date_value
    except ValueError:
        pass

    # noinspection PyBroadException
    try:
        return email.utils.parsedate_to_datetime(datetime_str)
    except Exception:
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

    if 'GMT' in datetime_str:
        dt_copy = datetime_str
        dt_split = dt_copy.split('GMT')
        tzinfo = dt_split[1].strip()
        if len(tzinfo) == 2:
            tzinfo = f'{tzinfo[0]}{tzinfo[1]:0>2}00'
        dt_copy = f'{dt_split[0]}GMT{tzinfo}'
        try:
            return datetime.strptime(dt_copy, '%Y%m%d_%H:%M:%S_%Z%z')
        except ValueError:
            pass

    try:
        return datetime.fromisoformat(datetime_str)
    except ValueError:
        pass

    raise TypeError(f'Invalid date literal {datetime_str}')


def normalize_datetime(value) -> Optional[datetime]:
    if value is None:
        return None

    if isinstance(value, str):
        value = __parse_datetime(value)

    elif isinstance(value, (int, float)):
        value = datetime.fromtimestamp(value)

    elif isinstance(value, pint.Quantity):
        value = datetime.fromtimestamp(value.magnitude)

    elif not isinstance(value, datetime) and isinstance(value, date):
        value = datetime.combine(value, datetime.min.time())

    if not isinstance(value, datetime):
        raise TypeError(f'{value} is not a datetime.')

    if value.tzinfo is None:
        value = value.replace(tzinfo=pytz.utc)
    else:
        value = value.astimezone(pytz.utc)

    return value


def camel_case_to_snake_case(obj: dict):
    for k, v in list(obj.items()):
        if k != k.lower() and k != k.upper() and "_" not in k:
            snake_case_key = re.sub(r'(?<!^)(?=[A-Z])', '_', k).lower()
            obj[snake_case_key] = v
            del obj[k]
            k = snake_case_key
        if isinstance(v, dict):
            obj[k] = camel_case_to_snake_case(v)
        if isinstance(v, list):
            for i, item in enumerate(v):
                if isinstance(item, dict):
                    obj[k][i] = camel_case_to_snake_case(item)
    return obj
