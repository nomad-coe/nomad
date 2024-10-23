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
from __future__ import annotations

import hashlib
import re
from difflib import SequenceMatcher
from typing import Any, Optional

import pint

from nomad.metainfo.data_type import Enum
from nomad.units import ureg

__hash_method = 'sha1'  # choose from hashlib.algorithms_guaranteed

MEnum = Enum  # type: ignore


class MQuantity:
    """
    A simple wrapper to represent complex quantities that may have multiple values,
    additional attributes, and more.
    """

    def __init__(
        self,
        in_name: str | None,
        in_value: Any,
        in_unit: pint.Unit | None = None,
        in_attributes: dict | None = None,
    ):
        """
        The validation of value/unit/attribute is performed at 'MSection' level.
        """
        self.name: str | None = in_name
        if self.name:
            assert isinstance(self.name, str), 'Name must be a string'

        self.unit: pint.Unit | None = None
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

        self.original_unit: pint.Unit | None = self.unit

        self.attributes: dict = {}
        if in_attributes is not None:
            self.attributes.update(**in_attributes)
            self.__dict__.update(**in_attributes)

    @staticmethod
    def from_dict(name: str, data: dict) -> MQuantity:
        m_quantity = MQuantity(
            name,
            data['m_value'],
            data.get('m_unit', None),
            data.get('m_attributes', None),
        )

        if 'm_original_unit' in data:
            m_quantity.original_unit = ureg.parse_units(data['m_original_unit'])

        return m_quantity

    @staticmethod
    def wrap(in_value: Any, in_name: str | None = None):
        """
        Syntax sugar to wrap a value into a MQuantity. The name is optional.

        This would be useful for non-variadic primitive quantities with additional attributes.
        """
        return MQuantity(in_name, in_value)

    def __repr__(self):
        return self.name if self.name else 'Unnamed quantity'

    def m_set_attribute(self, name, value):
        """
        Validation is done outside this container
        """
        self.attributes[name] = value

    def get(self):
        if self.unit:
            return ureg.Quantity(self.value, self.unit)

        return self.value


class MSubSectionList(list):
    def __init__(self, section, sub_section_def):
        self.section = section
        self.sub_section_def = sub_section_def
        super().__init__()

    def __setitem__(self, key, value):
        old_value = self[key]
        super().__setitem__(key, value)
        # noinspection PyProtectedMember
        self.section._on_add_sub_section(self.sub_section_def, value, key)
        # noinspection PyProtectedMember
        self.section._on_remove_sub_section(self.sub_section_def, old_value)

    def __getitem__(self, item):
        if isinstance(item, str):
            for sub_section in self:
                if sub_section.m_key == item:
                    return sub_section
            raise KeyError(f'No subsection keyed {item} found.')

        return super().__getitem__(item)

    def __delitem__(self, key):
        self.pop(key)

    def __setslice__(self, i, j, sequence):
        raise NotImplementedError('You can only append subsections.')

    def __delslice__(self, i, j):
        raise NotImplementedError('You can only append subsections.')

    def append(self, value):
        super().append(value)
        # noinspection PyProtectedMember
        self.section._on_add_sub_section(self.sub_section_def, value, len(self) - 1)

    def pop(self, key=-1):
        if key < 0:
            key = key + len(self)
        old_value = super().pop(key)
        for index in range(key, len(self)):
            # noinspection PyProtectedMember
            self.section._on_add_sub_section(self.sub_section_def, self[index], index)

        # noinspection PyProtectedMember
        self.section._on_remove_sub_section(self.sub_section_def, old_value)

        return old_value

    def extend(self, new_value):
        start_index = len(self)
        super().extend(new_value)
        for index, value in enumerate(new_value):
            # noinspection PyProtectedMember
            self.section._on_add_sub_section(
                self.sub_section_def, value, start_index + index
            )

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
        super().clear()
        for old_value in old_values:
            # noinspection PyProtectedMember
            self.section._on_remove_sub_section(self.sub_section_def, old_value)

    def has_duplicated_key(self) -> bool:
        """
        Check if each section has a unique key.
        The key is from `.m_key` or the index of the section.
        """
        unique_keys: set = set()
        for index, section in enumerate(self):
            if section is not None:
                if section.m_key in unique_keys:
                    return True
                unique_keys.add(section.m_key)

        return False


class Annotation:
    """Base class for annotations."""

    def m_to_dict(self):
        """
        Returns a JSON serializable representation that is used for exporting the
        annotation to JSON.
        """
        return str(self.__class__.__name__)


class DefinitionAnnotation(Annotation):
    """Base class for annotations for definitions."""

    def __init__(self):
        self.definition = None

    def init_annotation(self, definition):
        self.definition = definition


class SectionAnnotation(DefinitionAnnotation):
    """
    Special annotation class for section definition that allows to auto add annotations
    to section instances.
    """

    def new(self, section) -> dict[str, Any]:
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


def convert_to(from_magnitude, from_unit: ureg.Unit | None, to_unit: ureg.Unit | None):
    """
    Convert a magnitude from one unit to another.

    Arguments:
        from_magnitude: the magnitude to be converted
        from_unit: the unit of the magnitude
        to_unit: the unit to convert to

    Return:
        the converted magnitude
    """

    if to_unit is None:
        return from_magnitude

    from_quantity: ureg.Quantity = from_magnitude * from_unit

    return from_quantity.to(to_unit).m


def get_namefit(name: str, concept_name: str, name_any: bool = False) -> int:
    """
    Checks if a given name corresponds to a specified concept name.

    The function evaluates whether the provided `name` matches the `concept_name`, allowing
    for certain variations. A group of uppercase letters in the `concept_name` is treated
    as a freely choosable part of this name.

    If a match is found, this function returns twice the length of the `concept_name` for an
    exact match. Otherwise, it returns the number of matching characters (case insensitive) or
    zero if `name_any` is set to True. If no match is found, it returns -1.

    The function treats uppercase groups independently and counts lowercase matches without
    regard to uppercase group lengths. For example, calling `get_namefit("my_fancy_yet_long_name", "my_SOME_name")`
    would yield a score of 8 for the lowercase matches `my_..._name`.

    All characters in `[a-zA-Z0-9_.]` are considered for matching against uppercase letters.
    Using any other character in the `name` will result in a return value of -1.
    Periods at the beginning or end of the `name` are not allowed, and only exact matches will be considered.

    Examples:

        * `get_namefit("test_name", "TEST_name")` returns 9
        * `get_namefit("te_name", "TEST_name")` returns 7
        * `get_namefit("my_other_name", "TEST_name")` returns 5
        * `get_namefit("test_name", "test_name")` returns 18
        * `get_namefit("test_other", "test_name")` returns -1
        * `get_namefit("something", "XXXX")` returns 0
        * `get_namefit("something", "OTHER")` returns 1

    Args:
        name (str): The name to check for matching.
        concept_name (str): The concept name to match against.
        name_any (bool, optional):
            If True, accepts any name and returns either 0 (match) or -1 (no match).
            Defaults to False.

    Returns:
        int: -1 if no match is found, the number of matching characters (case insensitive),
             or twice the length of `concept_name` for an exact match.
    """
    if concept_name == name:
        return len(concept_name) * 2
    if name.startswith('.') or name.endswith('.'):
        # Don't match anything with a dot at the beginning or end
        return -1

    uppercase_parts_pattern = re.compile(r'[A-Z]+(?:_[A-Z]+)*')
    uppercase_parts = uppercase_parts_pattern.findall(concept_name)

    path_regex = r'([a-zA-Z0-9_.]+)'
    regex_name = concept_name
    uppercase_count = sum(len(part) for part in uppercase_parts)

    for up in uppercase_parts:
        regex_name = regex_name.replace(up, path_regex)

    # Compile the full regex for matching the HDF name
    name_regex = re.compile(rf'^{regex_name}$')
    name_match = name_regex.fullmatch(name)

    if name_match is None:
        return 0 if name_any else -1

    match_count = sum(
        1
        for up, match in zip(uppercase_parts, name_match.groups())
        for s1, s2 in zip(up.upper(), match.upper())
        if s1 == s2
    )

    return len(concept_name) + match_count - uppercase_count


def resolve_variadic_name(definitions: dict, name: str, hint: Optional[str] = None):
    """
    Resolves a property name with variadic patterns to its corresponding definition in the schema.

    For properties with variadic names, it is necessary to check all possible definitions
    in the schema to find the unique and correct definition that matches the naming pattern.

    In the schema, definitions may include properties with placeholders that can match any
    naming segment. For example, the definition name 'FOO_bar' indicates that 'FOO' is a
    placeholder, meaning the actual name in the data could be anything like 'a_bar' or 'b_bar'.

    This function checks each definition to find matches based on the following criteria:

    1. **Exact Match**: If the provided name exactly matches a definition, that definition is returned.
    2. **Pattern Matching**: Definitions that match the naming pattern derived from the name are collected.
    3. **Candidate Collection**: All candidates matching the pattern are gathered.
    4. **Hint Prioritization**: Candidates containing the hint attribute are prioritized.
    5. **Similarity Assessment**: If multiple candidates exist, name similarity is used to determine the best match.

    In case of multiple quantities with identical template/variadic patterns, the following strategy
    is used:
        1. Check all quantities and collect all qualified quantities that match the naming pattern
           in a candidate dictionary, with the values being an integer indicating how big the match is.
        2. Use the optionally provided hint string, which shall be one of attribute names of the desired
            quantity. Check all candidates if this attribute exists. The existence of a hint attribute
            prioritize this quantity, and it will be put into a prioritized dictionary.
        3. If the prioritized candidate dictionary contains multiple matches, the strongest match is selected
           (i.e., the one with the highest integer returned by `get_namefit`.
        4. If no hint is provided, or no candidate has the hint attribute, check all quantities in the
            first candidate dictionary and select the strongest match.

    Args:
        definitions (dict): A dictionary of definitions where keys are property names and values are their corresponding definitions.
        name (str): The property name to resolve against the definitions.
        hint (Optional[str]): An optional hint that, if present in the candidate attributes, prioritizes that candidate.

    Returns:
        The best-matching definition based on the criteria.

    Raises:
        ValueError: If the definitions dictionary is empty or if no proper definition can be found for the given name.
    """
    if len(definitions) == 0:
        raise ValueError('The definitions dictionary cannot be empty.')

    # Check for an exact name match
    if name in definitions:
        return definitions[name]

    candidates = {}
    hint_candidates = {}

    for definition in definitions:
        match_score = get_namefit(name, definition)
        if match_score > 0:
            candidates[definition] = match_score
            # Check if the hint exists in the definition
            if hint and hint in definition.all_attributes:
                hint_candidates[definition] = match_score

    # Prioritize hint candidates if any
    if hint_candidates:
        return definitions[max(hint_candidates, key=hint_candidates.get)]

    # If no hint candidates, find the best candidate from all
    if candidates:
        return definitions[max(candidates, key=candidates.get)]

    raise ValueError(f'Cannot find a proper definition for name "{name}".')


def validate_allowable_unit(
    dimensionality: str | None,
    allowable_list: str | list | pint.Unit | pint.Quantity,
) -> bool:
    """
    For a given list of units, e.g., ['m', 'cm', 'mm'], and a target NX unit token such as 'NX_LENGTH',
    this function checks the compatibility of the target unit with the list of units.

    Returns:
        True if ALL units are compatible with the unit token (dimensionality).
        False if at least one unit cannot be represented by the unit token (dimensionality).
    """
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
    """
    Returns a hash object using the designated hash algorithm.
    """
    return hashlib.new(__hash_method)


def split_python_definition(definition_with_id: str) -> tuple[list, str | None]:
    """
    Split a Python type name into names and an optional id.

    Example:
        my_package.my_section@my_id  ==> (['my_package', 'my_section'], 'my_id')

        my_package.my_section       ==> (['my_package', 'my_section'], None)

        my_package/section_definitions/0 ==> (['my_package', 'section_definitions/0'], None)
    """

    def __split(name: str):
        # The definition name must contain at least one dot which comes from the module name.
        # The actual definition could be either a path (e.g., my_package/section_definitions/0)
        # or a name (e.g., my_section).
        # If it is a path (e.g., a.b.c/section_definitions/0), after split at '.', the last segment
        # (c/section_definitions/0) contains the package name (c). It needs to be relocated.
        segments: list = name.split('.')
        if '/' in segments[-1]:
            segments.extend(segments.pop().split('/', 1))
        return segments

    if '@' not in definition_with_id:
        return __split(definition_with_id), None

    definition_names, definition_id = definition_with_id.split('@')
    return __split(definition_names), definition_id


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


def camel_case_to_snake_case(obj: dict):
    for k, v in list(obj.items()):
        if k != k.lower() and k != k.upper() and '_' not in k:
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
