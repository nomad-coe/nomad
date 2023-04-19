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
import base64
from copy import deepcopy
import importlib
import inspect
import itertools
import json
import re
import sys
from collections.abc import Iterable as IterableABC
from functools import reduce
from pydantic import parse_obj_as, ValidationError, BaseModel, Field
from typing import (
    Any, Callable as TypingCallable, Dict, Iterable, List, Optional, Set, Tuple, Type, TypeVar, Union, cast, ClassVar)
import docstring_parser
import jmespath
import numpy as np
import pandas as pd
import pint

from nomad.config import process
from nomad.metainfo.util import (
    Annotation, DefinitionAnnotation, MEnum, MQuantity, MRegEx, MSubSectionList, MTypes, ReferenceURL,
    SectionAnnotation, _delta_symbols, check_dimensionality, check_unit, convert_to, default_hash, dict_to_named_list,
    normalize_complex, normalize_datetime, resolve_variadic_name, retrieve_attribute, serialize_complex,
    split_python_definition, to_dict, to_numpy, to_section_def, validate_shape, validate_url)
from nomad.units import ureg as units

# todo: remove magic comment after upgrading pylint
# https://github.com/PyCQA/pylint/issues/5342
m_package: Optional[Package] = None  # pylint: disable=used-before-assignment

is_bootstrapping = True
Elasticsearch = TypeVar('Elasticsearch')
MSectionBound = TypeVar('MSectionBound', bound='MSection')
SectionDefOrCls = Union['Section', 'SectionProxy', Type['MSection']]
T = TypeVar('T')

_unset_value = '__UNSET__'
_HASH_OBJ = Type['hashlib._Hash']  # type: ignore


def _check_definition_id(target_id, tgt_section: MSectionBound) -> MSectionBound:
    '''
    Ensure section definition id matches the target id.
    '''
    if target_id is not None and tgt_section.definition_id != target_id:
        raise MetainfoReferenceError(f'Could not resolve {target_id}, id mismatch')
    return tgt_section


# Make pylint believe all bootstrap quantities are actual properties even though
# we have to initialize them to None due to bootstrapping
_placeholder_quantity: Quantity = property()  # type: ignore # pylint: disable=used-before-assignment
if True:
    _placeholder_quantity: Quantity = None  # type: ignore


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


class MProxy:
    '''
    A placeholder object that acts as reference to a value that is not yet resolved.

    A proxy is a replacement for an actual section (or quantity) that the proxy represents.
    The replaced section (or quantity) is identified by a reference. References
    are URL strings that identify a section (or quantity).

    If a proxy is accessed (i.e. like its proxies' counterpart would be accessed), it
    tries to resolve its reference and access the proxied element. If the reference
    cannot be resolved an exception is raised.

    There are different kinds of reference urls. Here are a few examples:

    .. code-block::
        /run/0/calculation/1  # same archive (legacy version)
        #/run/0/calculation/1  # same archive
        ../upload/archive/mainfile/{mainfile}#/run/0  # same upload
        /entries/{entry_id}/archive#/run/0/calculation/1  # same NOMAD
        /uploads/{upload_id}/archive/{entry_id}#/run/0/calculation/1  # same NOMAD
        https://my-oasis.de/api/v1/uploads/{upload_id}/archive/{entry_id}#/run/0/calculation/1  # global

    A URL can have 3 relevant path. The archive path as anchor string (after the `#`).
    The API part (more or less the URL's path). The oasis part (more or less the domain).

    The archive path starts at the `EntryArchive` root and uses property names and indices
    to navigate. The API and oasis parts correspond to NOMAD's v1 api. A path starting
    with `../upload` replaces `.../api/v1/uploads/{upload_id}` and allows to access the
    upload of the archive that contains the reference.

    The actual algorithm for resolving proxies is in `MSection.m_resolve()`.

    Attributes:
        m_proxy_value: The reference represented as a URL string.
        m_proxy_section:
            The section context, i.e. the section that this proxy is contained in. This
            section will provide the context for resolving the reference. For example,
            if the reference only has an archive part, this archive part is resolved
            starting with the root of `m_proxy_section`.
        m_proxy_context:
            Optional Context instance. Default is None and the m_context of the m_proxy_section
            is used.
        m_proxy_type:
            The quantity definition. Typically, MProxy is used for proxy-ing sections. With
            this set, the proxy will still act as a normal section proxy, but it will
            be used by quantities of type `QuantityReference` to resolve and return
            a quantity value.
    '''

    def __init__(
            self,
            m_proxy_value: Union[str, int, dict],
            m_proxy_section: MSection = None,
            m_proxy_context: Context = None,
            m_proxy_type: Reference = None):
        self.m_proxy_value = m_proxy_value
        self.m_proxy_section = m_proxy_section
        self.m_proxy_resolved = None
        self.m_proxy_type = m_proxy_type
        self.m_proxy_context = m_proxy_context

    def _set_resolved(self, resolved):
        self.m_proxy_resolved = resolved

        if self.m_proxy_resolved is not None and isinstance(self, MProxy):
            setattr(self, '__class__', self.m_proxy_resolved.__class__)
            self.__dict__.update(**self.m_proxy_resolved.__dict__)

    def m_proxy_resolve(self):
        if not self.m_proxy_resolved:
            if self.m_proxy_type and (self.m_proxy_context or self.m_proxy_section):
                self._set_resolved(self.m_proxy_type.resolve(self))

        return self.m_proxy_resolved

    def __getattr__(self, key):
        if self.m_proxy_resolve() is not None:
            return getattr(self.m_proxy_resolved, key)

        raise MetainfoReferenceError(f'could not resolve {self.m_proxy_value}')

    def __repr__(self):
        return f'{self.__class__.__name__}({self.m_proxy_value})'


class SectionProxy(MProxy):
    def __init__(self, m_proxy_value, **kwargs):
        if 'm_proxy_type' in kwargs:
            super().__init__(m_proxy_value=m_proxy_value, **kwargs)
        else:
            super().__init__(m_proxy_value=m_proxy_value, m_proxy_type=SectionReference, **kwargs)

    # TODO recursive proxy stuff
    def _resolve_name(self, name: str, context: Definition) -> Definition:
        if context is None:
            return None

        if context.name == name and context != self.m_proxy_section:
            return context

        if isinstance(context, Section):
            resolved = context.all_aliases.get(name)
            if resolved and resolved != self.m_proxy_section:
                return resolved

            resolved = context.all_inner_section_definitions.get(name)
            if resolved and resolved != self.m_proxy_section:
                return resolved

        if isinstance(context, Package):
            resolved = context.all_definitions.get(name)
            if resolved and resolved != self.m_proxy_section:
                return resolved

        parent = context.m_parent
        if isinstance(parent, Definition):
            return self._resolve_name(name, cast(Definition, parent))

        return None

    def m_proxy_resolve(self):
        if '#' in self.m_proxy_value or '/' in self.m_proxy_value:
            # This is not a python reference, use the usual mechanism
            return super().m_proxy_resolve()

        if '.' in self.m_proxy_value:
            # Try to interpret as python class name
            python_name, definition_id = split_python_definition(self.m_proxy_value)
            package_name = '.'.join(python_name[:-1])
            section_name = python_name[-1]

            try:
                module = importlib.import_module(package_name)
                cls = getattr(module, section_name)
                if cls.m_def:
                    if not definition_id or cls.m_def.definition_id == definition_id:
                        # matches, happy ending
                        self._set_resolved(cls.m_def)
                        return self.m_proxy_resolved

                    # mismatches, use the usual mechanism
                    return super().m_proxy_resolve()
            except Exception:
                pass

        # Try relative name
        if not self.m_proxy_section or self.m_proxy_resolved:
            return self.m_proxy_resolved

        python_name, definition_id = split_python_definition(self.m_proxy_value)
        current = self.m_proxy_section
        for name in python_name:
            current = self._resolve_name(name, current)

        if current is None:
            raise MetainfoReferenceError(
                f'could not resolve {self.m_proxy_value} from scope {self.m_proxy_section}')
        if not definition_id or current.m_def.definition_id == definition_id:
            # matches, happy ending
            self._set_resolved(current)
            return self.m_proxy_resolved

        # mismatches, use the usual mechanism
        return super().m_proxy_resolve()


class DataType:
    '''
    Allows to define custom data types that can be used in the meta-info.

    The metainfo supports the most types out of the box. These include the python build-in
    primitive types (int, bool, str, float, ...), references to sections, and enums.
    However, in some occasions you need to add custom data types.

    This base class lets you customize various aspects of value treatment. This includes
    type checks and various value transformations. This allows to store values in the
    section differently from how users might set/get them, and it allows to have
    non-serializable values that are transformed on de-/serialization.
    '''

    def set_normalize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        ''' Transforms the given value before it is set and checks its type. '''
        return value

    def get_normalize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        ''' Transforms the given value when it is get. '''
        return value

    def serialize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        ''' Transforms the given value when making the section serializable. '''
        return value

    def deserialize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        ''' Transforms the given value from its serializable form. '''
        return value


class _Dimension(DataType):
    def set_normalize(self, section, quantity_def: Quantity, value):
        if isinstance(value, int):
            return value

        if isinstance(value, str):
            if value.isidentifier():
                return value
            if re.match(MRegEx.index_range, value):
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

        raise TypeError(f'{str(value)} is not a valid dimension')


class _Unit(DataType):
    def set_normalize(self, section, quantity_def: Quantity, value):
        check_unit(value)

        if isinstance(value, str):
            value = units.parse_units(value)

        elif isinstance(value, pint.Quantity):
            value = value.units

        elif not isinstance(value, pint.Unit):
            raise TypeError('Units must be given as str or pint Unit instances.')

        check_dimensionality(quantity_def, value)

        return value

    def serialize(self, section, quantity_def: Quantity, value):
        if quantity_def.flexible_unit:
            return None

        value = value.__str__()
        # The delta prefixes are not serialized: only implicit deltas are
        # allowed currently.
        return reduce(lambda a, b: a.replace(b, ''), _delta_symbols, value)

    def deserialize(self, section, quantity_def: Quantity, value):
        check_unit(value)
        value = units.parse_units(value)
        check_dimensionality(quantity_def, value)
        return value


class _Callable(DataType):
    def serialize(self, section, quantity_def: Quantity, value):
        raise MetainfoError('Callables cannot be serialized')

    def deserialize(self, section, quantity_def: Quantity, value):
        raise MetainfoError('Callables cannot be serialized')


class _QuantityType(DataType):
    ''' Data type for defining the type of metainfo quantity.

    A metainfo quantity type can be one of

    - python build-in primitives: int, float, bool, str
    - numpy dtypes, e.g. np.int32
    - a section definition to define references
    - an MEnum instance to use its values as possible str values
    - a custom datatype, i.e. instance of :class:`DataType`
    - Any
    '''

    def set_normalize(self, section, quantity_def, value):
        if value in MTypes.primitive:
            return value

        if isinstance(value, MEnum):
            return value

        # we normalise all np.dtype to basic np.number types
        if isinstance(value, np.dtype):
            value = value.type

        if value in MTypes.numpy:
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

        raise MetainfoError(f'Type {value} of {quantity_def} is not a valid metainfo quantity type')

    def serialize(self, section, quantity_def, value):
        if value in MTypes.primitive:
            return dict(type_kind='python', type_data=value.__name__)

        if isinstance(value, MEnum):
            result = dict(type_kind='Enum', type_data=list(value))
            if len(value.get_all_descriptions()) > 0:
                result['type_descriptions'] = value.get_all_descriptions()
            return result

        if isinstance(value, np.dtype):
            value = value.type
        # serialise follows the same logic to use basic np.number only
        if value in MTypes.numpy:
            return dict(type_kind='numpy', type_data=str(value.__name__))

        if isinstance(value, Reference):
            if isinstance(value, QuantityReference):
                type_data = value.target_quantity_def.m_path()
                from nomad import config
                if config.process.store_package_definition_in_mongo:
                    type_data += f'@{value.target_quantity_def.definition_id}'
                return dict(type_kind='quantity_reference', type_data=type_data)

            section_root = section.m_root()
            context = cast(MSection, section_root).m_context
            if context is not None:
                try:
                    type_data = context.create_reference(section, quantity_def, value.target_section_def)
                except AssertionError:
                    pass

                if type_data is None:
                    # If no reference could be created from the context, we assume that
                    # the reference is "external" and only available as a Python-based definition.
                    type_data = value.target_section_def.qualified_name()
            else:
                type_data = value.target_section_def.m_path()

            from nomad import config
            if config.process.store_package_definition_in_mongo:
                type_data += f'@{value.target_section_def.definition_id}'

            return value.serialize_type(type_data)

        if isinstance(value, DataType):
            module = value.__class__.__module__
            if module is None or module == str.__class__.__module__:
                type_data = value.__class__.__name__
            else:
                type_data = f'{module}.{value.__class__.__name__}'

            return dict(type_kind='custom', type_data=type_data)

        if value == Any:
            return dict(type_kind='Any')

        raise MetainfoError(f'Type {value} of {quantity_def} is not a valid metainfo quantity type')

    def deserialize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        if isinstance(value, dict):
            if 'type_kind' not in value:
                raise MetainfoError(f'{value} is not a valid quantity type specification.')

            type_kind, type_data = value['type_kind'], value.get('type_data')
            if type_kind == 'python':
                return MTypes.primitive_name[type_data]
            if type_kind == 'Enum':
                return MEnum(*type_data)
            reference = Reference.deserialize_type(type_kind, type_data, section)
            if reference:
                return reference

            if type_kind == 'quantity_reference':
                return QuantityReference(cast(Quantity, MProxy(
                    type_data, m_proxy_section=section, m_proxy_type=Reference(Quantity.m_def))))
            if type_kind == 'Any':
                return Any
            if type_kind == 'custom':
                try:
                    module_name, impl_name = type_data.rsplit('.', 1)
                    module = importlib.import_module(module_name)
                    return getattr(module, impl_name.replace('_', ''))
                except Exception:
                    raise MetainfoError(
                        f'Could not load python implementation of custom datatype {type_data}')
            if type_kind == 'numpy':
                try:
                    return np.dtype(type_data).type
                except Exception:
                    raise MetainfoError(f'{type_data} is not a valid numpy type.')

            raise MetainfoError(f'{type_kind} is not a valid quantity type kind.')

        if value in MTypes.primitive_name:
            return MTypes.primitive_name[value]

        if isinstance(value, str):
            if value.startswith('np.') or value.startswith('numpy.'):
                try:
                    resolved = getattr(np, value.split('.', 1)[1])
                except Exception:
                    raise MetainfoError(f'{value.split(".", 1)[1]} is not a valid numpy type.')
                if resolved:
                    return resolved

            if value in predefined_datatypes:
                return predefined_datatypes[value]

            return Reference(SectionProxy(value, m_proxy_section=section))

        if isinstance(value, list):
            return MEnum(*value)

        return super().deserialize(section, quantity_def, value)


class Reference(DataType):
    '''
    Datatype used for quantities that use other sections as values.

    The target section definition is required to instantiate a Reference type.

    The behavior in this DataType class uses URLs to serialize references. In memory, the
    actual referenced section instance (or respective MProxy instances) are used as values.
    During de-serialization, MProxy instances that auto-resolve on usage, will be used.
    The reference datatype will also accept MProxy instances or URL strings as values
    when set in Python and replace the value with the resolved section instance.

    Subclasses might exchange URLs with a different string serialization, e.g. Python
    qualified names.

    Arguments:
        - section_def: A section definition (Python class or Section instance) that
            determines the type of the referenced sections. Support polymorphism and
            sections that inherit from the given section can also be used as values.
    '''

    def __init__(self, section_def: Optional[SectionDefOrCls]):
        self._target_section_def = to_section_def(section_def)

    @property
    def target_section_def(self):
        return self._target_section_def

    def resolve(self, proxy) -> MSection:
        '''
        Resolve the given proxy. The proxy is guaranteed to have a context and
        will to be not yet resolved.
        '''
        url = ReferenceURL(proxy.m_proxy_value)
        context_section = proxy.m_proxy_section
        if context_section is not None:
            context_section = context_section.m_root()
        if url.archive_url or '@' in url.fragment:
            context = proxy.m_proxy_context
            if context is None:
                context = context_section.m_context
            if not context:
                raise MetainfoReferenceError('Proxy with archive url, but no context to resolve it.')
            if '@' in url.fragment:
                # It's a reference to a section definition
                definition, definition_id = f'{url.archive_url}#{url.fragment}'.split('@')
                return context.resolve_section_definition(definition, definition_id).m_def

            context_section = context.resolve_archive_url(url.archive_url)

        return self.resolve_fragment(context_section, url.fragment)

    def resolve_fragment(self, context_section: MSection, fragment: str) -> MSection:
        return context_section.m_resolve(fragment)

    def serialize_proxy_value(self, proxy):
        return proxy.m_proxy_value

    def set_normalize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        if isinstance(self.target_section_def, MProxy):
            # TODO? This assumes that the type Reference is only used for Quantity.type
            proxy = self.target_section_def
            proxy.m_proxy_section = quantity_def
            proxy.m_proxy_type = Quantity.type.type
            self._target_section_def = proxy.m_proxy_resolve()

        if self.target_section_def.m_follows(Definition.m_def):
            # special case used in metainfo definitions, where we reference metainfo definitions
            # using their Python class. E.g. referencing a section definition using its
            # class instead of the object: Run vs. Run.m_def
            if isinstance(value, type):
                definition = getattr(value, 'm_def', None)
                if definition is not None and definition.m_follows(self.target_section_def):
                    return definition

        if isinstance(value, (str, int, dict)):
            if isinstance(value, str):
                value = self.normalize_reference(section, quantity_def, value)
            return self.deserialize(section, quantity_def, value)

        if isinstance(value, MProxy):
            value.m_proxy_section = section
            value.m_proxy_type = quantity_def.type
            return value

        if not isinstance(value, MSection):
            raise TypeError(f'The value {value} is not a section and can not be used as a reference.')

        if not value.m_follows(self.target_section_def.m_resolved()):  # type: ignore
            raise TypeError(
                f'{value} is not a {self.target_section_def} and therefore an invalid value of {quantity_def}.')

        return value

    @staticmethod
    def normalize_reference(section: MSection, quantity_def: Quantity, value: str):
        context = cast(MSection, section.m_root()).m_context
        return context.normalize_reference(section, value) if context else value

    def serialize_type(self, type_data):
        return dict(type_kind='reference', type_data=type_data)

    def serialize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        context = cast(MSection, section.m_root()).m_context
        return context.create_reference(section, quantity_def, value) if context else value.m_path()

    @classmethod
    def deserialize_type(cls, type_kind, type_data, section):
        if type_kind == 'reference':
            return Reference(SectionProxy(type_data, m_proxy_section=section))
        return None

    def deserialize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        return MProxy(value, m_proxy_section=section, m_proxy_type=quantity_def.type)


# TODO has to deal with URLs, Python qualified names, and Metainfo references
class _SectionReference(Reference):
    def __init__(self):
        super().__init__(None)

    @property
    def target_section_def(self):
        return Section.m_def

    def resolve_fragment(self, context_section: MSection, fragment_with_id: str) -> MSection:
        # First, we try to resolve based on definition names
        if '@' in fragment_with_id:
            fragment, definition_id = fragment_with_id.split('@')
        else:
            definition_id = None
            fragment = fragment_with_id

        definitions = None
        if isinstance(getattr(context_section, 'definitions', None), Definition):
            definitions = getattr(context_section, 'definitions')

        if isinstance(context_section, Definition):
            definitions = context_section

        if definitions:
            split_fragment = fragment.lstrip('/').split('/', 1)
            if len(split_fragment) == 2:
                first_segment, remaining_fragment = split_fragment
            else:
                first_segment, remaining_fragment = split_fragment[0], None

            resolved: Optional[MSection] = None
            for content in definitions.m_contents():
                if isinstance(content, Definition) and content.name == first_segment:
                    if remaining_fragment:
                        resolved = self.resolve_fragment(content, remaining_fragment)
                    else:
                        return _check_definition_id(definition_id, content)

            if resolved:
                return _check_definition_id(definition_id, resolved)

        # Resolve regularly as a fallback
        return super().resolve_fragment(context_section, fragment_with_id)

    def set_normalize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        if isinstance(value, str) and MRegEx.python_definition.match(value):
            return SectionProxy(value, m_proxy_section=section, m_proxy_type=quantity_def.type)

        return super().set_normalize(section, quantity_def, value)

    def serialize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        from nomad import config

        def _append_definition_id(section_name) -> str:
            if config.process.store_package_definition_in_mongo:
                return f'{section_name}@{value.definition_id}'
            return section_name

        # First we try to use a potentially available Python name to serialize
        if isinstance(value, Section):
            pkg: MSection = value.m_root()
            if isinstance(pkg, Package) and pkg.name not in [None, '*']:
                return f'{pkg.name}.{value.name}'

        # Default back to URL
        return _append_definition_id(super().serialize(section, quantity_def, value))

    def deserialize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        proxy_type = quantity_def.type if quantity_def else SectionReference
        if isinstance(value, str) and MRegEx.python_definition.match(value):
            # First assume it's a python name and try to resolve it.
            if '.' in value:
                python_name, definition_id = split_python_definition(value)
                package_name = '.'.join(python_name[:-1])
                section_name = python_name[-1]

                try:
                    module = importlib.import_module(package_name)
                    cls = getattr(module, section_name)
                    if cls:
                        m_def = getattr(cls, 'm_def')
                        if m_def and (definition_id is None or m_def.definition_id == definition_id):
                            # id matches, happy ending
                            return m_def
                except ModuleNotFoundError:
                    pass

            # If it's not a python name or definition id mismatches
            # we assume its referring to a local metainfo definition.
            return SectionProxy(value, m_proxy_section=section, m_proxy_type=proxy_type)

        # Default back to value being a URL
        return MProxy(value, m_proxy_section=section, m_proxy_type=proxy_type)


SectionReference = _SectionReference()


class QuantityReference(Reference):
    ''' Datatype used for reference quantities that reference other quantities. '''

    def __init__(self, quantity_def: Quantity):
        super().__init__(None)
        self.target_quantity_def = quantity_def

    @property
    def target_section_def(self):
        return cast(Section, self.target_quantity_def.m_parent)

    def serialize_proxy_value(self, proxy):
        return f'{proxy.m_proxy_value}/{self.target_quantity_def.name}'

    def set_normalize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        if not value.m_is_set(self.target_quantity_def):
            return _unset_value

        return super().set_normalize(section, quantity_def, value)

    def get_normalize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        section = super().get_normalize(section, quantity_def, value)
        return getattr(section, self.target_quantity_def.name)

    def serialize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        target_path = super().serialize(section, quantity_def, value)
        return f'{target_path}/{self.target_quantity_def.name}'

    def deserialize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        target_path = value.rsplit('/', 1)[0]
        return MProxy(target_path, m_proxy_section=section, m_proxy_type=quantity_def.type)


class _File(DataType):
    def set_normalize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        if not isinstance(value, str):
            raise TypeError('Files need to be given as URL strings')

        root_section: MSection = section.m_root()
        context = cast(MSection, root_section).m_context
        if context:
            return context.normalize_reference(root_section, value)

        return value


class _URL(DataType):
    def set_normalize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        return validate_url(value)

    def serialize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        return validate_url(value)

    def deserialize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        return validate_url(value)


class _Datetime(DataType):
    def set_normalize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        return normalize_datetime(value)

    def serialize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        return None if value is None else value.isoformat()

    def deserialize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        return normalize_datetime(value)


class _JSON(DataType):
    pass


class _Capitalized(DataType):
    def set_normalize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        if value is not None and len(value) >= 1:
            return value[0].capitalize() + value[1:]

        return value


class _Bytes(DataType):
    def serialize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        return base64.b64encode(value).decode('ascii')

    def deserialize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        return base64.b64decode(value)


class _HDF5Reference(DataType):
    pass


Dimension = _Dimension()
Unit = _Unit()
QuantityType = _QuantityType()
Callable = _Callable()
Datetime = _Datetime()
JSON = _JSON()
Capitalized = _Capitalized()
Bytes = _Bytes()
File = _File()
URL = _URL()
HDF5Reference = _HDF5Reference()

predefined_datatypes = {
    'Dimension': Dimension, 'Unit': Unit, 'Datetime': Datetime,
    'JSON': JSON, 'Capitalized': Capitalized, 'bytes': Bytes, 'File': File, 'URL': URL,
    'HDF5Reference': HDF5Reference}


# Metainfo data storage and reflection interface

class MObjectMeta(type):

    def __new__(self, cls_name, bases, dct):
        do_init = dct.get('do_init', None)
        if do_init is not None:
            del dct['do_init']
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
- the definition of a subsection
- or the section definition Python class
'''


def constraint(warning):
    ''' A decorator for methods implementing constraints. '''
    f = None
    if not isinstance(warning, bool):
        f = warning
        warning = False

    def decorator(_f):
        setattr(_f, 'm_constraint', True)
        setattr(_f, 'm_warning', warning)
        return _f

    return decorator if f is None else decorator(f)


class Context:
    '''
    The root of a metainfo section hierarchy can have a Context. Contexts allow to customize
    the resolution of references based on how and in what context a metainfo-based
    archive (or otherwise top-level section is used). This allows to logically combine
    multiple hierarchies (e.g. archives) with references.
    '''

    def warning(self, event, **kwargs):
        '''
        Used to log (or otherwise handle) warning that are issued, e.g. while serialization,
        reference resolution, etc.
        '''
        pass

    def create_reference(
        self, section: MSection, quantity_def: Quantity, value: MSection,
        global_reference: bool = False
    ) -> str:
        '''
        Returns a reference for the given target section (value) based on the given context.
        Allows subclasses to build references across resources, if necessary.

        Arguments:
            global_reference: A boolean flag that forces references with upload_ids.
                Should be used if the reference needs to be used outside the context
                of the own upload.

        Raises: MetainfoReferenceError
        '''
        return value.m_path()

    def normalize_reference(self, source: MSection, url: str) -> str:
        '''
        Rewrites the url into a normalized form. E.g., it replaces `..` with absolute paths,
        or replaces mainfiles with ids, etc.

        Arguments:
            source: The source section or root section of the source of the reference.
            url: The reference to normalize.

        Raises: MetainfoReferenceError
        '''
        return url

    def resolve_archive_url(self, url: str) -> MSection:
        '''
        Resolves the archive part of the given URL and returns the root section of the
        referenced archive.

        Raises: MetainfoReferenceError
        '''
        raise NotImplementedError()

    def resolve_archive(self, *args, **kwargs):
        return self.resolve_archive_url(*args, **kwargs)

    def cache_archive(self, url: str, archive):
        raise NotImplementedError()

    def retrieve_package_by_section_definition_id(self, definition_reference: str, definition_id: str) -> dict:
        raise NotImplementedError()

    def resolve_section_definition(self, definition_reference: str, definition_id: str) -> Type[MSectionBound]:
        pkg_definition = self.retrieve_package_by_section_definition_id(definition_reference, definition_id)

        entry_id_based_name = pkg_definition['entry_id_based_name']
        del pkg_definition['entry_id_based_name']
        pkg = Package.m_from_dict(pkg_definition)
        if entry_id_based_name != '*':
            pkg.entry_id_based_name = entry_id_based_name

        pkg.m_context = self
        pkg.init_metainfo()

        for section in pkg.section_definitions:
            if section.definition_id == definition_id:
                return section.section_cls
        return None


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
    without a priori knowledge of the `section definition`.

    .. automethod:: m_set
    .. automethod:: m_get
    .. automethod:: m_add_values
    .. automethod:: m_get_sub_section
    .. automethod:: m_get_sub_sections
    .. automethod:: m_create
    .. automethod:: m_add_sub_section
    .. automethod:: m_remove_sub_section

    There are some specific attributes for section instances that are subsections of
    another section. While subsections are directly accessible from the containing
    section by using the Python property that represents the subsection (e.g.
    `run.section_system`), there is also a way to navigate from the subsection to
    the containing section (`parent section`) using these Python properties:

    Attributes:
        m_parent:
            If this section is a subsection, this references the parent section instance.

        m_parent_sub_section:
            If this section is a subsection, this is the :class:`SubSection` that defines
            this relationship.

        m_parent_index:
            For repeatable sections, parent keep a list of subsections. This is the index
            of this section in the respective parent subsection list.

        m_context: The :class:`MContext` that manages this (root-)section.

    Often some general tasks have to be performed on a whole tree of sections without
    knowing about the definitions in advance. The following methods allow to access
    subsections reflectively.

    .. automethod:: m_traverse
    .. automethod:: m_all_contents
    .. automethod:: m_contents
    .. automethod:: m_xpath

    Each section and all its quantities and contents can be transformed into a general
    JSON-serializable Python dictionary. Similarly, a section can be instantiated from
    such a Python dictionary. This allows to save and load sections to JSON-files or
    by other compatible means (e.g. document databases, binary JSON flavours).

    .. automethod:: m_to_dict
    .. automethod:: m_from_dict
    .. automethod:: m_update_from_dict
    .. automethod:: m_to_json
    '''

    m_def: Section = None

    def __init__(
            self, m_def: Section = None, m_context: Context = None, **kwargs):

        self.m_def: Section = m_def
        self.m_parent: MSection = None
        self.m_parent_sub_section: SubSection = None
        self.m_parent_index = -1
        self.m_context = m_context
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

        elif not is_bootstrapping:
            MetainfoError('Section has no m_def.')

        # get annotations from kwargs
        self.m_annotations: Dict[str, Any] = kwargs.get('m_annotations', {})
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

        self.m_parse_annotations()

        # set remaining kwargs
        if is_bootstrapping:
            self.__dict__.update(**other_kwargs)  # type: ignore
        else:
            self.m_update(**other_kwargs)

    @classmethod
    def __init_cls__(cls):
        # ensure that the m_def is defined
        m_def = cls.__dict__.get('m_def')  # do not accidentally get the m_def from a potential base section
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
                        'Section defining classes must have MSection or a descendant of base classes.')
                base_sections.append(base_section)

        if len(base_sections) > 0:
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
                    description = re.sub(r'\(https?://[^)]*\)', lambda m: re.sub(r'\n', '', m.group(0)), description)
                    attr.description = description
                    attr.__doc__ = attr.description

                if isinstance(attr, Quantity):
                    m_def.m_add_sub_section(Section.quantities, attr)
                elif isinstance(attr, SubSection):
                    m_def.m_add_sub_section(Section.sub_sections, attr)
                else:
                    raise NotImplementedError('Unknown property kind.')
            elif isinstance(attr, Attribute):
                attr.name = name
                m_def.m_add_sub_section(Section.attributes, attr)

            if inspect.isclass(attr):
                inner_section_def = getattr(attr, 'm_def', None)
                if isinstance(inner_section_def, Section):
                    m_def.m_add_sub_section(Section.inner_section_definitions, inner_section_def)

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

        if len(constraints) > 0:
            m_def.constraints = list(sorted(constraints))
        if len(event_handlers) > 0:
            m_def.event_handlers = list(sorted(event_handlers))

        # add section cls' section to the module's package
        module_name = cls.__module__
        pkg = Package.from_module(module_name)
        pkg.m_add_sub_section(Package.section_definitions, cls.m_def)

        # apply_google_docstrings
        # Parses the Google doc string of the given class and properly updates the
        # definition descriptions.

        # This allows to document quantities and subsections with 'Args:' in the section
        # class. It will remove the 'Args' section from the section definition and will
        # set the respective pieces to the quantity and subsection descriptions.
        docstring = cls.__doc__
        if docstring is not None:
            parsed_docstring = docstring_parser.parse(docstring)
            short = parsed_docstring.short_description
            dsc = parsed_docstring.long_description

            if short and dsc:
                description = f'{short.strip()} {dsc.strip()}'
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
        if self.m_def is None:
            return super().__setattr__(name, value)

        alias_pool = self.m_def.all_aliases

        if alias_pool is not None and name in alias_pool:
            name = alias_pool[name].name
        elif self.m_def.has_variable_names and not MRegEx.reserved_name.match(name):
            resolved_name = resolve_variadic_name(self.m_def.all_properties, name)
            if resolved_name:
                name = resolved_name.name

        return super().__setattr__(name, value)

    def __getattr__(self, name):
        # The existence of __getattr__ will make mypy and pylint ignore 'missing' dynamic
        # attributes and functions and wrong types of those.
        # Ideally we have a plugin for both that add the correct type info

        if name in self.m_def.all_aliases:
            return getattr(self, self.m_def.all_aliases[name].name)

        if self.m_def.has_variable_names:
            m_definition: Definition = resolve_variadic_name(self.m_def.all_properties, name)
            if m_definition:
                if not isinstance(m_definition, Quantity) or not m_definition.use_full_storage:
                    return getattr(self, m_definition.name)

                m_storage: dict = self.__dict__.get(m_definition.name, None)
                if m_storage is None:
                    return None

                m_quantity = m_storage.get(name, None)
                if m_quantity is None:
                    return None

                if m_quantity.value is not None:
                    if m_quantity.unit is not None:
                        return units.Quantity(m_quantity.value, m_quantity.unit)

                    return m_quantity.value

        if name.startswith('a_'):
            annotation_name = name[2:]
            if annotation_name not in self.m_annotations:
                raise AttributeError(name)

            return self.m_get_annotations(annotation_name)

        raise AttributeError(name)

    def __set_normalize(self, quantity_def: Quantity, value: Any) -> Any:
        target_type = quantity_def.type

        if isinstance(target_type, DataType):
            return target_type.set_normalize(self, quantity_def, value)

        if isinstance(target_type, Section):
            if isinstance(value, MProxy):
                return value

            if not isinstance(value, MSection):
                raise TypeError(
                    f'The value {value} for reference quantity {quantity_def} is not a section instance.')

            if not value.m_follows(target_type):
                raise TypeError(
                    f'The value {value} for quantity {quantity_def} does not follow {target_type}')

            return value

        if isinstance(target_type, MEnum):
            if value not in cast(MEnum, target_type).get_all_values():
                raise TypeError(f'The value {value} is not an enum value for {quantity_def}.')
            return value

        if target_type == Any:
            return value

        if target_type == str and type(value) == np.str_:
            return str(value)

        if target_type == bool and type(value) == np.bool_:
            return bool(value)

        if target_type == int and type(value) == np.float_:
            return int(value)

        if target_type in MTypes.complex:
            return normalize_complex(value, target_type, quantity_def.unit)

        if type(value) != target_type:
            if target_type in MTypes.primitive:
                try:
                    return MTypes.primitive[target_type](value)  # type: ignore
                except ValueError as e:
                    raise TypeError(e)

            if value is not None:
                raise TypeError(f'The value {value} for {quantity_def} is not of type {target_type}.')

        return value

    def m_parse_annotations(self):
        for annotation_name, annotation in self.m_annotations.items():
            annotation_model = AnnotationModel.m_registry.get(annotation_name)
            if not annotation_model:
                continue

            def to_model(annotation_model, annotation_data):
                if annotation_data is None:
                    return None

                try:
                    if isinstance(annotation_data, AnnotationModel):
                        annotation = annotation_data
                    else:
                        annotation = parse_obj_as(annotation_model, annotation_data)

                    if isinstance(self, Definition):
                        annotation.m_definition = self

                except ValidationError as e:
                    annotation = AnnotationModel(m_error=str(e))

                return annotation

            if isinstance(annotation, list):
                for index, item in enumerate(annotation):
                    annotation[index] = to_model(annotation_model, item)
            else:
                annotation = to_model(annotation_model, annotation)

            if annotation:
                self.m_annotations[annotation_name] = annotation

    def m_set(self, quantity_def: Quantity, value: Any) -> None:
        ''' Set the given value for the given quantity. '''
        self.m_mod_count += 1

        if quantity_def.derived is not None:
            raise MetainfoError(f'The quantity {quantity_def} is derived and cannot be set.')

        item_name: str = quantity_def.name

        if value is None:
            # This implements the implicit "unset" semantics of assigned None as a value
            to_remove = self.__dict__.pop(item_name, None)
            # if full storage is used, also need to clear quantities created for convenient access
            if quantity_def.use_full_storage and to_remove:
                # self.__dict__[full_name] is guaranteed to be a 'dict[str, MQuantity]'
                for key in to_remove.keys():
                    self.__dict__.pop(key, None)
            return

        if not quantity_def.use_full_storage:
            # handles the non-repeating and no attribute case, store the value directly under the name
            if quantity_def.type in MTypes.numpy or isinstance(quantity_def.type, pd.DataFrame):
                value = to_numpy(quantity_def.type, quantity_def.shape, quantity_def.unit, quantity_def, value)
            else:
                dimensions = len(quantity_def.shape)
                if dimensions == 0:
                    value = self.__set_normalize(quantity_def, value)
                    if value == _unset_value:
                        return

                elif dimensions == 1:
                    if type(value) == str or not isinstance(value, IterableABC):
                        raise TypeError(
                            f'The shape of {quantity_def} requires an iterable value, but {value} is not iterable.')

                    if quantity_def.type == complex:
                        value = normalize_complex(value, complex, quantity_def.unit)
                    else:
                        value = [v for v in list(
                            self.__set_normalize(quantity_def, item) for item in value) if v != _unset_value]

                else:
                    raise MetainfoError(
                        f'Only numpy arrays and dtypes can be used for higher dimensional quantities: {quantity_def}')

            self.__dict__[item_name] = value
        else:
            # it is a repeating quantity w/o attributes
            # the actual value/name/unit would be wrapped into 'MQuantity'
            # check if there is an existing item
            m_quantity: MQuantity
            m_attribute: dict = {}
            if isinstance(value, MQuantity):
                m_quantity = value
                if not quantity_def.variable:
                    if not m_quantity.name:
                        m_quantity.name = item_name
                    elif m_quantity.name != item_name:
                        raise MetainfoError(f"The name of {value} must match definition name {item_name}")
                else:
                    if not m_quantity.name:
                        raise MetainfoError(f"The name must be provided for variadic quantity {item_name}")

                # swap to add attributes via the setter to allow validation
                m_attribute = m_quantity.attributes
                m_quantity.attributes = {}
            elif not quantity_def.variable:
                try:
                    m_quantity = self.__dict__[item_name][item_name]
                    if isinstance(value, pint.Quantity):
                        m_quantity.value = value.m
                        m_quantity.unit = value.u
                    else:
                        m_quantity.value = value
                except KeyError:
                    m_quantity = MQuantity(item_name, value)
            else:
                raise MetainfoError("Variadic quantities only accept raw values wrapped in 'MQuantity'")

            if not validate_shape(self, quantity_def, m_quantity.value):
                raise MetainfoError(f"The shape of {m_quantity} does not match {quantity_def.shape}")

            if quantity_def.unit is None:
                # no prescribed unit, need to check dimensionality, no need to convert
                check_dimensionality(quantity_def, m_quantity.unit)
            else:
                try:
                    m_quantity.value = convert_to(m_quantity.value, m_quantity.unit, quantity_def.unit)
                except (ValueError, TypeError):
                    raise MetainfoError(f'Could not convert {m_quantity.unit} to {quantity_def.unit}')
                m_quantity.unit = quantity_def.unit

            if quantity_def.type in MTypes.numpy or isinstance(quantity_def.type, pd.DataFrame):
                m_quantity.value = to_numpy(
                    quantity_def.type, quantity_def.shape, quantity_def.unit, quantity_def, m_quantity.value)
            else:
                dimensions = len(quantity_def.shape)
                if dimensions == 0:
                    m_quantity.value = self.__set_normalize(quantity_def, m_quantity.value)
                    if m_quantity.value == _unset_value:
                        return

                elif dimensions == 1:
                    if type(m_quantity.value) == str or not isinstance(m_quantity.value, IterableABC):
                        raise TypeError(
                            f'The shape of {quantity_def} requires an iterable value, '
                            f'but {m_quantity.value} is not iterable.')

                    if quantity_def.type == complex:
                        m_quantity.value = normalize_complex(m_quantity.value, complex, quantity_def.unit)
                    else:
                        m_quantity.value = [v for v in list(
                            self.__set_normalize(quantity_def, item) for item in m_quantity.value) if v != _unset_value]

                else:
                    raise MetainfoError(
                        f'Only numpy arrays and dtypes can be used for higher dimensional quantities: {quantity_def}')

            # store under variable name with suffix
            if item_name in self.__dict__:
                self.__dict__[item_name][m_quantity.name] = m_quantity
            else:
                self.__dict__[item_name] = {m_quantity.name: m_quantity}

            for k, v in m_attribute.items():
                self.m_set_quantity_attribute(m_quantity.name, k, v)

        for handler in self.m_def.event_handlers:
            if handler.__name__.startswith('on_set'):
                handler(self, quantity_def, value)

    def m_get(self, quantity_def: Quantity, full: bool = False) -> Any:
        ''' Retrieve the given value for the given quantity. '''
        if not full:
            return quantity_def.__get__(self, Quantity)

        return self.__dict__[quantity_def.name]

    def m_get_quantity_definition(self, quantity_name: str, hint: Optional[str] = None):
        '''
        Get the definition of the quantity with the target name.

        An optional hint string can be provided. The hint should be the name of one of attributes
        defined in the target quantity.
        '''
        return resolve_variadic_name(self.m_def.all_quantities, quantity_name, hint)

    def m_is_set(self, quantity_def: Quantity) -> bool:
        ''' True if the given quantity is set. '''
        if quantity_def.derived is not None:
            return True

        return quantity_def.name in self.__dict__

    def m_add_values(self, quantity_def: Quantity, values: Any, offset: int) -> None:
        ''' Add (partial) values for the given quantity of higher dimensionality. '''
        # TODO
        raise NotImplementedError()

    def _get_sub_sections(self, sub_section_def: SubSection):
        sub_section_lst = self.__dict__.get(sub_section_def.name)
        if sub_section_lst is None:
            sub_section_lst = MSubSectionList(self, sub_section_def)
            self.__dict__[sub_section_def.name] = sub_section_lst
        return sub_section_lst

    def _on_add_sub_section(
            self, sub_section_def: SubSection, sub_section: MSection, parent_index: int) -> None:
        self.m_mod_count += 1

        sub_section.m_parent = self
        sub_section.m_parent_sub_section = sub_section_def
        sub_section.m_parent_index = parent_index

        for handler in self.m_def.event_handlers:
            if handler.__name__.startswith('on_add_sub_section'):
                handler(self, sub_section_def, sub_section)

    def _on_remove_sub_section(self, sub_section_def: SubSection, sub_section: MSection) -> None:
        self.m_mod_count += 1

        sub_section.m_parent = None
        sub_section.m_parent_index = -1

    def m_add_sub_section(self, sub_section_def: SubSection, sub_section: MSection, index: int = -1) -> None:
        ''' Adds the given section instance as a subsection of the given subsection definition. '''

        sub_section_name = sub_section_def.name
        if sub_section_def.repeats:
            if index != -1:
                raise NotImplementedError('You can only append subsections.')

            sub_section_lst = self._get_sub_sections(sub_section_def)
            sub_section_lst.append(sub_section)
            if sub_section_lst.__class__ != MSubSectionList:
                self._on_add_sub_section(sub_section_def, sub_section, len(sub_section_lst) - 1)
        else:
            if (old_sub_section := self.__dict__.get(sub_section_name, None)) is sub_section:
                return
            self.__dict__[sub_section_name] = sub_section
            if sub_section is not None:
                self._on_add_sub_section(sub_section_def, sub_section, -1)
            if old_sub_section is not None:
                self._on_remove_sub_section(sub_section_def, old_sub_section)

    def m_remove_sub_section(self, sub_section_def: SubSection, index: int) -> None:
        ''' Removes the exiting section for a non-repeatable subsection '''
        self.m_mod_count += 1

        if sub_section_def.name not in self.__dict__:
            return

        if sub_section_def.repeats:
            sub_section = self.__dict__[sub_section_def.name][index]
            del self.__dict__[sub_section_def.name][index]
        else:
            sub_section = self.__dict__[sub_section_def.name]
            del self.__dict__[sub_section_def.name]

        self._on_remove_sub_section(sub_section_def, sub_section)

    def m_get_sub_section(self, sub_section_def: SubSection, index: Any) -> Optional[MSection]:
        ''' Retrieves a single subsection of the given subsection definition. '''
        if not sub_section_def.repeats:
            return self.__dict__.get(sub_section_def.name, None)

        if isinstance(index, int):
            return self.__dict__[sub_section_def.name][index]

        if isinstance(index, str):
            try:
                sub_sections: List[MSection] = [
                    section for section in self.__dict__[sub_section_def.name] if index == section.name]
                if len(sub_sections) > 1:
                    raise MetainfoReferenceError(f'multiple sections with this section id were found.')
                if len(sub_sections) == 1:
                    return sub_sections[0]
            except KeyError:
                raise MetainfoReferenceError(f'{index} is not a valid subsection.')

        return None

    def m_get_sub_sections(self, sub_section_def: SubSection) -> List[MSection]:
        ''' Retrieves  all subsections of the given subsection definition. '''
        if sub_section_def.repeats:
            return self._get_sub_sections(sub_section_def)

        try:
            return [self.__dict__[sub_section_def.name]]
        except KeyError:
            return []

    def m_sub_section_count(self, sub_section_def: SubSection) -> int:
        ''' Returns the number of subsections for the given subsection definition. '''
        try:
            value = self.__dict__[sub_section_def.name]
            return len(value) if sub_section_def.repeats else 1
        except KeyError:
            return 0

    def m_set_section_attribute(self, name: str, value: Any) -> None:
        '''
        Set attribute for the current section.
        '''
        self.__set_attribute(None, name, value)

    def m_set_quantity_attribute(self, quantity_def: Union[str, Quantity], name: str, value: Any, quantity: Quantity = None) -> None:
        '''
        Set attribute for the given quantity.
        '''
        self.__set_attribute(quantity_def, name, value, quantity=quantity)

    def __set_attribute(self, tgt_property: Union[Optional[str], Definition], attr_name: str, attr_value: Any, quantity: Quantity = None):
        '''
        Set attribute for current section for a quantity of the current section.

        For attributes of the current section, use None as the target property.
        For attributes of a quantity, use the quantity name/definition as the target property.

        Both the quantity name and the attribute name can be variadic.

        Arguments:
            tgt_property: The name or definition of the quantity to set the attribute for, can be None.
            attr_name: The name of the attribute to set.
            attr_value: The value of the attribute to set.
        '''
        tgt_name: Optional[str] = tgt_property.name if isinstance(tgt_property, Definition) else tgt_property

        tgt_def, tgt_attr = retrieve_attribute(self.m_def, tgt_property, attr_name)

        if tgt_attr.type in MTypes.numpy:
            attr_value = to_numpy(tgt_attr.type, [], None, tgt_attr, attr_value)
        else:
            dimension = len(tgt_attr.shape)
            if dimension == 0:
                attr_value = self.__set_normalize(tgt_attr, attr_value)
            elif dimension == 1:
                if type(attr_value) == str or not isinstance(attr_value, IterableABC):
                    raise TypeError(f'The shape requires an iterable value, but {attr_value} is not.')

                if tgt_attr.type == complex:
                    attr_value = normalize_complex(attr_value, complex, None)
                else:
                    attr_value = list(self.__set_normalize(tgt_attr, item) for item in attr_value)
            else:
                raise MetainfoError(f'Only numpy arrays can be used for higher dimensional quantities: {tgt_attr}.')

        if not validate_shape(self, tgt_attr, attr_value):
            raise MetainfoError(f'Invalid shape for attribute: {tgt_attr}.')

        if isinstance(tgt_def, Quantity) and tgt_def.use_full_storage:
            m_storage: Optional[dict] = self.__dict__.get(tgt_def.name, None)
            m_quantity: Optional[MQuantity] = m_storage.get(tgt_property if quantity is None else quantity.name, None) if m_storage else None
            if m_quantity is None:
                m_quantity = MQuantity(tgt_name, None)
                self.m_set(tgt_def, m_quantity)
            m_quantity.m_set_attribute(attr_name, attr_value)
        elif tgt_property is None:
            # indicating that the attribute is for the current section
            if 'm_attributes' not in self.__dict__:
                self.__dict__['m_attributes'] = {}
            self.__dict__['m_attributes'][attr_name] = attr_value

    def m_get_section_attribute(self, name: str) -> Any:
        '''
        Get attribute for the current section.
        '''
        return self.__get_attribute(None, name)

    def m_get_quantity_attribute(self, quantity_def: str, name: str) -> Any:
        '''
        Get attribute for the given quantity.
        '''
        return self.__get_attribute(quantity_def, name)

    def __get_attribute(self, tgt_property: Optional[str], attr_name: str):
        '''
        Get the attribute of a quantity of the current section, or of the current section itself.
        '''
        tgt_def: Definition = tgt_property if tgt_property is None else retrieve_attribute(
            self.m_def, tgt_property, attr_name)[0]

        # section attributes
        if tgt_def is None:
            if 'm_attributes' not in self.__dict__:
                return None

            return self.__dict__['m_attributes'].get(attr_name, None)

        # quantity attributes
        m_storage: Optional[dict] = self.__dict__.get(tgt_def.name, None)
        if m_storage is None:
            return None

        m_quantity: Optional[MQuantity] = m_storage.get(tgt_property, None)
        if m_quantity is None:
            return None

        return m_quantity.attributes.get(attr_name, None)

    def m_create(
            self, section_cls: Type[MSectionBound], sub_section_def: SubSection = None,
            **kwargs) -> MSectionBound:
        ''' Creates a section instance and adds it to this section provided there is a
        corresponding subsection.

        Args:
            section_cls: The section class for the subsection to create
            sub_section_def: If there are multiple subsections for the given class,
                this must be used to explicitly state the subsection definition.
        '''

        section_def = section_cls.m_def
        sub_section_defs = self.m_def.all_sub_sections_by_section.get(section_def, [])
        n_sub_section_defs = len(sub_section_defs)
        if n_sub_section_defs == 0:
            raise TypeError(f'There is no subsection to hold a {section_def} in {self.m_def}.')

        if n_sub_section_defs > 1 and sub_section_def is None:
            raise MetainfoError(
                f'There are multiple subsection to hold a {section_def} in {self.m_def}, '
                f'but no subsection was explicitly given.')

        if sub_section_def is not None and sub_section_def not in sub_section_defs:
            raise MetainfoError(
                f'The given subsection class {section_cls} does not '
                f'match the given subsection definition {sub_section_def}.')

        if sub_section_def is None:
            sub_section_def = sub_section_defs[0]

        sub_section = section_cls(**kwargs)
        self.m_add_sub_section(sub_section_def, sub_section)

        return cast(MSectionBound, sub_section)

    def m_update(self, m_ignore_additional_keys: bool = False, **kwargs):
        ''' Updates all quantities and subsections with the given arguments. '''
        self.m_mod_count += 1

        for name, value in kwargs.items():
            prop = self.m_def.all_aliases.get(name, None)
            if prop is None:
                if m_ignore_additional_keys:
                    continue
                raise KeyError(f'{name} is not an attribute of this section {self}')

            if isinstance(prop, SubSection):
                if prop.repeats:
                    if isinstance(value, List):
                        for item in value:
                            self.m_add_sub_section(prop, item)
                    else:
                        raise TypeError(f'Subsection {prop.name} repeats, but no list was given')
                else:
                    self.m_add_sub_section(prop, value)

            else:
                self.m_set(prop, value)

    def m_as(self, section_cls: Type[MSectionBound]) -> MSectionBound:
        ''' 'Casts' this section to the given extending sections. '''
        return cast(MSectionBound, self)

    def m_follows(self, definition: Section) -> bool:
        ''' Determines if this section's definition is or is derived from the given definition. '''
        if not isinstance(definition, Section):
            raise TypeError(f'{definition} is not an instance of class Section')
        return self.m_def == definition or definition in self.m_def.all_base_sections

    def m_to_dict(
            self, with_meta: bool = False, with_root_def: bool = False,
            with_out_meta: bool = False, with_def_id: bool = False,
            include_defaults: bool = False,
            include_derived: bool = False,
            resolve_references: bool = False,
            categories: List[Union[Category, Type['MCategory']]] = None,
            include: TypingCallable[[Definition, MSection], bool] = None,
            exclude: TypingCallable[[Definition, MSection], bool] = None,
            transform: TypingCallable[[Definition, MSection, Any, str], Any] = None) -> Dict[str, Any]:
        '''
        Returns the data of this section as a (json serializable) dictionary.

        With its default configuration, it is the opposite to :func:`MSection.m_from_dict`.

        There are a lot of ways to customize the behavior, e.g. to generate JSON for
        databases, search engines, etc.

        Arguments:
            with_meta: Include information about the section definition, the sections
                position in its parent, and annotations. For Definition instances this
                information will be included regardless; the section definition will
                always be included if the subsection definition references a base section
                and the concrete subsection is derived from this base section.
            with_out_meta: Exclude information `with_meta` information, even from
                Definition instances.
            with_root_def: Include the m_def for the top-level section. This allows to
                identify the used "schema" based on the root section definition.
            with_def_id: Include the definition id for the top-level section. This
                allows detection different versions of section definition used.
            include_defaults: Include default values of unset quantities.
            include_derived: Include values of derived quantities.
            resolve_references:
                Treat references as the sections and values they represent. References
                must not create circles; there is no check and danger of endless looping.
            categories: A list of category classes or category definitions that is used
                to filter the included quantities and subsections. Only applied to
                properties of this section, not on subsections. Is overwritten
                by partial.
            include: A function that determines if a property (quantity or subsection) will
                be included in the results. It takes the property definition and the current
                section as arguments. The function returns true for including and false for
                excluding the property. Include is applied recursively on subsections.
                Overrides categories.
            exclude: A function that determines if a property (quantity or subsection) will
                be excluded from the results. It takes the property definition and the current
                section as arguments. The function returns true for excluding and false for
                including the property. Exclude is applied recursively on subsections.
                Overrides categories.
            transform: A function that determines serialized quantity values.
                It takes the quantity definition, current section, the default
                serialized value and the metainfo path with respect to the
                document root as arguments. Depending on where this is used, you
                might have to ensure that the result is JSON-serializable.  By
                default, values are serialized to JSON according to the quantity
                type.
        '''
        if isinstance(self, Definition) and not with_out_meta:
            with_meta = True

        kwargs: Dict[str, Any] = dict(
            with_meta=with_meta, with_out_meta=with_out_meta, with_def_id=with_def_id,
            include_defaults=include_defaults,
            include_derived=include_derived,
            resolve_references=resolve_references,
            exclude=exclude,
            transform=transform)

        assert not (include is not None and exclude is not None), 'You can only include or exclude, not both.'

        if include is not None:
            def exclude(*args, **kwargs):  # pylint: disable=function-redefined
                return not include(*args, **kwargs)

            kwargs['exclude'] = exclude

        elif exclude is None:
            if categories is None:
                def exclude(prop, section):  # pylint: disable=function-redefined
                    return False

                kwargs['exclude'] = exclude
            else:
                category_defs: List[Category] = []
                for category in categories:
                    if issubclass(category, MCategory):  # type: ignore
                        category_defs.append(category.m_def)  # type: ignore
                    elif isinstance(category, Category):
                        category_defs.append(category)
                    else:
                        raise TypeError(f'{category} is not a category')

                def exclude(prop, section):  # pylint: disable=function-redefined
                    return not any(prop in v.get_all_definitions() for v in category_defs)

        def serialize_quantity(quantity, is_set, is_derived, path, target_value=None):
            quantity_type = quantity.type

            if resolve_references and isinstance(quantity_type, QuantityReference):
                quantity_type = quantity_type.target_quantity_def.type

            serialize: TypingCallable[[Any], Any]

            # define serialization functions for all valid data types
            is_reference = False
            if isinstance(quantity_type, Reference):
                is_reference = True

                def serialize_reference(value, path_override):
                    if resolve_references:
                        assert not isinstance(quantity_type, QuantityReference)
                        value = value.m_resolved()
                        ref_kwargs = dict(kwargs)
                        if kwargs["transform"]:
                            ref_kwargs["transform"] = lambda q, s, v, p: kwargs["transform"](q, s, v, path_override)
                        return value.m_to_dict(**ref_kwargs)

                    if isinstance(value, MProxy):
                        if value.m_proxy_resolved is not None:
                            return quantity_type.serialize(self, quantity, value)

                        return quantity_type.serialize_proxy_value(value)

                    return quantity_type.serialize(self, quantity, value)

                serialize = serialize_reference

            elif isinstance(quantity_type, DataType):

                def serialize_data_type(value):
                    return quantity_type.serialize(self, quantity, value)

                serialize = serialize_data_type

            elif quantity_type in MTypes.complex:

                serialize = serialize_complex

            elif quantity_type in MTypes.primitive:

                serialize = MTypes.primitive[quantity_type]

            elif quantity_type in MTypes.numpy:

                def serialize_dtype(value):
                    if not (isinstance(value, np.ndarray) ^ quantity.is_scalar):
                        self.m_warning('numpy quantity has wrong shape', quantity=str(quantity))

                    return value.tolist() if isinstance(value, np.ndarray) else value.item()

                serialize = serialize_dtype

            elif isinstance(quantity_type, MEnum):
                def serialize_enum(value):
                    return None if value is None else str(value)

                serialize = serialize_enum

            elif quantity_type == Any:
                def serialize_any(value: Any):
                    if type(value) not in [str, int, float, bool, np.bool_, list, dict, type(None)]:
                        raise MetainfoError(
                            f'Only python primitives are allowed for Any typed non-virtual '
                            f'quantities: {value} of quantity {quantity} in section {self}')

                    return value

                serialize = serialize_any

            else:
                raise MetainfoError(
                    f'Do not know how to serialize data with type {quantity_type} for quantity {quantity}')

            quantity_type = quantity.type
            if resolve_references and isinstance(quantity_type, QuantityReference):
                serialize_before_reference_resolution = serialize

                def serialize_reference_v2(value: Any):
                    value = getattr(value.m_resolved(), quantity_type.target_quantity_def.name)

                    return serialize_before_reference_resolution(value)

                serialize = serialize_reference_v2

            # get the value to be serialized
            # explicitly assigning the target value overrides the value from the section
            if target_value is None:
                if is_set:
                    target_value = self.__dict__[quantity.name]
                elif is_derived:
                    try:
                        target_value = quantity.derived(self)
                    except Exception:
                        target_value = quantity.default
                else:
                    target_value = quantity.default

            if transform is not None:
                serialize_before_transform = serialize

                def serialize_and_transform(value: Any, path_override=None):
                    if not is_reference:
                        return transform(quantity, self, serialize_before_transform(value), path_override)

                    return transform(quantity, self, serialize_before_transform(value, path_override), path_override)

                serialize = serialize_and_transform

            # serialization starts here
            if quantity_type in MTypes.numpy or quantity_type in MTypes.complex:
                return serialize(target_value)

            if len(quantity.shape) == 0:
                return serialize(target_value, path) if is_reference else serialize(target_value)

            if len(quantity.shape) == 1:
                if not is_reference:
                    return [serialize(item) for item in target_value]

                return [serialize(item, f"{path}/{index}") for index, item in enumerate(target_value)]

            raise NotImplementedError(f'Higher shapes ({quantity.shape}) not supported: {quantity}')

        def serialize_attribute(attribute: 'Attribute', value: Any) -> Any:
            if isinstance(attribute.type, DataType):
                return attribute.type.serialize(self, None, value)

            if attribute.type in MTypes.complex:
                return serialize_complex(value)

            if attribute.type in MTypes.primitive:
                if len(attribute.shape) == 0:
                    return MTypes.primitive[attribute.type](value)  # type: ignore

                return [MTypes.primitive[attribute.type](v) for v in value]  # type: ignore

            if isinstance(attribute.type, MEnum):
                return str(value)

            if isinstance(attribute.type, np.dtype):
                return value.item()

            return value

        def collect_attributes(attr_map: dict, all_attr: dict):
            result: dict = {}
            for attr_key, attr_value in attr_map.items():
                attr_def = resolve_variadic_name(all_attr, attr_key)
                result[attr_key] = serialize_attribute(attr_def, attr_value)
            return result

        def serialize_full_quantity(quantity_def: Quantity, values: Dict[str, MQuantity]):
            result: dict = {}
            for m_quantity in values.values():
                m_result: dict = {
                    'm_value': serialize_quantity(quantity_def, True, False, None, m_quantity.value)}
                if m_quantity.unit:
                    m_result['m_unit'] = str(m_quantity.unit)
                if m_quantity.original_unit:
                    m_result['m_original_unit'] = str(m_quantity.original_unit)
                if m_quantity.attributes:
                    a_result: dict = collect_attributes(m_quantity.attributes, quantity_def.all_attributes)
                    if a_result:
                        m_result['m_attributes'] = a_result
                result[m_quantity.name] = m_result

            return result

        def serialize_annotation(annotation):
            if isinstance(annotation, Annotation):
                return annotation.m_to_dict()
            elif isinstance(annotation, Dict):
                try:
                    json.dumps(annotation)
                    return annotation
                except Exception:
                    return str(annotation)
            else:
                return str(annotation)

        def items() -> Iterable[Tuple[str, Any]]:
            # metadata
            if with_meta:
                yield 'm_def', self.m_def.definition_reference(self)
                if with_def_id:
                    yield 'm_def_id', self.m_def.definition_id
                if self.m_parent_index != -1:
                    yield 'm_parent_index', self.m_parent_index
                if self.m_parent_sub_section is not None:
                    yield 'm_parent_sub_section', self.m_parent_sub_section.name

                annotations = {}
                for annotation_name, annotation in self.m_annotations.items():
                    if isinstance(annotation, list):
                        annotation_value = [serialize_annotation(item) for item in annotation]
                    else:
                        annotation_value = [serialize_annotation(annotation)]
                    annotations[annotation_name] = annotation_value
                if len(annotations) > 0:
                    yield 'm_annotations', annotations
            elif with_root_def:
                yield 'm_def', self.m_def.definition_reference(self)
                if with_def_id:
                    yield 'm_def_id', self.m_def.definition_id
            elif self.m_parent and self.m_parent_sub_section.sub_section != self.m_def:
                # The subsection definition's section def is different from our
                # own section def. We are probably a specialized derived section
                # from the base section that was used in the subsection def. To allow
                # clients to recognize the concrete section def, we force the export
                # of the section def.
                yield 'm_def', self.m_def.definition_reference(self)
                if with_def_id:
                    yield 'm_def_id', self.m_def.definition_id

            # quantities
            sec_path = self.m_path()
            for name, quantity in self.m_def.all_quantities.items():
                path = f"{sec_path}/{name}"
                if exclude(quantity, self):
                    continue

                try:
                    if quantity.virtual:
                        if include_derived and quantity.derived is not None:
                            yield name, serialize_quantity(quantity, False, True, path)
                        continue

                    is_set = self.m_is_set(quantity)
                    if not is_set:
                        if not include_defaults or not quantity.m_is_set(Quantity.default):
                            continue

                    if not quantity.use_full_storage:
                        yield name, serialize_quantity(quantity, is_set, False, path)
                    else:
                        yield name, serialize_full_quantity(quantity, self.__dict__[quantity.name])

                except ValueError as e:
                    raise ValueError(f'Value error ({str(e)}) for {quantity}')

            # section attributes
            if 'm_attributes' in self.__dict__:
                yield 'm_attributes', collect_attributes(
                    self.__dict__['m_attributes'], self.m_def.all_attributes)

            # subsections
            for name, sub_section_def in self.m_def.all_sub_sections.items():
                if exclude(sub_section_def, self):
                    continue

                is_set = False
                if sub_section_def.repeats:
                    if self.m_sub_section_count(sub_section_def) > 0:
                        is_set = True
                        yield name, [
                            None if item is None else item.m_to_dict(**kwargs)
                            for item in self.m_get_sub_sections(sub_section_def)]
                else:
                    sub_section = self.m_get_sub_section(sub_section_def, -1)
                    if sub_section is not None:
                        is_set = True
                        yield name, sub_section.m_to_dict(**kwargs)

                # attributes are disabled for subsections
                # if is_set:
                #     yield from collect_attributes(sub_section_def.all_attributes)

        return {key: value for key, value in items()}

    @staticmethod
    def __deserialize(section: MSection, quantity_def: Quantity, quantity_value: Any):
        tgt_type = quantity_def.type

        if tgt_type in MTypes.complex:
            return normalize_complex(quantity_value, tgt_type, quantity_def.unit)

        if tgt_type in MTypes.numpy:
            if not isinstance(quantity_value, list):
                return tgt_type(quantity_value)

            return np.asarray(quantity_value).astype(tgt_type)

        if isinstance(tgt_type, DataType):
            def __type_specific_deserialize(v):
                return tgt_type.deserialize(section, quantity_def, v)

            dimensions = len(quantity_def.shape)

            if dimensions == 0:
                return __type_specific_deserialize(quantity_value)
            if dimensions == 1:
                return list(__type_specific_deserialize(item) for item in quantity_value)

            raise MetainfoError('Only numpy quantities can have more than 1 dimension.')

        return quantity_value

    def m_update_from_dict(self, dct: Dict[str, Any]) -> None:
        '''
        Updates this section with the serialized data from the given dict, e.g. data
        produced by :func:`m_to_dict`.
        '''
        section_def = self.m_def
        section = self
        m_context = self.m_context if self.m_context else self

        if 'definitions' in dct:
            definition_def = section_def.all_aliases['definitions']
            definition_cls = definition_def.sub_section.section_cls
            definition_section = definition_cls.m_from_dict(
                dct['definitions'], m_parent=self, m_context=m_context)
            section.m_add_sub_section(definition_def, definition_section)

        for name, property_def in section_def.all_aliases.items():
            if name not in dct or name == 'definitions':
                continue

            if isinstance(property_def, SubSection):
                sub_section_def = property_def
                sub_section_value = dct.get(name)
                sub_section_cls = sub_section_def.sub_section.section_cls
                if sub_section_def.repeats:
                    for sub_section_dct in sub_section_value:
                        sub_section = None if sub_section_dct is None else sub_section_cls.m_from_dict(
                            sub_section_dct, m_parent=self, m_context=m_context)
                        section.m_add_sub_section(sub_section_def, sub_section)
                else:
                    sub_section = sub_section_cls.m_from_dict(
                        sub_section_value, m_parent=self, m_context=m_context)
                    section.m_add_sub_section(sub_section_def, sub_section)

            if isinstance(property_def, Quantity):
                quantity_def = property_def
                quantity_value = dct[name]

                if quantity_def.virtual:
                    # We silently ignore this, similar to how we ignore additional values.
                    continue

                if quantity_def.use_full_storage:
                    if not isinstance(quantity_value, dict):
                        raise MetainfoError('Full storage quantity must be a dict')

                    for each_name, each_quantity in quantity_value.items():
                        try:
                            m_value = self.__deserialize(section, quantity_def, each_quantity['m_value'])
                        except KeyError:
                            raise MetainfoError(f'Set full storage quantity {property_def} must have a value')
                        m_quantity = MQuantity(each_name, m_value)
                        if 'm_unit' in each_quantity:
                            m_quantity.unit = units.parse_units(each_quantity['m_unit'])
                        if 'm_original_unit' in each_quantity:
                            m_quantity.original_unit = units.parse_units(each_quantity['m_original_unit'])
                        if 'm_attributes' in each_quantity:
                            m_quantity.attributes = each_quantity['m_attributes']

                        section.m_set(quantity_def, m_quantity)
                else:
                    section.__dict__[property_def.name] = self.__deserialize(section, quantity_def, quantity_value)

        if 'm_attributes' in dct:
            for attr_key, attr_value in dct['m_attributes'].items():
                section.m_set_section_attribute(attr_key, attr_value)

    @classmethod
    def m_from_dict(cls: Type[MSectionBound], data: Dict[str, Any], **kwargs) -> MSectionBound:
        ''' Creates a section from the given serializable data dictionary.

        This is the 'opposite' of :func:`m_to_dict`. It takes a deserialized dict, e.g.
        loaded from JSON, and turns it into a proper section, i.e. instance of the given
        section class.
        '''
        return MSection.from_dict(data, cls=cls, **kwargs)

    @staticmethod
    def from_dict(
            dct: Dict[str, Any],
            cls: Type[MSectionBound] = None,
            m_parent: MSection = None,
            m_context: 'Context' = None,
            **kwargs
    ) -> MSectionBound:
        ''' Creates a section from the given serializable data dictionary.

        This is the 'opposite' of :func:`m_to_dict`. It is similar to the classmethod
        `m_from_dict`, but does not require a specific class. You can provide a class
        through the optional parameter. Otherwise, the section definition is read from
        the `m_def` key in the section data.
        '''
        if 'm_ref_archives' in dct and isinstance(m_context, Context):
            # dct['m_ref_archives'] guarantees that 'm_def' exists
            for entry_url, archive_json in dct['m_ref_archives'].items():
                m_context.cache_archive(
                    entry_url, MSection.from_dict(archive_json, m_parent=m_parent, m_context=m_context))
            del dct['m_ref_archives']

        # first try to find a m_def in the data
        if 'm_def' in dct:
            # We re-use the _SectionReference implementation for m_def
            m_def = dct['m_def']
            context_section = m_parent
            archive_root: MSection = m_parent.m_root() if m_parent else None
            if archive_root:
                definitions = getattr(archive_root, 'definitions', None)
                if isinstance(definitions, Package):
                    context_section = definitions
            m_def_proxy = SectionReference.deserialize(context_section, None, m_def)
            m_def_proxy.m_proxy_context = m_context
            cls = m_def_proxy.section_cls

        # if 'm_def_id' exists, check if id matches
        # in case of mismatch, retrieve the Package and use the corresponding section definition
        if 'm_def_id' in dct:
            if cls is None or cls.m_def is None or dct['m_def_id'] != cls.m_def.definition_id:
                if not isinstance(m_context, Context):
                    raise MetainfoError(f"A context object is needed to resolve definition {dct['m_def_id']}")
                cls = m_context.resolve_section_definition(dct.get('m_def', None), dct['m_def_id'])

        assert cls is not None, 'Section definition or cls needs to be known'

        if isinstance(m_context, Context):
            section = cls(m_context=m_context, **kwargs)
        else:
            section = cls(**kwargs)

        # We have to set this prematurely. It would be set later on, but to late to get
        # the proper root for subsequent resolution of m_def references.
        section.m_parent = m_parent

        if 'm_annotations' in dct:
            m_annotations = dct['m_annotations']
            if not isinstance(m_annotations, dict):
                raise MetainfoError(
                    f'The provided m_annotations is of a wrong type. {type(m_annotations).__name__} was provided.')
            section.m_annotations.update(m_annotations)
            section.m_parse_annotations()

        section.m_update_from_dict(dct)
        return section

    def m_to_json(self, **kwargs):
        ''' Returns the data of this section as a json string. '''
        return json.dumps(self.m_to_dict(), **kwargs)

    def m_all_contents(
            self, depth_first: bool = False, include_self: bool = False,
            stop: TypingCallable[[MSection], bool] = None) -> Iterable[MSection]:
        '''
        Returns an iterable over all sub and sub subsections.

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
        parent index for all set properties. If the section has no property the empty
        section is returned.
        '''

        empty = True
        for key in self.__dict__:
            property_def = self.m_def.all_properties.get(key)
            if property_def is None:
                continue

            empty = False

            if isinstance(property_def, SubSection):
                for sub_section in self.m_get_sub_sections(property_def):
                    for i in sub_section.m_traverse():
                        yield i
                    yield self, property_def, sub_section.m_parent_index

            else:
                yield self, property_def, -1

        if empty:
            yield self, None, -1

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

    def m_contents(self) -> Iterable[MSection]:
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

    def m_path(self, quantity_def: Quantity = None) -> str:
        ''' Returns the path of this section or the given quantity within the section hierarchy. '''
        if self.m_parent is None:
            return '/'

        if self.m_parent_index == -1:
            segment = self.m_parent_sub_section.name
        else:
            segment = f'{self.m_parent_sub_section.name}/{self.m_parent_index:d}'

        if quantity_def is not None:
            segment = f'{segment}/{quantity_def.name}'

        return f'{self.m_parent.m_path().rstrip("/")}/{segment}'

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

    def m_resolve(self, path_with_id: str, cls: Type[MSectionBound] = None) -> MSectionBound:
        '''
        Resolves the given path or dotted quantity name using this section as context and
        returns the sub_section or value.

        Arguments:
            path_with_id: The reference URL. See `MProxy` for details on reference URLs.
        '''
        section: MSection = self

        if '@' in path_with_id:
            path, target_id = path_with_id.split('@')
        else:
            target_id = None
            path = path_with_id

        if path.startswith('/'):
            section = section.m_root()

        path_stack = path.strip('/').split('/')
        path_stack.reverse()
        while len(path_stack) > 0:
            prop_name = path_stack.pop()
            prop_def = section.m_def.all_properties.get(prop_name, None)

            if prop_def is None:
                prop_def = section.m_def.all_aliases.get(prop_name, None)

            if prop_def is None:
                raise MetainfoReferenceError(
                    f'Could not resolve {path}, property {prop_name} does not exist')

            if isinstance(prop_def, SubSection):
                if prop_def.repeats:
                    if len(path_stack) == 0:
                        return _check_definition_id(target_id, section.m_get_sub_sections(prop_def))  # type: ignore

                    try:
                        sub_section: str = path_stack.pop()
                        index: Any = int(sub_section) if sub_section.isnumeric() else sub_section
                    except ValueError:
                        raise MetainfoReferenceError(
                            f'Could not resolve {path}, {prop_name} repeats but there is no '
                            'index in the path')

                    try:
                        section = section.m_get_sub_section(prop_def, index)
                    except Exception:
                        raise MetainfoReferenceError(
                            f'Could not resolve {path}, there is no subsection for '
                            f'{prop_name} at {index}')

                else:
                    section = section.m_get_sub_section(prop_def, -1)
                    if section is None:
                        raise MetainfoReferenceError(
                            f'Could not resolve {path}, there is no subsection {prop_name}')

            elif isinstance(prop_def, Quantity):
                if not section.m_is_set(prop_def):
                    raise MetainfoReferenceError(
                        f'Could not resolve {path}, {prop_name} is not set in {section}')

                quantity = section.m_get(prop_def)
                while len(path_stack) > 0:
                    if not isinstance(quantity, list):
                        raise MetainfoReferenceError(f'Could not resolve {path}, no property {prop_name}')
                    quantity = quantity[int(path_stack.pop())]

                return _check_definition_id(target_id, quantity)

        return _check_definition_id(target_id, cast(MSectionBound, section))

    def m_get_annotations(self, key: Union[str, type], default=None, as_list: bool = False):
        '''
        Convenience method to get annotations

        Arguments:
            key: Either the optional annotation name or an annotation class. In the first
                case the annotation is returned, regardless of its type. In the second
                case, all names and list for names are iterated and all annotations of the
                given class are returned.
            default: The default, if no annotation is found. None is the default `default`.
            as_list: Returns a list, no matter how many annotations have been found.
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

    def m_validate(self) -> Tuple[List[str], List[str]]:
        ''' Evaluates all constraints and shapes of this section and returns a list of errors. '''
        errors: List[str] = []
        warnings: List[str] = []
        if self.m_parent and hasattr(self.m_parent, 'all_base_sections'):
            for base_sections in self.m_parent.all_base_sections:
                for constraint_name in base_sections.constraints:
                    constraint = getattr(base_sections.section_cls, constraint_name, None)
                    if constraint is None:
                        raise MetainfoError(
                            f'Could not find implementation for constraint {constraint_name} of section {self.m_def}.')
                    try:
                        constraint(self)
                    except AssertionError as e:
                        error_str = str(e).strip()
                        if error_str == '':
                            error_str = f'Constraint {constraint_name} violated.'
                        if getattr(constraint, 'm_warning', False):
                            warnings.append(error_str)
                        else:
                            errors.append(error_str)

        for constraint_name in self.m_def.constraints:
            constraint = getattr(self, constraint_name, None)
            if constraint is None:
                raise MetainfoError(
                    f'Could not find implementation for constraint {constraint_name} of section {self.m_def}.')

            try:
                constraint()
            except AssertionError as e:
                error_str = str(e).strip()
                if error_str == '':
                    error_str = f'Constraint {constraint_name} violated.'
                if getattr(constraint, 'm_warning', False):
                    warnings.append(error_str)
                else:
                    errors.append(error_str)

        for quantity in self.m_def.all_quantities.values():
            if self.m_is_set(quantity) and not quantity.derived and quantity != Quantity.default:
                if not validate_shape(self, quantity, self.m_get(quantity)):
                    errors.append(f'The shape of quantity {quantity} does not match its value.')

        for annotation in self.m_annotations.values():
            def validate_annotation(annotation):
                if isinstance(annotation, AnnotationModel):
                    if annotation.m_error:
                        # This annotation could not be parsed and only contains the error.
                        errors.append(annotation.m_error)
                        return

                    try:
                        # Trigger model validation by re-assigning the definition
                        annotation.m_definition = self
                    except ValidationError as e:
                        errors.append(f'Annotation validation error for {self}: {str(e)}')

            if isinstance(annotation, list):
                for item in annotation:
                    validate_annotation(item)
            else:
                validate_annotation(annotation)

        return errors, warnings

    def m_copy(
            self: MSectionBound,
            deep=False,
            parent=None,
            es_annotation: List[Elasticsearch] = None) -> MSectionBound:
        '''
            es_annotation: Optional annotation for ElasticSearch. Will override
                any existing annotation.
        '''
        # TODO this a shallow copy, but should be a deep copy
        copy = self.m_def.section_cls()
        copy.__dict__.update(**self.__dict__)
        copy.m_parent = parent
        copy.m_parent_index = -1 if parent is None else self.m_parent_index
        copy.m_parent_sub_section = None if parent is None else self.m_parent_sub_section
        copy.m_context = self.m_context

        if deep:
            if isinstance(copy, Definition):
                copy.more = deepcopy(self.more)
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

        if es_annotation:
            copy.m_annotations["elasticsearch"] = es_annotation
        return cast(MSectionBound, copy)

    def m_all_validate(self):
        ''' Evaluates all constraints in the whole section hierarchy, incl. this section. '''
        errors: List[str] = []
        warnings: List[str] = []
        for section in self.m_all_contents(include_self=True):
            more_errors, more_warnings = section.m_validate()
            errors.extend(more_errors)
            warnings.extend(more_warnings)

        return errors, warnings

    def m_warning(self, *args, **kwargs):
        if self.m_context is not None:
            self.m_context.warning(*args, **kwargs)

    def __repr__(self):
        m_section_name = self.m_def.name
        # name_quantity_def = self.m_def.all_quantities.get('name', None)
        # if name_quantity_def is not None:
        #     name = self.m_get(name_quantity_def)
        try:
            name = self.__dict__['name']
            main = f'{name}:{m_section_name}'
        except KeyError:
            main = m_section_name

        more = ''
        props = [prop for prop in self.m_def.all_properties if prop in self.__dict__]

        if len(props) > 10:
            more = f', +{len(props) - 10:d} more properties'

        return f'{main}({", ".join(props[0:10])}{more})'

    def __getitem__(self, key):
        try:
            key = key.replace('.', '/')
            return self.m_resolve(key)
        except MetainfoReferenceError:
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

        return to_dict(jmespath.search(expression, self))


class MCategory(metaclass=MObjectMeta):
    m_def: Category = None

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

    All metainfo `definitions` (sections, quantities, subsections, packages, ...) share
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
            case `snake_case` identifier. Subsections are often prefixed with ``section_``
            to clearly separate subsections from quantities.

            Generally, you do not have to set this attribute manually, it will be derived
            from Python identifiers automatically.

        label: Each `definition` can have an optional label. Label are like names, but
            do not have to adhere to the Python identfier syntax.

        description: The description can be an arbitrary human-readable text that explains
            what a definition is about. For section definitions you do not have to set
            this manually as it will be derived from the classes doc string. Quantity and
            subsection descriptions can also be taken from the containing section class'
            doc-string ``Attributes:`` section.

        links: Each definition can be accompanied by a list of URLs. These should point
            to resources that further explain the definition.

        aliases: A list of alternative names. For quantities and subsections these
            can be used to access the respective property with a different name from
            its containing section.

        variable:
            A boolean that indicates this property as variable parts in its name.
            If this is set to true, all capital letters in the name can be
            replaced with arbitrary strings. However, variable names work similar to
            aliases and can be considered on-demand aliases. Other aliases and the
            defined name will work as well. Thus, variable names are only resolved
            at runtime by the Python interface and are not directly serialized.
            However, the variable name is set in a meta attribute `m_source_name`
            automatically for properties (but not attributes).
            Variable names are only implemented for Quantity, SubSection,
            Attribute.

        deprecated: If set this definition is marked deprecated. The value should be a
            string that describes how to replace the deprecated definition.

        categories: All metainfo definitions can be put into one or more `categories`.
            Categories allow to organize the definitions themselves. It is different from
            sections, which organize the data (e.g. quantity values) and not the definitions
            of data (e.g. quantities definitions). See :ref:`metainfo-categories` for more
            details.

        more: A dictionary that contains additional definition properties that are not
            part of the metainfo. Those can be passed as additional kwargs to definition
            constructors. The values must be JSON serializable.

        attributes:
            The attributes that can further qualify property values.

        all_attributes:
            A virtual convenient property that provides all attributes as a dictionary
            from attribute name to attribute. This includes meta attributes (starting with m_)
            that are defined for all properties of the same kind (sub_section or quantity).
    '''

    name: Quantity = _placeholder_quantity
    label: Quantity = _placeholder_quantity
    description: Quantity = _placeholder_quantity
    links: Quantity = _placeholder_quantity
    categories: Quantity = _placeholder_quantity
    deprecated: Quantity = _placeholder_quantity
    aliases: Quantity = _placeholder_quantity
    variable: Quantity = _placeholder_quantity
    more: Quantity = _placeholder_quantity

    attributes: SubSection = None  # type: ignore

    all_attributes: Quantity = _placeholder_quantity

    # store the hash object generated
    _cached_hash: _HASH_OBJ = None  # type: ignore

    def __init__(self, *args, **kwargs):
        if is_bootstrapping:
            super().__init__(*args, **kwargs)
            return

        # We add all kwargs that are not meta props, annotations, or metainfo properties
        # to the more property.
        more = {}
        new_kwargs = {}
        for key, value in kwargs.items():
            if key.startswith('m_') or key.startswith('a_') or key in self.__class__.m_def.all_aliases:
                new_kwargs[key] = value
            else:
                more[key] = value

        super().__init__(*args, **new_kwargs)
        self.more = more  # type: ignore

        self._cached_hash = None  # type: ignore

    def __init_metainfo__(self):
        '''
        An initialization method that is called after the class context of the definition
        has been initialized. For example, it is called on all quantities of a section
        class after the class was created. If metainfo definitions are created without
        a class context, this method must be called manually on all definitions.
        '''

        # for base_section in self.all_base_sections:
        #     for constraint in base_section.constraints:
        #         constraints.add(constraint)
        #     for event_handler in base_section.event_handlers:
        #         event_handlers.add(event_handler)

        # initialize definition annotations
        for annotation in self.m_get_annotations(DefinitionAnnotation, as_list=True):
            annotation.init_annotation(self)

    def init_metainfo(self):
        '''
        Calls __init_metainfo__ on all its children. This is necessary if the
        package, section was created without corresponding python classes
        packages, etc.
        '''
        self.__init_metainfo__()
        for content in self.m_all_contents(depth_first=True):
            content.__init_metainfo__()

    def __getattr__(self, name):
        if self.more and name in self.more:
            return self.more[name]

        return super().__getattr__(name)

    def m_is_set(self, quantity_def: Quantity) -> bool:
        if quantity_def == Definition.more:
            return len(self.more) > 0

        return super().m_is_set(quantity_def)

    def qualified_name(self):
        name = self.name if self.name else '*'
        if self.m_parent and self.m_parent.m_follows(Definition.m_def):
            return f'{self.m_parent.qualified_name()}.{name}'

        return name

    def on_set(self, quantity_def, value):
        if quantity_def == Definition.categories:
            for category in value:
                category.definitions.add(self)

    def m_to_dict(self, **kwargs) -> Dict[str, Any]:  # type: ignore
        if kwargs.get('with_def_id', False):
            value: dict = dict(definition_id=self.definition_id)
            value.update(super(Definition, self).m_to_dict(**kwargs))
            return value

        return super(Definition, self).m_to_dict(**kwargs)

    def __repr__(self):
        return f'{self.qualified_name()}:{self.m_def.name}'

    def _hash_seed(self) -> str:
        '''
        Generates a unique representation for computing the hash of a definition.

        The order of aliases is not important.
        '''
        seed: str = str(self.name)
        seed += 'T' if self.variable else 'F'
        if len(self.aliases) > 0:
            seed += ''.join([str(i) for i in sorted(self.aliases)])

        return seed

    def _hash(self, regenerate=False) -> _HASH_OBJ:
        '''
        Generates a hash object based on the unique representation of the definition.
        '''
        if self._cached_hash is None or regenerate:
            self._cached_hash = default_hash()
            self._cached_hash.update(self._hash_seed().encode('utf-8'))

        if self.attributes:
            for item in self.attributes:  # pylint: disable=not-an-iterable
                if id(self) != id(item):
                    self._cached_hash.update(item._hash(regenerate).digest())

        return self._cached_hash

    @property
    def definition_id(self) -> str:
        '''
        Syntax sugar.

        Returns the hash digest.
        '''
        return self._hash().hexdigest()

    def definition_reference(self, source, **kwargs):
        '''
        Creates a reference string that points to this definition from the
        given source section.
        '''
        definition_reference = self.qualified_name()
        if definition_reference.startswith('entry_id:'):
            # This is not from a python module, use archive reference instead
            # two cases:
            # 1. loaded from file so archive.definitions.archive is set by parser
            # 2. loaded from versioned mongo so entry_id_based_name is set by mongo
            # second one has no metadata, so do not create reference
            context = self.m_root().m_context
            if context:
                relative_name = context.create_reference(source, None, self, **kwargs)
                if relative_name:
                    definition_reference = relative_name

        if process.add_definition_id_to_reference and '@' not in definition_reference:
            definition_reference += '@' + self.definition_id

        return definition_reference


class Attribute(Definition):
    '''
    Attributes can be used to qualify all properties (subsections and quantities)
    with simple scalar values.

    Attributes:
        type: The type of the attribute. Can be any primitive type, including
            numpy types, Datetime and enums.
    '''

    type: Quantity = _placeholder_quantity
    shape: Quantity = _placeholder_quantity

    def __init__(self, *args, **kwargs):
        super(Attribute, self).__init__(*args, **kwargs)

    @constraint(warning=False)
    def is_primitive(self):
        if self.type in MTypes.primitive or self.type in MTypes.num:
            return

        if isinstance(self.type, (MEnum, np.dtype, _Datetime)):
            return

        assert False, 'Attributes must have primitive type.'

    def _hash_seed(self) -> str:
        seed = super(Attribute, self)._hash_seed()
        type_id = QuantityType.serialize(self, Quantity.type, self.type)
        if 'type_data' in type_id and isinstance(type_id['type_data'], list):
            type_id['type_data'].sort()
        seed += json.dumps(type_id)
        for dim in self.shape:
            seed += str(dim)
        return seed


class Property(Definition):
    '''
    A common base-class for section properties: subsections and quantities.
    '''

    def get_from_dict(self, data: Dict[str, Any], default_value: Any = None) -> Tuple[Optional[str], Any]:
        '''
        Attempts to read the property from a dict. Returns the used alias and value as
        tuple.
        '''
        for name in self.aliases + [self.name]:
            if name in data:
                return name, data[name]
        return None, default_value

    def get_base_property(self) -> Optional['Property']:
        '''
        Retrieve a potential overwritten property from a base-class.
        '''
        assert self.m_parent and isinstance(self.m_parent, Section), 'Property must be property of a section.'
        section = cast(Section, self.m_parent)
        for base_section in section.all_base_sections.values():
            base_property = base_section.all_properties.get(self.name)
            if base_property:
                if base_property.m_def != self.m_def:
                    raise MetainfoError('Cannot overwrite a property of different metainfo type.')
                return base_property

        return None


class Quantity(Property):
    '''
    To define quantities, instantiate :class:`Quantity` as a class attribute values in
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
            - a numpy `dtype`, e.g. ``np.float32``
            - ``typing.Any`` to support any value

            If set to `dtype`, this quantity will use a numpy array or scalar to store values
            internally. If a regular (nested) Python list or Python scalar is given, it will
            be automatically converted. The given `dtype` will be used in the numpy value.

            To define a reference, either a `section class` or instance of :class:`Section`
            can be given. See :ref:`metainfo-sections` for details. Instances of the given section
            constitute valid values for this type. Upon serialization, references section
            instance will represent with metainfo URLs. See :ref:`metainfo-urls`.

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
            numpy arrays. Their type must be a `dtype`.

        is_scalar:
            Derived quantity that is True, iff this quantity has shape of length 0

        unit:
            The physics unit for this quantity. It is optional.

            Units are represented with the Pint Python package. Pint defines units and
            their algebra. You can either use *pint* units directly, e.g. ``units.m / units.s``.
            The metainfo provides a preconfigured *pint* unit registry :py:data:`ureg`.
            You can also provide the unit as *pint* parsable string, e.g. ``'meter / seconds'`` or
            ``'m/s'``.

        dimensionality:
            The dimensionality of the quantity. It is optional.

            If set, it will be used to validate the compatibility between chosen unit and the target
            dimensionality.

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
            be got/set like regular quantities, but their values are not (de-)serialized,
            hence never permanently stored.
    '''

    type: Quantity = _placeholder_quantity
    shape: Quantity = _placeholder_quantity
    unit: Quantity = _placeholder_quantity
    dimensionality: Quantity = _placeholder_quantity
    default: Quantity = _placeholder_quantity
    derived: Quantity = _placeholder_quantity
    cached: Quantity = _placeholder_quantity
    virtual: Quantity = _placeholder_quantity

    is_scalar: Quantity = _placeholder_quantity
    repeats: Quantity = _placeholder_quantity
    use_full_storage: Quantity = _placeholder_quantity
    flexible_unit: Quantity = _placeholder_quantity

    # TODO derived_from = Quantity(type=Quantity, shape=['0..*'])
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __init_metainfo__(self):
        super().__init_metainfo__()

        if self.derived is not None:
            self.virtual = True  # type: ignore

        # replace the quantity implementation with an optimized version for the most
        # primitive quantities if applicable
        is_primitive = not self.derived and not self.use_full_storage
        is_primitive = is_primitive and len(self.shape) <= 1
        is_primitive = is_primitive and self.type in [str, bool, float, int]
        is_primitive = is_primitive and self.type not in MTypes.num_numpy
        if is_primitive:
            self._default = self.default
            self._name = self.name
            self._type = self.type
            self._list = len(self.shape) == 1
            self.__class__ = PrimitiveQuantity

        check_dimensionality(self, self.unit)

    def __get__(self, obj, cls):
        try:
            value = obj.__dict__[self.name]
            # appears to be a quantity using full storage
            # cannot use .use_full_storage as this is not set yet
            if isinstance(value, dict) and self.name in value:
                m_quantity = value[self.name]
                if m_quantity.unit:
                    value = units.Quantity(m_quantity.value, m_quantity.unit)
                else:
                    value = m_quantity.value
        except KeyError:
            if self.derived is not None:
                try:
                    if self.cached:
                        cached = obj.__dict__.setdefault(self.name + '_cached', [-1, None])
                        if cached[0] != obj.m_mod_count:
                            cached[0] = obj.m_mod_count
                            cached[1] = self.derived(obj)  # pylint: disable=not-callable
                        return cached[1]

                    return self.derived(obj)  # pylint: disable=not-callable
                except Exception as e:
                    raise DeriveError(f'Could not derive value for {self}: {str(e)}')

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
                value = list(self.type.get_normalize(obj, self, item) for item in value)

            else:
                raise MetainfoError(
                    'Only numpy arrays and dtypes can be used for higher dimensional quantities.')

        # no need to append unit if it is already a quantity from full storage
        if isinstance(value, units.Quantity):
            return value

        if self.unit is not None and self.type in MTypes.num:
            return value * self.unit

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

    @constraint(warning=False)
    def has_type(self):
        assert self.type is not None, f'The quantity {self.qualified_name()} must define a type'
        if isinstance(self.type, Reference):
            try:
                self.type.target_section_def.m_resolved()
            except MetainfoReferenceError as e:
                assert False, f'Cannot resolve "type" of {self.qualified_name()}: {str(e)}'

    @constraint(warning=True)
    def dimensions(self):
        for dimension in self.shape:
            if not isinstance(dimension, str) or not dimension.isidentifier():
                continue

            dim_quantity = self.m_parent.all_quantities.get(dimension, None)

            assert dim_quantity is not None, \
                f'Dimensions ({dimension}) must be quantities of the same section ({self.m_parent}).'

            assert dim_quantity.type in [int, np.int16, np.int32, np.int8, np.uint8] and len(
                dim_quantity.shape) == 0, \
                f'Dimensions ({dimension}) must be shapeless ({dim_quantity.shape}) ' \
                f'and int ({dim_quantity.type}) typed.'

    @constraint(warning=True)
    def higher_shapes_require_dtype(self):
        if len(self.shape) > 1:
            assert self.type in MTypes.numpy, \
                f'Higher dimensional quantities ({self}) need a dtype and will be treated as numpy arrays.'

    def _hash_seed(self) -> str:
        '''
        Generate a unique representation for this quantity.

        The custom type MUST have the method `_hash_seed()` to be used to generate a proper id.

        Returns:
            str: The generated unique representation.
        '''
        new_id = super(Quantity, self)._hash_seed()

        if isinstance(self.type, Reference):
            new_id += f'Ref->{self.type.target_section_def.qualified_name()}'
        else:
            try:
                type_id = QuantityType.serialize(self, Quantity.type, self.type)
                if 'type_data' in type_id:
                    if isinstance(type_id['type_data'], list):
                        type_id['type_data'].sort()
                new_id += json.dumps(type_id)
            except MetainfoError as e:
                # unlikely but for completeness
                if hasattr(self.type, '_hash_seed'):
                    new_id += self.type._hash_seed()
                else:
                    raise e

        for dim in self.shape:
            new_id += dim if isinstance(dim, str) else str(dim)

        new_id += f'{str(self.unit)}{"N" if self.default is None else self.default}'
        new_id += str(self.dimensionality)
        new_id += "T" if self.virtual else "F"

        return new_id


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
            value = obj.__dict__[self._name]
        except KeyError:
            value = self._default
        except AttributeError:
            return self
        if value is not None and self.unit is not None and self.type in MTypes.num:
            return value * self.unit  # type: ignore
        return value

    def __set__(self, obj, value):
        obj.m_mod_count += 1

        if value is None:
            obj.__dict__.pop(self.name, None)
            return

        # Handle pint quantities. Conversion is done automatically between
        # units. Notice that currently converting from float to int or vice
        # versa is not allowed for primitive types.
        if isinstance(value, pint.Quantity):
            if self.unit is None:
                raise TypeError(f'The quantity {self} does not have a unit, but value {value} has.')
            if self.type in MTypes.int:
                raise TypeError(
                    f'Cannot save data with unit conversion into the quantity {self} '
                    'with integer data type due to possible precision loss.')
            value = value.to(self.unit).magnitude

        if self._list:
            if not isinstance(value, list):
                if hasattr(value, 'tolist'):
                    value = value.tolist()
                else:
                    raise TypeError(f'The value {value} for quantity {self} has no shape {self.shape}')

            if any(v is not None and type(v) != self._type for v in value):
                raise TypeError(
                    f'The value {value} with type {type(value)} for quantity {self} is not of type {self.type}')

        elif type(value) != self._type:
            raise TypeError(
                f'The value {value} with type {type(value)} for quantity {self} is not of type {self.type}')

        try:
            obj.__dict__[self._name] = value
        except AttributeError:
            raise KeyError('Cannot overwrite quantity definition. Only values can be set.')


class SubSection(Property):
    '''
    Like quantities, subsections are defined in a `section class` as attributes
    of this class. Unlike quantities, each subsection definition becomes a property of
    the corresponding `section definition` (parent). A subsection definition references
    another `section definition` as the subsection (child). As a consequence, parent
    `section instances` can contain child `section instances` as subsections.

    Contrary to the old NOMAD metainfo, we distinguish between subsection the section
    and subsection the property. This allows to use on child `section definition` as
    subsection of many parent `section definitions`.

    Attributes:
        sub_section: A :class:`Section` or Python class object for a `section class`. This
            will be the child `section definition`. The defining section the child
            `section definition`.

        repeats: A boolean that determines whether this subsection can appear multiple
            times in the parent section.
    '''

    _used_sections: Dict[Section, Set[SubSection]] = {}

    sub_section: Quantity = _placeholder_quantity
    repeats: Quantity = _placeholder_quantity

    def __get__(self, obj, type=None):
        # the class attribute case
        if obj is None:
            return self

        # the object attribute case
        if self.repeats:
            return obj.m_get_sub_sections(self)

        return obj.m_get_sub_section(self, -1)

    def __set__(self, obj, value):
        if obj is None:
            raise NotImplementedError()

        if self.repeats:
            if isinstance(value, (list, set)):
                obj.m_get_sub_sections(self).clear()
                for item in value:
                    obj.m_add_sub_section(self, item)
                return

            if value is not None:
                raise NotImplementedError(
                    'Cannot set a repeating subsection directly, modify the list, e.a. via append.')

            obj.m_get_sub_sections(self).clear()

        elif value is None:
            obj.m_remove_sub_section(self, -1)

        else:
            obj.m_add_sub_section(self, value)

    def __delete__(self, obj):
        raise NotImplementedError('Deleting subsections is not supported.')

    @constraint(warning=False)
    def has_sub_section(self):
        assert self.sub_section is not None, \
            'Each subsection must define the section that is used as subsection via the "sub_section" quantity'
        try:
            assert not isinstance(self.sub_section.m_resolved(), MProxy), 'Cannot resolve "sub_section"'
        except MetainfoReferenceError as e:
            assert False, f'Cannot resolve "sub_section": {str(e)}'

    def _hash(self, regenerate=False) -> _HASH_OBJ:
        if self._cached_hash is not None and not regenerate:
            return self._cached_hash

        self._cached_hash = super(SubSection, self)._hash(regenerate)

        self._cached_hash.update(('T' if self.repeats else 'F').encode('utf-8'))

        for item in itertools.chain(
                self.sub_section.quantities,
                self.sub_section.base_sections,
                self.sub_section.sub_sections,
                self.sub_section.inner_section_definitions):
            if id(self) != id(item):
                self._cached_hash.update(item._hash(regenerate).digest())

        return self._cached_hash


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
            The subsection definitions of this section definition as list of :class:`SubSection`.
            Will be automatically set from the `section class`.

        base_sections:
            A list of `section definitions` (:class:`Section`). By default this definition will
            inherit all quantity and subsection definitions from the given section definitions.
            This behavior might be altered with ``extends_base_section``.

            If there are no base sections to define, you have to use :class:`MSection`.

    The Metainfo supports two inheritance mechanism. By default, it behaves like regular
    Python inheritance and the class inherits all its base
    classes' properties. The other mode (enabled via ``extends_base_section=True``), will
    add all subclass properties to the base-class. This is used throughout the NOMAD metainfo
    to add code-specific metadata to common section definitions. Here is an example:

    .. code-block:: python

        class Method(MSection):
            code_name = Quantity(str)

        class VASPMethod(Method):
            m_def = Section(extends_base_section=True)
            x_vasp_some_incar_parameter = Quantity(str)

        method = Method()
        method.x_vasp_same_incar_parameter = 'value'

    In this example, the section class ``VASPMethod`` defines a section definition that inherits
    from section definition ``Method``. The quantity `x_vasp_some_incar_parameter` will
    be added to `Method` and can be used in regular `Method` instances.

    The following :class:`Section` attributes manipulate the inheritance semantics:

    Attributes:
        extends_base_section:
            If True, this definition must have exactly one ``base_sections``.
            Instead of inheriting properties, the quantity and subsection definitions
            of this section will be added to the base section.

            This allows to add further properties to an existing section definition.
            To use such extension on section instances in a type-safe manner
            :py:func:`MSection.m_as` can be used to cast the base section to the extending
            section.

        extending_sections:
            A list of `section definitions` (:class:`Section`). These are those sections
            that add their properties to this section via :attr:`extends_base_section`.
            This quantity will be set automatically.

        inheriting_sections:
            A list of `section definitions` (:class:`Section`). These are those sections
            that inherit (i.e. are subclasses) of this section.


    Besides defining quantities and subsections, a section definition can also provide
    constraints that are used to validate a section and its quantities and subsections.
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
            with the :func:`constraint` decorator. They can raise :class:`ConstraintViolated`
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

        all_inheriting_section:
            A helper attribute that gives direct and indirect inheriting sections.

        all_properties:
            A helper attribute that gives all properties (subsection and quantity) definitions
            including inherited properties and properties from extending sections as a
            dictionary with names and definitions.

        all_quantities:
            A helper attribute that gives all quantity definition including inherited ones
            and ones from extending sections as a dictionary that maps names (strings)
            to :class:`Quantity`.

        all_sub_sections:
            A helper attribute that gives all subsection definition including inherited ones
            and ones from extending sections as a dictionary that maps names (strings)
            to :class:`SubSection`.

        all_sub_sections_by_section:
            A helper attribute that gives all subsection definition including inherited ones
            and ones from extending sections as a dictionary that maps section classes
            (i.e. Python class objects) to lists of :class:`SubSection`.

        all_aliases:
            A helper attribute that gives all aliases for all properties including
            inherited properties and properties form extending sections as a
            dictionary with aliases and the definitions.

        all_inner_section_definitions:
            A helper attribute that gives all inner_section_definitions including
            their aliases by name.

        path: Shortest path from a root section to this section. This is not the path
            in the metainfo schema (`m_path`) but an archive path in potential data.

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

    quantities: SubSection = None
    sub_sections: SubSection = None
    inner_section_definitions: SubSection = None

    base_sections: Quantity = _placeholder_quantity
    extending_sections: Quantity = _placeholder_quantity
    extends_base_section: Quantity = _placeholder_quantity
    inheriting_sections: Quantity = _placeholder_quantity
    constraints: Quantity = _placeholder_quantity
    event_handlers: Quantity = _placeholder_quantity

    inherited_sections: Quantity = _placeholder_quantity
    all_base_sections: Quantity = _placeholder_quantity
    all_inheriting_sections: Quantity = _placeholder_quantity
    all_properties: Quantity = _placeholder_quantity
    all_quantities: Quantity = _placeholder_quantity
    all_sub_sections: Quantity = _placeholder_quantity
    all_sub_sections_by_section: Quantity = _placeholder_quantity
    all_aliases: Quantity = _placeholder_quantity
    all_inner_section_definitions: Quantity = _placeholder_quantity
    has_variable_names: Quantity = _placeholder_quantity
    path: Quantity = _placeholder_quantity

    def __init__(self, *args, validate: bool = True, **kwargs):
        self._section_cls: Type[MSection] = None  # type: ignore

        super().__init__(*args, **kwargs)
        self.validate = validate

    def m_validate(self) -> Tuple[List[str], List[str]]:
        if not self.validate:
            return [], []

        return super().m_validate()

    @property
    def section_cls(self) -> Type[MSection]:
        if self._section_cls is None:
            # set a temporary to avoid endless recursion
            self._section_cls = type(self.name, (MSection,), dict(do_init=False))

            # Create a section class if this does not exist. This happens if the section
            # is not created through a class definition.
            attrs = {prop.name: prop for prop in itertools.chain(self.quantities, self.sub_sections)}

            for name, inner_section_def in self.all_inner_section_definitions.items():
                attrs[name] = inner_section_def.section_cls

            attrs.update(m_def=self, do_init=False)
            if len(self.base_sections) > 0:
                bases = tuple(base_section.section_cls for base_section in self.base_sections)
            else:
                bases = (MSection,)
            self._section_cls = type(self.name, bases, attrs)

        return self._section_cls

    def __init_metainfo__(self):
        super().__init_metainfo__()

        # Init extending_sections
        if self.extends_base_section:
            base_sections_count = len(self.base_sections)
            if base_sections_count == 0:
                raise MetainfoError(
                    f'Section {self} extend the base section, but has no base section.')

            if base_sections_count > 1:
                raise MetainfoError(
                    f'Section {self} extend the base section, but has more than one base section.')

            base_section = self.base_sections[0]
            for name, attr in self.section_cls.__dict__.items():
                if isinstance(attr, Property):
                    setattr(base_section.section_cls, name, attr)

            base_section.extending_sections = base_section.extending_sections + [self]

        # Init inheriting_sections
        if not self.extends_base_section:
            for base_section in self.base_sections:
                base_section.inheriting_sections = base_section.inheriting_sections + [self]

        # Transfer properties of inherited and overwritten property definitions that
        # have not been overwritten
        inherited_properties: Dict[str, Property] = dict()
        for base_section in self.all_base_sections:
            inherited_properties.update(**base_section.all_properties)

        for property in itertools.chain(self.quantities, self.sub_sections):
            inherited_property = inherited_properties.get(property.name)
            if inherited_property is None:
                continue

            for m_quantity in property.m_def.all_quantities.values():
                if not property.m_is_set(m_quantity) and inherited_property.m_is_set(m_quantity):
                    property.m_set(m_quantity, inherited_property.m_get(m_quantity))

    @constraint
    def unique_names(self):
        names: Set[str] = set()
        for base in self.extending_sections:
            for quantity in itertools.chain(base.quantities, base.sub_sections):
                for alias in quantity.aliases:
                    names.add(alias)
                names.add(quantity.name)

        for definition in itertools.chain(self.quantities, self.sub_sections):
            assert definition.name not in names, \
                f'All names in a section must be unique. ' \
                f'Name {definition.name} of {definition} in {definition.m_parent} already exists in {self}.'
            names.add(definition.name)
            for alias in definition.aliases:
                assert alias not in names, \
                    f'All names (incl. aliases) in a section must be unique. ' \
                    f'Alias {alias} of {definition} in {definition.m_parent} already exists in {self}.'
                names.add(alias)

    @constraint
    def resolved_base_sections(self):
        for base_section in self.base_sections:
            try:
                base_section.m_resolved()
            except MetainfoReferenceError as e:
                assert False, f'Cannot resolve base_section: {str(e)}'

    @classmethod
    def m_from_dict(cls, data: Dict[str, Any], **kwargs):
        for list_quantity in [Section.quantities, Section.sub_sections, Section.inner_section_definitions]:
            alias, potential_dict_value = list_quantity.get_from_dict(data)
            if alias:
                data[alias] = dict_to_named_list(potential_dict_value)

        _, sub_sections = Section.sub_sections.get_from_dict(data, [])
        for sub_section in sub_sections:
            alias, section_def = SubSection.sub_section.get_from_dict(sub_section)
            if not section_def or isinstance(section_def, str):
                continue

            inner_sections_alias, inner_sections = Section.inner_section_definitions.get_from_dict(data)
            if not inner_sections_alias:
                inner_sections = []
                data['inner_section_definitions'] = inner_sections

            inner_sections.append(section_def)
            name = sub_section.get('name')
            if name:
                name = name.title().replace('_', '')
            section_def['name'] = name
            sub_section[alias] = name

        if 'base_section' in data:
            data['base_sections'] = [data.pop('base_section')]

        return super(Section, cls).m_from_dict(data, **kwargs)

    def _hash(self, regenerate=False) -> _HASH_OBJ:
        if self._cached_hash is not None and not regenerate:
            return self._cached_hash

        self._cached_hash = super(Section, self)._hash(regenerate)
        self._cached_hash.update(('T' if self.extends_base_section else 'F').encode('utf-8'))

        for item in itertools.chain(
                self.quantities,
                self.base_sections,
                self.sub_sections,
                self.inner_section_definitions):
            if id(self) != id(item):
                self._cached_hash.update(item._hash(regenerate).digest())

        return self._cached_hash


class Package(Definition):
    ''' Packages organize metainfo definitions alongside Python modules

    Each Python module with metainfo Definition (explicitly or implicitly) has a member
    ``m_package`` with an instance of this class. Definitions (categories, sections) in
    Python modules are automatically added to the module's :class:`Package`.
    Packages are not nested and rather have the fully qualified Python module name as
    name.

    This allows to inspect all definitions in a Python module and automatically puts
    module name and docstring as :class:`Package` name and description.

    Besides the regular :class:`Definition` attributes, packages can have the following
    attributes:

    Attributes:
        section_definitions: All `section definitions` in this package as :class:`Section`
            objects.

        category_definitions: All `category definitions` in this package as :class:`Category`
            objects.

        all_definitions: A helper attribute that provides all section and category definitions
            by name and aliases.
    '''

    section_definitions: SubSection = None
    category_definitions: SubSection = None

    all_definitions: Quantity = _placeholder_quantity
    dependencies: Quantity = _placeholder_quantity

    registry: Dict[str, Package] = {}
    ''' A static member that holds all currently known packages. '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.errors, self.warnings = [], []
        self.archive = None
        self.entry_id_based_name = None

    def __init_metainfo__(self):
        super().__init_metainfo__()

        Package.registry[self.name] = self

        # access potential SectionProxies to resolve them
        for content in self.m_all_contents():
            if isinstance(content, Quantity):
                if isinstance(content.type, MProxy):
                    content.type.m_proxy_resolve()
            elif isinstance(content, SubSection):
                target = content.sub_section
                if isinstance(target, MProxy):
                    target = target.m_proxy_resolve()
                SubSection._used_sections.setdefault(target, []).append(content)
            elif isinstance(content, Section):
                for base_section in content.base_sections:
                    if isinstance(base_section, MProxy):
                        base_section.m_proxy_resolve()

        # validate
        if is_bootstrapping:
            return

        self.errors, self.warnings = self.m_all_validate()
        if len(self.errors) > 0:
            raise MetainfoError(
                f'One constraint was violated: {str(self.errors[0]).strip(".")} '
                f'(there are {len(self.errors) - 1:d} more violations)')

    @staticmethod
    def from_module(module_name: str):
        module = sys.modules[module_name]

        pkg: Package = getattr(module, 'm_package', None)
        if pkg is None:
            pkg = Package()
            setattr(module, 'm_package', pkg)

        pkg.name = module_name
        if pkg.description is None and module.__doc__ is not None:
            pkg.description = inspect.cleandoc(module.__doc__).strip()

        return pkg

    @classmethod
    def m_from_dict(cls, data: Dict[str, Any], **kwargs):
        for list_quantity in [Package.section_definitions, Package.category_definitions]:
            alias, potential_dict_value = list_quantity.get_from_dict(data)
            if alias:
                data[alias] = dict_to_named_list(potential_dict_value)

        return super(Package, cls).m_from_dict(data, **kwargs)

    def qualified_name(self):
        # packages loaded from files have a hot qualified name based on entry id
        # this name is not serialized which causes '*' name when reloaded from cold
        # we store this name in a `str` and it will be reloaded from cold
        # see Context.resolve_section_definition()
        if self.entry_id_based_name:
            return self.entry_id_based_name

        if self.archive:
            # If the package was defined within a regular uploaded archive file, we
            # use its id, which is a globally unique identifier for the package.
            if self.archive.metadata and self.archive.metadata.entry_id:
                self.entry_id_based_name = f'entry_id:{self.archive.metadata.entry_id}'
            else:
                self.entry_id_based_name = f'entry_id:*'

            return self.entry_id_based_name

        return super().qualified_name()

    def normalize(self, archive, logger=None):
        if archive.definitions == self and archive.metadata:
            if archive.metadata.entry_name is None and self.name and self.name != '*':
                archive.metadata.entry_name = self.name

    def _hash(self, regenerate=False) -> _HASH_OBJ:
        if self._cached_hash is not None and not regenerate:
            return self._cached_hash

        self._cached_hash = super(Package, self)._hash(regenerate)

        for item in self.section_definitions:  # pylint: disable=not-an-iterable
            if id(self) != id(item):
                self._cached_hash.update(item._hash(regenerate).digest())

        return self._cached_hash


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


Section.m_def = Section(name='Section')
Section.m_def.m_def = Section.m_def
Section.m_def._section_cls = Section

Definition.m_def = Section(name='Definition')
Attribute.m_def = Section(name='Attribute')
Property.m_def = Section(name='Property')
Quantity.m_def = Section(name='Quantity')
SubSection.m_def = Section(name='SubSection')
Category.m_def = Section(name='Category')
Package.m_def = Section(name='Package')

Attribute.type = DirectQuantity(type=QuantityType, name='type')
Attribute.shape = DirectQuantity(type=Dimension, shape=['0..*'], name='shape', default=[])

Definition.name = DirectQuantity(type=str, name='name')
Definition.label = DirectQuantity(type=str, name='name')
Definition.description = Quantity(type=str, name='description')
Definition.links = Quantity(type=str, shape=['0..*'], name='links')
Definition.categories = Quantity(
    type=Reference(Category.m_def), shape=['0..*'], default=[], name='categories')
Definition.deprecated = Quantity(type=str, name='deprecated')
Definition.aliases = Quantity(type=str, shape=['0..*'], default=[], name='aliases')
Definition.variable = Quantity(type=bool, name='variable', default=False)
Definition.more = Quantity(type=JSON, name='more', default={})
Definition.attributes = SubSection(sub_section=Attribute.m_def, name='attributes', repeats=True)


@derived(cached=True, virtual=True)  # Virtual has to be set manually, due to bootstrapping hen-egg problems
def all_attributes(self: Property) -> Dict[str, Attribute]:
    result: Dict[str, Attribute] = {}
    for definition in self.attributes:
        result[definition.name] = definition

    return result


Definition.all_attributes = all_attributes

Section.quantities = SubSection(
    sub_section=Quantity.m_def, name='quantities', repeats=True)
Section.sub_sections = SubSection(
    sub_section=SubSection.m_def, name='sub_sections', repeats=True)
Section.inner_section_definitions = SubSection(
    sub_section=Section.m_def, name='inner_section_definitions', repeats=True,
    aliases=['inner_section_defs', 'section_defs', 'inner_sections', 'sections'])

Section.base_sections = Quantity(
    type=SectionReference, shape=['0..*'], default=[], name='base_sections')
Section.extending_sections = Quantity(
    type=SectionReference, shape=['0..*'], default=[], name='extending_sections')
Section.extends_base_section = Quantity(type=bool, default=False, name='extends_base_section')
Section.inheriting_sections = Quantity(
    type=SectionReference, shape=['0..*'], default=[], name='inheriting_sections', virtual=True)
Section.constraints = Quantity(type=str, shape=['0..*'], default=[], name='constraints')
Section.event_handlers = Quantity(
    type=Callable, shape=['0..*'], name='event_handlers', virtual=True, default=[])


@derived(cached=True)
def inherited_sections(self) -> List[Section]:
    result: List[Section] = []
    for base_section in self.all_base_sections:
        result.append(base_section)
    for extending_section in self.extending_sections:
        result.append(extending_section)
    result.append(self)
    return result


@derived(cached=False)
def all_base_sections(self) -> List[Section]:
    result: List[Section] = []
    for base_section in self.base_sections:
        if isinstance(base_section, SectionProxy):
            # In some reference resolution contexts, it is important to reevaluate later
            self.m_mod_count += 1
            continue
        for base_base_section in base_section.all_base_sections:
            if isinstance(base_base_section, SectionProxy):
                # In some reference resolution contexts, it is important to reevaluate later
                self.m_mod_count += 1
                continue
            result.append(base_base_section)
        result.append(base_section)
    return result


@derived(cached=True)
def all_inheriting_sections(self) -> List[Section]:
    result: Set[Section] = set()
    for inheriting_section in self.inheriting_sections:
        if isinstance(inheriting_section, SectionProxy):
            # In some reference resolution contexts, it is important to reevaluate later
            self.m_mod_count += 1
            continue
        for inheriting_inheriting_section in inheriting_section.all_inheriting_sections:
            if isinstance(inheriting_inheriting_section, SectionProxy):
                # In some reference resolution contexts, it is important to reevaluate later
                self.m_mod_count += 1
                continue
            result.add(inheriting_inheriting_section)
        result.add(inheriting_section)
    return list(result)


@derived(cached=True)
def all_properties(self) -> Dict[str, Union[SubSection, Quantity]]:
    result: Dict[str, Union[SubSection, Quantity]] = dict()
    for section in self.inherited_sections:
        for definition in itertools.chain(section.quantities, section.sub_sections):
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
        for definition in itertools.chain(section.quantities, section.sub_sections):
            for alias in definition.aliases:
                result[alias] = definition
            result[definition.name] = definition

    return result


@derived(cached=True)
def all_inner_section_definitions(self) -> Dict[str, Section]:
    result: Dict[str, Section] = dict()
    for base_section_or_self in self.all_base_sections + [self]:  # pylint: disable=not-an-iterable
        for section in base_section_or_self.inner_section_definitions:  # pylint: disable=not-an-iterable
            if section.name:
                result[section.name] = section
            for alias in section.aliases:
                result[alias] = section

    return result


@derived(cached=True)
def has_variable_names(self) -> bool:
    return any(value.variable for value in self.all_properties.values())


@derived(cached=True)
def section_path(self) -> str:
    used_in_sub_sections: List[SubSection] = SubSection._used_sections.get(self, [])  # type: ignore
    if len(used_in_sub_sections) == 0:
        return None if self.name == 'EntryArchive' else '__no_archive_path__'

    if len(used_in_sub_sections) > 1:
        return '__ambiguous__'

    parent_section = used_in_sub_sections[0].m_parent
    parent_path = parent_section.path

    if parent_path is None:
        return used_in_sub_sections[0].name

    if parent_path.startswith('__'):
        return parent_path

    return f'{parent_path}.{used_in_sub_sections[0].name}'


Section.inherited_sections = inherited_sections
Section.all_base_sections = all_base_sections
Section.all_inheriting_sections = all_inheriting_sections
Section.all_properties = all_properties
Section.all_quantities = all_quantities
Section.all_sub_sections = all_sub_sections
Section.all_sub_sections_by_section = all_sub_sections_by_section
Section.all_aliases = all_aliases
Section.all_inner_section_definitions = all_inner_section_definitions
Section.has_variable_names = has_variable_names
Section.path = section_path

SubSection.repeats = Quantity(type=bool, name='repeats', default=False)

SubSection.sub_section = Quantity(
    type=SectionReference, name='sub_section',
    aliases=['section_definition', 'section_def', 'section'])

Quantity.m_def._section_cls = Quantity
Quantity.type = DirectQuantity(type=QuantityType, name='type')
Quantity.shape = DirectQuantity(type=Dimension, shape=['0..*'], name='shape', default=[])
Quantity.unit = Quantity(type=Unit, name='unit')
Quantity.dimensionality = DirectQuantity(type=str, name='dimensionality')
Quantity.default = DirectQuantity(type=Any, default=None, name='default')
Quantity.derived = DirectQuantity(type=Callable, default=None, name='derived', virtual=True)
Quantity.virtual = DirectQuantity(type=bool, default=False, name='virtual')
Quantity.is_scalar = Quantity(
    type=bool, name='is_scalar', derived=lambda quantity: len(quantity.shape) == 0)
Quantity.use_full_storage = Quantity(
    type=bool, name='use_full_storage',
    derived=lambda quantity: quantity.repeats or quantity.variable or len(quantity.attributes) > 0)
Quantity.flexible_unit = Quantity(type=bool, name='flexible_unit', default=False)
Quantity.repeats = Quantity(type=bool, name='repeats', default=False)
Quantity.cached = Quantity(type=bool, name='cached', default=False)

Package.section_definitions = SubSection(
    sub_section=Section.m_def, name='section_definitions', repeats=True,
    aliases=['section_defs', 'sections'])

Package.category_definitions = SubSection(
    sub_section=Category.m_def, name='category_definitions', repeats=True,
    aliases=['category_defs'])


@derived(cached=True)
def all_definitions(self):
    result: Dict[str, Definition] = dict()
    for sub_section_def in [Package.section_definitions, Package.category_definitions]:
        for definition in self.m_get_sub_sections(sub_section_def):
            result[definition.name] = definition
            for alias in definition.aliases:
                result[alias] = definition
    return result


@derived(cached=True)
def dependencies(self):
    '''
    All packages which have definitions that definitions from this package need. Being
    'needed' includes categories, base sections, and referenced definitions.
    '''
    result = set()
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
            if dependency not in result:
                result.add(dependency)
                more_dependencies.extend(dependency.dependencies)

    return result


Package.all_definitions = all_definitions
Package.dependencies = dependencies

is_bootstrapping = False

Definition.__init_cls__()
Attribute.__init_cls__()
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
            for definition in self.all_definitions_by_name.get(name, [])  # pylint: disable=no-member
            if isinstance(definition, section_cls)
            if not (isinstance(definition, Section) and definition.extends_base_section)
            if filter is None or filter(definition)]

    def resolve_definition(
            self, name, section_cls: Type[MSectionBound],
            filter: TypingCallable[[MSection], bool] = None) -> MSectionBound:

        defs = self.resolve_definitions(name, section_cls, filter=filter)
        if len(defs) == 1:
            return defs[0]

        if len(defs) > 1:
            raise KeyError(f'Could not uniquely identify {name}, candidates are {defs}')

        raise KeyError(f'Could not resolve {name}')


class AnnotationModel(Annotation, BaseModel):
    '''
    Base class for defining annotation models. Annotations used with simple dict-based
    values, can be validated by defining and registering a formal pydantic-based
    model.
    '''

    m_definition: Definition = Field(
        None, description='The definition that this annotation is annotating.')

    m_error: str = Field(None, description='Holds a potential validation error.')

    m_registry: ClassVar[Dict[str, Type['AnnotationModel']]] = {}
    ''' A static member that holds all currently known annotations with pydantic model. '''

    def m_to_dict(self, *args, **kwargs):
        return self.dict(exclude_unset=True)

    class Config:
        fields = {
            'm_definition': {
                'exclude': True,
            },
            'm_error': {
                'exclude': True
            }
        }

        validate_assignment = True
        arbitrary_types_allowed = True


AnnotationModel.update_forward_refs()
