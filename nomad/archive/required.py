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

import re
from typing import Any, Dict, List, Tuple, Union, cast

from nomad import utils
from nomad.metainfo import Definition, Quantity, QuantityReference, Reference, Section, SubSection
from .query import ArchiveQueryError
from .storage import ArchiveError, ArchiveList, ArchiveObject, ArchiveReader


def _fix_index(index, length):
    if index is None:
        return index

    if index < 0:
        return max(-length, index)

    return min(length, index)


class RequiredValidationError(Exception):
    def __init__(self, msg, loc):
        super().__init__(msg)
        self.msg = msg
        self.loc = loc


# noinspection PyTypeChecker
class RequiredReader:
    '''
    Clients can read only the required parts of an archive. They specify the required
    parts with a required specification like this.

    .. code-block:: python
        {
            "results": "include-resolved",  // contains all the results and all its references
                                            // resolved
            "workflow": {
                "calculation_result_ref": {  // implicitly resolves the reference
                    "calculation_to_system_ref": "*"  // shorthand for include
                    "calculation_to_method_ref": "include"  // resolves and includes the target,
                                                            // but no references in the target
                                                            // are resolved.
                    "*": "*",  // includes everything else ...
                    "section_eigenvalues": "exclude"  // ... but explicitly excluded parts
                }
            }
        }

    The structure has to adhere to metainfo definitions of an archive's subsections and
    quantities. At each point in the specification, children can be replaced with certain
    directives.

    The different directives are:
        * include ('*' is an alias), includes whole subtree, resolves on a ref
        * include-resolved, includes whole subtree, resolves all refs recursively
        * exclude (in combination with wildcard keys), omits it from the result

    This class allows to keep a requirement specification and use it to read with it
    from given upload files and entry ids.

    Attributes:
        - required: The requirement specification as a python dictionary or directive string.
    '''

    _query_archive_key_pattern = re.compile(r'^([\s\w\-]+)(\[([-?0-9]*)(:([-?0-9]*))?])?$')

    def __init__(self, required: Union[dict, str], root_section_def: Section = None,
                 resolve_inplace: bool = False):
        if root_section_def is None:
            from nomad import datamodel
            self._root_section_def = datamodel.EntryArchive.m_def
        else:
            self._root_section_def = root_section_def

        self._result_root: dict = None
        self._archive_root: dict = None  # it is actually an ArchiveReader, but we use it as dict

        self._resolve_inplace = resolve_inplace
        self._required = self.validate(required, is_root=True)

    @staticmethod
    def __parse_required_key(key: str) -> Tuple[str, Union[Tuple[int, int], int]]:
        '''
        Returns:
            key, None
            key, n
            key, (0, n)
            key, (0, None)
        '''
        key = key.strip()

        match = RequiredReader._query_archive_key_pattern.match(key)
        if not match:
            raise Exception(f'invalid key format: {key}')

        key = match.group(1)

        index: Union[Tuple[int, int], int] = None

        # check if we have indices
        if match.group(2) is not None:
            group = match.group(3)
            first_index = None if group == '' else int(group)

            if match.group(4) is not None:
                group = match.group(5)
                last_index = None if group == '' else int(group)

                index = (0 if first_index is None else first_index, last_index)
            else:
                index = first_index  # one item

        return key, index

    def __to_son(self, data):
        if isinstance(data, (ArchiveList, List)):
            return [self.__to_son(item) for item in data]

        if isinstance(data, ArchiveObject):
            return data.to_dict()

        return data

    def __resolve_ref(self, required: dict, path: str) -> Union[dict, str]:
        # The archive item is a reference, the required is still a dict, the references
        # This is a simplified version of the metainfo implementation (m_resolve).
        # It implements the same semantics, but does not apply checks.
        # TODO the metainfo should also provide this implementation

        # resolve to archive_item
        if not path.startswith('/'):
            # TODO support custom reference resolution, e.g. user_id based
            return path

        resolved = self._archive_root
        path_stack = path.strip('/').split('/')
        path_stack.reverse()

        try:
            while len(path_stack) > 0:
                prop = path_stack.pop()
                resolved = resolved[prop]
                if len(path_stack) > 0 and path_stack[-1].isdigit():
                    resolved = resolved[int(path_stack.pop())]
        except Exception:
            raise ArchiveError('could not resolve reference')

        # apply required to resolved archive_item
        if isinstance(required['_def'], Quantity):
            resolved_result = self.__to_son(resolved)
        else:
            resolved_result = self.__apply_required(required, resolved)

        # return or add to root depending on self.resolve_inplace
        if self._resolve_inplace:
            return resolved_result

        def _setdefault(target: Union[dict, list], local_prop, value_type: type):
            if isinstance(target, list):
                if target[local_prop] is None:
                    target[local_prop] = value_type()
                return target[local_prop]

            if local_prop not in target:
                target[local_prop] = value_type()

            return target[local_prop]

        path_stack = path.strip('/').split('/')
        path_stack.reverse()
        target_container: Union[dict, list] = self._result_root
        prop_or_index: Union[str, int] = None
        while len(path_stack) > 0:
            if prop_or_index is not None:
                target_container = _setdefault(target_container, prop_or_index, dict)

            prop_or_index = path_stack.pop()
            if len(path_stack) > 0 and path_stack[-1].isdigit():
                target_list = _setdefault(target_container, prop_or_index, list)
                index = int(path_stack.pop())
                for _ in range(len(target_list), index + 1):
                    target_list.append(None)
                target_container = target_list
                prop_or_index = index

        target_container[prop_or_index] = resolved_result  # type: ignore
        return path

    def __resolve_refs(self, archive: dict, definition: Definition) -> dict:
        ''' Resolves all references in archive. '''
        if isinstance(definition, Quantity):
            # it's a quantity ref, the archive is already resolved
            return self.__to_son(archive[definition.name])

        # it's a section ref
        section_def = cast(Section, definition)
        archive = self.__to_son(archive)
        result = {}
        for prop in archive:
            value = archive[prop]
            prop_def = section_def.all_properties[prop]
            if isinstance(prop_def, SubSection):
                def handle_item(v):
                    return self.__resolve_refs(v, prop_def.sub_section.m_resolved())
            elif isinstance(prop_def.type, Reference):
                if isinstance(prop_def.type, QuantityReference):
                    target_def = prop_def.type.target_quantity_def.m_resolved()
                else:
                    target_def = prop_def.type.target_section_def.m_resolved()

                required = dict(_directive='include-resolved', _def=target_def, _ref=prop_def.type)

                def handle_item(v):
                    return self.__resolve_ref(required, v)
            else:
                def handle_item(v):
                    return v

            if isinstance(value, (list, ArchiveList)):
                result[prop] = [handle_item(item) for item in value]
            else:
                result[prop] = handle_item(value)

        return result

    def __apply_required(self, required: dict, archive_item: Union[dict, str]) -> Union[Dict, str]:
        if archive_item is None:
            return None

        directive = required.get('_directive')
        if directive is not None:
            if directive == 'include-resolved':
                if isinstance(archive_item, str):
                    return self.__resolve_ref(required, archive_item)

                return self.__resolve_refs(archive_item, required['_def'])

            if directive in ['*', 'include']:
                return self.__to_son(archive_item)

            raise ArchiveQueryError(f'unknown directive {required}')

        if isinstance(archive_item, str):
            return self.__resolve_ref(required, archive_item)

        result: dict = {}
        for key, val in required.items():
            if key.startswith('_'):
                continue

            prop, index = val['_prop'], val['_index']

            try:
                archive_child = archive_item[prop]
                is_list = isinstance(archive_child, (ArchiveList, list))

                if index is not None and not is_list:
                    raise ArchiveQueryError(f'cannot use list key on none list {prop}')

                if index is None and is_list:
                    index = (0, None)

                if index is not None:
                    # is_list === True
                    # index = n
                    # index = (0, n)
                    # index = (0, None)
                    length = len(archive_child)

                    if isinstance(index, tuple):
                        index = (_fix_index(index[0], length), _fix_index(index[1], length))
                        if index[0] == index[1]:
                            archive_child = [archive_child[index[0]]]
                        else:
                            archive_child = archive_child[index[0]: index[1]]
                    else:
                        archive_child = [archive_child[_fix_index(index, length)]]

                if isinstance(archive_child, (ArchiveList, list)):
                    result[prop] = [self.__apply_required(val, item) for item in archive_child]
                else:
                    result[prop] = self.__apply_required(val, archive_child)

            except (KeyError, IndexError):
                continue

        return result

    def validate(self, required: Union[str, dict], definition: Definition = None, loc: list = None,
                 is_root=False) -> dict:
        '''
        Validates the required specification of this instance. It will replace all
        string directives with dicts. Those will have keys `_def` and `_directive`. It
        will add a key `_def` to all dicts. The `_def` will be the respective metainfo
        definition. It will also add key `_ref`. It will be None or contain a
        Reference instance, if the definition is a reference target.

        This method will raise an exception (:class:`RequiredValidationError`) to denote
        any mismatches between the required specification and the metainfo. Also, to denote
        misused directives, bad structure, etc.

        raises:
            - RequiredValidationError
        '''
        if is_root and isinstance(required, dict):
            resolve_inplace = required.get('resolve-inplace', None)
            if isinstance(resolve_inplace, bool):
                self._resolve_inplace = resolve_inplace
            elif resolve_inplace is not None:
                raise RequiredValidationError('resolve-inplace is not a bool', ['resolve-inplace'])

        if definition is None:
            definition = self._root_section_def
        if loc is None:
            loc = []

        # replace definition with the target definition, if its reference or subsection
        reference = None
        if isinstance(definition, Quantity):
            if isinstance(definition.type, Reference):
                reference = definition.type
                if isinstance(definition.type, QuantityReference):
                    definition = definition.type.target_quantity_def.m_resolved()
                else:
                    definition = definition.type.target_section_def.m_resolved()
        elif isinstance(definition, SubSection):
            definition = definition.sub_section.m_resolved()

        is_directive = isinstance(required, str)
        if not isinstance(definition, Section) and not is_directive:
            raise RequiredValidationError(f'{definition.name} is not a section or reference', loc)

        if is_directive:
            if required not in ['*', 'include', 'exclude', 'include-resolved']:
                raise RequiredValidationError(f'{required} is not a valid directive', loc)

            # TODO support 'exclude'
            if required == 'exclude':
                raise RequiredValidationError('exclude is not supported yet', loc)

            return dict(_def=definition, _directive=required, _ref=reference)

        result: Dict[str, Any] = dict(_def=definition, _ref=reference)
        for key, value in cast(dict, required).items():
            if key == 'resolve-inplace':
                continue

            loc.append(key)
            try:
                prop, index = self.__parse_required_key(key)
            except Exception:
                raise RequiredValidationError(f'invalid key format {key}', loc)

            if prop == '*':
                # TODO support wildcards
                raise RequiredValidationError('wildcard (*) keys are not supported yet', loc)

            try:
                prop_def = cast(Section, definition).all_properties[prop]
            except KeyError:
                raise RequiredValidationError(f'{definition.name} has not property {prop}', loc)

            result[key] = self.validate(value, prop_def, loc)
            result[key].update(_prop=prop, _index=index)
            loc.pop()

        return result

    def read(self, archive_reader: ArchiveReader, entry_id: str) -> dict:
        '''
        Reads the archive of the given entry id from the given archive reader and applies
        the instance's requirement specification.
        '''
        assert self._result_root is None and self._archive_root is None, \
            'instance cannot be used concurrently for multiple reads at the same time'

        result: dict = {}

        self._archive_root = archive_reader[utils.adjust_uuid_size(entry_id)]
        self._result_root = result

        result.update(**cast(dict, self.__apply_required(self._required, self._archive_root)))

        self._result_root = None
        self._archive_root = None

        return result
