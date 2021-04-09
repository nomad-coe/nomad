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

from typing import Any, Tuple, Dict, Callable, Union, List
from io import BytesIO
import re

from nomad import utils

from .storage import ArchiveReader, ArchiveList, ArchiveObject


class ArchiveQueryError(Exception):
    '''
    An error that indicates that an archive query is either not valid or does not fit to
    the queried archive.
    '''
    pass


__query_archive_key_pattern = re.compile(r'^([\s\w\-]+)(\[([-?0-9]*)(:([-?0-9]*))?\])?$')


def query_archive(
        f_or_archive_reader: Union[str, ArchiveReader, BytesIO], query_dict: dict,
        **kwargs) -> Dict:
    '''
    Takes an open msg-pack based archive (either as str, reader, or BytesIO) and returns
    the archive as JSON serializable dictionary filtered based on the given required
    (query_dict) specification.

    Required example with some "required" (query_dict) features:

    .. code-block:: json
        {
            "results": "recursive-resolve-inplace",
            "workflow": {
                "final": "resolve"
            },
            "simulation": {
                "scc[-1]": {
                    "method": "recursive-resolve"
                    "system": {
                        "symmetry": "exclude",
                        "*": "include"
                    }
                    "dos": "recursive-include"
                }
            }
        }

    The different directives are:
        * include ('*' is an alias), includes whole subtree, resolves on a ref
        * recursive-resolve, includes whole subtree, resolves all refs recursively
        * recursive-resolve-inplace, includes whole subtree, resolves all refs recursively
            and replaces the ref with the resolved data
        * exclude (in combination with wildcard keys), replaces the value with null
    '''

    def _to_son(data):
        if isinstance(data, (ArchiveList, List)):
            data = [_to_son(item) for item in data]

        elif isinstance(data, ArchiveObject):
            data = data.to_dict()

        return data

    def _load_data(query_dict: Dict[str, Any], archive_item: ArchiveObject) -> Dict:
        query_dict_with_fixed_ids = {
            utils.adjust_uuid_size(key): value for key, value in query_dict.items()}
        return filter_archive(query_dict_with_fixed_ids, archive_item, transform=_to_son)

    if isinstance(f_or_archive_reader, ArchiveReader):
        return _load_data(query_dict, f_or_archive_reader)

    elif isinstance(f_or_archive_reader, (BytesIO, str)):
        with ArchiveReader(f_or_archive_reader, **kwargs) as archive:
            return _load_data(query_dict, archive)

    else:
        raise TypeError('%s is neither a file-like nor ArchiveReader' % f_or_archive_reader)


def filter_archive(
        required: Union[str, Dict[str, Any]], archive_item: Union[Dict, ArchiveObject, str],
        transform: Callable, result_root: Dict = None, resolve_inplace: bool = False) -> Dict:

    def _fix_index(index, length):
        if index is None:
            return index
        if index < 0:
            return max(-(length), index)
        else:
            return min(length, index)

    if archive_item is None:
        return None

    if isinstance(required, str):
        if required == 'exclude':
            return None

        if required == 'resolve':
            # TODO this requires to reflect on the definition to determine what are references!
            pass

        elif required in ['*', 'include']:
            pass

        else:
            raise ArchiveQueryError(f'unknown directive {required}')

        return transform(archive_item)

    elif not isinstance(required, dict):
        raise ArchiveQueryError('a value in required is neither dict not string directive')

    if isinstance(archive_item, str):
        # The archive item is a reference, the required is still a dict, the references
        # needs to be resolved
        # TODO
        raise ArchiveQueryError(f'resolving references in non partial archives is not yet implemented')

    result: Dict[str, Any] = {}
    for key, val in required.items():
        key = key.strip()

        # process array indices
        match = __query_archive_key_pattern.match(key)
        index: Union[Tuple[int, int], int] = None
        if match:
            key = match.group(1)

            # check if we have indices
            if match.group(2) is not None:
                first_index, last_index = None, None
                group = match.group(3)
                first_index = None if group == '' else int(group)

                if match.group(4) is not None:
                    group = match.group(5)
                    last_index = None if group == '' else int(group)
                    index = (0 if first_index is None else first_index, last_index)

                else:
                    index = first_index  # one item

            else:
                index = None

        elif key == '*':
            # TODO
            raise ArchiveQueryError('key wildcards not yet implemented')
        else:
            raise ArchiveQueryError('invalid key format: %s' % key)

        try:
            archive_child = archive_item[key]
            is_list = isinstance(archive_child, (ArchiveList, list))

            if index is None and is_list:
                index = (0, None)

            elif index is not None and not is_list:
                raise ArchiveQueryError('cannot use list key on none list %s' % key)

            if index is None:
                pass
            else:
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
                result[key] = [
                    filter_archive(val, item, transform=transform)
                    for item in archive_child]
            else:
                result[key] = filter_archive(val, archive_child, transform=transform)

        except (KeyError, IndexError):
            continue

    return result
