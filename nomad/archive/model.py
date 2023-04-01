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

import functools
import re
from enum import Enum
from hashlib import sha1
from typing import FrozenSet, Optional, Union, Tuple, List, Dict

from pydantic import BaseModel, Field, Extra, ValidationError, validator

from nomad.app.v1.models import MetadataPagination, WithQuery
from nomad.app.v1.routers.datasets import DatasetPagination
from nomad.app.v1.routers.uploads import UploadProcDataQuery, UploadProcDataPagination


class DatasetQuery(BaseModel):
    dataset_id: str = Field(None)
    dataset_name: str = Field(None)
    user_id: List[str] = Field(None)
    dataset_type: str = Field(None)
    doi: str = Field(None)
    prefix: str = Field(None)


class DirectiveType(Enum):
    plain = 'plain'
    resolved = 'resolved'

    def __repr__(self):
        return self.value


class DefinitionType(Enum):
    specialized_only = 'specialized_only'
    both = 'both'
    none = 'none'

    def __repr__(self):
        return self.value


class ResolveType(Enum):
    upload = 'upload'
    user = 'user'
    dataset = 'dataset'
    entry = 'entry'

    def __repr__(self):
        return self.value


class RequestConfig(BaseModel):
    '''
    A class to represent the query configuration.
    An instance of `RequestConfig` shall be attached to each required field.
    The `RequestConfig` is used to determine the following.
        1. Whether the field should be included/excluded.
        2. For reference, whether the reference should be resolved, and how to resolve it.
    Each field can be handled differently.
    '''
    property_name: str = Field(None, description='''
        The name of the current field, either a quantity or a subsection.
        This may differ from the key in the query as the key may contain indices.
        This is for internal use and is automatically set by the `RequestConfig` class.
        User should NOT set this field.
        ''')
    directive: DirectiveType = Field(DirectiveType.plain, description='''
        Indicate whether to include or exclude the current quantity/section.
        References can be resolved using `resolved`.
        The `*` is a shortcut of `plain`.
        ''')
    include: Optional[FrozenSet[str]] = Field(None, regex=r'^[*?+a-zA-z_]+$', description='''
        A list of patterns to match the quantities and subsections of the current section.
        The quantities/sections that match the include patterns AND do not match the include patterns are included.
        Only one of `include` and `exclude` can be set.
        ''')
    exclude: Optional[FrozenSet[str]] = Field(None, regex=r'^[*?+a-zA-z_]+$', description='''
        A list of patterns to match the quantities and subsections of the current section.
        The quantities/sections that match the include patterns AND do not match the include patterns are included.
        Only one of `include` and `exclude` can be set.
        ''')
    depth: int = Field(None, ge=0, description='''
        Indicate the maximum depth to be retrieved for the current section.
        If `None`, the depth is unlimited.
        This option does not apply to primitive quantities, which are always retrieved.
        ''')
    resolve_depth: int = Field(None, ge=0, description='''
        Indicate the maximum depth to be resolved for references.
        If `None`, the depth is unlimited.
        ''')
    resolve_type: ResolveType = Field(None, description='''
        Indicate how the current data should be interpreted.
        This option does not affect normal quantities/sections and should be left unassigned in most cases.
        If a value is assigned, for example `upload`, the target data will be treated as an upload id and
        the corresponding information will be retrieved.
        The original data will be left as it is if the assigned resolve type cannot find additional information.
        ''')
    max_list_size: int = Field(None, ge=0, description='''
        Indicate the size limit of lists. If assigned, lists longer than this limit will be ignored.
        ''')
    max_dict_size: int = Field(None, ge=0, description='''
        Indicate the size limit of dictionaries. If assigned, dictionaries larger than this limit will be ignored.
        ''')
    resolve_inplace: bool = Field(False, description='''
        Indicate whether to resolve references in-place.
        Deprecated, always set to `False`.
        If `false`, the reference string will be kept unchanged.
        The resolved quantity/section will be placed in the same archive.
        ''')
    include_definition: DefinitionType = Field(DefinitionType.none, description='''
        Indicate whether to include the definition of the current section.
        If `default`, the default original standard definition bundled with the NOMAD will be included.
        If `custom`, the custom definition will be included.
        If `both`, both original and custom definitions will be included, with custom definitions taking precedence.
        If `none`, no definition will be included.
        ''')
    index: Union[Tuple[int], Tuple[Optional[int], Optional[int]]] = Field(None, description='''
        The start and end index of the current field if it is a list.
        Can be a tuple of one index: (index).
        Or a tuple of two indices: (start, end), in which one of two can be `None`.
        This index field can be optionally used to slice the list, but the indices in key name has a higher priority.
        ''')
    inherit_from_parent: bool = Field(True, description='''
        Indicate whether to inherit the configuration from the parent section.
        This field only applies to the target section only, i.e., it does not propagate to its children.
        ''')
    pagination: Union[dict, DatasetPagination, UploadProcDataPagination, MetadataPagination] = Field(
        None, description='''
        The pagination configuration used for MongoDB search.
        This setting does not propagate to its children.
        For 'm_entries', 'm_uploads' and 'm_datasets', different validation rules apply.
        Please refer to `DatasetPagination`, `UploadProcDataPagination`, `MetadataPagination` for details.
        ''')
    query: Union[dict, DatasetQuery, UploadProcDataQuery, WithQuery] = Field(None, description='''
        The query configuration used for either mongo or elastic search.
        This setting does not propagate to its children.
        It can only be defined at the root levels including 'm_entries', 'm_uploads' and 'm_datasets'.
        For 'm_entries', the query is used in elastic search. It must comply with `WithQuery`.
        For 'm_uploads', the query is used in mongo search. It must comply with `UploadProcDataQuery`.
        For 'm_datasets', the query is used in mongo search. It must comply with `DatasetQuery`.
        ''')

    class Config:
        # do NOT allow extra fields
        extra = Extra.forbid
        keep_untouched = (functools.cached_property,)

    def __init__(self, **data):
        super().__init__(**data)
        if self.include and self.exclude:
            raise ValueError('One and only one of include and exclude can be set.')

        # do not set default using `Field`
        if self.include is None and self.exclude is None:
            self.include = frozenset({'*'})

    @validator('resolve_inplace')
    def _validate_directive(cls, _v):  # pylint: disable=no-self-argument
        return False

    def new(self, query: dict, *, retain_pattern: bool = False) -> RequestConfig:
        '''
        Create a new `RequestConfig` instance from a dictionary if possible.
        '''
        query_copy: dict = {k.replace('-', '_'): v for k, v in query.items()}

        try:
            if query_copy.pop('inherit_from_parent', None) is False:
                return RequestConfig.parse_obj(query_copy)

            if not retain_pattern:
                # override parent's include pattern
                if self.include is not None:
                    query_copy.setdefault('include', None)

                # override parent's exclude pattern
                if self.exclude is not None:
                    query_copy.setdefault('exclude', None)

            return RequestConfig.parse_obj(dict(self.dict(exclude_defaults=True, exclude_none=True), **query_copy))
        except ValidationError:
            raise ValueError(f'Invalid query config: {query}.')

    @staticmethod
    @functools.lru_cache(maxsize=1024)
    def _normalise_pattern(pattern: FrozenSet[str]) -> FrozenSet[str]:  # pylint: disable=no-self-argument
        '''
        Normalise the patterns.
        The received pattern can be regular expression and glob pattern, such as `quantity` and `quantity*`.
        In order to use it with regular expressions, we need to convert the glob patterns to regular expressions.
        1. Remove consecutive `*`s such that `**` becomes `*`.
        2. Replace `*` with `[a-zA-z_]*` such that `quantity*` becomes `quantity[a-zA-z_]*`.
        3. Add `^` and `$` to the beginning and end of the pattern such that `quantity` becomes `^quantity$`.
        '''
        pattern = frozenset(re.sub(r'\*+', '*', v) for v in pattern)
        # replace wildcard with regex
        pattern = frozenset(v.replace('*', r'\w*').replace('?', r'\w') for v in pattern)
        # add ^ and $ to the pattern if not present
        pattern = frozenset('^' + v if not v.startswith('^') else v for v in pattern)
        pattern = frozenset(v + '$' if not v.endswith('$') else v for v in pattern)
        return pattern

    def if_include(self, key: str) -> bool:
        '''
        For a given key, check whether it should be included.
        '''
        if self.include:
            return any(re.match(pattern, key) for pattern in self._normalise_pattern(self.include))
        if self.exclude:
            return not any(re.match(pattern, key) for pattern in self._normalise_pattern(self.exclude))

        # should not reach here
        raise ValueError('Invalid config: neither include nor exclude is set.')

    @functools.cached_property
    def hash(self) -> str:
        return sha1(self.json(exclude_defaults=True, exclude_none=True).encode('utf-8')).hexdigest()


class RequestQuery(Dict[str, Union["RequestQuery", Dict[str, RequestConfig]]]):
    pass
