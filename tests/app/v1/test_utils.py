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

import pytest

from nomad.app.v1.utils import get_query_keys


@pytest.mark.parametrize(
    'source, exclude_keys, expected_keys',
    [
        pytest.param({'a': 1, 'b': 2}, None, {'a', 'b'}, id='simple_dict'),
        pytest.param([{'a': 1}, {'b': 2}], None, {'a', 'b'}, id='list_of_dicts'),
        pytest.param({'a': {'b': 1, 'c': 2}}, None, {'a.b', 'a.c'}, id='nested_dict'),
        pytest.param(
            {'a': [{'b': 1}, {'c': 2}]}, None, {'a.b', 'a.c'}, id='dict_with_list'
        ),
        pytest.param({'a': 1, 'b': 2}, ['a'], {'b'}, id='with_exclude'),
        pytest.param(
            {'a': {'b': 1, 'c': 2}}, ['a'], {'b', 'c'}, id='nested_with_exclude'
        ),
        pytest.param({}, None, set(), id='empty_dict'),
        pytest.param([], None, set(), id='empty_list'),
        pytest.param({'a': []}, None, {'a'}, id='dict_with_empty_list'),
        pytest.param({'a': {}}, None, {'a'}, id='dict_with_empty_dict'),
    ],
)
def test_get_query_keys(source, exclude_keys, expected_keys):
    assert get_query_keys(source, exclude_keys) == expected_keys
