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

from nomad.config.models.config import Options


@pytest.mark.parametrize(
    'include, exclude, expected_keys',
    [
        pytest.param(['B', 'A'], None, ['B', 'A'], id='custom order'),
        pytest.param(['A', 'B'], ['A'], ['B'], id='exclude takes precedence'),
        pytest.param(['*'], None, ['A', 'B'], id='include all'),
        pytest.param(['A', 'B'], ['*'], [], id='exclude all'),
    ],
)
def test_options(include, exclude, expected_keys):
    options = Options(options={'A': 'A', 'B': 'B'}, include=include, exclude=exclude)
    assert options.filtered_keys() == expected_keys
