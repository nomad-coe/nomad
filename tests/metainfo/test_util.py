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

from nomad.metainfo.util import MDefNotFound, MDefWithoutMetainfo, resolve_m_def
from tests.metainfo.test_metainfo import SectionWithBoth


@pytest.mark.parametrize(
    'identifier, expected',
    [
        pytest.param(
            'tests.metainfo.test_metainfo.SectionWithBoth',
            SectionWithBoth.m_def,
            id='section',
        ),
        pytest.param(
            'tests.metainfo.test_metainfo.SectionWithBoth.quantity',
            SectionWithBoth.quantity,
            id='quantity',
        ),
        pytest.param(
            'tests.metainfo.test_metainfo.SectionWithBoth.subsection_norepeat',
            SectionWithBoth.subsection_norepeat,
            id='subsection',
        ),
    ],
)
def test_resolve_m_def(identifier, expected):
    assert resolve_m_def(identifier) == expected


@pytest.mark.parametrize(
    'identifier, expected_error',
    [
        pytest.param('nonexistent', MDefNotFound, id='non-existent-module'),
        pytest.param('non.existent', MDefNotFound, id='non-existent-class'),
        pytest.param(
            'tests.metainfo.test_metainfo.SectionWithBoth.nonexistent',
            MDefNotFound,
            id='non-existent-property',
        ),
        pytest.param(
            'nomad.metainfo.data_type.Datatype',
            MDefWithoutMetainfo,
            id='valid-class-no-m_def',
        ),
    ],
)
def test_resolve_m_def_errors(identifier, expected_error):
    with pytest.raises(expected_error):
        resolve_m_def(identifier)
