# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import pytest

from nomad.parsing import LocalBackend
from nomad.normalizing import normalizers

from tests.test_parsing import parsed_vasp_example  # pylint: disable=unused-import


@pytest.fixture
def normalized_vasp_example(parsed_vasp_example: LocalBackend) -> LocalBackend:
    status, _ = parsed_vasp_example.status

    assert status == 'ParseSuccess'

    for normalizer_class in normalizers:
        normalizer = normalizer_class(parsed_vasp_example)
        normalizer.normalize()

    return parsed_vasp_example


def test_normalizer(normalized_vasp_example: LocalBackend):
    assert normalized_vasp_example.get_value('atom_species', 0) is not None
    assert normalized_vasp_example.get_value('system_type', 0) is not None
    assert normalized_vasp_example.get_value('crystal_system', 0) is not None
    assert normalized_vasp_example.get_value('space_group_number', 0) is not None
    assert normalized_vasp_example.get_value('XC_functional_name', 0) is not None
    assert normalized_vasp_example.get_value('chemical_composition', 0) is not None
    assert normalized_vasp_example.get_value('chemical_composition_bulk_reduced', 0) is not None
