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
from tests.test_parsing import parsed_template_example  # pylint: disable=unused-import
from tests.test_parsing import parsed_example  # pylint: disable=unused-import
from tests.test_parsing import parsed_faulty_unknown_matid_example  # pylint: disable=unused-import
from tests.utils import assert_log

def run_normalize(backend: LocalBackend) -> LocalBackend:
    status, _ = backend.status

    assert status == 'ParseSuccess'

    for normalizer_class in normalizers:
        normalizer = normalizer_class(backend)
        normalizer.normalize()

    return backend


@pytest.fixture
def normalized_vasp_example(parsed_vasp_example: LocalBackend) -> LocalBackend:
    return run_normalize(parsed_vasp_example)


@pytest.fixture
def normalized_example(parsed_example: LocalBackend) -> LocalBackend:
    return run_normalize(parsed_example)


@pytest.fixture
def normalized_template_example(parsed_template_example) -> LocalBackend:
    return run_normalize(parsed_template_example)


def assert_normalized(backend):
    # The assertions are based on the quanitites need for the repository.
    assert backend.get_value('atom_species', 0) is not None
    assert backend.get_value('system_type', 0) is not None
    assert backend.get_value('chemical_composition', 0) is not None
    assert backend.get_value('chemical_composition_bulk_reduced', 0) is not None
    # The below tests are not always present for non periodic
    # cells that don't have a simulation_cell or lattice_vectors.
    if backend.get_value('system_type', 0) not in ['Atom', 'Molecule / Cluster']:
        assert backend.get_value('crystal_system', 0) is not None
        assert backend.get_value('space_group_number', 0) is not None
    # The NWChem example for MD does not have functional information in its output.
    if backend.get_value('program_name', 0) != 'NWChem':
        assert backend.get_value('XC_functional_name', 0) is not None


def test_normalizer(normalized_example: LocalBackend, no_warn):
    assert_normalized(normalized_example)


def test_normalizer_faulty_matid(
        parsed_faulty_unknown_matid_example: LocalBackend, caplog):
    """ Runs normalizer on an example w/ bools for atom pos. Should force matid error."""
    run_normalize(parsed_faulty_unknown_matid_example)
    unknown_class_error = (
        'Matid classfication has given us an unexpected type')

    wrong_class_for_no_sim_cell = (
        'Matid classified more than 1D despite having no simulation_cell')

    assert_log(caplog, 'ERROR', unknown_class_error)
    assert_log(caplog, 'ERROR', wrong_class_for_no_sim_cell)

