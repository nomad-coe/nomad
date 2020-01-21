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

from ase import Atoms

from nomad.parsing import LocalBackend
from nomad.normalizing import normalizers

from tests.test_parsing import parsed_vasp_example  # pylint: disable=unused-import
from tests.test_parsing import parsed_template_example  # pylint: disable=unused-import
from tests.test_parsing import parsed_example  # pylint: disable=unused-import
from tests.test_parsing import parse_file
from tests.test_parsing import parsed_template_no_system  # pylint: disable=unused-import


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


def run_normalize_for_structure(atoms: Atoms) -> LocalBackend:
    template = parsed_template_no_system()

    # Fill structural information
    gid = template.openSection("section_system")
    template.addArrayValues("atom_positions", atoms.get_positions() * 1E-10)
    template.addArrayValues("atom_labels", atoms.get_chemical_symbols())
    template.addArrayValues("simulation_cell", atoms.get_cell() * 1E-10)
    template.addArrayValues("configuration_periodic_dimensions", atoms.get_pbc())
    template.closeSection("section_system", gid)

    return run_normalize(template)


@pytest.fixture(scope='session')
def geometry_optimization() -> LocalBackend:
    parser_name = "parsers/template"
    filepath = "tests/data/normalizers/fcc_crystal_structure.json"
    backend = parse_file((parser_name, filepath))
    backend = run_normalize(backend)
    return backend


@pytest.fixture(scope='session')
def molecular_dynamics(bulk) -> LocalBackend:
    return bulk


@pytest.fixture(scope='session')
def phonon() -> LocalBackend:
    parser_name = "parsers/phonopy"
    filepath = "tests/data/parsers/phonopy/phonopy-FHI-aims-displacement-01/control.in"
    backend = parse_file((parser_name, filepath))
    backend = run_normalize(backend)
    return backend


@pytest.fixture(scope='session')
def bulk() -> LocalBackend:
    parser_name = "parsers/cp2k"
    filepath = "tests/data/normalizers/cp2k_bulk_md/si_md.out"
    backend = parse_file((parser_name, filepath))
    backend = run_normalize(backend)
    return backend


@pytest.fixture(scope='session')
def two_d() -> LocalBackend:
    parser_name = "parsers/fhi-aims"
    filepath = "tests/data/normalizers/fhiaims_2d_singlepoint/aims.out"
    backend = parse_file((parser_name, filepath))
    backend = run_normalize(backend)
    return backend


@pytest.fixture(scope='session')
def surface() -> LocalBackend:
    parser_name = "parsers/fhi-aims"
    filepath = "tests/data/normalizers/fhiaims_surface_singlepoint/PBE-light+tight-rho2.out"
    backend = parse_file((parser_name, filepath))
    backend = run_normalize(backend)
    return backend


@pytest.fixture(scope='session')
def molecule() -> LocalBackend:
    parser_name = "parsers/fhi-aims"
    filepath = "tests/data/normalizers/fhiaims_molecule_singlepoint/aims.out"
    backend = parse_file((parser_name, filepath))
    backend = run_normalize(backend)
    return backend


@pytest.fixture(scope='session')
def atom() -> LocalBackend:
    parser_name = "parsers/gaussian"
    filepath = "tests/data/normalizers/gaussian_atom_singlepoint/m9b7.out"
    backend = parse_file((parser_name, filepath))
    backend = run_normalize(backend)
    return backend


@pytest.fixture(scope='session')
def one_d() -> LocalBackend:
    parser_name = "parsers/exciting"
    filepath = "tests/data/normalizers/exciting_1d_singlepoint/INFO.OUT"
    backend = parse_file((parser_name, filepath))
    backend = run_normalize(backend)
    return backend
