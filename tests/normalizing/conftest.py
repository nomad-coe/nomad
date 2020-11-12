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

from ase import Atoms

from nomad.normalizing import normalizers
from nomad.datamodel import EntryArchive
from nomad.datamodel.metainfo.public import section_system as System

from tests.test_parsing import parsed_vasp_example  # pylint: disable=unused-import
from tests.test_parsing import parsed_template_example  # pylint: disable=unused-import
from tests.test_parsing import parsed_example  # pylint: disable=unused-import
from tests.test_parsing import parse_file
from tests.test_parsing import parsed_template_no_system  # pylint: disable=unused-import


def run_normalize(entry_archive: EntryArchive) -> EntryArchive:
    for normalizer_class in normalizers:
        normalizer = normalizer_class(entry_archive)
        normalizer.normalize()
    return entry_archive


@pytest.fixture
def normalized_vasp_example(parsed_vasp_example: EntryArchive) -> EntryArchive:
    return run_normalize(parsed_vasp_example)


@pytest.fixture
def normalized_example(parsed_example: EntryArchive) -> EntryArchive:
    return run_normalize(parsed_example)


@pytest.fixture
def normalized_template_example(parsed_template_example) -> EntryArchive:
    return run_normalize(parsed_template_example)


def run_normalize_for_structure(atoms: Atoms) -> EntryArchive:
    template = parsed_template_no_system()

    # Fill structural information
    system = template.section_run[0].m_create(System)
    system.atom_positions = atoms.get_positions() * 1E-10
    system.atom_labels = atoms.get_chemical_symbols()
    system.simulation_cell = atoms.get_cell() * 1E-10
    system.configuration_periodic_dimensions = atoms.get_pbc()

    return run_normalize(template)


@pytest.fixture(scope='session')
def single_point(two_d) -> EntryArchive:
    return two_d


@pytest.fixture(scope='session')
def gw(two_d) -> EntryArchive:
    parser_name = "parsers/template"
    filepath = "tests/data/normalizers/gw.json"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def geometry_optimization() -> EntryArchive:
    parser_name = "parsers/template"
    filepath = "tests/data/normalizers/fcc_crystal_structure.json"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def molecular_dynamics(bulk) -> EntryArchive:
    return bulk


@pytest.fixture(scope='session')
def phonon() -> EntryArchive:
    parser_name = "parsers/phonopy"
    filepath = "tests/data/parsers/phonopy/phonopy-FHI-aims-displacement-01/control.in"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def elastic() -> EntryArchive:
    parser_name = "parsers/elastic"
    filepath = "tests/data/parsers/elastic/diamond/INFO_ElaStic"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def bulk() -> EntryArchive:
    parser_name = "parsers/cp2k"
    filepath = "tests/data/normalizers/cp2k_bulk_md/si_md.out"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def two_d() -> EntryArchive:
    parser_name = "parsers/fhi-aims"
    filepath = "tests/data/normalizers/fhiaims_2d_singlepoint/aims.out"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def surface() -> EntryArchive:
    parser_name = "parsers/fhi-aims"
    filepath = "tests/data/normalizers/fhiaims_surface_singlepoint/PBE-light+tight-rho2.out"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def molecule() -> EntryArchive:
    parser_name = "parsers/fhi-aims"
    filepath = "tests/data/normalizers/fhiaims_molecule_singlepoint/aims.out"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def atom() -> EntryArchive:
    parser_name = "parsers/gaussian"
    filepath = "tests/data/normalizers/gaussian_atom_singlepoint/m9b7.out"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def one_d() -> EntryArchive:
    parser_name = "parsers/exciting"
    filepath = "tests/data/normalizers/exciting_1d_singlepoint/INFO.OUT"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def bands_unpolarized_gap_indirect() -> EntryArchive:
    parser_name = "parsers/vasp"
    filepath = "tests/data/normalizers/band_structure/unpolarized_gap/vasprun.xml.bands.xz"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def bands_polarized_no_gap() -> EntryArchive:
    parser_name = "parsers/vasp"
    filepath = "tests/data/normalizers/band_structure/polarized_no_gap/vasprun.xml.bands.xz"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def bands_unpolarized_no_gap() -> EntryArchive:
    parser_name = "parsers/vasp"
    filepath = "tests/data/normalizers/band_structure/unpolarized_no_gap/vasprun.xml.bands.xz"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def bands_polarized_gap_indirect() -> EntryArchive:
    parser_name = "parsers/vasp"
    filepath = "tests/data/normalizers/band_structure/polarized_gap/vasprun.xml.bands.xz"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def dos_si_vasp() -> EntryArchive:
    parser_name = "parsers/vasp"
    filepath = "tests/data/normalizers/dos/dos_si_vasp/vasprun.xml.relax2.xz"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def dos_si_exciting() -> EntryArchive:
    parser_name = "parsers/exciting"
    filepath = "tests/data/normalizers/dos/dos_si_exciting/INFO.OUT"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def dos_si_fhiaims() -> EntryArchive:
    parser_name = "parsers/fhi-aims"
    filepath = "tests/data/normalizers/dos/dos_si_fhiaims/aims.log"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def dos_polarized_vasp() -> EntryArchive:
    parser_name = "parsers/vasp"
    filepath = "tests/data/normalizers/dos/polarized_vasp/vasprun.xml.relax2.xz"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def dos_unpolarized_vasp() -> EntryArchive:
    parser_name = "parsers/vasp"
    filepath = "tests/data/normalizers/dos/unpolarized_vasp/vasprun.xml.xz"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def hash_exciting() -> EntryArchive:
    parser_name = "parsers/exciting"
    filepath = "tests/data/normalizers/hashes/exciting/INFO.OUT"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def hash_vasp(bands_unpolarized_gap_indirect) -> EntryArchive:
    return bands_unpolarized_gap_indirect


@pytest.fixture(scope='session')
def band_path_cF(bands_unpolarized_gap_indirect) -> EntryArchive:
    """Band structure calculation for a cP Bravais lattice.
    """
    return bands_unpolarized_gap_indirect


@pytest.fixture(scope='session')
def band_path_tP() -> EntryArchive:
    """Band structure calculation for a tP Bravais lattice.
    """
    parser_name = "parsers/vasp"
    filepath = "tests/data/normalizers/band_structure/tP/vasprun.xml.bands.xz"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def band_path_hP() -> EntryArchive:
    """Band structure calculation for a hP Bravais lattice.
    """
    parser_name = "parsers/vasp"
    filepath = "tests/data/normalizers/band_structure/hP/vasprun.xml.bands.xz"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def band_path_mP_nonstandard() -> EntryArchive:
    """Band structure calculation for a mP Bravais lattice with a non-standard
    lattice ordering.
    """
    parser_name = "parsers/vasp"
    filepath = "tests/data/normalizers/band_structure/mP_nonstandard/vasprun.xml.bands.xz"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def band_path_cF_nonstandard() -> EntryArchive:
    """Band structure calculation for a mP Bravais lattice with a non-standard
    lattice ordering.
    """
    parser_name = "parsers/exciting"
    filepath = "tests/data/normalizers/band_structure/cF_nonstandard/INFO.OUT"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)
