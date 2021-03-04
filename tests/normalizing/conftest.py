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

import numpy as np

import pytest

from ase import Atoms
import ase.build

from nomad.normalizing import normalizers
from nomad.datamodel import EntryArchive
from nomad.datamodel.metainfo.public import (
    section_system as System,
    Run,
    SingleConfigurationCalculation,
    XCFunctionals
)
from nomad.datamodel.metainfo.common_dft import (
    Method,
    MethodToMethodRefs,
    XCFunctionals,
    FrameSequence,
    SamplingMethod
)

from tests.parsing.test_parsing import parsed_vasp_example  # pylint: disable=unused-import
from tests.parsing.test_parsing import parsed_template_example  # pylint: disable=unused-import
from tests.parsing.test_parsing import parsed_example  # pylint: disable=unused-import
from tests.parsing.test_parsing import parse_file


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


def get_template() -> EntryArchive:
    """Returns a basic archive template.
    """
    template = EntryArchive()
    run = template.m_create(Run)
    run.program_name = "VASP"
    run.program_version = "4.6.35  3Apr08 complex  parallel LinuxIFC"
    run.program_basis_set_type = "plane waves"
    method = run.m_create(Method)
    method.electronic_structure_method = "DFT"
    functional = method.m_create(XCFunctionals)
    functional.XC_functional_name = "GGA_X_PBE"
    system = run.m_create(System)
    system.simulation_cell = [
        [5.76372622e-10, 0.0, 0.0],
        [0.0, 5.76372622e-10, 0.0],
        [0.0, 0.0, 4.0755698899999997e-10]
    ]
    system.atom_positions = [
        [2.88186311e-10, 0.0, 2.0377849449999999e-10],
        [0.0, 2.88186311e-10, 2.0377849449999999e-10],
        [0.0, 0.0, 0.0],
        [2.88186311e-10, 2.88186311e-10, 0.0],
    ]
    system.atom_labels = ["Br", "K", "Si", "Si"]
    system.configuration_periodic_dimensions = [True, True, True]
    scc = run.m_create(SingleConfigurationCalculation)
    scc.single_configuration_calculation_to_system_ref = system
    scc.single_configuration_to_calculation_method_ref = method
    scc.energy_free = -1.5936767191492225e-18
    scc.energy_total = -1.5935696296699573e-18
    scc.energy_total_T0 = -3.2126683561907e-22
    sampling_method = run.m_create(SamplingMethod)
    sampling_method.sampling_method = "geometry_optimization"
    frame_sequence = run.m_create(FrameSequence)
    frame_sequence.frame_sequence_to_sampling_ref = sampling_method
    frame_sequence.frame_sequence_local_frames_ref = [scc]
    return template


def get_template_for_structure(atoms: Atoms) -> EntryArchive:
    template = get_template()
    template.section_run[0].section_single_configuration_calculation[0].single_configuration_calculation_to_system_ref = None
    template.section_run[0].section_system = None

    # Fill structural information
    system = template.section_run[0].m_create(System)
    system.atom_positions = atoms.get_positions() * 1E-10
    system.atom_labels = atoms.get_chemical_symbols()
    system.simulation_cell = atoms.get_cell() * 1E-10
    system.configuration_periodic_dimensions = atoms.get_pbc()

    return run_normalize(template)


@pytest.fixture(scope='session')
def dft() -> EntryArchive:
    """DFT calculation."""
    template = get_template()
    template.section_run[0].section_method = None
    run = template.section_run[0]
    method_dft = run.m_create(Method)
    method_dft.electronic_structure_method = "DFT"
    method_dft.smearing_kind = "gaussian"
    method_dft.smearing_width = 1e-20
    method_dft.scf_threshold_energy_change = 1e-24
    method_dft.number_of_spin_channels = 2
    method_dft.van_der_Waals_method = "G06"
    method_dft.relativity_method = "scalar_relativistic"
    C = method_dft.m_create(XCFunctionals)
    C.XC_functional_name = "GGA_C_PBE"
    C.XC_functional_weight = 1.0
    X = method_dft.m_create(XCFunctionals)
    X.XC_functional_name = "GGA_X_PBE"
    X.XC_functional_weight = 1.0
    return run_normalize(template)


@pytest.fixture(scope='session')
def dft_plus_u() -> EntryArchive:
    """DFT+U calculation."""
    template = get_template()
    template.section_run[0].section_method = None
    run = template.section_run[0]
    method_dft = run.m_create(Method)
    method_dft.electronic_structure_method = "DFT+U"
    method_dft.smearing_kind = "gaussian"
    method_dft.smearing_width = 1e-20
    method_dft.scf_threshold_energy_change = 1e-24
    method_dft.number_of_spin_channels = 2
    method_dft.van_der_Waals_method = "G06"
    method_dft.relativity_method = "scalar_relativistic"
    C = method_dft.m_create(XCFunctionals)
    C.XC_functional_name = "GGA_C_PBE"
    C.XC_functional_weight = 1.0
    X = method_dft.m_create(XCFunctionals)
    X.XC_functional_name = "GGA_X_PBE"
    X.XC_functional_weight = 1.0
    return run_normalize(template)


@pytest.fixture(scope='session')
def gw() -> EntryArchive:
    """GW calculation."""
    template = get_template()
    template.section_run[0].section_method = None
    run = template.section_run[0]
    method_dft = run.m_create(Method)
    method_dft.electronic_structure_method = "DFT"
    method_dft.smearing_kind = "gaussian"
    method_dft.smearing_width = 1e-20
    method_dft.scf_threshold_energy_change = 1e-24
    method_dft.number_of_spin_channels = 2
    method_dft.van_der_Waals_method = "G06"
    method_dft.relativity_method = "scalar_relativistic"
    C = method_dft.m_create(XCFunctionals)
    C.XC_functional_name = "GGA_C_PBE"
    C.XC_functional_weight = 1.0
    X = method_dft.m_create(XCFunctionals)
    X.XC_functional_name = "GGA_X_PBE"
    X.XC_functional_weight = 1.0
    method_gw = run.m_create(Method)
    method_gw.electronic_structure_method = "G0W0"
    method_gw.gw_type = "G0W0"
    method_gw.gw_starting_point = "GGA_C_PBE GGA_X_PBE"
    ref = method_gw.m_create(MethodToMethodRefs)
    ref.method_to_method_kind = "starting_point"
    ref.method_to_method_ref = run.section_method[0]
    return run_normalize(template)


@pytest.fixture(scope='session')
def bulk() -> EntryArchive:
    atoms = ase.build.bulk('Si', 'diamond', cubic=True, a=5.431)
    return get_template_for_structure(atoms)


@pytest.fixture(scope='session')
def two_d() -> EntryArchive:
    atoms = Atoms(
        symbols=["C", "C"],
        scaled_positions=[
            [0, 0, 0.5],
            [1 / 3, 1 / 3, 0.5],
        ],
        cell=[
            [2.461, 0, 0],
            [np.cos(np.pi / 3) * 2.461, np.sin(np.pi / 3) * 2.461, 0],
            [0, 0, 20]
        ],
        pbc=True
    )
    return get_template_for_structure(atoms)


@pytest.fixture(scope='session')
def surface() -> EntryArchive:
    atoms = ase.build.fcc111('Al', size=(2, 2, 3), vacuum=10.0)
    return get_template_for_structure(atoms)


@pytest.fixture(scope='session')
def molecule() -> EntryArchive:
    atoms = ase.build.molecule("CO2")
    return get_template_for_structure(atoms)


@pytest.fixture(scope='session')
def atom() -> EntryArchive:
    atoms = Atoms(
        symbols=["H"],
        scaled_positions=[[0.5, 0.5, 0.5]],
        cell=[10, 10, 10],
        pbc=True,
    )
    return get_template_for_structure(atoms)


@pytest.fixture(scope='session')
def one_d() -> EntryArchive:
    atoms = ase.build.graphene_nanoribbon(1, 1, type='zigzag', vacuum=10, saturated=True)
    return get_template_for_structure(atoms)


@pytest.fixture(scope='session')
def single_point() -> EntryArchive:
    """Single point calculation."""
    template = get_template()
    template.section_run[0].section_frame_sequence = None
    template.section_run[0].section_sampling_method = None

    return run_normalize(template)


@pytest.fixture(scope='session')
def geometry_optimization() -> EntryArchive:
    parser_name = "parsers/template"
    filepath = "tests/data/normalizers/fcc_crystal_structure.json"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def molecular_dynamics() -> EntryArchive:
    """Molecular dynamics calculation."""
    template = get_template()
    template.section_run[0].section_frame_sequence = None
    template.section_run[0].section_sampling_method = None
    run = template.section_run[0]
    sampling_method = run.m_create(SamplingMethod)
    sampling_method.sampling_method = "molecular_dynamics"
    frame_sequence = run.m_create(FrameSequence)
    frame_sequence.frame_sequence_local_frames_ref = [run.section_single_configuration_calculation[0]]
    frame_sequence.frame_sequence_to_sampling_ref = sampling_method
    frame_sequence.number_of_frames_in_sequence = 1

    return run_normalize(template)


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
