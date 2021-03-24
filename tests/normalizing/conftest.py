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

from typing import List

import pytest

from ase import Atoms
import ase.build

from nomad.units import ureg
from nomad.normalizing import normalizers
from nomad.datamodel import EntryArchive
from nomad.datamodel.metainfo.common_dft import (
    Method,
    MethodToMethodRefs,
    XCFunctionals,
    FrameSequence,
    SamplingMethod,
    SingleConfigurationCalculation,
    Run,
    System,
    Dos,
    KBand,
    KBandSegment,
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
    run.program_version = "4.6.35"
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
    # system = template.section_run[0].m_create(System)
    # system.atom_positions = atoms.get_positions() * 1E-10
    # system.atom_labels = atoms.get_chemical_symbols()
    # system.simulation_cell = atoms.get_cell() * 1E-10
    # system.configuration_periodic_dimensions = atoms.get_pbc()
    system = get_section_system(atoms)
    template.section_run[0].m_add_sub_section(Run.section_system, system)

    return run_normalize(template)


def get_section_system(atoms: Atoms):
    system = System()
    system.atom_positions = atoms.get_positions() * 1E-10
    system.atom_labels = atoms.get_chemical_symbols()
    system.simulation_cell = atoms.get_cell() * 1E-10
    system.configuration_periodic_dimensions = atoms.get_pbc()
    return system


def get_template_dos(
        spin_polarized: bool = False,
        type: str = "electronic",
        has_references: bool = True,
        normalize: bool = True) -> EntryArchive:
    """Used to create a test data for DOS.

    Args:
        spin_polarized: Whether the DOS is spin-polarized or not.
        type: "electronic" or "vibrational"
        has_references: Whether the DOS has energy references or not.
        normalize: Whether the returned value is already normalized or not.
    """
    if spin_polarized and type != "electronic":
        raise ValueError("Cannot create spin polarized DOS for non-electronic data.")
    template = get_template()
    scc = template.section_run[0].section_single_configuration_calculation[0]
    idx_valence = 20
    idx_conduction = 150
    n_values = 200
    dos = scc.m_create(Dos)
    dos.dos_kind = type
    energies = np.linspace(-5, 20, n_values)
    if spin_polarized:
        dos.dos_values = np.zeros((2, n_values))
    else:
        dos.dos_values = np.zeros((1, n_values))
    dos.dos_values[:, 0:idx_valence] = 1
    dos.dos_values[:, idx_conduction:] = 1
    dos.dos_energies = energies

    if has_references:
        n_spin_channels = 2 if spin_polarized else 1
        scc.energy_reference_highest_occupied = [energies[idx_valence]] * n_spin_channels
        scc.energy_reference_lowest_unoccupied = [energies[idx_conduction]] * n_spin_channels
        scc.energy_reference_fermi = [energies[idx_valence + 1]] * n_spin_channels

    # import matplotlib.pyplot as mpl
    # if spin_polarized:
        # mpl.plot(dos.dos_energies, dos.dos_values[0, :])
        # mpl.plot(dos.dos_energies, dos.dos_values[1, :])
    # else:
        # mpl.plot(dos.dos_energies, dos.dos_values[0, :])
    # mpl.show()
    # raise

    if normalize:
        template = run_normalize(template)
    return template


def get_template_band_structure(
        band_gaps: List = None,
        type: str = "electronic",
        has_references: bool = True,
        normalize: bool = True) -> EntryArchive:
    """Used to create a test data for band structures.

    Args:
        band_gaps: List containing the band gap value and band gap type as a
            tuple, e.g. [(1, "direct"), (0.5, "indirect)]. Band gap values are
            in eV. Use a list of Nones if you don't want a gap for a specific
            channel.
        type: "electronic" or "vibrational"
        has_references: Whether the band structure has energy references or not.
        normalize: Whether the returned value is already normalized or not.
    """
    if band_gaps is None:
        band_gaps = [None]
    template = get_template()
    scc = template.section_run[0].section_single_configuration_calculation[0]
    bs = scc.m_create(KBand)
    bs.band_structure_kind = type
    if type == "electronic":
        n_spin_channels = len(band_gaps)
        fermi: List[float] = []
        highest: List[float] = []
        lowest: List[float] = []
        for gap in band_gaps:
            if gap is None:
                highest.append(0)
                lowest.append(0)
                fermi.append(0)
            else:
                fermi.append(1 * 1.60218e-19)

        if has_references:
            scc.energy_reference_fermi = fermi
            if len(highest) > 0:
                scc.energy_reference_highest_occupied = highest
            if len(lowest) > 0:
                scc.energy_reference_lowest_unoccupied = lowest
    else:
        n_spin_channels = 1
    n_segments = 2
    full_space = np.linspace(0, 2 * np.pi, 200)
    k, m = divmod(len(full_space), n_segments)
    space = list((full_space[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n_segments)))
    for i_seg in range(n_segments):
        krange = space[i_seg]
        n_points = len(krange)
        seg = bs.m_create(KBandSegment)
        energies = np.zeros((n_spin_channels, n_points, 2))
        k_points = np.zeros((n_points, 3))
        k_points[:, 0] = np.linspace(0, 1, n_points)
        if type == "electronic":
            for i_spin in range(n_spin_channels):
                if band_gaps[i_spin] is not None:
                    if band_gaps[i_spin][1] == "direct":
                        energies[i_spin, :, 0] = -np.cos(krange)
                        energies[i_spin, :, 1] = np.cos(krange)
                    elif band_gaps[i_spin][1] == "indirect":
                        energies[i_spin, :, 0] = -np.cos(krange)
                        energies[i_spin, :, 1] = np.sin(krange)
                    else:
                        raise ValueError("Invalid band gap type")
                    energies[i_spin, :, 1] += 2 + band_gaps[i_spin][0]
                else:
                    energies[i_spin, :, 0] = -np.cos(krange)
                    energies[i_spin, :, 1] = np.cos(krange)
        else:
            energies[0, :, 0] = -np.cos(krange)
            energies[0, :, 1] = np.cos(krange)
        seg.band_energies = energies * 1.60218e-19
        seg.band_k_points = k_points

    # For plotting
    # e = []
    # for i_seg in range(n_segments):
        # seg = bs.section_k_band_segment[i_seg]
        # e.append(seg.band_energies)
    # e = np.concatenate(e, axis=1)
    # import matplotlib.pyplot as mpl
    # mpl.plot(e[0], color="blue")
    # if n_spin_channels == 2:
        # mpl.plot(e[1], color="orange")
    # mpl.show()

    if normalize:
        template = run_normalize(template)
    return template


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
    template = get_template()
    template.section_run[0].section_frame_sequence = None
    template.section_run[0].section_sampling_method = None
    template.section_run[0].section_system = None
    template.section_run[0].section_single_configuration_calculation = None
    run = template.section_run[0]
    atoms1 = ase.build.bulk('Si', 'diamond', cubic=True, a=5.431)
    atoms2 = ase.build.bulk('Si', 'diamond', cubic=True, a=5.431)
    atoms2.translate([0.01, 0, 0])
    sys1 = get_section_system(atoms1)
    sys2 = get_section_system(atoms2)
    scc1 = run.m_create(SingleConfigurationCalculation)
    scc2 = run.m_create(SingleConfigurationCalculation)
    scc1.energy_total_T0 = 1e-19
    scc2.energy_total_T0 = 0.5e-19
    scc1.single_configuration_calculation_to_system_ref = sys1
    scc2.single_configuration_calculation_to_system_ref = sys2
    scc1.single_configuration_to_calculation_method_ref = run.section_method[0]
    scc2.single_configuration_to_calculation_method_ref = run.section_method[0]
    run.m_add_sub_section(Run.section_system, sys1)
    run.m_add_sub_section(Run.section_system, sys2)
    sampling_method = run.m_create(SamplingMethod)
    sampling_method.sampling_method = "geometry_optimization"
    sampling_method.geometry_optimization_energy_change = 1e-3 * ureg.electron_volt
    sampling_method.geometry_optimization_threshold_force = 1e-11 * ureg.newton
    sampling_method.geometry_optimization_geometry_change = 1e-3 * ureg.angstrom
    sampling_method.geometry_optimization_method = "bfgs"
    frame_sequence = run.m_create(FrameSequence)
    frame_sequence.frame_sequence_local_frames_ref = [scc1, scc2]
    frame_sequence.frame_sequence_to_sampling_ref = sampling_method
    frame_sequence.number_of_frames_in_sequence = 1
    return run_normalize(template)


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
def dos_electronic() -> EntryArchive:
    template = get_template_dos()
    return run_normalize(template)


@pytest.fixture(scope='session')
def bands_unpolarized_gap_indirect() -> EntryArchive:
    return get_template_band_structure([(1, "indirect")])


@pytest.fixture(scope='session')
def bands_polarized_no_gap() -> EntryArchive:
    return get_template_band_structure([None, None])


@pytest.fixture(scope='session')
def bands_unpolarized_no_gap() -> EntryArchive:
    return get_template_band_structure([None])


@pytest.fixture(scope='session')
def bands_polarized_gap_indirect() -> EntryArchive:
    return get_template_band_structure([(1, "indirect"), (0.8, "indirect")])


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
def band_path_cF() -> EntryArchive:
    """Band structure calculation for a cP Bravais lattice.
    """
    parser_name = "parsers/vasp"
    filepath = "tests/data/normalizers/band_structure/cF/vasprun.xml.bands.xz"
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


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
