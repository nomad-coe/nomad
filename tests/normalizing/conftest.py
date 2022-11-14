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

# from atomisticparsers.utils.mdanalysis import mean_squared_displacement
import numpy as np
from typing import List, Union
import pytest
from ase import Atoms
import ase.build
import re
import yaml
from warnings import warn
# from nomad.datamodel.results import MeanSquaredDisplacement

from nomad.utils import strip
from nomad.units import ureg
from nomad.normalizing import normalizers
from nomad.datamodel import EntryArchive
from nomad.datamodel.results import Relation, Symmetry, Prototype, Cell, Structure, LatticeParameters, WyckoffSet, System as resSystem
from nomad.datamodel.optimade import Species
from nomad.datamodel.metainfo.simulation.run import Run, Program
from nomad.datamodel.metainfo.simulation.method import (
    Method, BasisSet, Electronic, DFT, XCFunctional, Functional,
    Electronic, Smearing, Scf, GW, AtomParameters, HubbardModel)
from nomad.datamodel.metainfo.simulation.system import (
    AtomsGroup, System, Atoms as AtomsSystem)
from nomad.datamodel.metainfo.simulation.calculation import (
    Calculation, Energy, EnergyEntry, Dos, DosValues, BandStructure, BandEnergies)
from nomad.datamodel.metainfo.workflow import (
    DiffusionConstantValues,
    IntegrationParameters,
    MeanSquaredDisplacement,
    MeanSquaredDisplacementValues,
    MolecularDynamicsResults,
    RadialDistributionFunction,
    RadialDistributionFunctionValues,
    Workflow,
    GeometryOptimization,
    Elastic,
    MolecularDynamics,
    EquationOfState,
    EOSFit
)

from tests.parsing.test_parsing import parsed_vasp_example  # pylint: disable=unused-import
from tests.parsing.test_parsing import parsed_template_example  # pylint: disable=unused-import
from tests.parsing.test_parsing import parsed_example  # pylint: disable=unused-import
from tests.parsing.test_parsing import parse_file
from nomad.datamodel.context import ServerContext
from nomad.datamodel.datamodel import EntryArchive, EntryMetadata
from nomad.parsing.parser import ArchiveParser
from nomad.processing.data import Upload
from tests.test_files import create_test_upload_files


def run_normalize(entry_archive: EntryArchive) -> EntryArchive:
    for normalizer_class in normalizers:
        normalizer = normalizer_class(entry_archive)
        normalizer.normalize()
    return entry_archive


def run_processing(directory, mainfile):
    # create upload with example files
    upload_files = create_test_upload_files('test_upload_id', published=False, raw_files=directory)
    upload = Upload(upload_id='test_upload_id')

    # parse
    parser = ArchiveParser()
    context = ServerContext(upload=upload)
    test_archive = EntryArchive(m_context=context, metadata=EntryMetadata())
    parser.parse(
        upload_files.raw_file_object(mainfile).os_path,
        test_archive)
    run_normalize(test_archive)
    return test_archive


@pytest.fixture
def normalized_vasp_example(parsed_vasp_example: EntryArchive) -> EntryArchive:
    return run_normalize(parsed_vasp_example)


@pytest.fixture
def normalized_example(parsed_example: EntryArchive) -> EntryArchive:
    return run_normalize(parsed_example)


@pytest.fixture
def normalized_template_example(parsed_template_example) -> EntryArchive:
    return run_normalize(parsed_template_example)


def get_template_dft() -> EntryArchive:
    '''Returns a basic archive template for a DFT calculation.
    '''
    template = EntryArchive()
    run = template.m_create(Run)
    run.program = Program(name='VASP', version='4.6.35')
    method = run.m_create(Method)
    method.basis_set.append(BasisSet(type='plane waves'))
    method.electronic = Electronic(method='DFT')
    xc_functional = XCFunctional(exchange=[Functional(name='GGA_X_PBE')])
    method.dft = DFT(xc_functional=xc_functional)
    system = run.m_create(System)
    system.atoms = AtomsSystem(
        lattice_vectors=[
            [5.76372622e-10, 0.0, 0.0],
            [0.0, 5.76372622e-10, 0.0],
            [0.0, 0.0, 4.0755698899999997e-10]
        ],
        positions=[
            [2.88186311e-10, 0.0, 2.0377849449999999e-10],
            [0.0, 2.88186311e-10, 2.0377849449999999e-10],
            [0.0, 0.0, 0.0],
            [2.88186311e-10, 2.88186311e-10, 0.0],
        ],
        labels=['Br', 'K', 'Si', 'Si'],
        periodic=[True, True, True])
    scc = run.m_create(Calculation)
    scc.system_ref = system
    scc.method_ref = method
    scc.energy = Energy(
        free=EnergyEntry(value=-1.5936767191492225e-18),
        total=EnergyEntry(value=-1.5935696296699573e-18),
        total_t0=EnergyEntry(value=-3.2126683561907e-22))
    workflow = template.m_create(Workflow)
    workflow.type = 'geometry_optimization'
    return template


def get_template_eels() -> EntryArchive:
    '''Returns a basic archive template for an EELS experiment.
    '''
    # Ensure that the eels schema is loaded
    from eelsdbparser import eelsdb_parser  # pylint: disable=unused-import
    dct_data = yaml.safe_load(strip(f'''
        results:
            properties:
                spectroscopy:
                    spectrum: '#/measurement/0/eels/spectrum'
                    eels:
                        detector_type: Quantum GIF
                        min_energy: {(100 * ureg.electron_volt).to(ureg.joule).m}
                        max_energy: {(200 * ureg.electron_volt).to(ureg.joule).m}
                        resolution: {(1 * ureg.electron_volt).to(ureg.joule).m}
        measurement:
            - method_name: electron energy loss spectroscopy
              method_abbreviation: EELS
              sample:
                  - elements:
                        - Si
                        - O
                    chemical_formula: SiO
              eels:
                  spectrum: {{}}
    '''))
    archive = EntryArchive.m_from_dict(dct_data)
    archive.measurement[0].eels.spectrum.count = np.linspace(0, 100, 1)
    archive.measurement[0].eels.spectrum.energy = np.linspace(100, 200, 1)
    return archive


def get_template_for_structure(atoms: Atoms) -> EntryArchive:
    template = get_template_dft()
    template.run[0].calculation[0].system_ref = None
    template.run[0].system = None

    # Fill structural information
    # system = template.run[0].m_create(System)
    # system.atom_positions = atoms.get_positions() * 1E-10
    # system.atom_labels = atoms.get_chemical_symbols()
    # system.simulation_cell = atoms.get_cell() * 1E-10
    # system.configuration_periodic_dimensions = atoms.get_pbc()
    system = get_section_system(atoms)
    template.run[0].m_add_sub_section(Run.system, system)

    return run_normalize(template)


def get_section_system(atoms: Atoms):
    system = System()
    system.atoms = AtomsSystem(
        positions=atoms.get_positions() * 1E-10,
        labels=atoms.get_chemical_symbols(),
        lattice_vectors=atoms.get_cell() * 1E-10,
        periodic=atoms.get_pbc())
    return system


def get_template_dos(
        fill: List = [[[0, 1], [2, 3]]],
        energy_reference_fermi: Union[float, None] = None,
        energy_reference_highest_occupied: Union[float, None] = None,
        energy_reference_lowest_unoccupied: Union[float, None] = None,
        n_values: int = 101,
        type: str = 'electronic',
        normalize: bool = True) -> EntryArchive:
    '''Used to create a test data for DOS.

    Args:
        fill: List containing the energy ranges (eV) that should be filled with
            non-zero values, e.g. [[[0, 1], [2, 5]]. Defaults to single channel DOS
            with a gap.
        energy_fermi: Fermi energy (eV)
        energy_reference_highest_occupied: Highest occupied energy (eV) as given by a parser.
        energy_reference_lowest_unoccupied: Lowest unoccupied energy (eV) as given by a parser.
        type: 'electronic' or 'vibrational'
        has_references: Whether the DOS has energy references or not.
        normalize: Whether the returned value is already normalized or not.
    '''
    if len(fill) > 1 and type != 'electronic':
        raise ValueError('Cannot create spin polarized DOS for non-electronic data.')
    template = get_template_dft()
    scc = template.run[0].calculation[0]
    dos_type = Calculation.dos_electronic if type == 'electronic' else Calculation.dos_phonon
    dos = scc.m_create(Dos, dos_type)
    energies = np.linspace(-5, 5, n_values)
    for i, range_list in enumerate(fill):
        dos_total = dos.m_create(DosValues, Dos.total)
        dos_total.spin = i
        dos_value = np.zeros(n_values)
        for r in range_list:
            idx_bottom = (np.abs(energies - r[0])).argmin()
            idx_top = (np.abs(energies - r[1])).argmin()
            dos_value[idx_bottom:idx_top] = 1
        dos_total.value = dos_value

    dos.energies = energies * ureg.electron_volt
    if energy_reference_fermi is not None:
        energy_reference_fermi *= ureg.electron_volt
    if energy_reference_highest_occupied is not None:
        energy_reference_highest_occupied *= ureg.electron_volt
    if energy_reference_lowest_unoccupied is not None:
        energy_reference_lowest_unoccupied *= ureg.electron_volt
    scc.energy = Energy(
        fermi=energy_reference_fermi,
        highest_occupied=energy_reference_highest_occupied,
        lowest_unoccupied=energy_reference_lowest_unoccupied)
    # import matplotlib.pyplot as mpl
    # if n_channels == 2:
    #     mpl.plot(energies, dos.dos_values[0, :])
    #     mpl.plot(energies, dos.dos_values[1, :])
    # else:
    #     mpl.plot(energies, dos.dos_values[0, :])
    # mpl.show()

    if normalize:
        template = run_normalize(template)
    return template


def get_template_band_structure(
        band_gaps: List = None,
        type: str = 'electronic',
        has_references: bool = True,
        normalize: bool = True,
        has_reciprocal_cell: bool = True) -> EntryArchive:
    '''Used to create a test data for band structures.

    Args:
        band_gaps: List containing the band gap value and band gap type as a
            tuple, e.g. [(1, 'direct'), (0.5, 'indirect)]. Band gap values are
            in eV. Use a list of Nones if you don't want a gap for a specific
            channel.
        type: 'electronic' or 'vibrational'
        has_references: Whether the band structure has energy references or not.
        normalize: Whether the returned value is already normalized or not.
        has_reciprocal_cell: Whether the reciprocal cell is available or not
    '''
    if band_gaps is None:
        band_gaps = [None]
    template = get_template_dft()
    if not has_reciprocal_cell:
        template.run[0].system[0].atoms = None
    scc = template.run[0].calculation[0]
    if type == 'electronic':
        bs = scc.m_create(BandStructure, Calculation.band_structure_electronic)
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
            scc.energy = Energy(fermi=fermi[0])
            if len(highest) > 0:
                scc.energy.highest_occupied = highest[0]
            if len(lowest) > 0:
                scc.energy.lowest_unoccupied = lowest[0]
    else:
        bs = scc.m_create(BandStructure, Calculation.band_structure_phonon)
        n_spin_channels = 1
    n_segments = 2
    full_space = np.linspace(0, 2 * np.pi, 200)
    k, m = divmod(len(full_space), n_segments)
    space = list((full_space[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n_segments)))
    for i_seg in range(n_segments):
        krange = space[i_seg]
        n_points = len(krange)
        seg = bs.m_create(BandEnergies)
        energies = np.zeros((n_spin_channels, n_points, 2))
        k_points = np.zeros((n_points, 3))
        k_points[:, 0] = np.linspace(0, 1, n_points)
        if type == 'electronic':
            for i_spin in range(n_spin_channels):
                if band_gaps[i_spin] is not None:
                    if band_gaps[i_spin][1] == 'direct':
                        energies[i_spin, :, 0] = -np.cos(krange)
                        energies[i_spin, :, 1] = np.cos(krange)
                    elif band_gaps[i_spin][1] == 'indirect':
                        energies[i_spin, :, 0] = -np.cos(krange)
                        energies[i_spin, :, 1] = np.sin(krange)
                    else:
                        raise ValueError('Invalid band gap type')
                    energies[i_spin, :, 1] += 2 + band_gaps[i_spin][0]
                else:
                    energies[i_spin, :, 0] = -np.cos(krange)
                    energies[i_spin, :, 1] = np.cos(krange)
        else:
            energies[0, :, 0] = -np.cos(krange)
            energies[0, :, 1] = np.cos(krange)
        seg.energies = energies * 1.60218e-19
        seg.kpoints = k_points

    # For plotting
    # e = []
    # for i_seg in range(n_segments):
        # seg = bs.section_k_band_segment[i_seg]
        # e.append(seg.band_energies)
    # e = np.concatenate(e, axis=1)
    # import matplotlib.pyplot as mpl
    # mpl.plot(e[0], color='blue')
    # if n_spin_channels == 2:
        # mpl.plot(e[1], color='orange')
    # mpl.show()

    if normalize:
        template = run_normalize(template)
    return template


def set_dft_values(xc_functional_names: list) -> EntryArchive:
    ''''''
    template = get_template_dft()
    template.run[0].method = None
    run = template.run[0]
    method_dft = run.m_create(Method)
    method_dft.basis_set.append(BasisSet(type='plane waves'))
    method_dft.dft = DFT()
    method_dft.electronic = Electronic(
        method='DFT',
        smearing=Smearing(kind='gaussian', width=1e-20),
        n_spin_channels=2,
        van_der_waals_method='G06',
        relativity_method='scalar_relativistic')
    method_dft.scf = Scf(threshold_energy_change=1e-24)
    method_dft.dft.xc_functional = XCFunctional()
    xc = method_dft.dft.xc_functional
    for xc_functional_name in xc_functional_names:
        if re.search('^HYB_', xc_functional_name):
            xc.hybrid.append(Functional(name=xc_functional_name, weight=1.0))
            continue
        if re.search('_X?C_', xc_functional_name):
            xc.correlation.append(Functional(name=xc_functional_name, weight=1.0))
        if re.search('_XC?_', xc_functional_name):
            xc.exchange.append(Functional(name=xc_functional_name, weight=1.0))
        xc.correlation.append(Functional(name=xc_functional_name, weight=1.0))
    return template


@pytest.fixture(scope='session')
def dft() -> EntryArchive:
    '''DFT calculation.'''
    template = set_dft_values(['GGA_C_PBE', 'GGA_X_PBE'])
    return run_normalize(template)


@pytest.fixture(scope='session')
def dft_method_referenced() -> EntryArchive:
    '''DFT calculation with two methods: one referencing the other.'''
    template = get_template_dft()
    template.run[0].method = None
    run = template.run[0]
    method_dft = run.m_create(Method)
    method_dft.basis_set.append(BasisSet(type='plane waves'))
    method_dft.electronic = Electronic(
        smearing=Smearing(kind='gaussian', width=1e-20),
        n_spin_channels=2, van_der_waals_method='G06',
        relativity_method='scalar_relativistic')
    method_dft.scf = Scf(threshold_energy_change=1e-24)
    method_dft.dft = DFT(xc_functional=XCFunctional())
    method_dft.dft.xc_functional.correlation.append(Functional(name='GGA_C_PBE', weight=1.0))
    method_dft.dft.xc_functional.exchange.append(Functional(name='GGA_X_PBE', weight=1.0))
    method_ref = run.m_create(Method)
    method_ref.basis_set.append(BasisSet(type='plane waves'))
    method_ref.electronic = Electronic(method='DFT')
    method_ref.core_method_ref = method_dft
    run.calculation[0].method_ref = method_ref
    return run_normalize(template)


@pytest.fixture(scope='session')
def dft_exact_exchange() -> EntryArchive:
    ''''''
    template = set_dft_values(['GGA_C_PBE', 'GGA_X_PBE'])
    template.run[0].method[0].dft.xc_functional.hybrid.append(Functional())
    template.run[0].method[0].dft.xc_functional.hybrid[0].parameters = {'exact_exchange_mixing_factor': .25}
    template.run[0].method[0].dft.xc_functional.hybrid[0].name = '+alpha'
    return run_normalize(template)


@pytest.fixture(scope='session')
def dft_b3lyp() -> EntryArchive:
    ''''''
    template = set_dft_values(['HYB_GGA_XC_B3LYP'])
    return run_normalize(template)


@pytest.fixture(scope='session')
def dft_pbeh() -> EntryArchive:
    ''''''
    template = set_dft_values(['HYB_GGA_XC_PBEH'])
    return run_normalize(template)


@pytest.fixture(scope='session')
def dft_m05() -> EntryArchive:
    ''''''
    template = set_dft_values(['MGGA_C_M05', 'HYB_MGGA_X_M05'])
    return run_normalize(template)


@pytest.fixture(scope='session')
def dft_pbe0_13() -> EntryArchive:
    ''''''
    template = set_dft_values(['HYB_GGA_XC_PBE0_13'])
    return run_normalize(template)


@pytest.fixture(scope='session')
def dft_pbe38() -> EntryArchive:
    ''''''
    template = set_dft_values(['HYB_GGA_XC_PBE38'])
    return run_normalize(template)


@pytest.fixture(scope='session')
def dft_pbe50() -> EntryArchive:
    ''''''
    template = set_dft_values(['HYB_GGA_XC_PBE50'])
    return run_normalize(template)


@pytest.fixture(scope='session')
def dft_m06_2x() -> EntryArchive:
    ''''''
    template = set_dft_values(['HYB_MGGA_X_M06_2X'])
    return run_normalize(template)


@pytest.fixture(scope='session')
def dft_m05_2x() -> EntryArchive:
    ''''''
    template = set_dft_values(['MGGA_C_M05_2X', 'HYB_MGGA_X_M05_2X'])
    return run_normalize(template)


@pytest.fixture(scope='session')
def dft_plus_u() -> EntryArchive:
    '''DFT+U calculation with a Hubbard model.'''
    template = get_template_dft()
    template.run[0].method = None
    run = template.run[0]
    method_dft = run.m_create(Method)
    method_dft.basis_set.append(BasisSet(type='plane waves'))
    method_dft.electronic = Electronic(
        method='DFT+U',
        smearing=Smearing(kind='gaussian', width=1e-20),
        n_spin_channels=2, van_der_waals_method='G06',
        relativity_method='scalar_relativistic')
    method_dft.scf = Scf(threshold_energy_change=1e-24)
    method_dft.dft = DFT(xc_functional=XCFunctional())
    method_dft.dft.xc_functional.correlation.append(Functional(name='GGA_C_PBE', weight=1.0))
    method_dft.dft.xc_functional.exchange.append(Functional(name='GGA_X_PBE', weight=1.0))
    method_dft.atom_parameters.append(AtomParameters(label='Ti'))
    method_dft.atom_parameters[0].hubbard_model = HubbardModel(orbital='3d', u=4.5e-19,
                                                               j=1e-19, method='Dudarev', projection_type='on-site')
    return run_normalize(template)


@pytest.fixture(scope='session')
def gw() -> EntryArchive:
    '''GW calculation.'''
    template = get_template_dft()
    template.run[0].method = None
    run = template.run[0]
    method_dft = run.m_create(Method)
    method_dft.electronic = Electronic(
        method='DFT',
        smearing=Smearing(kind='gaussian', width=1e-20),
        n_spin_channels=2, van_der_waals_method='G06',
        relativity_method='scalar_relativistic')
    method_dft.scf = Scf(threshold_energy_change=1e-24)
    method_dft.dft = DFT(xc_functional=XCFunctional())
    method_dft.dft.xc_functional.correlation.append(Functional(name='GGA_C_PBE', weight=1.0))
    method_dft.dft.xc_functional.exchange.append(Functional(name='GGA_X_PBE', weight=1.0))

    method_gw = run.m_create(Method)
    method_gw.gw = GW(type='G0W0', starting_point='GGA_C_PBE GGA_X_PBE')
    method_gw.starting_method_ref = run.method[0]
    return run_normalize(template)


@pytest.fixture(scope='session')
def eels() -> EntryArchive:
    '''EELS experiment.'''
    template = get_template_eels()
    return run_normalize(template)


@pytest.fixture(scope='session')
def mechanical() -> EntryArchive:
    '''Entry with mechanical properties.'''
    template = get_template_dft()

    # Elastic workflow
    workflow = template.m_create(Workflow)
    workflow.type = 'elastic'
    workflow.elastic = Elastic(
        shear_modulus_hill=10000,
        shear_modulus_reuss=10000,
        shear_modulus_voigt=10000,
    )

    # EOS workflow
    workflow = template.m_create(Workflow)
    workflow.type = 'equation_of_state'
    equation_of_state = EquationOfState(
        volumes=np.linspace(0, 10, 10) * ureg.angstrom ** 3,
        energies=np.linspace(0, 10, 10) * ureg.electron_volt,
    )
    eos_fit = equation_of_state.m_create(EOSFit)
    eos_fit.function_name = 'murnaghan'
    eos_fit.fitted_energies = np.linspace(0, 10, 10) * ureg.electron_volt
    eos_fit.bulk_modulus = 10000
    workflow.equation_of_state = equation_of_state

    return run_normalize(template)


@pytest.fixture(scope='session')
def bulk() -> EntryArchive:
    atoms = ase.build.bulk('Si', 'diamond', cubic=True, a=5.431)
    return get_template_for_structure(atoms)


@pytest.fixture(scope='session')
def two_d() -> EntryArchive:
    atoms = Atoms(
        symbols=['C', 'C'],
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
    atoms = ase.build.molecule('CO2')
    return get_template_for_structure(atoms)


@pytest.fixture(scope='session')
def atom() -> EntryArchive:
    atoms = Atoms(
        symbols=['H'],
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
    '''Single point calculation.'''
    template = get_template_dft()

    return run_normalize(template)


@pytest.fixture(scope='session')
def geometry_optimization() -> EntryArchive:
    template = get_template_dft()
    template.run[0].system = None
    template.run[0].calculation = None
    run = template.run[0]
    atoms1 = ase.build.bulk('Si', 'diamond', cubic=True, a=5.431)
    atoms2 = ase.build.bulk('Si', 'diamond', cubic=True, a=5.431)
    atoms2.translate([0.01, 0, 0])
    sys1 = get_section_system(atoms1)
    sys2 = get_section_system(atoms2)
    scc1 = run.m_create(Calculation)
    scc2 = run.m_create(Calculation)
    scc1.energy = Energy(total=EnergyEntry(value=1e-19), total_t0=EnergyEntry(value=1e-19))
    scc2.energy = Energy(total=EnergyEntry(value=0.5e-19), total_t0=EnergyEntry(value=0.5e-19))
    scc1.system_ref = sys1
    scc2.system_ref = sys2
    scc1.method_ref = run.method[0]
    scc2.method_ref = run.method[0]
    run.m_add_sub_section(Run.system, sys1)
    run.m_add_sub_section(Run.system, sys2)
    workflow = template.m_create(Workflow)
    workflow.type = 'geometry_optimization'
    workflow.geometry_optimization = GeometryOptimization(
        convergence_tolerance_energy_difference=1e-3 * ureg.electron_volt,
        convergence_tolerance_force_maximum=1e-11 * ureg.newton,
        convergence_tolerance_displacement_maximum=1e-3 * ureg.angstrom,
        method='bfgs')
    return run_normalize(template)


@pytest.fixture(scope='session')
def molecular_dynamics() -> EntryArchive:
    '''Molecular dynamics calculation.'''
    template = get_template_dft()
    run = template.run[0]

    # Create calculations
    n_steps = 10
    calcs = []
    for step in range(n_steps):
        system = System()
        run.m_add_sub_section(Run.system, system)
        calc = Calculation()
        calc.system_ref = system
        calc.time = step
        calc.step = step
        calc.volume = step
        calc.pressure = step
        calc.temperature = step
        calc.energy = Energy(
            potential=EnergyEntry(value=step),
        )
        calcs.append(calc)
        run.m_add_sub_section(Run.calculation, calc)

    # Create workflow
    workflow = template.m_create(Workflow)
    workflow.type = 'molecular_dynamics'
    workflow.calculation_result_ref = calcs[-1]
    workflow.calculations_ref = calcs
    diff_values = DiffusionConstantValues(
        value=2.1,
        error_type='Pearson correlation coefficient',
        errors=0.98,
    )
    msd_values = MeanSquaredDisplacementValues(
        times=[0, 1, 2],
        n_times=3,
        value=[0, 1, 2],
        label='MOL',
        errors=[0, 1, 2],
        diffusion_constant=diff_values,
    )
    msd = MeanSquaredDisplacement(
        type='molecular',
        direction='xyz',
        error_type='bootstrapping',
        mean_squared_displacement_values=[msd_values],
    )
    rdf_values = RadialDistributionFunctionValues(
        bins=[0, 1, 2],
        n_bins=3,
        value=[0, 1, 2],
        frame_start=0,
        frame_end=100,
        label='MOL-MOL',
    )
    rdf = RadialDistributionFunction(
        type='molecular',
        radial_distribution_function_values=[rdf_values],
    )
    results = MolecularDynamicsResults(
        radial_distribution_functions=[rdf],
        mean_squared_displacements=[msd],
    )
    md = MolecularDynamics(
        thermodynamic_ensemble='NVT',
        integration_parameters=IntegrationParameters(
            integration_timestep=0.5 * ureg('fs'),
        ),
        results=results
    )
    workflow.molecular_dynamics = md

    return run_normalize(template)


def get_template_topology(pbc=False) -> EntryArchive:
    template = get_template_dft()
    run = template.run[0]
    del run.system[0]

    # System
    water1 = ase.build.molecule('H2O')
    water2 = ase.build.molecule('H2O')
    water2.translate([5, 0, 0])
    sys = water1 + water2
    sys.set_cell([10, 10, 10])
    sys.set_pbc(pbc)
    system = get_section_system(sys)
    run.m_add_sub_section(Run.system, system)

    # Topology
    molecule_group = AtomsGroup(
        label='MOL_GROUP',
        type='molecule_group',
        index=0,
        composition_formula='H(4)O(2)',
        n_atoms=6,
        atom_indices=[0, 1, 2, 3, 4, 5]
    )
    system.m_add_sub_section(System.atoms_group, molecule_group)
    molecule1 = AtomsGroup(
        label='MOL',
        type='molecule',
        index=0,
        composition_formula='H(2)O(1)',
        n_atoms=3,
        atom_indices=[0, 1, 2]
    )
    molecule_group.m_add_sub_section(AtomsGroup.atoms_group, molecule1)
    molecule2 = AtomsGroup(
        label='MOL',
        type='molecule',
        index=0,
        composition_formula='H(2)O(1)',
        n_atoms=3,
        atom_indices=[3, 4, 5]
    )
    molecule_group.m_add_sub_section(AtomsGroup.atoms_group, molecule2)
    monomer_group = AtomsGroup(
        label='MON_GROUP',
        type='monomer_group',
        index=0,
        composition_formula='H(2)',
        n_atoms=2,
        atom_indices=[0, 1]
    )
    molecule1.m_add_sub_section(AtomsGroup.atoms_group, monomer_group)
    monomer = AtomsGroup(
        label='MON',
        type='monomer',
        index=0,
        composition_formula='H(2)',
        n_atoms=2,
        atom_indices=[0, 1]
    )
    monomer_group.m_add_sub_section(AtomsGroup.atoms_group, monomer)

    # Calculation
    calc = Calculation()
    calc.system_ref = system
    run.m_add_sub_section(Run.calculation, calc)

    return run_normalize(template)


@pytest.fixture(scope='session')
def phonon() -> EntryArchive:
    parser_name = 'parsers/phonopy'
    filepath = 'tests/data/parsers/phonopy/phonopy-FHI-aims-displacement-01/control.in'
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def elastic() -> EntryArchive:
    parser_name = 'parsers/elastic'
    filepath = 'tests/data/parsers/elastic/diamond/INFO_ElaStic'
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def dos_electronic() -> EntryArchive:
    template = get_template_dos()
    return run_normalize(template)


@pytest.fixture(scope='session')
def dos_si_vasp() -> EntryArchive:
    parser_name = 'parsers/vasp'
    filepath = 'tests/data/normalizers/dos/dos_si_vasp/vasprun.xml.relax2.xz'
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def dos_si_exciting() -> EntryArchive:
    parser_name = 'parsers/exciting'
    filepath = 'tests/data/normalizers/dos/dos_si_exciting/INFO.OUT'
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def dos_si_fhiaims() -> EntryArchive:
    parser_name = 'parsers/fhi-aims'
    filepath = 'tests/data/normalizers/dos/dos_si_fhiaims/aims.log'
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def dos_polarized_vasp() -> EntryArchive:
    parser_name = 'parsers/vasp'
    filepath = 'tests/data/normalizers/dos/polarized_vasp/vasprun.xml.relax2.xz'
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def dos_unpolarized_vasp() -> EntryArchive:
    parser_name = 'parsers/vasp'
    filepath = 'tests/data/normalizers/dos/unpolarized_vasp/vasprun.xml.xz'
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def hash_exciting() -> EntryArchive:
    parser_name = 'parsers/exciting'
    filepath = 'tests/data/normalizers/hashes/exciting/INFO.OUT'
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def hash_vasp() -> EntryArchive:
    return get_template_band_structure([(1, 'indirect')])


@pytest.fixture(scope='session')
def band_path_cF() -> EntryArchive:
    '''Band structure calculation for a cP Bravais lattice.
    '''
    parser_name = 'parsers/vasp'
    filepath = 'tests/data/normalizers/band_structure/cF/vasprun.xml.bands.xz'
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def band_path_tP() -> EntryArchive:
    '''Band structure calculation for a tP Bravais lattice.
    '''
    parser_name = 'parsers/vasp'
    filepath = 'tests/data/normalizers/band_structure/tP/vasprun.xml.bands.xz'
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def band_path_hP() -> EntryArchive:
    '''Band structure calculation for a hP Bravais lattice.
    '''
    parser_name = 'parsers/vasp'
    filepath = 'tests/data/normalizers/band_structure/hP/vasprun.xml.bands.xz'
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def band_path_mP_nonstandard() -> EntryArchive:
    '''Band structure calculation for a mP Bravais lattice with a non-standard
    lattice ordering.
    '''
    parser_name = 'parsers/vasp'
    filepath = 'tests/data/normalizers/band_structure/mP_nonstandard/vasprun.xml.bands.xz'
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


@pytest.fixture(scope='session')
def band_path_cF_nonstandard() -> EntryArchive:
    '''Band structure calculation for a mP Bravais lattice with a non-standard
    lattice ordering.
    '''
    parser_name = 'parsers/exciting'
    filepath = 'tests/data/normalizers/band_structure/cF_nonstandard/INFO.OUT'
    archive = parse_file((parser_name, filepath))
    return run_normalize(archive)


def create_system(label: str,
                  structural_type: str,
                  elements: List[str],
                  formula_hill: str,
                  formula_reduced: str,
                  formula_anonymous: str,
                  system_relation: Relation,
                  indices: List[int] = None,
                  material_id: str = None,
                  atoms: Structure = None,
                  cell: Cell = None,
                  symmetry: Symmetry = None,
                  prototype: Prototype = None) -> resSystem:
    system = resSystem()
    system.label = label
    system.structural_type = structural_type
    system.elements = elements
    system.formula_hill = formula_hill
    system.formula_reduced = formula_reduced
    system.formula_anonymous = formula_anonymous
    system.system_relation = system_relation
    if label == 'subsystem':
        system.indices = indices
    elif label == 'conventional cell':
        system.material_id = material_id
        system.atoms = atoms
        system.cell = cell
        system.symmetry = Symmetry()
        system.symmetry = symmetry
        system.prototype = prototype
    else:
        warn('Warning: subsystem label is missing')
    return system


def crystal_structure_properties(crystal_structure):
    if crystal_structure == 'fcc':
        crystal_structure_property = Symmetry(bravais_lattice="cF", crystal_system="cubic", hall_number=523, hall_symbol="-F 4 2 3", point_group="m-3m", space_group_number=225, space_group_symbol="Fm-3m")
    elif crystal_structure == 'bcc':
        crystal_structure_property = Symmetry(bravais_lattice="cI", crystal_system="cubic", hall_number=529, hall_symbol="-I 4 2 3", point_group="m-3m", space_group_number=229, space_group_symbol="Im-3m")
    elif crystal_structure is None:
        warn('Warning: crystal structure is missing')
        return None
    return crystal_structure_property


def single_Cu_surface_atoms():
    # create simple Cu surface
    Cu_fcc_100 = ase.build.fcc100('Cu', (3, 5, 5), vacuum=10, periodic=True)
    Cu_fcc_100.rattle(stdev=0.001, seed=None, rng=None)
    Cu_fcc_110 = ase.build.fcc110('Cu', (3, 5, 5), vacuum=10, periodic=True)
    Cu_fcc_110.rattle(stdev=0.001, seed=None, rng=None)
    return Cu_fcc_100, Cu_fcc_110


def single_Cu_surface_prototype():
    prototype_Cu_fcc = Prototype()
    prototype_Cu_fcc.aflow_id = "A_cF4_225_a"
    prototype_Cu_fcc.assignment_method = "normalized-wyckoff"
    prototype_Cu_fcc.label = "225-Cu-cF4"
    prototype_Cu_fcc.formula = "Cu4"
    return prototype_Cu_fcc


def single_Cu_surface_topology():
    Cell_Cu_fcc = [Cell(a=3.610000000000001 * ureg.angstrom,
                        b=3.610000000000001 * ureg.angstrom,
                        c=3.610000000000001 * ureg.angstrom,
                        alpha=1.5707963267948966 * ureg.rad,
                        beta=1.5707963267948966 * ureg.rad,
                        gamma=1.5707963267948966 * ureg.rad)]
    prototype_Cu_fcc = single_Cu_surface_prototype()

    # create Cu topology system
    label_Cu = 'subsystem'
    structural_type_Cu = 'surface'
    elements_Cu = ['Cu']
    formula_hill_Cu = 'Cu75'
    formula_reduced_Cu = 'Cu75'
    formula_anonymous_Cu = 'A75'
    system_relation = Relation()
    system_relation.type = "subsystem"
    indices_Cu = [i for i in range(75)]

    subsystem_Cu = create_system(label_Cu, structural_type_Cu, elements_Cu, formula_hill_Cu, formula_reduced_Cu, formula_anonymous_Cu, system_relation, indices=indices_Cu)
    topologies_Cu = [subsystem_Cu]

    label_Cu_conv = 'conventional cell'
    structural_type_Cu_conv = 'bulk'
    formula_hill_Cu_conv = 'Cu4'
    formula_reduced_Cu_conv = 'Cu4'
    formula_anonymous_Cu_conv = 'A4'
    material_id_Cu_conv = "3M6onRRrQbutydx916-Y15I79Z_X"
    atoms_Cu_conv = Structure()
    atoms_Cu_conv.dimension_types = [0, 0, 0]
    atoms_Cu_conv.lattice_vectors = [[3.609999999999999, 0.0, 0.0], [0.0, 3.609999999999999e-10, 0.0], [0.0, 0.0, 3.609999999999999]] * ureg.angstrom
    atoms_Cu_conv.cartesian_site_positions = [[0.0, 0.0, 0.0], [0.0, 1.8049999999999995, 1.8049999999999995], [1.8049999999999995, 0.0, 1.8049999999999995], [1.8049999999999995, 1.8049999999999995, 0.0]] * ureg.angstrom
    atoms_Cu_conv.species_at_sites = ["Cu", "Cu", "Cu", "Cu"]
    atoms_Cu_conv.cell_volume = 4.704588099999991e-29
    species = Species()
    species.name = "Cu"
    species.chemical_symbols = ["Cu"]
    species.concentration = [1.0]
    atoms_Cu_conv.species = [species]
    atoms_Cu_conv.lattice_parameters = LatticeParameters()
    atoms_Cu_conv.lattice_parameters.a = 3.609999999999999 * ureg.angstrom
    atoms_Cu_conv.lattice_parameters.b = 3.609999999999999 * ureg.angstrom
    atoms_Cu_conv.lattice_parameters.c = 3.609999999999999 * ureg.angstrom
    atoms_Cu_conv.lattice_parameters.alpha = 1.5707963267948966 * ureg.rad
    atoms_Cu_conv.lattice_parameters.beta = 1.5707963267948966 * ureg.rad
    atoms_Cu_conv.lattice_parameters.gamma = 1.5707963267948966 * ureg.rad
    wyckoff_sets = WyckoffSet()
    wyckoff_sets.wyckoff_letter = "a"
    wyckoff_sets.indices = [0, 1, 2, 3]
    wyckoff_sets.element = "Cu"
    atoms_Cu_conv.wyckoff_sets = [wyckoff_sets]

    Symmetry_fcc = crystal_structure_properties('fcc')
    convsystem_Cu = create_system(label_Cu_conv, structural_type_Cu_conv, elements_Cu, formula_hill_Cu_conv, formula_reduced_Cu_conv, formula_anonymous_Cu_conv, system_relation, material_id=material_id_Cu_conv, atoms=atoms_Cu_conv, cell=Cell_Cu_fcc[0], symmetry=Symmetry_fcc, prototype=prototype_Cu_fcc)
    topologies_Cu.append(convsystem_Cu)
    return topologies_Cu


def single_Cr_surface_atoms():
    Cr_bcc_100 = ase.build.bcc100('Cr', (5, 5, 5), vacuum=10, periodic=True)
    Cr_bcc_100.rattle(stdev=0.001, seed=None, rng=None)
    Cr_bcc_110 = ase.build.bcc110('Cr', (5, 5, 5), vacuum=10, periodic=True)
    Cr_bcc_110.rattle(stdev=0.001, seed=None, rng=None)
    return Cr_bcc_100, Cr_bcc_110


def single_Cr_surface_topology():
    Cell_Cr_bcc = [Cell(a=2.8800000000000014 * ureg.angstrom,
                        b=2.8800000000000014 * ureg.angstrom,
                        c=2.8800000000000014 * ureg.angstrom,
                        alpha=1.5707963267948966 * ureg.rad,
                        beta=1.5707963267948966 * ureg.rad,
                        gamma=1.5707963267948966 * ureg.rad)]
    prototype_Cr_bcc = Prototype()
    prototype_Cr_bcc.aflow_id = "A_cI2_229_a"
    prototype_Cr_bcc.assignment_method = "normalized-wyckoff"
    prototype_Cr_bcc.label = "229-W-cI2"
    prototype_Cr_bcc.formula = "W2"

    # create Cr topology system
    label = 'subsystem'
    structural_type = 'surface'
    elements = ['Cr']
    formula_hill = 'Cr125'
    formula_reduced = 'Cr125'
    formula_anonymous = 'A125'
    system_relation = Relation()
    system_relation.type = "subsystem"
    indices = [i for i in range(75)]

    subsystem_Cr = create_system(label, structural_type, elements, formula_hill, formula_reduced, formula_anonymous, system_relation, indices=indices)
    topologies_Cr = [subsystem_Cr]

    label = 'conventional cell'
    structural_type = 'bulk'
    elements = ['Cr']
    formula_hill = 'Cr2'
    formula_reduced = 'Cr2'
    formula_anonymous = 'A2'
    material_id = "MDlo8h4C2Ppy-kLY9fHRovgnTN9T"
    atoms = Structure()
    atoms.dimension_types = [0, 0, 0]
    atoms.lattice_vectors = [[2.8800000000000014, 0.0, 0.0], [0.0, 2.8800000000000014e-10, 0.0], [0.0, 0.0, 2.8800000000000014]] * ureg.angstrom
    atoms.cartesian_site_positions = [[0.0, 0.0, 0.0], [1.4400000000000007, 1.4400000000000007, 1.4400000000000007]] * ureg.angstrom
    atoms.species_at_sites = ["Cr", "Cr"]
    atoms.cell_volume = 2.388787200000006e-29
    species = Species()
    species.name = "Cr"
    species.chemical_symbols = ["Cr"]
    species.concentration = [1.0]
    atoms.species = [species]
    atoms.lattice_parameters = LatticeParameters()
    atoms.lattice_parameters.a = 2.8800000000000014 * ureg.angstrom
    atoms.lattice_parameters.b = 2.8800000000000014 * ureg.angstrom
    atoms.lattice_parameters.c = 2.8800000000000014 * ureg.angstrom
    atoms.lattice_parameters.alpha = 1.5707963267948966 * ureg.rad
    atoms.lattice_parameters.beta = 1.5707963267948966 * ureg.rad
    atoms.lattice_parameters.gamma = 1.5707963267948966 * ureg.rad
    wyckoff_sets = WyckoffSet()
    wyckoff_sets.wyckoff_letter = "a"
    wyckoff_sets.indices = [0, 1]
    wyckoff_sets.element = "Cr"
    atoms.wyckoff_sets = [wyckoff_sets]

    Symmetry_bcc = crystal_structure_properties('bcc')
    convsystem_Cr = create_system(label, structural_type, elements, formula_hill, formula_reduced, formula_anonymous, system_relation, material_id=material_id, atoms=atoms, cell=Cell_Cr_bcc[0], symmetry=Symmetry_bcc, prototype=prototype_Cr_bcc)
    topologies_Cr.append(convsystem_Cr)
    return topologies_Cr


def single_Ni_surface_topology():
    # create Ni topology
    Cell_Ni_fcc = [Cell(a=3.52 * ureg.angstrom,
                        b=3.52 * ureg.angstrom,
                        c=3.52 * ureg.angstrom,
                        alpha=1.5707963267948966 * ureg.rad,
                        beta=1.5707963267948966 * ureg.rad,
                        gamma=1.5707963267948966 * ureg.rad)]
    label = 'subsystem'
    structural_type = 'surface'
    elements = ['Ni']
    formula_hill = 'Ni75'
    formula_reduced = 'Ni75'
    formula_anonymous = 'A75'
    system_relation = Relation()
    system_relation.type = "subsystem"
    indices = [i for i in range(75, 150)]

    subsystem_Ni = create_system(label, structural_type, elements, formula_hill, formula_reduced, formula_anonymous, system_relation, indices=indices)
    topologies_Ni = [subsystem_Ni]

    label = 'conventional cell'
    structural_type = 'bulk'
    elements = ['Ni']
    formula_hill = 'Ni4'
    formula_reduced = 'Ni4'
    formula_anonymous = 'A4'
    material_id = "NdIWxnQzlp-aeP1IM2d8YJ04h6T0"
    atoms = Structure()
    atoms.dimension_types = [0, 0, 0]
    atoms.lattice_vectors = [[3.52, 0.0, 0.0], [0.0, 3.52, 0.0], [0.0, 0.0, 3.52]] * ureg.angstrom
    atoms.cartesian_site_positions = [[0.0, 0.0, 0.0], [0.0, 1.76, 1.76], [1.76, 0.0, 1.76], [1.76, 1.76, 0.0]] * ureg.angstrom
    atoms.species_at_sites = ["Ni", "Ni", "Ni", "Ni"]
    atoms.cell_volume = 4.3614208000000044e-29
    species = Species()
    species.name = "Ni"
    species.chemical_symbols = ["Ni"]
    species.concentration = [1.0]
    atoms.species = [species]
    atoms.lattice_parameters = LatticeParameters()
    atoms.lattice_parameters.a = 3.52 * ureg.angstrom
    atoms.lattice_parameters.b = 3.52 * ureg.angstrom
    atoms.lattice_parameters.c = 3.52 * ureg.angstrom
    atoms.lattice_parameters.alpha = 1.5707963267948966 * ureg.rad
    atoms.lattice_parameters.beta = 1.5707963267948966 * ureg.rad
    atoms.lattice_parameters.gamma = 1.5707963267948966 * ureg.rad
    wyckoff_sets = WyckoffSet()
    wyckoff_sets.wyckoff_letter = "a"
    wyckoff_sets.indices = [0, 1, 2, 3]
    wyckoff_sets.element = "Ni"
    atoms.wyckoff_sets = [wyckoff_sets]

    Symmetry_fcc = crystal_structure_properties('fcc')
    prototype_Cu_fcc = single_Cu_surface_prototype()
    convsystem_Ni = create_system(label, structural_type, elements, formula_hill, formula_reduced, formula_anonymous, system_relation, material_id=material_id, atoms=atoms, cell=Cell_Ni_fcc[0], symmetry=Symmetry_fcc, prototype=prototype_Cu_fcc)
    topologies_Ni.append(convsystem_Ni)
    return topologies_Ni


def stacked_Cu_Ni_surface():
    # stacked Cu and Ni surface
    Cu_fcc_111 = ase.build.fcc111('Cu', (3, 5, 5), vacuum=None, periodic=False)
    topologies_Cu = single_Cu_surface_topology()

    Ni_fcc_111 = ase.build.fcc111('Ni', (3, 5, 5), vacuum=None, periodic=False)
    topologies_Ni = single_Ni_surface_topology()

    CuNi_fcc_111 = ase.build.stack(Cu_fcc_111, Ni_fcc_111, axis=2, distance=2, maxstrain=2.4)
    CuNi_fcc_111.rattle(stdev=0.001, seed=None, rng=None)

    # create stacked Cu and Ni surface topology
    topologies_Cu_Ni = topologies_Cu + topologies_Ni
    return CuNi_fcc_111, topologies_Cu_Ni


def single_2D_graphene_layer_atoms():
    # Graphene
    symbols_C = ['C', 'C']
    positions_C = [[0.0, 0.0, 2.1712595], [1.2338620706831436, -0.712370598651782, 2.1712595]] * ureg.angstrom
    cell_C = [[1.2338620706831436, -2.137111795955346, 0.0], [1.2338620706831436, 2.137111795955346, 0.0], [0.0, 0.0, 8.685038]] * ureg.angstrom
    system_C = Atoms(
        symbols=symbols_C,
        positions=positions_C,
        cell=cell_C,
        pbc=True
    ) * [4, 4, 1]
    system_C.rattle(stdev=0.001, seed=None, rng=None)
    return system_C


def single_2D_graphene_layer_topology():
    # create graphene topology
    label_C = 'subsystem'
    structural_type_C = '2D'
    elements_C = ['C']
    formula_hill_C = 'C32'
    formula_reduced_C = 'C32'
    formula_anonymous_C = 'A32'
    system_relation = Relation()
    system_relation.type = "subsystem"
    indices_C = [i for i in range(32)]

    subsystem_C = create_system(label_C, structural_type_C, elements_C, formula_hill_C, formula_reduced_C, formula_anonymous_C, system_relation, indices=indices_C)
    topologies_C = [subsystem_C]

    label_C_conv = 'conventional cell'
    structural_type_C_conv = '2D'
    elements_C_conv = ['C']
    formula_hill_C_conv = 'C2'
    formula_reduced_C_conv = 'C2'
    formula_anonymous_C_conv = 'A2'
    material_id_C_conv = "jdP9AhZIFuYhubLWkm2FPtEV5IZA"
    atoms_C_conv = Structure()
    atoms_C_conv.dimension_types = [1, 1, 0]

    atoms_C_conv.lattice_vectors = [[2.4677241413662866, 0.0, 0.0], [-1.2338620706831433, 2.1371117959553457, 0.0], [0.0, 0.0, 1]] * ureg.angstrom
    atoms_C_conv.cartesian_site_positions = [[1.2338620706831433, 0.712370598651782, 0.5], [-2.7636130944313266e-16, 1.4247411973035641, 0.5]] * ureg.angstrom
    atoms_C_conv.species_at_sites = ["C", "C"]
    species = Species()
    species.name = "C"
    species.chemical_symbols = ["C"]
    species.concentration = [1.0]
    atoms_C_conv.species = [species]
    atoms_C_conv.lattice_parameters = LatticeParameters()
    atoms_C_conv.lattice_parameters.a = 2.4677241413662866 * ureg.angstrom
    atoms_C_conv.lattice_parameters.b = 2.4677241413662866 * ureg.angstrom
    atoms_C_conv.lattice_parameters.c = 0
    atoms_C_conv.lattice_parameters.gamma = 2.0943951023931957 * ureg.rad
    atoms_C_conv.wyckoff_sets = None

    Cell_C_conv = [Cell(a=atoms_C_conv.lattice_parameters.a, b=atoms_C_conv.lattice_parameters.b, c=atoms_C_conv.lattice_parameters.c, gamma=atoms_C_conv.lattice_parameters.gamma)]

    convsystem_C = create_system(label_C_conv, structural_type_C_conv, elements_C_conv, formula_hill_C_conv, formula_reduced_C_conv, formula_anonymous_C_conv, system_relation, material_id=material_id_C_conv, atoms=atoms_C_conv, cell=Cell_C_conv[0], symmetry=None, prototype=None)
    topologies_C.append(convsystem_C)
    return topologies_C


def single_2D_BN_layer_atoms():
    # boron nitride
    symbols_BN = ['B', 'N']
    positions_BN = [[1.2557999125000436, -0.7250364175302085, 6.200847], [0.0, 0.0, 6.200847]] * ureg.angstrom
    cell_BN = [[1.2557999125000436, -2.1751092525906257, 0.0], [1.2557999125000436, 2.1751092525906257, 0.0], [0.0, 0.0, 8.267796]] * ureg.angstrom
    system_BN = Atoms(
        symbols=symbols_BN,
        positions=positions_BN,
        cell=cell_BN,
        pbc=True
    )
    BN_2D = ase.build.surface(system_BN, (0, 0, 1), layers=1, periodic=True)
    BN_2 = ase.build.stack(BN_2D, BN_2D, axis=0)
    BN_4 = ase.build.stack(BN_2, BN_2, axis=1)
    BN_8 = ase.build.stack(BN_4, BN_4, axis=0)
    BN_16 = ase.build.stack(BN_8, BN_8, axis=1)
    BN_16.rattle(stdev=0.001, seed=None, rng=None)
    return BN_16


def single_2D_BN_layer_topology():
    # create boron nitrid topology
    label = 'subsystem'
    structural_type = '2D'
    elements = ['B', 'N']
    formula_hill = 'B16N16'
    formula_reduced = 'B16N16'
    formula_anonymous = 'A16B16'
    system_relation = Relation()
    system_relation.type = "subsystem"
    indices_BN = [i for i in range(32)]

    subsystem_BN = create_system(label, structural_type, elements, formula_hill, formula_reduced, formula_anonymous, system_relation, indices=indices_BN)
    topologies_BN = [subsystem_BN]

    label = 'conventional cell'
    structural_type = '2D'
    elements = ['B', 'N']
    formula_hill = 'BN'
    formula_reduced = 'BN'
    formula_anonymous = 'AB'
    material_id = "RxRsol0dp1vDkU7-pE3v2exglkpM"
    atoms = Structure()
    atoms.dimension_types = [1, 1, 0]
    atoms.lattice_vectors = [[2.510266994011973, 0.0, 0.0], [-1.2551334970059864, 2.1739549870959678, 0.0], [0.0, 0.0, 1]] * ureg.angstrom
    atoms.species_at_sites = ["B", "N"]
    species_B = Species()
    species_B.name = "B"
    species_B.chemical_symbols = ["B"]
    species_B.concentration = [1.0]
    species_N = Species()
    species_N.name = "N"
    species_N.chemical_symbols = ["N"]
    species_N.concentration = [1.0]
    atoms.species = [species_B, species_N]
    atoms.lattice_parameters = LatticeParameters()
    atoms.lattice_parameters.a = 2.510266994011973 * ureg.angstrom
    atoms.lattice_parameters.b = 2.510266994011973 * ureg.angstrom
    atoms.lattice_parameters.c = 0
    atoms.lattice_parameters.gamma = 2.0943951023931957 * ureg.rad
    atoms.wyckoff_sets = None

    Cell_BN = [Cell(a=atoms.lattice_parameters.a, b=atoms.lattice_parameters.b, c=atoms.lattice_parameters.c, alpha=atoms.lattice_parameters.alpha, beta=atoms.lattice_parameters.beta, gamma=atoms.lattice_parameters.gamma)]
    convsystem_BN = create_system(label, structural_type, elements, formula_hill, formula_reduced, formula_anonymous, system_relation, material_id=material_id, atoms=atoms, cell=Cell_BN[0], symmetry=None, prototype=None)
    topologies_BN.append(convsystem_BN)
    return topologies_BN


def single_2D_MoS2_layer_atoms():
    # MoS2
    symbols_MoS2 = ['Mo', 'S', 'S']
    positions_MoS2 = [[0.0, 0.0, 9.063556323175761], [1.5920332323422965, 0.9191608152516547, 10.62711264635152], [1.5920332323422965, 0.9191608152516547, 7.5]] * ureg.angstrom
    cell_MoS2 = [[3.184066464684593, 0.0, 0.0], [-1.5920332323422965, 2.7574824457549643, 0.0], [0.0, 0.0, 18.127112646351521]] * ureg.angstrom
    system_MoS2 = Atoms(
        symbols=symbols_MoS2,
        positions=positions_MoS2,
        cell=cell_MoS2,
        pbc=True
    )
    MoS2_2D = ase.build.surface(system_MoS2, (1, 1, 0), layers=4, vacuum=None, periodic=True)
    stacked_2D_MoS2 = ase.build.stack(MoS2_2D, MoS2_2D, axis=2, distance=2.5)
    stacked_2D_MoS2_2 = ase.build.stack(stacked_2D_MoS2, stacked_2D_MoS2, axis=2)
    stacked_2D_MoS2_2.rattle(stdev=0.001, seed=None, rng=None)
    return stacked_2D_MoS2_2


def single_2D_MoS2_layer_topology():
    # create MoS2 topology
    label = 'subsystem'
    structural_type = '2D'
    elements = ['Mo', 'S']
    formula_hill = 'Mo16S32'
    formula_reduced = 'Mo16S32'
    formula_anonymous = 'A32B16'
    system_relation = Relation()
    system_relation.type = "subsystem"
    indices_MoS2 = [i for i in range(48)]
    subsystem_MoS2 = create_system(label, structural_type, elements, formula_hill, formula_reduced, formula_anonymous, system_relation, indices=indices_MoS2)
    topologies_MoS2 = [subsystem_MoS2]

    label = 'conventional cell'
    structural_type = '2D'
    elements = ['Mo', 'S']
    formula_hill = 'MoS2'
    formula_reduced = 'MoS2'
    formula_anonymous = 'A2B'
    material_id = "KV4aYm-S1VJOH-SKeXXuG8JkTiGF"
    atoms = Structure()
    atoms.dimension_types = [1, 1, 0]
    atoms.lattice_vectors = [[3.253646631826119, 0.0, 0.0], [-1.6268233159130596, 2.8177406380990937, 0.0], [0.0, 0.0, 3.124912396241947]] * ureg.angstrom
    atoms.cartesian_site_positions = [[0.0, 0.0, 1.562456198120974], [1.626823332181293, 0.9392468699738958, 3.124912396241947], [1.626823332181293, 0.9392468699738958, 0.0]] * ureg.angstrom
    atoms.species_at_sites = ["Mo", "S", "S"]
    species_Mo = Species()
    species_Mo.name = "Mo"
    species_Mo.chemical_symbols = ["Mo"]
    species_Mo.concentration = [1.0]
    species_S = Species()
    species_S.name = "S"
    species_S.chemical_symbols = ["S"]
    species_S.concentration = [1.0]
    atoms.species = [species_Mo, species_S]
    atoms.lattice_parameters = LatticeParameters()
    atoms.lattice_parameters.a = 3.253646631826119 * ureg.angstrom
    atoms.lattice_parameters.b = 3.253646631826119 * ureg.angstrom
    atoms.lattice_parameters.c = 0.
    atoms.lattice_parameters.gamma = 2.0943951023931957 * ureg.rad
    atoms.wyckoff_sets = None
    Cell_MoS2 = [Cell(a=atoms.lattice_parameters.a, b=atoms.lattice_parameters.b, c=atoms.lattice_parameters.c, gamma=atoms.lattice_parameters.gamma)]
    convsystem_MoS2 = create_system(label, structural_type, elements, formula_hill, formula_reduced, formula_anonymous, system_relation, material_id=material_id, atoms=atoms, cell=Cell_MoS2[0], symmetry=None, prototype=None)
    topologies_MoS2.append(convsystem_MoS2)
    return topologies_MoS2


def stacked_C_BN_2D_layers():
    # stacked 2D of C and BN
    C_32 = single_2D_graphene_layer_atoms()
    topologies_C = single_2D_graphene_layer_topology()
    # create stacked C and BN topologies
    topologies_C[0].indices = [i for i in range(32, 64)]
    BN_16 = single_2D_BN_layer_atoms()
    stacked_C_BN = ase.build.stack(BN_16, C_32, axis=2, maxstrain=6.7)
    stacked_C_BN = ase.build.surface(stacked_C_BN, (0, 0, 1), layers=1, vacuum=10, periodic=True)
    stacked_C_BN.rattle(stdev=0.001, seed=None, rng=None)
    topologies_BN = single_2D_BN_layer_topology()
    topologies_C_BN = topologies_C + topologies_BN
    return stacked_C_BN, topologies_C_BN
