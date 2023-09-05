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
from nomad.units import ureg
import pytest


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


def test_method_k_mesh(bulk):
    """K-point mesh from a DFT calculation."""
    k_mesh = bulk.run[-1].method[-1].k_mesh
    assert k_mesh.n_points == 1
    assert k_mesh.grid.tolist() == [1, 1, 1]
    assert k_mesh.points.tolist() == [[0, 0, 0]]


def test_precision(bulk):
    """Precision from a bulk DFT calculation."""
    precision = bulk.results.method.simulation.precision
    lattice_length = bulk.run[0].system[0].atoms.lattice_vectors.magnitude[0, 0]
    assert precision.k_line_density.magnitude == lattice_length / (2 * np.pi)


def test_method_dft(dft):
    """Methodology from a DFT calculation."""
    method = dft.results.method
    assert method.method_name == "DFT"
    assert method.simulation.program_name == "VASP"
    assert method.simulation.program_version == "4.6.35"
    assert method.simulation.dft.basis_set_type == "plane waves"
    assert method.simulation.dft.core_electron_treatment == "pseudopotential"
    assert method.simulation.dft.xc_functional_names == ["GGA_C_PBE", "GGA_X_PBE"]
    assert method.simulation.dft.xc_functional_type == "GGA"
    assert method.simulation.dft.smearing_kind == "gaussian"
    assert method.simulation.dft.smearing_width == 1e-20
    assert method.simulation.dft.spin_polarized is True
    assert method.simulation.dft.scf_threshold_energy_change == 1e-24 * ureg.joule
    assert method.simulation.dft.van_der_Waals_method == "G06"
    assert method.simulation.dft.relativity_method == "scalar_relativistic"


def test_method_referenced(dft_method_referenced):
    """Methodology from a calculation which uses references to tie together
    several methodologies.
    """
    method = dft_method_referenced.results.method
    assert method.method_name == "DFT"
    assert method.simulation.program_name == "VASP"
    assert method.simulation.program_version == "4.6.35"
    assert method.simulation.dft.basis_set_type == "plane waves"
    assert method.simulation.dft.core_electron_treatment == "pseudopotential"
    assert method.simulation.dft.xc_functional_names == ["GGA_C_PBE", "GGA_X_PBE"]
    assert method.simulation.dft.xc_functional_type == "GGA"
    assert method.simulation.dft.smearing_kind == "gaussian"
    assert method.simulation.dft.smearing_width == 1e-20
    assert method.simulation.dft.spin_polarized is True
    assert method.simulation.dft.scf_threshold_energy_change == 1e-24 * ureg.joule
    assert method.simulation.dft.van_der_Waals_method == "G06"
    assert method.simulation.dft.relativity_method == "scalar_relativistic"


@pytest.mark.parametrize('entry, expected', [
    ('dft_exact_exchange', .25),
    ('dft_b3lyp', .2),
    ('dft_pbeh', .25),
    ('dft_m05', .28),
    ('dft_pbe0_13', 1 / 3),
    ('dft_pbe38', 3 / 8),
    ('dft_pbe50', .5),
    ('dft_m06_2x', .54),
    ('dft_m05_2x', .56),
])
def test_exact_exchange_mixing_factor(entry, expected, request):
    """Exact exchange mixing factor in hybrid functionals"""
    entry = request.getfixturevalue(entry)
    assert entry.results.method.simulation.dft.exact_exchange_mixing_factor == expected


@pytest.mark.parametrize('entry, expected', [
    pytest.param('dft_empty', 'unavailable',
                 id='DFT with no XC functional specification'),
    pytest.param('dft_wrong', 'unavailable',
                 id='DFT with non-sensical XC functional specification'),
    pytest.param('dft_pw', 'LDA', id='PW functional'),
    pytest.param('dft', 'GGA', id='PBE functional'),
    pytest.param('dft_m06', 'meta-GGA', id='M06 functional'),
    pytest.param('dft_m05', 'hyper-GGA', id='M05 functional'),
    pytest.param('dft_b3lyp', 'hybrid', id='B3LYP functional'),
])
def test_jacobs_ladder_value(entry, expected, request):
    """Test assignment of the rung on Jacob's ladder."""
    entry = request.getfixturevalue(entry)
    assert entry.results.method.simulation.dft.jacobs_ladder == expected
    assert entry.results.method.simulation.dft.xc_functional_type == expected


def test_method_dft_plus_u(dft_plus_u):
    """Methodology from a DFT+U calculation with a Hubbard model."""
    method = dft_plus_u.results.method
    assert method.method_name == "DFT"
    assert method.simulation.program_name == "VASP"
    assert method.simulation.program_version == "4.6.35"
    assert method.simulation.dft.basis_set_type == "plane waves"
    assert method.simulation.dft.core_electron_treatment == "pseudopotential"
    assert method.simulation.dft.xc_functional_names == ["GGA_C_PBE", "GGA_X_PBE"]
    assert method.simulation.dft.xc_functional_type == "GGA"
    assert method.simulation.dft.smearing_kind == "gaussian"
    assert method.simulation.dft.smearing_width == 1e-20
    assert method.simulation.dft.spin_polarized is True
    assert method.simulation.dft.scf_threshold_energy_change == 1e-24 * ureg.joule
    assert method.simulation.dft.van_der_Waals_method == "G06"
    assert method.simulation.dft.relativity_method == "scalar_relativistic"
    assert method.simulation.dft.hubbard_kanamori_model[0].atom_label == 'Ti'
    assert method.simulation.dft.hubbard_kanamori_model[0].orbital == '3d'
    assert method.simulation.dft.hubbard_kanamori_model[0].u_effective.magnitude == 3.5e-19
    assert method.simulation.dft.hubbard_kanamori_model[0].u.magnitude == 4.5e-19
    assert method.simulation.dft.hubbard_kanamori_model[0].j.magnitude == 1e-19
    assert method.simulation.dft.hubbard_kanamori_model[0].double_counting_correction == 'Dudarev'


def test_method_projection(projection):
    """Methodology from a Projection calculation"""
    method = projection.results.method
    assert method.method_name == "Projection"
    assert method.simulation.program_name == "Wannier90"
    assert method.simulation.program_version == "3.1.0"
    assert method.simulation.projection.type == 'wannier'
    assert method.simulation.projection.localization_type == "single_shot"


def test_method_gw(gw):
    """Methodology from a SinglePoint GW calculation."""
    method = gw.results.method
    assert method.method_name == "GW"
    assert method.simulation.program_name == "VASP"
    assert method.simulation.program_version == "4.6.35"
    assert method.simulation.gw.type == "G0W0"


def test_method_bse(bse):
    """Methodology from a SinglePoint GW calculation."""
    method = bse.results.method
    assert method.method_name == "BSE"
    assert method.simulation.program_name == "VASP"
    assert method.simulation.program_version == "4.6.35"
    assert method.simulation.bse.type == "Singlet"
    assert method.simulation.bse.solver == "Lanczos-Haydock"


def test_method_dmft(dmft):
    """Methodology from a SinglePoint DMFT calculation"""
    method = dmft.results.method
    assert method.method_name == "DMFT"
    assert method.simulation.program_name == "w2dynamics"
    assert method.simulation.dmft.impurity_solver_type == "CT-HYB"
    assert method.simulation.dmft.inverse_temperature.magnitude == 60.0
    assert method.simulation.dmft.magnetic_state == "paramagnetic"
    assert method.simulation.dmft.u.magnitude == 4.0e-19
    assert method.simulation.dmft.jh.magnitude == 0.6e-19


def test_method_eels(eels):
    method = eels.results.method
    assert method.method_name == "EELS"


def test_basis_set(normalized_example):
    ref_basis_sets = {
        'exciting': ['APW', 'LAPW', 'APW+lo', 'LAPW+lo', '(L)APW', '(L)APW+lo'],
        'fleur': ['APW', 'LAPW', 'APW+lo', 'LAPW+lo', '(L)APW', '(L)APW+lo'],  # TODO: does not match the metainfo name
        'WIEN2k': ['APW', 'LAPW', 'APW+lo', 'LAPW+lo', '(L)APW', '(L)APW+lo'],
        'Elk': ['APW', 'LAPW', 'APW+lo', 'LAPW+lo', '(L)APW', '(L)APW+lo'],
        'ABINIT': ['plane waves'],
        'AFlow': ['plane waves'],
        'atomate': ['plane waves'],
        'CASTEP': ['plane waves'],
        'CPMD': ['plane waves'],
        'MaterialsProject': ['plane waves'],
        'Quantum Espresso': ['plane waves'],
        'VASP': ['plane waves'],
        'AMS': ['atom-centered orbitals'],
        'Crystal': ['atom-centered orbitals'],
        'DMol3': ['atom-centered orbitals'],
        'FHI-aims': ['atom-centered orbitals'],
        'FHI-vibes': ['atom-centered orbitals'],
        'GAMESS': ['atom-centered orbitals'],
        'Gaussian': ['atom-centered orbitals'],
        'LOBSTER': ['atom-centered orbitals'],
        'NWChem': ['atom-centered orbitals'],
        'ORCA': ['atom-centered orbitals'],
        'Psi4': ['atom-centered orbitals'],
        'turbomole': ['atom-centered orbitals'],
        'CP2K': ['gaussians + plane waves', 'plane waves'],
        'BigDFT': ['real-space grid'],
        'Octopus': ['real-space grid'],
        'GPAW': ['real-space grid', 'plane waves', 'atom-centered orbitals'],
        'Siesta': ['real-space grid', 'atom-centered orbitals'],
        'ONETEP': ['support functions'],
        'DL_POLY_4': [],
        'Phonopy': [],
        'ATK': [],
        'gulp': [],
        'elastic': [],
        'qbox': [],
        'YAMBO': [],
    }
    codes_withouth_methods = ['libAtoms', 'MOLCAS', 'ASR']
    try:
        base = normalized_example.results.method.simulation
    except AttributeError:
        return
    if base is not None:
        program_name = base.program_name
        if program_name in codes_withouth_methods:
            assert not base.precision  # precision is generated only if run.method exists
        elif program_name != 'not processed':
            reference = ref_basis_sets[program_name]
            if basis_set := base.precision.basis_set:
                assert basis_set in reference + ['unavailable']


@pytest.mark.parametrize('entry, method_identified', [
    ('hash_exciting', True),
    ('hash_vasp', False)
])
def test_method_id(entry, method_identified, request):
    """Test that method_id can be detected or is left undetected from certain
    calculations.
    """
    entry = request.getfixturevalue(entry)
    assert (entry.results.method.method_id is not None) == method_identified
    assert (entry.results.method.equation_of_state_id is not None) == method_identified
    assert (entry.results.method.parameter_variation_id is not None) == method_identified
