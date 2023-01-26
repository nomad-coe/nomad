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

from nomad.units import ureg


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
    ('dft_m05_2x', .56)
])
def test_exact_exchange_mixing_factor(entry, expected, request):
    """Exact exchange mixing factor in hybrid functionals"""
    entry = request.getfixturevalue(entry)
    assert entry.results.method.simulation.dft.exact_exchange_mixing_factor == expected


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
    """Methodology from a GW calculation."""
    method = gw.results.method
    assert method.method_name == "GW"
    assert method.simulation.program_name == "VASP"
    assert method.simulation.program_version == "4.6.35"
    assert method.simulation.gw.type == "G0W0"
    assert method.simulation.gw.starting_point_names == ["GGA_C_PBE", "GGA_X_PBE"]
    assert method.simulation.gw.starting_point_type == "GGA"
    assert method.simulation.gw.basis_set_type == "plane waves"


def test_method_dmft(dmft):
    """Methodology from a DMFT calculation"""
    method = dmft.results.method
    assert method.method_name == "DMFT"
    assert method.simulation.program_name == "w2dynamics"
    assert method.simulation.dmft.impurity_solver_type == "CT-HYB"
    assert method.simulation.dmft.total_filling == 0.5 / 3.0
    assert method.simulation.dmft.inverse_temperature.magnitude == 60.0
    assert method.simulation.dmft.magnetic_state == "paramagnetic"
    assert method.simulation.dmft.u.magnitude == 4.0e-19
    assert method.simulation.dmft.hunds_hubbard_ratio == 0.6e-19 / 4.0e-19


def test_method_eels(eels):
    method = eels.results.method
    assert method.method_name == "EELS"


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
