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
import numpy as np

from nomad.parsing import LocalBackend
from nomad.metainfo.encyclopedia import Encyclopedia
from tests.normalizing.conftest import geometry_optimization, molecular_dynamics, phonon, two_d, bulk   # pylint: disable=unused-import


def test_geometry_optimization(geometry_optimization: LocalBackend):
    """Tests that geometry optimizations are correctly processed."
    """
    enc = geometry_optimization.get_mi2_section(Encyclopedia.m_def)
    run_type = enc.calculation.run_type
    assert run_type == "geometry optimization"


def test_molecular_dynamics(molecular_dynamics: LocalBackend):
    """Tests that geometry optimizations are correctly processed."
    """
    enc = molecular_dynamics.get_mi2_section(Encyclopedia.m_def)
    run_type = enc.calculation.run_type
    assert run_type == "molecular dynamics"


def test_phonon(phonon: LocalBackend):
    """Tests that geometry optimizations are correctly processed."
    """
    enc = phonon.get_mi2_section(Encyclopedia.m_def)
    run_type = enc.calculation.run_type
    assert run_type == "phonon calculation"


def test_system_type(geometry_optimization: LocalBackend):
    """Tests that geometry optimizations are correctly processed.
    """
    enc = geometry_optimization.get_mi2_section(Encyclopedia.m_def)
    system_type = enc.material.system_type
    assert system_type == "bulk"


def test_bulk_metainfo(bulk: LocalBackend):
    """Tests that information for bulk systems is correctly processed.
    """
    enc = bulk.get_mi2_section(Encyclopedia.m_def)
    assert enc.material.system_type == "bulk"
    assert enc.material.number_of_atoms == 4
    assert enc.material.atom_labels == ["Na", "Na", "Na", "Na"]
    assert enc.material.atom_positions is not None
    assert enc.material.crystal_system == "cubic"
    assert enc.material.bravais_lattice == "cF"
    assert enc.material.formula == "Na"
    assert enc.material.formula_reduced == "Na"
    assert enc.material.has_free_wyckoff_parameters is False
    assert enc.material.material_name == "Sodium"
    assert enc.material.point_group == "m-3m"
    assert enc.material.cell_normalized is not None
    assert enc.material.cell_primitive is not None
    assert np.array_equal(enc.material.periodicity, [0, 1, 2])
    assert enc.material.wyckoff_groups is not None

    assert enc.calculation.atomic_density == pytest.approx(4.0e+30, rel=0.000001, abs=None)
    assert enc.calculation.lattice_parameters is not None
    assert enc.calculation.mass_density == 4 * 22.98976928 * 1.6605389e-27 / 1e-30  # Atomic mass in kg / cell volume
    assert enc.calculation.cell_volume == 1e-30


def test_2d_metainfo(two_d: LocalBackend):
    """Tests that information for 2D systems is correctly processed.
    """
    enc = two_d.get_mi2_section(Encyclopedia.m_def)
    assert enc.material.system_type == "2D"
    assert enc.material.number_of_atoms == 2
    assert enc.material.atom_labels == ["C", "C"]
    assert enc.material.atom_positions is not None
    assert enc.material.cell_normalized is not None
    assert enc.material.cell_primitive is not None
    assert enc.material.formula == "C2"
    assert enc.material.formula_reduced == "C"
    assert np.allclose(enc.calculation.lattice_parameters, [2.46559821e-10, 2.46559821e-10, 0, 120 / 180 * np.pi, 0, 0])
    assert np.array_equal(enc.material.periodicity, [0, 1])
