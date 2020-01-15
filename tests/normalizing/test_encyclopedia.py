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

from nomad.metainfo.encyclopedia import Encyclopedia
from tests.normalizing.conftest import geometry_optimization, molecular_dynamics, phonon   # pylint: disable=unused-import


def test_geometry_optimization(geometry_optimization: Encyclopedia):
    """Tests that geometry optimizations are correctly processed."
    """
    run_type = geometry_optimization.calculation.run_type
    assert run_type == "geometry optimization"


def test_molecular_dynamics(molecular_dynamics: Encyclopedia):
    """Tests that geometry optimizations are correctly processed."
    """
    run_type = molecular_dynamics.calculation.run_type
    assert run_type == "molecular dynamics"


def test_phonon(phonon: Encyclopedia):
    """Tests that geometry optimizations are correctly processed."
    """
    run_type = phonon.calculation.run_type
    assert run_type == "phonon calculation"


def test_system_type(geometry_optimization: Encyclopedia):
    """Tests that geometry optimizations are correctly processed."
    """
    system_type = geometry_optimization.material.system_type
    assert system_type == "bulk"


def test_bulk_information(geometry_optimization: Encyclopedia):
    """Tests that information for bulk systems is correctly processed."
    """
    go = geometry_optimization

    assert go.material.system_type == "bulk"
    assert go.material.number_of_atoms == 4
    assert go.material.atom_labels == ["Na", "Na", "Na", "Na"]
    assert go.material.bravais_lattice == "cF"
    assert go.material.point_group == "m-3m"
    assert go.material.cell_normalized is not None
    assert go.material.cell_primitive is not None
    assert np.array_equal(go.material.periodicity, [0, 1, 2])

    assert go.calculation.atomic_density == pytest.approx(4.0e+30, rel=0.000001, abs=None)
    assert go.calculation.lattice_parameters is not None
    assert go.calculation.cell_angles_string is not None
    assert go.calculation.mass_density == 4 * 22.98976928 * 1.6605389e-27 / 1e-30  # Atomic mass in kg / cell volume
