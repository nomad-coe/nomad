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
import ase

from tests.normalizing.conftest import get_template_for_structure
from tests.normalizing.test_topology import assert_topology


@pytest.mark.skip(reason='Porosity normalizer is not run by default.')
@pytest.mark.parametrize(
    'filepath',
    [
        'tests/data/normalizers/porous_systems/COF-1.cif',
        'tests/data/normalizers/porous_systems/IRR.cif',
        'tests/data/normalizers/porous_systems/SARSUC.cif',
    ],
)
def test_porosity(filepath):
    atoms = ase.io.read(filepath)
    archive = get_template_for_structure(atoms)
    topology = archive.results.material.topology
    assert_topology(topology)
    topology_porous = [
        top for top in archive.results.material.topology if top.method == 'porosity'
    ]
    porous_system = topology_porous[0]
    assert porous_system.largest_cavity_diameter != 0
    assert porous_system.pore_limiting_diameter != 0
    assert porous_system.largest_included_sphere_along_free_sphere_path != 0
    assert porous_system.accessible_surface_area != 0
    assert porous_system.accessible_volume != 0
    assert porous_system.void_fraction != 0
    assert porous_system.n_channels > 0


@pytest.mark.skip(reason='Porosity normalizer is not run by default.')
@pytest.mark.parametrize(
    'filepath',
    [
        'tests/data/normalizers/mofs/EDUSIF.cif',
        'tests/data/normalizers/mofs/RUBTAK01.cif',
        'tests/data/normalizers/mofs/SARSUC.cif',
    ],
)
def test_mof(filepath):
    atoms = ase.io.read(filepath)
    archive = get_template_for_structure(atoms)
    topology = archive.results.material.topology
    assert_topology(topology)
    topology_porous = [
        top for top in archive.results.material.topology if top.method == 'porosity'
    ]
    porous_system = topology_porous[0]
    metal_sbus = [
        m_sbu.label for m_sbu in topology_porous[1:] if 'metal_sbu' in m_sbu.label
    ]
    organic_sbus = [
        m_sbu.label for m_sbu in topology_porous[1:] if 'organic_sbu' in m_sbu.label
    ]
    ligands = [
        m_sbu.label for m_sbu in topology_porous[1:] if 'organic_ligand' in m_sbu.label
    ]
    assert porous_system.largest_cavity_diameter != 0
    assert porous_system.pore_limiting_diameter != 0
    assert porous_system.largest_included_sphere_along_free_sphere_path != 0
    assert porous_system.accessible_surface_area != 0
    assert porous_system.accessible_volume != 0
    assert porous_system.void_fraction != 0
    assert porous_system.n_channels != 0
    assert len(metal_sbus) > 0
    assert len(organic_sbus) > 0
    assert len(ligands) > 0
