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
from collections import defaultdict
import pytest

from tests.normalizing.conftest import (
    get_template_for_structure,
    get_template_topology,
    single_Cu_surface_atoms,
    single_Cu_surface_topology,
    single_Cr_surface_atoms,
    single_Cr_surface_topology,
    stacked_Cu_Ni_surface,
    single_2D_graphene_layer_atoms,
    single_2D_graphene_layer_topology,
    single_2D_BN_layer_atoms,
    single_2D_BN_layer_topology,
    single_2D_MoS2_layer_atoms,
    single_2D_MoS2_layer_topology,
    stacked_C_BN_2D_layers
)
from tests.normalizing.conftest import projection  # pylint: disable=unused-import


def assert_topology(topology):
    '''Checks that the given topology has a valid structure.'''
    child_map = {}
    child_map_determined = defaultdict(list)
    for top in topology:
        assert top.system_id is not None
        assert top.parent_system is not None or top.label == 'original'
        assert top.atoms is not None or top.atoms_ref is not None or top.indices is not None
        assert top.n_atoms is not None
        assert top.elements is not None
        assert top.n_elements is not None
        assert top.chemical_formula_hill is not None
        assert top.chemical_formula_reduced is not None
        assert top.chemical_formula_anonymous is not None
        if top.parent_system:
            child_map_determined[top.parent_system].append(top.system_id)
        if top.child_systems:
            child_map[top.system_id] = top.child_systems

    assert len(child_map) == len(child_map_determined)
    for key in child_map.keys():
        assert child_map[key] == child_map_determined[key]


@pytest.mark.parametrize(
    'pbc',
    [
        pytest.param(True, id='fully periodic'),
        pytest.param(False, id='non-periodic'),
    ]
)
def test_topology_calculation(pbc):
    '''Tests that a topology that originates from the calculation itself is
    correctly extracted.
    '''
    topology_calculation = get_template_topology(pbc)
    topology = topology_calculation.results.material.topology
    assert_topology(topology)
    assert len(topology) == 5

    # Test the original structure
    original = topology[0]
    assert original.structural_type == 'unavailable'
    assert original.atoms_ref.positions.shape == (6, 3)
    assert len(original.atoms_ref.labels) == 6
    assert original.atoms_ref.lattice_vectors.shape == (3, 3)
    expected_pbc = np.zeros(3, bool)
    expected_pbc[:] = pbc
    assert original.atoms_ref.periodic == expected_pbc.tolist()
    assert original.chemical_formula_hill == 'H4O2'
    assert original.chemical_formula_reduced == 'H2O'
    assert original.chemical_formula_anonymous == 'A2B'
    assert original.elements == ['H', 'O']
    assert original.n_elements == 2
    assert original.n_atoms == 6
    assert original.parent_system is None
    assert original.child_systems == ['results/material/topology/1']

    # Test molecule group
    mol_group = topology[1]
    assert mol_group.structural_type == 'group'
    assert np.array_equal(mol_group.indices, [[0, 1, 2, 3, 4, 5]])
    assert original.chemical_formula_hill == 'H4O2'
    assert original.chemical_formula_reduced == 'H2O'
    assert original.chemical_formula_anonymous == 'A2B'
    assert mol_group.elements == ['H', 'O']
    assert mol_group.n_elements == 2
    assert mol_group.n_atoms == 6
    assert mol_group.parent_system == 'results/material/topology/0'
    assert mol_group.child_systems == ['results/material/topology/2']

    # Test molecule
    mol = topology[2]
    assert mol.structural_type == 'molecule'
    assert np.array_equal(mol.indices, [[0, 1, 2], [3, 4, 5]])
    assert mol.chemical_formula_hill == 'H2O'
    assert mol.chemical_formula_reduced == 'H2O'
    assert mol.chemical_formula_anonymous == 'A2B'
    assert mol.elements == ['H', 'O']
    assert mol.n_elements == 2
    assert mol.n_atoms == 3
    assert mol.parent_system == 'results/material/topology/1'
    assert mol.child_systems == ['results/material/topology/3']

    # Test monomer group
    mon_group = topology[3]
    assert mon_group.structural_type == 'group'
    assert np.array_equal(mon_group.indices, [[1, 2]])
    assert mon_group.chemical_formula_hill == 'H2'
    assert mon_group.chemical_formula_reduced == 'H'
    assert mon_group.chemical_formula_anonymous == 'A'
    assert mon_group.elements == ['H']
    assert mon_group.n_elements == 1
    assert mon_group.n_atoms == 2
    assert mon_group.parent_system == 'results/material/topology/2'
    assert mon_group.child_systems == ['results/material/topology/4']

    # Test monomer
    mon = topology[4]
    assert mon.structural_type == 'monomer'
    assert np.array_equal(mon.indices, [[1, 2]])
    assert mon.chemical_formula_hill == 'H2'
    assert mon.chemical_formula_reduced == 'H'
    assert mon.chemical_formula_anonymous == 'A'
    assert mon.elements == ['H']
    assert mon.n_elements == 1
    assert mon.n_atoms == 2
    assert mon.parent_system == 'results/material/topology/3'
    assert mon.child_systems is None


@pytest.mark.parametrize('fixture', [
    pytest.param('atom', id='atom'),
    pytest.param('molecule', id='molecule'),
    pytest.param('one_d', id='1D'),
    pytest.param('bulk', id='bulk'),
])
def test_no_topology(fixture, request):
    # Test that some entries don't get a topology. This will changed later, but
    # for now we only create topologies for a subset of systems.
    entry = request.getfixturevalue(fixture)
    assert not entry.results.material.topology


@pytest.mark.skip
@pytest.mark.parametrize('surface, ref_topologies', [
    pytest.param(single_Cu_surface_atoms()[0], single_Cu_surface_topology(),
                 id='single surface Cu FCC 100'),
    pytest.param(single_Cu_surface_atoms()[1], single_Cu_surface_topology(),
                 id='single surface Cu FCC 110'),
    pytest.param(single_Cr_surface_atoms()[0], single_Cr_surface_topology(),
                 id='single surface Cr BCC 100'),
    pytest.param(single_Cr_surface_atoms()[1], single_Cr_surface_topology(),
                 id='single surface Cr BCC 110'),
    pytest.param(stacked_Cu_Ni_surface()[0], stacked_Cu_Ni_surface()[1],
                 id='stacked surfaces of Cu and Ni'),
    pytest.param(single_2D_graphene_layer_atoms(), single_2D_graphene_layer_topology(),
                 id='single 2D layer of graphene'),
    pytest.param(single_2D_BN_layer_atoms(), single_2D_BN_layer_topology(),
                 id='single 2D layer of BN'),
    pytest.param(single_2D_MoS2_layer_atoms(), single_2D_MoS2_layer_topology(),
                 id='single 2D layer of MoS2'),
    pytest.param(stacked_C_BN_2D_layers()[0], stacked_C_BN_2D_layers()[1],
                 id='stacked layer of BN and C')
])
def test_surface_2D_topology(surface, ref_topologies):
    entry_archive = get_template_for_structure(surface)
    topology = entry_archive.results.material.topology
    assert_topology(topology)
    subsystem_topologies = topology[1:]
    # Compare topology with reference system topology. topology[0] is the original system
    assert len(subsystem_topologies) == len(ref_topologies)
    for subsystem_topology in subsystem_topologies:
        formula_hill = subsystem_topology.chemical_formula_hill
        for ref_top_counter, ref_topology in enumerate(ref_topologies):
            if ref_topology.chemical_formula_hill == formula_hill:
                ref_formula_hill = ref_topology.chemical_formula_hill
                ref_index = ref_top_counter
                break
        ref_elements = ref_topologies[ref_index].elements
        elements = subsystem_topology.elements
        assert elements == ref_elements
        assert formula_hill == ref_formula_hill

        ref_structural_type = ref_topologies[ref_index].structural_type
        structural_type = subsystem_topology.structural_type
        assert ref_structural_type == structural_type

        if subsystem_topology.label == 'conventional cell':
            # Cell
            ref_cell = ref_topologies[ref_index].cell
            cell = subsystem_topology.cell
            if ref_structural_type == '2D':
                assert np.allclose(list(cell.values())[:4], list(ref_cell.values()), rtol=1e-05, atol=1e-9)
            else:
                assert np.allclose(list(cell.values())[:6], list(ref_cell.values()), rtol=1e-05, atol=1e-9)

            # Symmetry
            if ref_topologies[ref_index].symmetry:
                symmetry = subsystem_topology.symmetry.m_to_dict()
                ref_symmetry = ref_topologies[ref_index].symmetry.m_to_dict()
                for ref_symmetry_property_key, ref_symmetry_property in ref_symmetry.items():
                    symmetry_property = symmetry[ref_symmetry_property_key]
                    assert ref_symmetry_property == symmetry_property
            else:
                assert subsystem_topology.symmetry == ref_topologies[ref_index].symmetry

            # Prototype
            if ref_topologies[ref_index].prototype:
                prototype = subsystem_topology.prototype.m_to_dict()
                ref_prototype = ref_topologies[ref_index].prototype.m_to_dict()
                for ref_prototype_property_key, ref_prototype_property in ref_prototype.items():
                    prototype_property = prototype[ref_prototype_property_key]
                    assert ref_prototype_property == prototype_property
            else:
                assert ref_topologies[ref_index].prototype == subsystem_topology.prototype

            # Atoms
            atoms = subsystem_topology.atoms.m_to_dict()
            ref_atoms = ref_topologies[ref_index].atoms.m_to_dict()
            for ref_atoms_property_key, ref_atoms_property in ref_atoms.items():
                if ref_atoms_property_key == 'm_def':
                    continue
                atoms_property = atoms[ref_atoms_property_key]
                if type(atoms_property) == list:
                    property = atoms_property[0]
                    if type(property) == list:
                        assert np.allclose(atoms_property, ref_atoms_property, rtol=1e-05, atol=1e-9)
                    elif type(property) == dict:
                        for property_keys, property_values in property.items():
                            ref_property = ref_atoms_property[0][property_keys]
                            assert property_values == ref_property
                elif type(atoms_property) == dict:
                    for property_keys, property_values in atoms_property.items():
                        ref_property_value = ref_atoms_property[property_keys]
                        if type(property_values) == float:
                            assert np.allclose(property_values, ref_property_value, rtol=1e-05, atol=1e-9)
                        else:
                            assert ref_atoms_property == property_values
                elif type(atoms_property) == float:
                    assert np.allclose(ref_atoms_property, atoms_property, rtol=1e-05, atol=1e-9)
                else:
                    assert ref_atoms_property == atoms_property

        elif subsystem_topology.label == 'subsystem':
            # Indices: passes if the index overlapp is large enough
            ref_indices = ref_topologies[ref_index].indices
            indices = subsystem_topology.indices[0]
            indices_overlap = set(ref_indices).intersection(set(indices))
            assert len(indices_overlap) / \
                len(ref_indices) > 0.85


def test_topology_projection(projection):
    system = projection.run[-1].system[-1]
    assert system.type == 'bulk'
    assert len(system.atoms_group) == 1
    assert system.atoms_group[-1].label == 'Br'
    assert system.atoms_group[-1].type == 'projection'
    assert system.atoms_group[-1].n_atoms == 1
    assert not system.atoms_group[-1].is_molecule
    assert system.atoms_group[-1].atom_indices[0] == 0
    assert not system.atoms_group[-1].atoms_group
    material = projection.results.material
    assert material.structural_type == system.type
    assert material.m_xpath('topology')
    topology = material.topology
    assert_topology(topology)
    assert len(topology) == 2
    assert topology[0].label == 'original'
    assert topology[1].label == system.atoms_group[-1].label
    assert topology[0].structural_type == 'bulk'
    assert topology[1].structural_type == 'group'
    assert len(topology[0].child_systems) == 1
    assert topology[0].child_systems[-1] == topology[1].system_id
    assert topology[0].elements == ['Br', 'K', 'Si']
    assert topology[1].elements == ['Br']
