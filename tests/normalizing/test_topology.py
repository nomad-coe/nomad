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

from nomad.units import ureg
from tests.normalizing.conftest import (  # pylint: disable=unused-import
    get_template_for_structure,
    get_template_topology,
    conv_bcc,
    conv_fcc,
    rattle,
    stack,
    surf,
    single_cu_surface_topology,
    single_cr_surface_topology,
    stacked_cu_ni_surface_topology,
    graphene,
    graphene_topology,
    boron_nitride,
    boron_nitride_topology,
    mos2,
    mos2_topology,
    stacked_graphene_boron_nitride_topology,
    projection
)


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
    assert mol_group.dimensionality is None
    assert mol_group.building_block == 'group'
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
    assert mol_group.dimensionality is None
    assert mol.building_block == 'molecule'
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
    assert mol_group.dimensionality is None
    assert mon_group.building_block == 'group'
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
    assert mol_group.dimensionality is None
    assert mon.building_block == 'monomer'
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


@pytest.mark.parametrize('surface, ref_topologies', [
    pytest.param(rattle(surf(conv_fcc('Cu'), [1, 0, 0])), single_cu_surface_topology(),
                 id='single surface Cu FCC 100'),
    pytest.param(rattle(surf(conv_fcc('Cu'), [1, 1, 0])), single_cu_surface_topology(),
                 id='single surface Cu FCC 110'),
    pytest.param(rattle(surf(conv_bcc('Cr'), [1, 0, 0])), single_cr_surface_topology(),
                 id='single surface Cr BCC 100'),
    pytest.param(rattle(surf(conv_bcc('Cr'), [1, 1, 0])), single_cr_surface_topology(),
                 id='single surface Cr BCC 110'),
    pytest.param(rattle(
        stack(
            surf(conv_fcc('Cu'), [1, 0, 0], vacuum=0),
            surf(conv_fcc('Ni'), [1, 0, 0], vacuum=0)
        )),
        stacked_cu_ni_surface_topology(),
        id='stacked surfaces of Cu and Ni'),
    pytest.param(rattle(graphene()), graphene_topology(),
                 id='single 2D layer of graphene'),
    pytest.param(rattle(boron_nitride()), boron_nitride_topology(),
                 id='single 2D layer of BN'),
    pytest.param(rattle(mos2()), mos2_topology(),
                 id='single 2D layer of MoS2'),
    pytest.param(rattle(
        stack(
            graphene(),
            boron_nitride()
        )),
        stacked_graphene_boron_nitride_topology(),
        id='stacked layer of BN and C'
    )
])
def test_2D_topology(surface, ref_topologies):

    def compare_section(real, ref):
        """Used to compare two metainfo sections."""
        for name in ref.m_def.all_quantities.keys():
            compare_quantity(name, real, ref)

    def compare_quantity(name, value, ref):
        """Used to compare two metainfo quantities."""
        real_value = getattr(value, name)
        ref_value = getattr(ref, name)
        if ref_value is None:
            assert real_value is None
        elif isinstance(ref_value, ureg.Quantity):
            assert real_value.magnitude == pytest.approx(ref_value.magnitude, rel=0.01, abs=0)
        elif isinstance(ref_value, (np.ndarray, list)):
            real_array = np.array(real_value)
            ref_array = np.array(ref_value)
            if ref_array.dtype == bool:
                assert np.array_equal(real_array, ref_array)
            else:
                assert np.isclose(real_array, ref_array, rtol=0.01, atol=0)
        else:
            raise ValueError("Could not compare values.")

    entry_archive = get_template_for_structure(surface)
    topology = entry_archive.results.material.topology
    assert_topology(topology)
    subsystem_topologies = topology[1:]
    # Compare topology with reference system topology. topology[0] is the original system
    assert len(subsystem_topologies) == len(ref_topologies)
    for res_system in subsystem_topologies:
        formula_hill = res_system.chemical_formula_hill
        for ref_top_counter, ref_topology in enumerate(ref_topologies):
            if ref_topology.chemical_formula_hill == formula_hill:
                ref_formula_hill = ref_topology.chemical_formula_hill
                ref_index = ref_top_counter
                break
        ref_elements = ref_topologies[ref_index].elements
        elements = res_system.elements
        assert elements == ref_elements
        assert formula_hill == ref_formula_hill

        ref_system = ref_topologies[ref_index]
        assert res_system.structural_type == ref_system.structural_type
        assert res_system.dimensionality == ref_system.dimensionality
        assert res_system.building_block == ref_system.building_block

        if res_system.label == 'conventional cell':
            # Cell
            compare_section(res_system.cell, ref_topologies[ref_index].cell)

            # Symmetry
            if ref_topologies[ref_index].symmetry:
                symmetry = res_system.symmetry.m_to_dict()
                ref_symmetry = ref_topologies[ref_index].symmetry.m_to_dict()
                for ref_symmetry_property_key, ref_symmetry_property in ref_symmetry.items():
                    symmetry_property = symmetry[ref_symmetry_property_key]
                    assert ref_symmetry_property == symmetry_property
            else:
                assert res_system.symmetry == ref_topologies[ref_index].symmetry

            # Prototype
            if ref_topologies[ref_index].prototype:
                prototype = res_system.prototype.m_to_dict()
                ref_prototype = ref_topologies[ref_index].prototype.m_to_dict()
                for ref_prototype_property_key, ref_prototype_property in ref_prototype.items():
                    prototype_property = prototype[ref_prototype_property_key]
                    assert ref_prototype_property == prototype_property
            else:
                assert ref_topologies[ref_index].prototype == res_system.prototype

            # Atoms
            atoms = res_system.atoms.m_to_dict()
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

        elif res_system.label == 'subsystem':
            # Indices: passes if the index overlap is large enough
            ref_indices = ref_topologies[ref_index].indices
            indices = res_system.indices[0]
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
    assert material.topology
    topology = material.topology
    assert_topology(topology)
    assert len(topology) == 2
    assert topology[0].label == 'original'
    assert topology[1].label == system.atoms_group[-1].label
    assert topology[1].structural_type == 'group'
    assert topology[1].building_block == 'group'
    assert len(topology[0].child_systems) == 1
    assert topology[0].child_systems[-1] == topology[1].system_id
    assert topology[0].elements == ['Br', 'K', 'Si']
    assert topology[1].elements == ['Br']
