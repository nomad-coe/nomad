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
from matid.symmetry.wyckoffset import WyckoffSet  # pylint: disable=import-error

from nomad.units import ureg
from nomad import atomutils
from nomad.utils import hash
from tests.normalizing.conftest import get_template_for_structure


def assert_material(material):
    assert material.elements
    assert material.n_elements
    assert material.chemical_formula_hill
    assert material.chemical_formula_iupac
    assert material.chemical_formula_reduced
    assert material.chemical_formula_anonymous
    assert material.chemical_formula_descriptive
    assert material.chemical_formula_reduced_fragments


def assert_symmetry(symmetry):
    assert symmetry.bravais_lattice
    assert symmetry.crystal_system
    assert symmetry.hall_number
    assert symmetry.point_group
    assert symmetry.space_group_number
    assert symmetry.space_group_symbol
    assert symmetry.structure_name
    assert symmetry.strukturbericht_designation
    assert symmetry.prototype_formula
    assert symmetry.prototype_aflow_id


def assert_structure(structure, has_cell=True, has_wyckoff=False):
    assert len(structure.cartesian_site_positions) == structure.n_sites
    assert len(structure.species_at_sites) == structure.n_sites
    assert len(structure.species) > 0
    assert structure.species[0].name
    assert structure.species[0].concentration
    assert structure.species[0].chemical_symbols
    if has_cell:
        assert len(structure.dimension_types) == 3
        assert np.sum(structure.dimension_types) == structure.nperiodic_dimensions
        assert structure.lattice_vectors.shape == (3, 3)
        a = structure.lattice_parameters.a
        b = structure.lattice_parameters.b
        c = structure.lattice_parameters.c
        assert a is not None
        assert b is not None
        assert c is not None
        if b != 0 and c != 0:
            assert structure.lattice_parameters.alpha is not None
        else:
            assert structure.lattice_parameters.alpha is None
        if a != 0 and c != 0:
            assert structure.lattice_parameters.beta is not None
        else:
            assert structure.lattice_parameters.beta is None
        if a != 0 and b != 0:
            assert structure.lattice_parameters.gamma is not None
        else:
            assert structure.lattice_parameters.gamma is None
        assert structure.cell_volume is not None
        if structure.cell_volume != 0 and structure.nperiodic_dimensions == 3:
            assert structure.mass_density
            assert structure.atomic_density
    if has_wyckoff:
        assert len(structure.wyckoff_sets) > 0
        for wset in structure.wyckoff_sets:
            assert len(wset.indices) > 0
            assert wset.wyckoff_letter
            assert wset.element


def test_material_atom(atom):
    material = atom.results.material
    assert_material(material)
    assert material.material_id is None
    assert material.structural_type == 'atom'
    assert material.functional_type is None
    assert material.compound_type is None
    assert material.material_name is None
    assert material.symmetry is None

    properties = atom.results.properties
    assert_structure(properties.structures.structure_original)


def test_material_molecule(molecule):
    material = molecule.results.material
    assert_material(material)
    assert material.material_id is None
    assert material.structural_type == 'molecule / cluster'
    assert material.functional_type is None
    assert material.compound_type is None
    assert material.material_name is None
    assert material.symmetry is None

    properties = molecule.results.properties
    assert_structure(properties.structures.structure_original, has_cell=False)


def test_material_1d(one_d):
    # Material
    material = one_d.results.material
    assert_material(material)
    assert isinstance(material.material_id, str)
    assert material.structural_type == '1D'
    assert material.functional_type is None
    assert material.compound_type is None
    assert material.material_name is None
    assert material.chemical_formula_hill == 'C2H2'
    assert material.chemical_formula_iupac == 'CH'
    assert material.chemical_formula_descriptive == 'C2H2'
    assert material.chemical_formula_reduced == 'CH'
    assert material.chemical_formula_anonymous == 'AB'
    assert material.elements == ['C', 'H']
    assert material.n_elements == 2
    assert material.symmetry is None

    # Conventional structure
    conv = one_d.results.properties.structures.structure_conventional
    assert_structure(conv)
    assert conv.n_sites == 4
    assert conv.species_at_sites == ['C', 'C', 'H', 'H']
    assert np.array_equal(conv.dimension_types, [1, 0, 0])
    assert conv.lattice_parameters.a.to(ureg.angstrom).magnitude == pytest.approx(2.459, abs=1e-3)
    assert conv.lattice_parameters.b.to(ureg.angstrom).magnitude == 0
    assert conv.lattice_parameters.c.to(ureg.angstrom).magnitude == pytest.approx(2.890, abs=1e-3)
    assert conv.lattice_parameters.alpha is None
    assert conv.lattice_parameters.beta.magnitude == pytest.approx(np.pi / 2)
    assert conv.lattice_parameters.gamma is None

    # Original structure
    assert_structure(one_d.results.properties.structures.structure_original)

    # Primitive structure
    assert_structure(one_d.results.properties.structures.structure_primitive)


def test_material_2d(two_d):
    # Material
    material = two_d.results.material
    assert_material(material)
    assert isinstance(material.material_id, str)
    assert material.structural_type == '2D'
    assert material.functional_type is None
    assert material.compound_type is None
    assert material.material_name is None
    assert material.chemical_formula_hill == 'C2'
    assert material.chemical_formula_iupac == 'C'
    assert material.chemical_formula_descriptive == 'C2'
    assert material.chemical_formula_reduced == 'C'
    assert material.chemical_formula_anonymous == 'A'
    assert material.elements == ['C']
    assert material.n_elements == 1
    assert material.symmetry is None

    # Conventional structure
    conv = two_d.results.properties.structures.structure_conventional
    assert_structure(conv)
    assert conv.n_sites == 2
    assert conv.species_at_sites == ['C', 'C']
    assert np.array_equal(conv.dimension_types, [1, 1, 0])
    assert conv.lattice_parameters.a.to(ureg.angstrom).magnitude == pytest.approx(2.461, abs=1e-3)
    assert conv.lattice_parameters.b.to(ureg.angstrom).magnitude == pytest.approx(2.461, abs=1e-3)
    assert conv.lattice_parameters.c.to(ureg.angstrom).magnitude == 0
    assert conv.lattice_parameters.alpha is None
    assert conv.lattice_parameters.beta is None
    assert conv.lattice_parameters.gamma.magnitude == pytest.approx(120 / 180 * np.pi)

    # Original structure
    assert_structure(two_d.results.properties.structures.structure_original)

    # Primitive structure
    assert_structure(two_d.results.properties.structures.structure_primitive)


def test_material_surface(surface):
    material = surface.results.material
    assert_material(material)
    assert material.material_id is None
    assert material.structural_type == 'surface'
    assert material.functional_type is None
    assert material.compound_type is None
    assert material.material_name is None
    assert material.symmetry is None

    properties = surface.results.properties
    assert_structure(properties.structures.structure_original)


def test_material_bulk(bulk):
    # Material
    material = bulk.results.material
    assert_material(material)
    assert isinstance(material.material_id, str)
    assert material.structural_type == 'bulk'
    assert material.functional_type
    assert material.compound_type
    assert material.material_name == 'Silicon'
    assert material.chemical_formula_hill == 'Si8'
    assert material.chemical_formula_iupac == 'Si'
    assert material.chemical_formula_descriptive == 'Si8'
    assert material.chemical_formula_reduced == 'Si'
    assert material.chemical_formula_anonymous == 'A'
    assert material.elements == ['Si']
    assert material.n_elements == 1
    assert_symmetry(material.symmetry)

    # Conventional structure
    conv = bulk.results.properties.structures.structure_conventional
    assert_structure(conv, has_wyckoff=True)
    assert conv.n_sites == 8
    assert conv.species_at_sites == ['Si', 'Si', 'Si', 'Si', 'Si', 'Si', 'Si', 'Si']
    assert np.array_equal(conv.dimension_types, [1, 1, 1])
    assert conv.lattice_parameters.a.to(ureg.angstrom).magnitude == pytest.approx(5.431, abs=1e-3)
    assert conv.lattice_parameters.b.to(ureg.angstrom).magnitude == pytest.approx(5.431, abs=1e-3)
    assert conv.lattice_parameters.c.to(ureg.angstrom).magnitude == pytest.approx(5.431, abs=1e-3)
    assert conv.lattice_parameters.alpha.magnitude == pytest.approx(np.pi / 2)
    assert conv.lattice_parameters.beta.magnitude == pytest.approx(np.pi / 2)
    assert conv.lattice_parameters.gamma.magnitude == pytest.approx(np.pi / 2)

    # Original structure
    assert_structure(bulk.results.properties.structures.structure_original)

    # Primitive structure
    assert_structure(bulk.results.properties.structures.structure_primitive)


def test_material_eels(eels):
    material = eels.results.material
    assert material.n_elements == 2
    assert material.elements == ['Si', 'O']
    assert material.chemical_formula_hill == 'OSi'
    assert material.chemical_formula_iupac == 'SiO'
    assert material.chemical_formula_reduced == 'OSi'
    assert material.chemical_formula_descriptive == 'OSi'


def test_1d_material_identification():
    # Original nanotube
    nanotube1 = ase.build.nanotube(4, 4, vacuum=4)
    hash1 = get_template_for_structure(nanotube1).results.material.material_id

    # Rotated copy
    nanotube2 = nanotube1.copy()
    nanotube2.rotate(90, 'z', rotate_cell=True)
    hash2 = get_template_for_structure(nanotube2).results.material.material_id
    assert hash2 == hash1

    # Longer copy
    nanotube3 = nanotube1.copy()
    nanotube3 *= [1, 1, 2]
    hash3 = get_template_for_structure(nanotube3).results.material.material_id
    assert hash3 == hash1

    # Slightly distorted copies should match
    np.random.seed(4)
    for _ in range(10):
        nanotube4 = nanotube1.copy()
        pos = nanotube4.get_positions()
        pos += 0.2 * np.random.rand(pos.shape[0], pos.shape[1])
        nanotube4.set_positions(pos)
        hash4 = get_template_for_structure(nanotube4).results.material.material_id
        assert hash4 == hash1

    # Too distorted copy should not match
    nanotube5 = nanotube1.copy()
    pos = nanotube5.get_positions()
    np.random.seed(4)
    pos += 1 * np.random.rand(pos.shape[0], pos.shape[1])
    nanotube5.set_positions(pos)
    hash5 = get_template_for_structure(nanotube5).results.material.material_id
    assert hash5 != hash1


def test_2d_material_identification():
    # Expected information for graphene. Graphene is an example of completely
    # flat 2D material.
    wyckoff_sets = [WyckoffSet(
        wyckoff_letter='c',
        element='C',
        indices=[0, 1]
    )]
    space_group_number = 191
    norm_hash_string = atomutils.get_symmetry_string(space_group_number, wyckoff_sets, is_2d=True)
    graphene_material_id = hash(norm_hash_string)

    # Graphene orthogonal cell
    graphene = Atoms(
        symbols=['C', 'C', 'C', 'C'],
        positions=[
            [2.84, 7.5, 6.148780366869514e-1],
            [3.55, 7.5, 1.8446341100608543],
            [7.1e-1, 7.5, 1.8446341100608543],
            [1.42, 7.5, 6.148780366869514e-1]
        ],
        cell=[
            [4.26, 0.0, 0.0],
            [0.0, 15, 0.0],
            [0.0, 0.0, 2.4595121467478055]
        ],
        pbc=True
    )
    material_id = get_template_for_structure(graphene).results.material.material_id
    assert material_id == graphene_material_id

    # Graphene orthogonal supercell
    graphene2 = graphene.copy()
    graphene2 *= [2, 1, 2]
    material_id = get_template_for_structure(graphene2).results.material.material_id
    assert material_id == graphene_material_id

    # Graphene primitive cell
    graphene3 = Atoms(
        symbols=['C', 'C'],
        positions=[
            [0, 1.42, 6],
            [1.2297560733739028, 7.100000000000001e-1, 6]
        ],
        cell=[
            [2.4595121467478055, 0.0, 0.0],
            [-1.2297560733739028, 2.13, 0.0],
            [0.0, 0.0, 12]
        ],
        pbc=True
    )
    material_id = get_template_for_structure(graphene3).results.material.material_id
    assert material_id == graphene_material_id

    # Slightly distorted system should match
    np.random.seed(4)
    for _ in range(10):
        graphene4 = graphene.copy()
        pos = graphene4.get_positions()
        pos += 0.05 * np.random.rand(pos.shape[0], pos.shape[1])
        graphene4.set_positions(pos)
        material_id = get_template_for_structure(graphene4).results.material.material_id
        assert material_id == graphene_material_id

    # Too distorted system should not match
    graphene5 = graphene.copy()
    pos = graphene5.get_positions()
    np.random.seed(4)
    pos += 1 * np.random.rand(pos.shape[0], pos.shape[1])
    graphene5.set_positions(pos)
    material_id = get_template_for_structure(graphene5).results.material.material_id
    assert material_id != graphene_material_id

    # Expected information for MoS2. MoS2 has finite thichkness unlike
    # graphene. The structure is thus treated differently and tested
    # separately.
    wyckoff_sets = [
        WyckoffSet(
            wyckoff_letter='e',
            element='S',
            indices=[2, 5]
        ),
        WyckoffSet(
            wyckoff_letter='e',
            element='S',
            indices=[3, 4]
        ),
        WyckoffSet(
            wyckoff_letter='e',
            element='Mo',
            indices=[0, 1]
        )
    ]
    space_group_number = 11
    norm_hash_string = atomutils.get_symmetry_string(space_group_number, wyckoff_sets, is_2d=True)
    mos2_material_id = hash(norm_hash_string)

    # MoS2 orthogonal cell
    atoms = Atoms(
        symbols=['Mo', 'Mo', 'S', 'S', 'S', 'S'],
        scaled_positions=[
            [0.000000, 0.022916, 0.630019],
            [0.500000, 0.418064, 0.635988],
            [0.500000, 0.795155, 0.68827],
            [0.000000, 0.299555, 0.70504],
            [0.500000, 0.141429, 0.56096],
            [0.000000, 0.645894, 0.57774],
        ],
        cell=[
            [3.193638, 0.000000, 0.000000],
            [0.000000, 5.738503, 0.110928],
            [0.000000, 0.021363, 24.194079],
        ],
        pbc=True
    )
    material_id = get_template_for_structure(atoms).results.material.material_id
    assert material_id == mos2_material_id

    # MoS2 orthogonal supercell
    atoms *= [2, 3, 1]
    material_id = get_template_for_structure(atoms).results.material.material_id
    assert material_id == mos2_material_id


def test_bulk_material_identification():
    # Original system
    wurtzite = ase.build.bulk('SiC', crystalstructure='wurtzite', a=3.086, c=10.053)
    material_id_wurtzite = get_template_for_structure(wurtzite).results.material.material_id

    # Rotated
    wurtzite2 = wurtzite.copy()
    wurtzite2.rotate(90, 'z', rotate_cell=True)
    material_id = get_template_for_structure(wurtzite2).results.material.material_id
    assert material_id == material_id_wurtzite

    # Supercell
    wurtzite3 = wurtzite.copy()
    wurtzite3 *= [2, 3, 1]
    materia_id = get_template_for_structure(wurtzite3).results.material.material_id
    assert materia_id == material_id_wurtzite

    # Slightly distorted system should match
    np.random.seed(4)
    for _ in range(10):
        wurtzite4 = wurtzite.copy()
        pos = wurtzite4.get_positions()
        pos += 0.05 * np.random.rand(pos.shape[0], pos.shape[1])
        wurtzite4.set_positions(pos)
        material_id = get_template_for_structure(wurtzite4).results.material.material_id
        assert material_id == material_id_wurtzite

    # Too distorted system should not match
    wurtzite5 = wurtzite.copy()
    pos = wurtzite5.get_positions()
    np.random.seed(4)
    pos += 1 * np.random.rand(pos.shape[0], pos.shape[1])
    wurtzite5.set_positions(pos)
    material_id = get_template_for_structure(wurtzite5).results.material.material_id
    assert material_id != material_id_wurtzite


one_d_split = Atoms(
    symbols=['H', 'C'],
    positions=[
        [0.0, 0.0, 0],
        [1.0, 0.0, 10.0]
    ],
    cell=[
        [0.0, 10, 0.0],
        [2, 0.0, 0.0],
        [0.0, 0.0, 10]
    ],
    pbc=True
)
one_d_split_expected = Atoms(
    symbols=['H', 'C'],
    positions=[
        [0, 0, 0],
        [0, 0, 1],
    ],
    cell=[
        [0, 0, 2],
        [0, 0, 0],
        [0, 0, 0],
    ],
    pbc=[True, False, False]
)
two_d_split = Atoms(
    symbols=['H', 'C'],
    positions=[
        [0.0, 0.0, 0],
        [0.0, 0.0, 13.800000000000002]
    ],
    cell=[
        [2, 0.0, 0.0],
        [0.0, 0.0, 15],
        [0.0, 2, 0.0],
    ],
    pbc=True
)
two_d_split_expected = Atoms(
    symbols=['H', 'C'],
    positions=[
        [0, 0, 1.2],
        [0, 0, 0],
    ],
    cell=[
        [2, 0, 0],
        [0, 2, 0],
        [0, 0, 1.2]
    ],
    pbc=[True, True, False]
)
two_d_swap = Atoms(
    symbols=['B', 'N'],
    positions=[
        [0, 0, 0],
        [-0.6, 0.3, 0],
    ],
    cell=[
        [1, 2, 0],
        [-2, 1, 0],
        [0, 0, 20]
    ],
    pbc=True
)
two_d_swap_expected = Atoms(
    symbols=['B', 'N'],
    positions=[
        [0, 0, 1.51589629],
        [0, 0, 0.84507589],
    ],
    cell=[
        [2.23607, 0, 0],
        [0, 0, 2.23607],
        [0, 0, 0]
    ],
    pbc=[True, True, False]
)


@pytest.mark.parametrize(
    'atoms, expected',
    [
        # 1D with cell boundary in the middle of the structure
        (one_d_split, one_d_split_expected),
        # 2D with cell boundary in the middle of the structure
        (two_d_split, two_d_split_expected),
        # 2D with cell where the nonperiodic axis is not last by default in the
        # conventional cell.
        (two_d_swap, two_d_swap_expected),
    ]
)
def test_conventional_structure(atoms, expected):
    '''Tests that the conventional structure has the correct form.
    '''
    entry = get_template_for_structure(atoms)
    structure_conventional = entry.results.properties.structures.structure_conventional
    pos = structure_conventional.cartesian_site_positions.to(ureg.angstrom).magnitude
    cell = structure_conventional.lattice_vectors.to(ureg.angstrom).magnitude
    pbc = np.array(structure_conventional.dimension_types, dtype=bool)

    assert np.array_equal(pbc, expected.get_pbc())
    assert np.allclose(pos, expected.get_positions())
    assert np.array_equal(structure_conventional.species_at_sites, expected.get_chemical_symbols())
    assert np.allclose(cell, expected.get_cell())
