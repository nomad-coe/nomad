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
from tests.normalizing.conftest import (  # pylint: disable=unused-import
    atom,
    molecule,
    one_d,
    two_d,
    surface,
    bulk,
)


def assert_material(material):
    assert material.elements
    assert material.nelements
    assert material.chemical_formula_descriptive
    assert material.chemical_formula_reduced
    assert material.chemical_formula_hill
    assert material.chemical_formula_anonymous
    assert material.chemical_formula_reduced_fragments


def assert_symmetry(symmetry):
    assert symmetry.bravais_lattice
    assert symmetry.crystal_system
    assert symmetry.hall_number
    assert symmetry.international_short_symbol
    assert symmetry.point_group
    assert symmetry.space_group_number
    assert symmetry.structure_name
    assert symmetry.strukturbericht_designation
    assert symmetry.prototype
    assert symmetry.aflow_prototype_label


def assert_structure(structure, has_cell=True, has_wyckoff=False):
    assert len(structure.cartesian_site_positions) == structure.nsites
    if has_cell:
        assert len(structure.dimension_types) == 3
        assert np.sum(structure.dimension_types) == structure.nperiodic_dimensions
        assert structure.lattice_vectors.shape == (3, 3)
        assert structure.lattice_parameters.a
        assert structure.lattice_parameters.b
        assert structure.lattice_parameters.c
        assert structure.lattice_parameters.alpha
        assert structure.lattice_parameters.beta
        assert structure.lattice_parameters.gamma
        assert structure.cell_volume
    if has_wyckoff:
        assert len(structure.wyckoff_sets) > 0
        for wset in structure.wyckoff_sets:
            assert len(wset.indices) > 0
            assert wset.wyckoff_letter
            assert wset.element


def test_material_atom(atom):
    material = atom.results.material
    assert_material(material)
    assert material.structural_type == "atom"
    assert material.material_id is None
    assert material.functional_type is None
    assert material.compound_type is None
    assert material.material_name is None
    assert material.symmetry is None
    assert_structure(material.structure_original, has_cell=False)


def test_material_molecule(molecule):
    material = molecule.results.material
    assert_material(material)
    assert material.structural_type == "molecule / cluster"
    assert material.material_id is None
    assert material.functional_type is None
    assert material.compound_type is None
    assert material.material_name is None
    assert material.symmetry is None
    assert_structure(material.structure_original, has_cell=False)


def test_material_1d(one_d):
    material = one_d.results.material
    assert_material(material)
    assert material.structural_type == "1D"
    assert isinstance(material.material_id, str)
    assert material.functional_type is None
    assert material.compound_type is None
    assert material.material_name is None
    assert material.symmetry is None
    assert_structure(material.structure_original)


def test_material_2d(two_d):
    material = two_d.results.material
    assert_material(material)
    assert material.structural_type == "2D"
    assert isinstance(material.material_id, str)
    assert material.functional_type is None
    assert material.compound_type is None
    assert material.material_name is None
    assert material.symmetry is None
    assert_structure(material.structure_original)


def test_material_surface(surface):
    material = surface.results.material
    assert_material(material)
    assert material.structural_type == "surface"
    assert material.material_id is None
    assert material.functional_type is None
    assert material.compound_type is None
    assert material.material_name is None
    assert material.symmetry is None
    assert_structure(material.structure_original)


def test_material_bulk(bulk):
    material = bulk.results.material
    assert_material(material)
    assert material.structural_type == "bulk"
    assert isinstance(material.material_id, str)
    assert material.functional_type
    assert material.compound_type
    assert isinstance(material.material_name, str)
    assert_symmetry(material.symmetry)
    assert_structure(material.structure_original)
    assert_structure(material.structure_primitive)
    assert_structure(material.structure_conventional, has_wyckoff=True)
