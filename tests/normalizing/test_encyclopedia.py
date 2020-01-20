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
from hashlib import sha512
import numpy as np
from ase import Atoms
from matid.symmetry.wyckoffset import WyckoffSet

from nomad.parsing import LocalBackend
from nomad.normalizing import structure
from nomad.metainfo.encyclopedia import Encyclopedia
from tests.normalizing.conftest import run_normalize_for_structure, geometry_optimization, molecular_dynamics, phonon, two_d, bulk   # pylint: disable=unused-import


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


def test_2d_flat():
    """Tests that a completely flat 2D material in different supercells is
    detected always as the same material.
    """
    # Expected information for graphene
    wyckoff_sets = [WyckoffSet(
        wyckoff_letter="c",
        element="C",
        indices=[0, 1]
    )]
    space_group_number = 191
    norm_hash_string = structure.get_symmetry_string(space_group_number, wyckoff_sets)
    graphene_material_hash = sha512(norm_hash_string.encode('utf-8')).hexdigest()

    # Graphene orthogonal cell
    atoms = Atoms(
        symbols=["C", "C", "C", "C"],
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
    backend = run_normalize_for_structure(atoms)
    enc = backend.get_mi2_section(Encyclopedia.m_def)
    assert enc.material.material_hash == graphene_material_hash

    # Graphene orthogonal cell, multiplied
    atoms = Atoms(
        symbols=["C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C"],
        positions=[
            [7.1, 6, 6.148780366869514e-1],
            [7.81, 6, 1.8446341100608543],
            [7.1, 6, 3.074390183434757],
            [7.81, 6, 4.30414625680866],
            [4.97, 6, 1.8446341100608543],
            [5.68, 6, 3.074390183434757],
            [4.97, 6, 4.30414625680866],
            [5.68, 6, 6.148780366869514e-1],
            [2.84, 6, 6.148780366869514e-1],
            [3.55, 6, 1.8446341100608543],
            [2.84, 6, 3.074390183434757],
            [3.55, 6, 4.30414625680866],
            [7.1e-1, 6, 1.8446341100608543],
            [1.42, 6, 3.074390183434757],
            [7.1e-1, 6, 4.30414625680866],
            [1.42, 6, 6.148780366869514e-1]
        ],
        cell=[
            [8.52, 0.0, 0.0],
            [0.0, 12, 0.0],
            [0.0, 0.0, 4.919024293495611]
        ],
        pbc=True
    )
    backend = run_normalize_for_structure(atoms)
    enc = backend.get_mi2_section(Encyclopedia.m_def)
    assert enc.material.material_hash == graphene_material_hash

    # Graphene primitive cell
    atoms = Atoms(
        symbols=["C", "C"],
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
    backend = run_normalize_for_structure(atoms)
    enc = backend.get_mi2_section(Encyclopedia.m_def)
    assert enc.material.material_hash == graphene_material_hash


def test_2d_finite_thickness():
    """Tests that a non-flat 2D material is detected properly.
    """
    # Expected information for MoS2
    wyckoff_sets = [
        WyckoffSet(
            wyckoff_letter="e",
            element="S",
            indices=[2, 5]
        ),
        WyckoffSet(
            wyckoff_letter="e",
            element="S",
            indices=[3, 4]
        ),
        WyckoffSet(
            wyckoff_letter="e",
            element="Mo",
            indices=[0, 1]
        )
    ]
    space_group_number = 11
    norm_hash_string = structure.get_symmetry_string(space_group_number, wyckoff_sets)
    graphene_material_hash = sha512(norm_hash_string.encode('utf-8')).hexdigest()

    # MoS2 orthogonal cell
    atoms = Atoms(
        symbols=["Mo", "Mo", "S", "S", "S", "S"],
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
    backend = run_normalize_for_structure(atoms)
    enc = backend.get_mi2_section(Encyclopedia.m_def)
    assert enc.material.material_hash == graphene_material_hash

    # MoS2 orthogonal supercell
    atoms *= [2, 3, 1]
    backend = run_normalize_for_structure(atoms)
    enc = backend.get_mi2_section(Encyclopedia.m_def)
    assert enc.material.material_hash == graphene_material_hash
