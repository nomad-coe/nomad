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
from ase import Atoms
import ase.build
from matid.symmetry.wyckoffset import WyckoffSet
from pint import UnitRegistry

from nomad.utils import hash
from nomad.parsing import LocalBackend
from nomad.normalizing import structure
from nomad.metainfo.encyclopedia import Encyclopedia
from tests.normalizing.conftest import (  # pylint: disable=unused-import
    run_normalize_for_structure,
    geometry_optimization,
    molecular_dynamics,
    phonon,
    two_d,
    bulk,
    bands_unpolarized_no_gap,
    bands_polarized_no_gap,
    bands_unpolarized_gap_indirect,
    bands_polarized_gap_indirect,
    dos_unpolarized_vasp,
    dos_polarized_vasp,
)

ureg = UnitRegistry()


def test_geometry_optimization(geometry_optimization: LocalBackend):
    """Tests that geometry optimizations are correctly processed."
    """
    enc = geometry_optimization.get_mi2_section(Encyclopedia.m_def)
    run_type = enc.run_type.run_type
    assert run_type == "geometry optimization"


def test_molecular_dynamics(molecular_dynamics: LocalBackend):
    """Tests that geometry optimizations are correctly processed."
    """
    enc = molecular_dynamics.get_mi2_section(Encyclopedia.m_def)
    run_type = enc.run_type.run_type
    assert run_type == "molecular dynamics"


def test_phonon(phonon: LocalBackend):
    """Tests that geometry optimizations are correctly processed."
    """
    enc = phonon.get_mi2_section(Encyclopedia.m_def)
    run_type = enc.run_type.run_type
    assert run_type == "phonon calculation"


def test_1d_metainfo(one_d: LocalBackend):
    """Tests that metainfo for 1D systems is correctly processed.
    """
    enc = one_d.get_mi2_section(Encyclopedia.m_def)
    # Material
    assert enc.material.system_type == "1D"
    assert enc.material.formula == "C6H4"
    assert enc.material.formula_reduced == "C3H2"
    assert np.array_equal(enc.material.periodicity, [True, False, False])

    # Representative system
    assert enc.material.number_of_atoms == 10
    assert enc.material.atom_labels == ["C", "C", "C", "C", "C", "C", "H", "H", "H", "H"]
    assert enc.material.atom_positions is not None
    assert enc.material.cell_normalized is not None
    assert np.allclose(enc.material.lattice_parameters, [4.33793652e-10, 0, 0, 0, 0, 0], atol=0)


def test_2d_metainfo(two_d: LocalBackend):
    """Tests that metainfo for 2D systems is correctly processed.
    """
    enc = two_d.get_mi2_section(Encyclopedia.m_def)
    # Material
    assert enc.material.system_type == "2D"
    assert enc.material.formula == "C2"
    assert enc.material.formula_reduced == "C"
    assert np.array_equal(enc.material.periodicity, [True, True, False])

    # Representative system
    assert enc.material.number_of_atoms == 2
    assert enc.material.atom_labels == ["C", "C"]
    assert enc.material.atom_positions is not None
    assert enc.material.cell_normalized is not None
    assert enc.material.cell_primitive is not None
    assert np.allclose(enc.material.lattice_parameters, [2.46559821e-10, 2.46559821e-10, 0, 120 / 180 * np.pi, 0, 0], atol=0)


def test_bulk_metainfo(bulk: LocalBackend):
    """Tests that metainfo for bulk systems is correctly processed.
    """
    enc = bulk.get_mi2_section(Encyclopedia.m_def)
    # Material
    assert enc.material.system_type == "bulk"
    assert enc.material.formula == "Si2"
    assert enc.material.formula_reduced == "Si"
    assert enc.material.material_name == "Silicon"
    assert enc.material.structure_type == "diamond"
    assert enc.material.structure_prototype == "C"
    assert enc.material.strukturbericht_designation == "A4"

    # Symmetry
    assert enc.material.crystal_system == "cubic"
    assert enc.material.bravais_lattice == "cF"
    assert enc.material.has_free_wyckoff_parameters is False
    assert enc.material.point_group == "m-3m"
    assert enc.material.wyckoff_sets is not None
    assert enc.material.space_group_number == 227
    assert enc.material.space_group_international_short_symbol == "Fd-3m"

    # Representative system
    assert enc.material.number_of_atoms == 8
    assert enc.material.atom_labels == ["Si", "Si", "Si", "Si", "Si", "Si", "Si", "Si"]
    assert enc.material.atom_positions is not None
    assert enc.material.cell_normalized is not None
    assert enc.material.cell_primitive is not None
    assert np.array_equal(enc.material.periodicity, [True, True, True])
    assert enc.material.lattice_parameters is not None
    assert enc.material.cell_volume == pytest.approx(5.431**3 * 1e-30)

    # Calculation
    assert enc.properties.atomic_density == pytest.approx(4.99402346512432e+28)
    assert enc.properties.mass_density == pytest.approx(8 * 28.0855 * 1.6605389e-27 / (5.431**3 * 1e-30))  # Atomic mass in kg/m^3


def test_1d_material_identification():
    # Original nanotube
    nanotube1 = ase.build.nanotube(4, 4, vacuum=4)
    backend = run_normalize_for_structure(nanotube1)
    enc = backend.get_mi2_section(Encyclopedia.m_def)
    hash1 = enc.material.material_hash

    # Rotated copy
    nanotube2 = nanotube1.copy()
    nanotube2.rotate(90, "z", rotate_cell=True)
    backend = run_normalize_for_structure(nanotube2)
    enc = backend.get_mi2_section(Encyclopedia.m_def)
    hash2 = enc.material.material_hash
    assert hash2 == hash1

    # Longer copy
    nanotube3 = nanotube1.copy()
    nanotube3 *= [1, 1, 2]
    backend = run_normalize_for_structure(nanotube3)
    enc = backend.get_mi2_section(Encyclopedia.m_def)
    hash3 = enc.material.material_hash
    assert hash3 == hash1

    # Slightly distorted copies should match
    np.random.seed(4)
    for _ in range(10):
        nanotube4 = nanotube1.copy()
        pos = nanotube4.get_positions()
        pos += 0.2 * np.random.rand(pos.shape[0], pos.shape[1])
        nanotube4.set_positions(pos)
        backend = run_normalize_for_structure(nanotube4)
        enc = backend.get_mi2_section(Encyclopedia.m_def)
        hash4 = enc.material.material_hash
        assert hash4 == hash1

    # Too distorted copy should not match
    nanotube5 = nanotube1.copy()
    pos = nanotube5.get_positions()
    np.random.seed(4)
    pos += 1 * np.random.rand(pos.shape[0], pos.shape[1])
    nanotube5.set_positions(pos)
    backend = run_normalize_for_structure(nanotube5)
    enc = backend.get_mi2_section(Encyclopedia.m_def)
    hash5 = enc.material.material_hash
    assert hash5 != hash1


def test_2d_material_identification():
    # Expected information for graphene. Graphene is an example of completely
    # flat 2D material.
    wyckoff_sets = [WyckoffSet(
        wyckoff_letter="c",
        element="C",
        indices=[0, 1]
    )]
    space_group_number = 191
    norm_hash_string = structure.get_symmetry_string(space_group_number, wyckoff_sets)
    graphene_material_hash = hash(norm_hash_string)

    # Graphene orthogonal cell
    graphene = Atoms(
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
    backend = run_normalize_for_structure(graphene)
    enc = backend.get_mi2_section(Encyclopedia.m_def)
    assert enc.material.material_hash == graphene_material_hash

    # Graphene orthogonal supercell
    graphene2 = graphene.copy()
    graphene2 *= [2, 1, 2]
    backend = run_normalize_for_structure(graphene2)
    enc = backend.get_mi2_section(Encyclopedia.m_def)
    assert enc.material.material_hash == graphene_material_hash

    # Graphene primitive cell
    graphene3 = Atoms(
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
    backend = run_normalize_for_structure(graphene3)
    enc = backend.get_mi2_section(Encyclopedia.m_def)
    assert enc.material.material_hash == graphene_material_hash

    # Slightly distorted system should match
    np.random.seed(4)
    for _ in range(10):
        graphene4 = graphene.copy()
        pos = graphene4.get_positions()
        pos += 0.05 * np.random.rand(pos.shape[0], pos.shape[1])
        graphene4.set_positions(pos)
        backend = run_normalize_for_structure(graphene4)
        enc = backend.get_mi2_section(Encyclopedia.m_def)
        hash4 = enc.material.material_hash
        assert hash4 == graphene_material_hash

    # Too distorted system should not match
    graphene5 = graphene.copy()
    pos = graphene5.get_positions()
    np.random.seed(4)
    pos += 1 * np.random.rand(pos.shape[0], pos.shape[1])
    graphene5.set_positions(pos)
    backend = run_normalize_for_structure(graphene5)
    enc = backend.get_mi2_section(Encyclopedia.m_def)
    hash5 = enc.material.material_hash
    assert hash5 != graphene_material_hash

    # Expected information for MoS2. MoS2 has finite thichkness unlike
    # graphene. The structure is thus treated differently and tested
    # separately.
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
    mos2_material_hash = hash(norm_hash_string)

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
    assert enc.material.material_hash == mos2_material_hash

    # MoS2 orthogonal supercell
    atoms *= [2, 3, 1]
    backend = run_normalize_for_structure(atoms)
    enc = backend.get_mi2_section(Encyclopedia.m_def)
    assert enc.material.material_hash == mos2_material_hash


def test_bulk_material_identification():
    # Original system
    wurtzite = ase.build.bulk("SiC", crystalstructure="wurtzite", a=3.086, c=10.053)
    backend = run_normalize_for_structure(wurtzite)
    enc = backend.get_mi2_section(Encyclopedia.m_def)
    hash1 = enc.material.material_hash

    # Rotated
    wurtzite2 = wurtzite.copy()
    wurtzite2.rotate(90, "z", rotate_cell=True)
    backend = run_normalize_for_structure(wurtzite2)
    enc = backend.get_mi2_section(Encyclopedia.m_def)
    hash2 = enc.material.material_hash
    assert hash2 == hash1

    # Supercell
    wurtzite3 = wurtzite.copy()
    wurtzite3 *= [2, 3, 1]
    backend = run_normalize_for_structure(wurtzite3)
    enc = backend.get_mi2_section(Encyclopedia.m_def)
    hash3 = enc.material.material_hash
    assert hash3 == hash1

    # Slightly distorted system should match
    np.random.seed(4)
    for _ in range(10):
        wurtzite4 = wurtzite.copy()
        pos = wurtzite4.get_positions()
        pos += 0.05 * np.random.rand(pos.shape[0], pos.shape[1])
        wurtzite4.set_positions(pos)
        backend = run_normalize_for_structure(wurtzite4)
        enc = backend.get_mi2_section(Encyclopedia.m_def)
        hash4 = enc.material.material_hash
        assert hash4 == hash1

    # Too distorted system should not match
    wurtzite5 = wurtzite.copy()
    pos = wurtzite5.get_positions()
    np.random.seed(4)
    pos += 1 * np.random.rand(pos.shape[0], pos.shape[1])
    wurtzite5.set_positions(pos)
    backend = run_normalize_for_structure(wurtzite5)
    enc = backend.get_mi2_section(Encyclopedia.m_def)
    hash5 = enc.material.material_hash
    assert hash5 != hash1


def test_1d_structure_structure_at_cell_boundary():
    """Tests that the visualization that is made for 1D systems has the
    correct form even if the cell boundary is at the middle of the
    structure.
    """
    atoms = Atoms(
        symbols=["H", "C"],
        positions=[
            [0.0, 0.0, 0],
            [1.0, 0.0, 10.0]
        ],
        cell=[
            [2, 0.0, 0.0],
            [0.0, 10, 0.0],
            [0.0, 0.0, 10]
        ],
        pbc=True
    )
    backend = run_normalize_for_structure(atoms)
    enc = backend.get_mi2_section(Encyclopedia.m_def)

    expected_cell = [
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 2e-10]
    ]
    expected_pos = [
        [0, 0, 0],
        [0, 0, 0.5],
    ]
    expected_labels = [
        "H",
        "C",
    ]

    assert np.allclose(enc.material.atom_positions, expected_pos)
    assert np.array_equal(enc.material.atom_labels, expected_labels)
    assert np.allclose(enc.material.cell_normalized, expected_cell)


def test_2d_structure_structure_at_cell_boundary():
    """Tests that the visualization that is made for 2D systems has the
    correct form even if the cell boundary is at the middle of the
    structure.
    """
    atoms = Atoms(
        symbols=["H", "C"],
        positions=[
            [0.0, 0.0, 0],
            [0.0, 0.0, 13.800000000000002]
        ],
        cell=[
            [2, 0.0, 0.0],
            [0.0, 2, 0.0],
            [0.0, 0.0, 15]
        ],
        pbc=True
    )
    backend = run_normalize_for_structure(atoms)
    enc = backend.get_mi2_section(Encyclopedia.m_def)

    expected_cell = [
        [2e-10, 0, 0],
        [0, 2e-10, 0],
        [0, 0, 1.2e-10]
    ]
    expected_pos = [
        [0, 0, 1],
        [0, 0, 0],
    ]
    expected_labels = [
        "H",
        "C",
    ]

    assert np.allclose(enc.material.atom_positions, expected_pos)
    assert np.array_equal(enc.material.atom_labels, expected_labels)
    assert np.allclose(enc.material.cell_normalized, expected_cell)


def test_method_dft_metainfo(single_point):
    enc = single_point.get_mi2_section(Encyclopedia.m_def)
    assert enc.method.basis_set_type == "Numeric AOs"
    assert enc.method.core_electron_treatment == "full all electron"
    assert enc.method.code_name == "FHI-aims"
    assert enc.method.code_version == "010314"
    assert enc.method.functional_long_name == "GGA_C_PBE+GGA_X_PBE"
    assert enc.method.functional_type == "GGA"


def test_method_gw_metainfo(gw):
    enc = gw.get_mi2_section(Encyclopedia.m_def)
    assert enc.method.code_name == "FHI-aims"
    assert enc.method.code_version == "180607"
    assert enc.method.gw_type == "G0W0"
    assert enc.method.gw_starting_point == "GGA_C_PBE+0.75*GGA_X_PBE+0.25*HF_X"


def test_band_structure(bands_unpolarized_no_gap, bands_polarized_no_gap, bands_unpolarized_gap_indirect, bands_polarized_gap_indirect):

    def test_generic(bs, n_channels):
        """Generic tests for band structure data."""
        for segment in bs.segments:
            assert segment.energies.shape[0] == n_channels
            assert len(segment.energies.shape) == 3
            assert segment.energies.shape[2] == segment.k_points.shape[0]
        assert bs.brillouin_zone is not None
        assert bs.reciprocal_cell.shape == (3, 3)

    # Unpolarized, no gaps
    enc = bands_unpolarized_no_gap.get_mi2_section(Encyclopedia.m_def)
    properties = enc.properties
    bs = properties.electronic_band_structure
    test_generic(bs, n_channels=1)
    assert bs.band_gap is None
    assert bs.band_gap_spin_up is None
    assert bs.band_gap_spin_down is None

    # Polarized, no gaps
    enc = bands_polarized_no_gap.get_mi2_section(Encyclopedia.m_def)
    properties = enc.properties
    bs = properties.electronic_band_structure
    test_generic(bs, n_channels=2)
    assert bs.band_gap is None
    assert bs.band_gap_spin_up is None
    assert bs.band_gap_spin_down is None

    # Unpolarized, finite gap, indirect
    enc = bands_unpolarized_gap_indirect.get_mi2_section(Encyclopedia.m_def)
    properties = enc.properties
    bs = properties.electronic_band_structure
    test_generic(bs, n_channels=1)
    gap_ev = (bs.band_gap.value * ureg.J).to(ureg.eV).magnitude
    assert gap_ev == pytest.approx(0.62, 0.01)
    assert bs.band_gap.type == "indirect"
    assert bs.band_gap_spin_up is None
    assert bs.band_gap_spin_down is None

    # Polarized, finite gap, indirect
    enc = bands_polarized_gap_indirect.get_mi2_section(Encyclopedia.m_def)
    properties = enc.properties
    bs = properties.electronic_band_structure
    test_generic(bs, n_channels=2)
    gap = bs.band_gap
    gap_up = bs.band_gap_spin_up
    gap_down = bs.band_gap_spin_down
    gap_ev = (gap.value * ureg.J).to(ureg.eV).magnitude
    gap_up_ev = (gap_up.value * ureg.J).to(ureg.eV).magnitude
    gap_down_ev = (gap_down.value * ureg.J).to(ureg.eV).magnitude
    assert gap_up.type == "indirect"
    assert gap_down.type == "indirect"
    assert gap_up_ev != gap_down_ev
    assert gap_up_ev == gap_ev
    assert gap_up_ev == pytest.approx(0.956, 0.01)
    assert gap_down_ev == pytest.approx(1.230, 0.01)


def test_dos(dos_unpolarized_vasp, dos_polarized_vasp):

    def test_generic(dos, n_channels):
        """Generic tests for DOS data."""
        assert len(dos.energies.shape) == 1
        assert len(dos.values.shape) == 2
        assert dos.values.shape[1] == dos.energies.shape[0]
        assert dos.values.shape[0] == n_channels

    # Unpolarized
    enc = dos_unpolarized_vasp.get_mi2_section(Encyclopedia.m_def)
    properties = enc.properties
    dos = properties.electronic_dos
    test_generic(dos, n_channels=1)

    # Polarized
    enc = dos_polarized_vasp.get_mi2_section(Encyclopedia.m_def)
    properties = enc.properties
    dos = properties.electronic_dos
    test_generic(dos, n_channels=2)
