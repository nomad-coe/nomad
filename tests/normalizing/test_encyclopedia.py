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
import numpy as np
from ase import Atoms
import ase.build
from matid.symmetry.wyckoffset import WyckoffSet

from nomad.utils import hash
from nomad.files import UploadFiles
from nomad import atomutils
from nomad.datamodel import EntryArchive
from nomad.datamodel.encyclopedia import (
    Calculation,
    EncyclopediaMetadata,
)
from tests.processing.test_data import run_processing
from tests.normalizing.conftest import get_template_for_structure


def test_geometry_optimization(geometry_optimization: EntryArchive):
    """Tests that geometry optimizations are correctly processed."
    """
    enc = geometry_optimization.metadata.encyclopedia
    calc_type = enc.calculation.calculation_type
    assert calc_type == "geometry optimization"


def test_molecular_dynamics(molecular_dynamics: EntryArchive):
    """Tests that molecular dynamics are correctly processed."
    """
    enc = molecular_dynamics.metadata.encyclopedia
    calc_type = enc.calculation.calculation_type
    assert calc_type == "molecular dynamics"


def test_1d_metainfo(one_d: EntryArchive):
    """Tests that metainfo for 1D systems is correctly processed.
    """
    enc = one_d.metadata.encyclopedia
    # Material
    material = enc.material
    assert material.material_type == "1D"
    assert material.formula == "C2H2"
    assert material.formula_reduced == "CH"

    # Idealized structure
    ideal = enc.material.idealized_structure
    assert ideal.number_of_atoms == 4
    assert ideal.atom_labels == ["C", "C", "H", "H"]
    assert ideal.atom_positions is not None
    assert ideal.lattice_vectors is not None
    assert np.array_equal(ideal.periodicity, [True, False, False])
    assert ideal.lattice_parameters.a == pytest.approx(2.45951214e-10)


def test_2d_metainfo(two_d: EntryArchive):
    """Tests that metainfo for 2D systems is correctly processed.
    """
    enc = two_d.metadata.encyclopedia
    # Material
    material = enc.material
    assert material.material_type == "2D"
    assert material.formula == "C2"
    assert material.formula_reduced == "C"

    # Idealized structure
    ideal = enc.material.idealized_structure
    assert ideal.number_of_atoms == 2
    assert ideal.atom_labels == ["C", "C"]
    assert ideal.atom_positions is not None
    assert ideal.lattice_vectors is not None
    assert ideal.lattice_vectors_primitive is not None
    assert np.array_equal(ideal.periodicity, [True, True, False])
    assert ideal.lattice_parameters.a == pytest.approx(2.461e-10)
    assert ideal.lattice_parameters.b == pytest.approx(2.461e-10)
    assert ideal.lattice_parameters.c is None
    assert ideal.lattice_parameters.alpha is None
    assert ideal.lattice_parameters.beta is None
    assert ideal.lattice_parameters.gamma == pytest.approx(120 / 180 * np.pi)


def test_bulk_metainfo(bulk: EntryArchive):
    """Tests that metainfo for bulk systems is correctly processed.
    """
    enc = bulk.metadata.encyclopedia
    # Material
    material = enc.material
    assert material.material_type == "bulk"
    assert material.formula == "Si2"
    assert material.formula_reduced == "Si"
    assert material.material_name == "Silicon"

    # Bulk
    bulk = enc.material.bulk
    assert bulk.crystal_system == "cubic"
    assert bulk.bravais_lattice == "cF"
    assert bulk.has_free_wyckoff_parameters is False
    assert bulk.point_group == "m-3m"
    assert bulk.space_group_number == 227
    assert bulk.structure_type == "diamond"
    assert bulk.structure_prototype == "C"
    assert bulk.strukturbericht_designation == "A4"
    assert bulk.space_group_international_short_symbol == "Fd-3m"

    # Idealized structure
    ideal = enc.material.idealized_structure
    assert ideal.wyckoff_sets is not None
    assert ideal.number_of_atoms == 8
    assert ideal.atom_labels == ["Si", "Si", "Si", "Si", "Si", "Si", "Si", "Si"]
    assert ideal.atom_positions is not None
    assert ideal.lattice_vectors is not None
    assert ideal.lattice_vectors_primitive is not None
    assert np.array_equal(ideal.periodicity, [True, True, True])
    assert ideal.lattice_parameters is not None
    assert ideal.cell_volume.magnitude == 5.431**3 * 1e-30

    # Properties
    prop = enc.properties
    assert prop.atomic_density.magnitude == pytest.approx(4.99402346512432e+28)
    assert prop.mass_density.magnitude == pytest.approx(8 * 28.0855 * 1.6605389e-27 / (5.431**3 * 1e-30))  # Atomic mass in kg/m^3


def test_1d_material_identification():
    # Original nanotube
    nanotube1 = ase.build.nanotube(4, 4, vacuum=4)
    enc = get_template_for_structure(nanotube1).metadata.encyclopedia
    hash1 = enc.material.material_id

    # Rotated copy
    nanotube2 = nanotube1.copy()
    nanotube2.rotate(90, "z", rotate_cell=True)
    enc = get_template_for_structure(nanotube2).metadata.encyclopedia
    hash2 = enc.material.material_id
    assert hash2 == hash1

    # Longer copy
    nanotube3 = nanotube1.copy()
    nanotube3 *= [1, 1, 2]
    enc = get_template_for_structure(nanotube3).metadata.encyclopedia
    hash3 = enc.material.material_id
    assert hash3 == hash1

    # Slightly distorted copies should match
    np.random.seed(4)
    for _ in range(10):
        nanotube4 = nanotube1.copy()
        pos = nanotube4.get_positions()
        pos += 0.2 * np.random.rand(pos.shape[0], pos.shape[1])
        nanotube4.set_positions(pos)
        enc = get_template_for_structure(nanotube4).metadata.encyclopedia
        hash4 = enc.material.material_id
        assert hash4 == hash1

    # Too distorted copy should not match
    nanotube5 = nanotube1.copy()
    pos = nanotube5.get_positions()
    np.random.seed(4)
    pos += 1 * np.random.rand(pos.shape[0], pos.shape[1])
    nanotube5.set_positions(pos)
    enc = get_template_for_structure(nanotube5).metadata.encyclopedia
    hash5 = enc.material.material_id
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
    norm_hash_string = atomutils.get_symmetry_string(space_group_number, wyckoff_sets, is_2d=True)
    graphene_material_id = hash(norm_hash_string)

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
    enc = get_template_for_structure(graphene).metadata.encyclopedia
    assert enc.material.material_id == graphene_material_id

    # Graphene orthogonal supercell
    graphene2 = graphene.copy()
    graphene2 *= [2, 1, 2]
    enc = get_template_for_structure(graphene2).metadata.encyclopedia
    assert enc.material.material_id == graphene_material_id

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
    enc = get_template_for_structure(graphene3).metadata.encyclopedia
    assert enc.material.material_id == graphene_material_id

    # Slightly distorted system should match
    np.random.seed(4)
    for _ in range(10):
        graphene4 = graphene.copy()
        pos = graphene4.get_positions()
        pos += 0.05 * np.random.rand(pos.shape[0], pos.shape[1])
        graphene4.set_positions(pos)
        entry_archive = get_template_for_structure(graphene4)
        enc = entry_archive.metadata.encyclopedia
        hash4 = enc.material.material_id
        assert hash4 == graphene_material_id

    # Too distorted system should not match
    graphene5 = graphene.copy()
    pos = graphene5.get_positions()
    np.random.seed(4)
    pos += 1 * np.random.rand(pos.shape[0], pos.shape[1])
    graphene5.set_positions(pos)
    enc = get_template_for_structure(graphene5).metadata.encyclopedia
    hash5 = enc.material.material_id
    assert hash5 != graphene_material_id

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
    norm_hash_string = atomutils.get_symmetry_string(space_group_number, wyckoff_sets, is_2d=True)
    mos2_material_id = hash(norm_hash_string)

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
    entry_archive = get_template_for_structure(atoms)
    enc = entry_archive.metadata.encyclopedia
    assert enc.material.material_id == mos2_material_id

    # MoS2 orthogonal supercell
    atoms *= [2, 3, 1]
    enc = get_template_for_structure(atoms).metadata.encyclopedia
    assert enc.material.material_id == mos2_material_id


def test_bulk_material_identification():
    # Original system
    wurtzite = ase.build.bulk("SiC", crystalstructure="wurtzite", a=3.086, c=10.053)
    enc = get_template_for_structure(wurtzite).metadata.encyclopedia
    hash1 = enc.material.material_id

    # Rotated
    wurtzite2 = wurtzite.copy()
    wurtzite2.rotate(90, "z", rotate_cell=True)
    enc = get_template_for_structure(wurtzite2).metadata.encyclopedia
    hash2 = enc.material.material_id
    assert hash2 == hash1

    # Supercell
    wurtzite3 = wurtzite.copy()
    wurtzite3 *= [2, 3, 1]
    enc = get_template_for_structure(wurtzite3).metadata.encyclopedia
    hash3 = enc.material.material_id
    assert hash3 == hash1

    # Slightly distorted system should match
    np.random.seed(4)
    for _ in range(10):
        wurtzite4 = wurtzite.copy()
        pos = wurtzite4.get_positions()
        pos += 0.05 * np.random.rand(pos.shape[0], pos.shape[1])
        wurtzite4.set_positions(pos)
        enc = get_template_for_structure(wurtzite4).metadata.encyclopedia
        hash4 = enc.material.material_id
        assert hash4 == hash1

    # Too distorted system should not match
    wurtzite5 = wurtzite.copy()
    pos = wurtzite5.get_positions()
    np.random.seed(4)
    pos += 1 * np.random.rand(pos.shape[0], pos.shape[1])
    wurtzite5.set_positions(pos)
    enc = get_template_for_structure(wurtzite5).metadata.encyclopedia
    hash5 = enc.material.material_id
    assert hash5 != hash1


def test_1d_idealized_structure():
    """Tests that the idealized structure for 1D systems has the correct form.
    """
    # Cell boundary in the middle of the structure.
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
    enc = get_template_for_structure(atoms).metadata.encyclopedia

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

    ideal = enc.material.idealized_structure
    assert np.allclose(ideal.atom_positions, expected_pos)
    assert np.array_equal(ideal.atom_labels, expected_labels)
    assert np.allclose(ideal.lattice_vectors, expected_cell)


def test_2d_idealized_structure():
    """Tests that the visualization that is made for 2D systems has the
    correct form.
    """
    # Cell boundary in the middle of the structure.
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
    enc = get_template_for_structure(atoms).metadata.encyclopedia

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

    ideal = enc.material.idealized_structure
    assert np.allclose(ideal.atom_positions, expected_pos)
    assert np.array_equal(ideal.atom_labels, expected_labels)
    assert np.allclose(ideal.lattice_vectors, expected_cell)


def test_method_dft_metainfo(dft):
    enc = dft.metadata.encyclopedia
    assert enc.method.core_electron_treatment == "pseudopotential"
    assert enc.method.functional_long_name == "1.0*GGA_C_PBE+1.0*GGA_X_PBE"
    assert enc.method.functional_type == "GGA"


def test_method_gw_metainfo(gw):
    enc = gw.metadata.encyclopedia
    assert enc.method.gw_type == "G0W0"
    assert enc.method.gw_starting_point == "1.0*GGA_C_PBE+1.0*GGA_X_PBE"


def test_hashes_exciting(hash_exciting):
    """Tests that the hashes has been successfully created for calculations
    from exciting.
    """
    enc = hash_exciting.metadata.encyclopedia
    method_id = enc.method.method_id
    group_eos_id = enc.method.group_eos_id
    group_parametervariation_id = enc.method.group_parametervariation_id
    assert method_id is not None
    assert group_eos_id is not None
    assert group_parametervariation_id is not None


def test_hashes_undefined(hash_vasp):
    """Tests that the hashes are not present when the method settings cannot be
    determined at a sufficient accuracy.
    """
    enc = hash_vasp.metadata.encyclopedia
    method_id = enc.method.method_id
    group_eos_id = enc.method.group_eos_id

    # If the method cannot be determined accurately, the method hash and group
    # hash cannot be set. Parametervariation has may still be valid, as it does
    # not really need the method to be accurately defined.
    assert method_id is None
    assert group_eos_id is None


def test_dos(dos_unpolarized_vasp, dos_polarized_vasp):
    """Tests that referenced DOS information is valid.
    """
    def generaltests(dos, n_channels):
        assert dos is not None
        assert len(dos.total) == n_channels
        assert dos.total[n_channels - 1].value.shape == (301,)
        assert dos.energies.shape == (301,)

    generaltests(dos_unpolarized_vasp.metadata.encyclopedia.properties.electronic_dos, n_channels=1)
    generaltests(dos_polarized_vasp.metadata.encyclopedia.properties.electronic_dos, n_channels=2)


def test_electronic_bands(bands_unpolarized_no_gap, bands_polarized_no_gap, band_path_cF_nonstandard):
    """Tests that referenced electronic band structure information is valid.
    """
    def generaltests(band):
        assert band is not None
        for segment in band.segment:
            assert segment.energies is not None
            assert segment.kpoints is not None
            assert segment.endpoints_labels is not None

    # VASP bands
    generaltests(bands_unpolarized_no_gap.metadata.encyclopedia.properties.electronic_band_structure)
    generaltests(bands_polarized_no_gap.metadata.encyclopedia.properties.electronic_band_structure)

    # Band structure from exciting calculation where there are multiple sccs
    # and multiple bands present for some reason...
    generaltests(band_path_cF_nonstandard.metadata.encyclopedia.properties.electronic_band_structure)


def test_phonon(test_user, proc_infra):
    """Tests that phonon calculations are correctly processed.
    """
    # Process a phonon calculation
    upload = run_processing(("phonon_id", "tests/data/normalizers/phonons.zip"), test_user)

    # Read the resulting archive
    upload_id = upload.upload_id
    calcs = upload.calcs()
    phonon_id = None
    for calc in calcs:
        if calc.parser == "parsers/phonopy":
            phonon_id = calc.calc_id
            break
    upload_file = UploadFiles.get(upload_id, is_authorized=lambda: True)
    archive_reader = upload_file.read_archive(phonon_id)
    phonon_archive = archive_reader[phonon_id].to_dict()
    phonon = EntryArchive.m_from_dict(phonon_archive)

    enc = phonon.metadata.encyclopedia
    calc_type = enc.calculation.calculation_type
    status = enc.status
    prop = enc.properties
    method = enc.method
    band = prop.phonon_band_structure
    dos = prop.phonon_dos
    thermo_props = prop.thermodynamical_properties
    assert calc_type == Calculation.calculation_type.type.phonon_calculation
    assert status == EncyclopediaMetadata.status.type.success

    # There should be a reference to the external calculation
    assert phonon.run[0].calculation[0].calculations_path[0] is not None

    # The method information should have been read from the referenced
    # calculation
    assert method.method_type == "DFT"
    assert method.core_electron_treatment == "full all electron"
    assert method.functional_type == "LDA"
    assert method.functional_long_name == "LDA_C_PW+LDA_X"

    # Check dos
    assert dos is not None
    assert dos.energies is not None
    assert len(dos.total) > 0

    # Check band structure
    assert band is not None
    for segment in band.segment:
        assert segment.energies is not None
        assert segment.kpoints is not None
        assert segment.endpoints_labels is not None

    # Check thermodynamical properties
    assert thermo_props is not None
    assert thermo_props.heat_capacity_c_v is not None
    assert thermo_props.specific_heat_capacity is not None
    assert thermo_props.temperature is not None
    assert thermo_props.vibrational_free_energy_at_constant_volume is not None


def test_elastic(elastic: EntryArchive):
    """Tests that elastic constants calculations are correctly processed. For
    now, the method information is not being processed, as it requires an
    additional processing step similar to phonon calculations.
    """
    enc = elastic.metadata.encyclopedia
    calc_type = enc.calculation.calculation_type
    status = enc.status
    assert calc_type == Calculation.calculation_type.type.elastic_constants
    assert status == EncyclopediaMetadata.status.type.unsupported_method_type
