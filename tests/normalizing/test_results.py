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

from nomad.units import ureg

from .conftest import (
    get_template_dos,
    get_template_band_structure,
    run_normalize
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
    assert symmetry.point_group
    assert symmetry.space_group_number
    assert symmetry.space_group_symbol
    assert symmetry.structure_name
    assert symmetry.strukturbericht_designation
    assert symmetry.prototype_formula
    assert symmetry.prototype_aflow_id


def assert_structure(structure, has_cell=True, has_wyckoff=False):
    assert len(structure.cartesian_site_positions) == structure.nsites
    assert len(structure.species_at_sites) == structure.nsites
    assert len(structure.species) > 0
    assert structure.species[0].name
    assert structure.species[0].concentration
    assert structure.species[0].chemical_elements
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
    assert material.material_id is None
    assert material.type_structural == "atom"
    assert material.type_functional is None
    assert material.type_compound is None
    assert material.material_name is None
    assert material.symmetry is None

    properties = atom.results.properties
    assert_structure(properties.structure_original)


def test_material_molecule(molecule):
    material = molecule.results.material
    assert_material(material)
    assert material.material_id is None
    assert material.type_structural == "molecule / cluster"
    assert material.type_functional is None
    assert material.type_compound is None
    assert material.material_name is None
    assert material.symmetry is None

    properties = molecule.results.properties
    assert_structure(properties.structure_original, has_cell=False)


def test_material_1d(one_d):
    material = one_d.results.material
    assert_material(material)
    assert isinstance(material.material_id, str)
    assert material.type_structural == "1D"
    assert material.type_functional is None
    assert material.type_compound is None
    assert material.material_name is None
    assert material.symmetry is None

    properties = one_d.results.properties
    assert_structure(properties.structure_original)


def test_material_2d(two_d):
    material = two_d.results.material
    assert_material(material)
    assert isinstance(material.material_id, str)
    assert material.type_structural == "2D"
    assert material.type_functional is None
    assert material.type_compound is None
    assert material.material_name is None
    assert material.symmetry is None

    properties = two_d.results.properties
    assert_structure(properties.structure_original)


def test_material_surface(surface):
    material = surface.results.material
    assert_material(material)
    assert material.material_id is None
    assert material.type_structural == "surface"
    assert material.type_functional is None
    assert material.type_compound is None
    assert material.material_name is None
    assert material.symmetry is None

    properties = surface.results.properties
    assert_structure(properties.structure_original)


def test_material_bulk(bulk):
    material = bulk.results.material
    assert_material(material)
    assert material.type_structural == "bulk"
    assert isinstance(material.material_id, str)
    assert material.type_functional
    assert material.type_compound
    assert isinstance(material.material_name, str)
    assert_symmetry(material.symmetry)

    properties = bulk.results.properties
    assert_structure(properties.structure_original)
    assert_structure(properties.structure_primitive)
    assert_structure(properties.structure_conventional, has_wyckoff=True)


def test_method_dft(dft):
    method = dft.results.method
    assert method.method_name == "DFT"
    assert method.simulation.dft.basis_set_type == "plane waves"
    assert method.simulation.dft.core_electron_treatment == "pseudopotential"
    assert method.simulation.dft.xc_functional_names == ["GGA_C_PBE", "GGA_X_PBE"]
    assert method.simulation.dft.xc_functional_type == "GGA"
    assert method.simulation.dft.smearing_kind == "gaussian"
    assert method.simulation.dft.smearing_width == 1e-20
    assert method.simulation.dft.spin_polarized is True
    assert method.simulation.dft.scf_threshold_energy_change == 1e-24 * ureg.joule
    assert method.simulation.dft.van_der_Waals_method == "G06"
    assert method.simulation.dft.relativity_method == "scalar_relativistic"


def test_method_dft_plus_u(dft_plus_u):
    method = dft_plus_u.results.method
    assert method.method_name == "DFT"
    assert method.simulation.dft.basis_set_type == "plane waves"
    assert method.simulation.dft.core_electron_treatment == "pseudopotential"
    assert method.simulation.dft.xc_functional_names == ["GGA_C_PBE", "GGA_X_PBE"]
    assert method.simulation.dft.xc_functional_type == "GGA"
    assert method.simulation.dft.smearing_kind == "gaussian"
    assert method.simulation.dft.smearing_width == 1e-20
    assert method.simulation.dft.spin_polarized is True
    assert method.simulation.dft.scf_threshold_energy_change == 1e-24 * ureg.joule
    assert method.simulation.dft.van_der_Waals_method == "G06"
    assert method.simulation.dft.relativity_method == "scalar_relativistic"


def test_method_gw(gw):
    method = gw.results.method
    assert method.method_name == "GW"
    assert method.simulation.gw.gw_type == "G0W0"
    assert method.simulation.gw.starting_point == ["GGA_C_PBE", "GGA_X_PBE"]


def test_dos_electronic():
    # DOS without energy references
    archive = get_template_dos(has_references=False)
    dos = archive.results.properties.dos_electronic
    assert dos.spin_polarized is False
    assert dos.densities.shape == (1, 200)
    assert dos.energies.shape == (200, )

    # Unpolarized DOS
    archive = get_template_dos(spin_polarized=False)
    dos = archive.results.properties.dos_electronic
    assert dos.spin_polarized is False
    assert dos.densities.shape == (1, 200)
    assert dos.energies.shape == (200, )

    # Polarized DOS
    archive = get_template_dos(spin_polarized=True)
    dos = archive.results.properties.dos_electronic
    assert dos.spin_polarized is True
    assert dos.densities.shape == (2, 200)
    assert dos.energies.shape == (200, )

    # Vibrational instead of electronic
    archive = get_template_dos(type="vibrational")
    dos = archive.results.properties.dos_electronic
    assert dos is None

    # Empty values
    archive = get_template_dos(normalize=False)
    archive.section_run[0].section_single_configuration_calculation[0].section_dos[0].dos_values = []
    archive = run_normalize(archive)
    dos = archive.results.properties.dos_electronic
    assert dos is None

    # Empty energies
    archive = get_template_dos(normalize=False)
    archive.section_run[0].section_single_configuration_calculation[0].section_dos[0].dos_energies = []
    archive = run_normalize(archive)
    dos = archive.results.properties.dos_electronic
    assert dos is None


def test_band_structure_electronic():
    # Band structure without energy reference
    archive = get_template_band_structure([(1, "direct")], has_references=False)
    bs = archive.results.properties.band_structure_electronic
    assert bs.reciprocal_cell.shape == (3, 3)
    assert bs.band_gap is None
    assert bs.band_gap_type is None
    assert bs.spin_polarized is False
    assert bs.energy_highest_occupied is None
    assert bs.energy_lowest_unoccupied is None
    assert bs.segments[0].band_energies.shape == (1, 100, 2)
    assert bs.segments[0].band_k_points.shape == (100, 3)

    # Unpolarized band structure with no gap
    archive = get_template_band_structure([None])
    bs = archive.results.properties.band_structure_electronic
    assert bs.reciprocal_cell.shape == (3, 3)
    assert bs.band_gap.magnitude == [0]
    assert bs.band_gap_type == ["no_gap"]
    assert bs.spin_polarized is False
    assert bs.energy_highest_occupied.shape == (1,)
    assert bs.energy_lowest_unoccupied.shape == (1,)
    assert bs.segments[0].band_energies.shape == (1, 100, 2)
    assert bs.segments[0].band_k_points.shape == (100, 3)

    # Polarized band structure with no gap
    archive = get_template_band_structure([None, None])
    bs = archive.results.properties.band_structure_electronic
    assert bs.reciprocal_cell.shape == (3, 3)
    assert np.allclose(bs.band_gap.magnitude, [0, 0])
    assert bs.band_gap_type == ["no_gap", "no_gap"]
    assert bs.spin_polarized is True
    assert bs.energy_highest_occupied.shape == (2,)
    assert bs.energy_lowest_unoccupied.shape == (2,)
    assert bs.segments[0].band_energies.shape == (2, 100, 2)
    assert bs.segments[0].band_k_points.shape == (100, 3)

    # Unpolarized band structure with direct gap
    gap = 1  # eV
    gap_type = "direct"
    archive = get_template_band_structure([(gap, gap_type)])
    bs = archive.results.properties.band_structure_electronic
    assert bs.reciprocal_cell.shape == (3, 3)
    assert bs.band_gap.magnitude == [pytest.approx((gap * ureg.electron_volt).to(ureg.joule).magnitude)]
    assert bs.band_gap_type == [gap_type]
    assert bs.spin_polarized is False
    assert bs.energy_fermi.shape == (1,)
    assert bs.energy_highest_occupied.shape == (1,)
    assert bs.energy_lowest_unoccupied.shape == (1,)
    assert bs.segments[0].band_energies.shape == (1, 100, 2)
    assert bs.segments[0].band_k_points.shape == (100, 3)

    # Unpolarized band structure with indirect gap
    gap = 1   # eV
    gap_type = "indirect"
    archive = get_template_band_structure([(gap, gap_type)])
    bs = archive.results.properties.band_structure_electronic
    assert bs.reciprocal_cell.shape == (3, 3)
    assert bs.band_gap.magnitude == [pytest.approx((gap * ureg.electron_volt).to(ureg.joule).magnitude)]
    assert bs.band_gap_type == [gap_type]
    assert bs.spin_polarized is False
    assert bs.energy_fermi.shape == (1,)
    assert bs.energy_highest_occupied.shape == (1,)
    assert bs.energy_lowest_unoccupied.shape == (1,)
    assert bs.segments[0].band_energies.shape == (1, 100, 2)
    assert bs.segments[0].band_k_points.shape == (100, 3)

    # Polarized band structure with direct gap
    gap1 = 1  # eV
    gap2 = 2  # eV
    gap_type = "direct"
    archive = get_template_band_structure([(gap1, gap_type), (gap2, gap_type)])
    bs = archive.results.properties.band_structure_electronic
    gaps_joule = (np.array([gap1, gap2]) * ureg.electron_volt).to(ureg.joule).magnitude
    assert bs.reciprocal_cell.shape == (3, 3)
    assert np.allclose(bs.band_gap.magnitude, gaps_joule)
    assert bs.band_gap_type == [gap_type, gap_type]
    assert bs.spin_polarized is True
    assert bs.energy_fermi.shape == (2,)
    assert bs.energy_highest_occupied.shape == (2,)
    assert bs.energy_lowest_unoccupied.shape == (2,)
    assert bs.segments[0].band_energies.shape == (2, 100, 2)
    assert bs.segments[0].band_k_points.shape == (100, 3)

    # Polarized band structure with indirect gap
    gap1 = 1  # eV
    gap2 = 2  # eV
    gap1_type = "indirect"
    gap2_type = "indirect"
    archive = get_template_band_structure([(gap1, gap1_type), (gap2, gap2_type)])
    bs = archive.results.properties.band_structure_electronic
    gaps_joule = (np.array([gap1, gap2]) * ureg.electron_volt).to(ureg.joule).magnitude
    assert bs.reciprocal_cell.shape == (3, 3)
    assert np.allclose(bs.band_gap.magnitude, gaps_joule)
    assert bs.band_gap_type == [gap1_type, gap2_type]
    assert bs.spin_polarized is True
    assert bs.energy_fermi.shape == (2,)
    assert bs.energy_highest_occupied.shape == (2,)
    assert bs.energy_lowest_unoccupied.shape == (2,)
    assert bs.segments[0].band_energies.shape == (2, 100, 2)
    assert bs.segments[0].band_k_points.shape == (100, 3)


def test_dos_phonon():
    # DOS with all correct metainfo
    archive = get_template_dos(type="vibrational")
    dos = archive.results.properties.dos_phonon
    assert dos.densities.shape == (1, 200)
    assert dos.energies.shape == (200, )

    # Electronic instead of vibrational
    archive = get_template_dos(type="electronic")
    dos = archive.results.properties.dos_phonon
    assert dos is None

    # Empty values
    archive = get_template_dos(type="vibrational", normalize=False)
    archive.section_run[0].section_single_configuration_calculation[0].section_dos[0].dos_values = []
    archive = run_normalize(archive)
    dos = archive.results.properties.dos_phonon
    assert dos is None

    # Empty energies
    archive = get_template_dos(type="vibrational", normalize=False)
    archive.section_run[0].section_single_configuration_calculation[0].section_dos[0].dos_energies = []
    archive = run_normalize(archive)
    dos = archive.results.properties.dos_phonon
    assert dos is None


def test_band_structure_phonon():
    # Valid phonon band structure
    archive = get_template_band_structure(type="vibrational")
    bs = archive.results.properties.band_structure_phonon
    assert bs.segments[0].band_energies.shape == (1, 100, 2)
    assert bs.segments[0].band_k_points.shape == (100, 3)


def pprint(root, indent=None):
    ''' Pretty prints the containment hierarchy '''
    contents = list(root.m_contents())
    quantities = root.m_def.quantities
    n_sub = len(contents)
    n_q = len(quantities)

    if indent is None:
        indent = []
    if len(indent) == 0:
        indent_str = ''
    else:
        indent_str = ''.join(['   ' if last else '   ' for last in indent[:-1]])
        indent_str += '   ' if indent[-1] else '   '
    msg = indent_str + str(root.m_def.name)
    if hasattr(root.m_def, "repeats") and root.m_def.repeats:
        msg += " (repeats)"
    print(msg)

    for i, q in enumerate(quantities):
        istr = ''.join(['   ' if last else '   ' for last in indent])
        istr += '   ' if i == n_q - 1 and n_sub == 0 else '   '
        try:
            dtype = q.type.__name__
        except Exception:
            if isinstance(q.type, np.dtype):
                dtype = str(q.type)
            else:
                dtype = type(q.type).__name__
        shape = q.shape
        if shape:
            tp = "{}{}".format(dtype, shape)
        else:
            tp = dtype
        print(istr + str(q.name) + ": " + tp)

    for i, content in enumerate(contents):
        pprint(content, indent + [i == n_sub - 1])
