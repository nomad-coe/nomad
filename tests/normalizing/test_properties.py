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

from nomad.units import ureg

from .test_material import assert_structure
from .conftest import (
    get_template_dos,
    get_template_band_structure,
    run_normalize
)


def test_eels(eels):
    assert eels.results.properties.spectroscopy.eels is not None
    spectroscopy_data = eels.results.properties.spectroscopy.spectrum
    eels_data = eels.results.properties.spectroscopy.eels
    assert eels_data.resolution.to(ureg.electron_volt).magnitude == pytest.approx(1)
    assert eels_data.min_energy.to(ureg.electron_volt).magnitude == pytest.approx(100)
    assert eels_data.max_energy.to(ureg.electron_volt).magnitude == pytest.approx(200)
    assert eels_data.detector_type == "Quantum GIF"
    assert spectroscopy_data.n_values == spectroscopy_data.count.shape[0] == spectroscopy_data.energy.shape[0]


def test_bulk_modulus(mechanical):
    bulk_modulus = mechanical.results.properties.mechanical.bulk_modulus
    assert len(bulk_modulus) == 1
    modulus = bulk_modulus[0]
    assert modulus.type == "murnaghan"
    assert modulus.value.magnitude == pytest.approx(10000)


def test_shear_modulus(mechanical):
    shear_modulus = mechanical.results.properties.mechanical.shear_modulus
    assert len(shear_modulus) == 3
    modulus_voigt_reuss_hill = shear_modulus[0]
    assert modulus_voigt_reuss_hill.type == "voigt_reuss_hill_average"
    assert modulus_voigt_reuss_hill.value.magnitude == 10000
    modulus_voigt = shear_modulus[1]
    assert modulus_voigt.type == "voigt_average"
    assert modulus_voigt.value.magnitude == 10000
    modulus_reuss_hill = shear_modulus[2]
    assert modulus_reuss_hill.type == "reuss_average"
    assert modulus_reuss_hill.value.magnitude == 10000


def test_energy_volume_curve(mechanical):
    ev = mechanical.results.properties.mechanical.energy_volume_curve
    assert len(ev) == 2
    ev_raw = ev[0]
    assert ev_raw.type == "raw"
    assert ev_raw.energies_raw.shape == (10,)
    assert ev_raw.volumes.shape == (10,)
    ev_murhagnan = ev[1]
    assert ev_murhagnan.type == "murnaghan"
    assert ev_murhagnan.energies_fit.shape == (10,)
    assert ev_murhagnan.volumes.shape == (10,)


def test_dos_electronic():
    gap_fill = [[0, 1], [2, 3]]

    # DOS without energy references
    archive = get_template_dos()
    dos = archive.results.properties.electronic.dos_electronic
    assert dos.spin_polarized is False
    assert dos.total[0].value.shape == (101,)
    assert dos.energies.shape == (101,)

    # Unpolarized DOS with gap:
    efermi = 1.5
    archive = get_template_dos(energy_reference_fermi=efermi)
    dos = archive.results.properties.electronic.dos_electronic
    assert len(dos.band_gap) == 1
    assert dos.band_gap[0].energy_highest_occupied is not None
    assert dos.band_gap[0].energy_lowest_unoccupied is not None
    assert dos.spin_polarized is False
    assert dos.total[0].value.shape == (101,)
    assert dos.energies.shape == (101,)

    # Polarized DOS
    efermi = 1.5
    archive = get_template_dos(fill=[gap_fill, gap_fill], energy_reference_fermi=efermi)
    dos = archive.results.properties.electronic.dos_electronic
    assert len(dos.band_gap) == 2
    assert dos.band_gap[0].energy_highest_occupied is not None
    assert dos.band_gap[0].energy_lowest_unoccupied is not None
    assert dos.spin_polarized is True
    assert dos.total[0].value.shape == (101,)
    assert dos.total[1].value.shape == (101,)
    assert dos.energies.shape == (101,)

    # Vibrational instead of electronic
    archive = get_template_dos(type="vibrational")
    electronic = archive.results.properties.electronic
    assert electronic is None

    # Empty values
    archive = get_template_dos(normalize=False)
    archive.run[0].calculation[0].dos_electronic[0].total = None
    archive = run_normalize(archive)
    electronic = archive.results.properties.electronic
    assert electronic is None

    # Empty energies
    archive = get_template_dos(normalize=False)
    archive.run[0].calculation[0].dos_electronic[0].energies = None
    archive = run_normalize(archive)
    electronic = archive.results.properties.electronic
    assert electronic is None


def test_band_structure_electronic():
    # Band structure without energy reference
    archive = get_template_band_structure([(1, "direct")], has_references=False)
    bs = archive.results.properties.electronic.band_structure_electronic
    band_gaps = bs.band_gap
    assert bs.reciprocal_cell.shape == (3, 3)
    assert bs.spin_polarized is False
    assert len(band_gaps) == 0
    assert bs.segment[0].energies.shape == (1, 100, 2)
    assert bs.segment[0].kpoints.shape == (100, 3)

    # Unpolarized band structure with no gap
    archive = get_template_band_structure([None])
    bs = archive.results.properties.electronic.band_structure_electronic
    band_gaps = bs.band_gap
    assert bs.reciprocal_cell.shape == (3, 3)
    assert bs.spin_polarized is False
    assert bs.energy_fermi is not None
    assert len(band_gaps) == 1
    assert band_gaps[0].value == 0
    assert band_gaps[0].type is None
    assert band_gaps[0].energy_highest_occupied is not None
    assert band_gaps[0].energy_lowest_unoccupied is not None
    assert bs.segment[0].energies.shape == (1, 100, 2)
    assert bs.segment[0].kpoints.shape == (100, 3)

    # Polarized band structure with no gap
    archive = get_template_band_structure([None, None])
    bs = archive.results.properties.electronic.band_structure_electronic
    band_gaps = bs.band_gap
    assert bs.reciprocal_cell.shape == (3, 3)
    assert bs.spin_polarized is True
    assert bs.energy_fermi is not None
    assert len(band_gaps) == 2
    assert band_gaps[0].energy_highest_occupied is not None
    assert band_gaps[0].energy_lowest_unoccupied is not None
    assert band_gaps[0].value == 0
    assert band_gaps[0].type is None
    assert band_gaps[1].energy_highest_occupied is not None
    assert band_gaps[1].energy_lowest_unoccupied is not None
    assert band_gaps[1].value == 0
    assert band_gaps[1].type is None
    assert bs.segment[0].energies.shape == (2, 100, 2)
    assert bs.segment[0].kpoints.shape == (100, 3)

    # Unpolarized band structure with direct gap
    gap = 1  # eV
    gap_type = "direct"
    archive = get_template_band_structure([(gap, gap_type)])
    bs = archive.results.properties.electronic.band_structure_electronic
    band_gaps = bs.band_gap
    assert bs.reciprocal_cell.shape == (3, 3)
    assert bs.spin_polarized is False
    assert bs.energy_fermi is not None
    assert len(band_gaps) == 1
    assert band_gaps[0].energy_highest_occupied is not None
    assert band_gaps[0].energy_lowest_unoccupied is not None
    assert band_gaps[0].value == pytest.approx((gap * ureg.electron_volt).to(ureg.joule).magnitude)
    assert band_gaps[0].type == gap_type
    assert bs.segment[0].energies.shape == (1, 100, 2)
    assert bs.segment[0].kpoints.shape == (100, 3)

    # Unpolarized band structure with indirect gap
    gap = 1   # eV
    gap_type = "indirect"
    archive = get_template_band_structure([(gap, gap_type)])
    bs = archive.results.properties.electronic.band_structure_electronic
    band_gaps = bs.band_gap
    assert bs.reciprocal_cell.shape == (3, 3)
    assert bs.spin_polarized is False
    assert bs.energy_fermi is not None
    assert len(band_gaps) == 1
    assert band_gaps[0].energy_highest_occupied is not None
    assert band_gaps[0].energy_lowest_unoccupied is not None
    assert band_gaps[0].value == pytest.approx((gap * ureg.electron_volt).to(ureg.joule).magnitude)
    assert band_gaps[0].type == gap_type
    assert bs.segment[0].energies.shape == (1, 100, 2)
    assert bs.segment[0].kpoints.shape == (100, 3)

    # Polarized band structure with direct gap
    gap1 = 1  # eV
    gap2 = 2  # eV
    gap_type = "direct"
    archive = get_template_band_structure([(gap1, gap_type), (gap2, gap_type)])
    bs = archive.results.properties.electronic.band_structure_electronic
    band_gaps = bs.band_gap
    assert bs.reciprocal_cell.shape == (3, 3)
    assert bs.spin_polarized is True
    assert bs.energy_fermi is not None
    assert len(band_gaps) == 2
    assert band_gaps[0].energy_highest_occupied is not None
    assert band_gaps[0].energy_lowest_unoccupied is not None
    assert band_gaps[0].value == pytest.approx((gap1 * ureg.electron_volt).to(ureg.joule).magnitude)
    assert band_gaps[0].type == gap_type
    assert band_gaps[1].energy_highest_occupied is not None
    assert band_gaps[1].energy_lowest_unoccupied is not None
    assert band_gaps[1].value == pytest.approx((gap2 * ureg.electron_volt).to(ureg.joule).magnitude)
    assert band_gaps[1].type == gap_type
    assert bs.segment[0].energies.shape == (2, 100, 2)
    assert bs.segment[0].kpoints.shape == (100, 3)

    # Polarized band structure with indirect gap
    gap1 = 1  # eV
    gap2 = 2  # eV
    gap_type = "indirect"
    archive = get_template_band_structure([(gap1, gap_type), (gap2, gap_type)])
    bs = archive.results.properties.electronic.band_structure_electronic
    band_gaps = bs.band_gap
    assert bs.reciprocal_cell.shape == (3, 3)
    assert bs.spin_polarized is True
    assert bs.energy_fermi is not None
    assert len(band_gaps) == 2
    assert band_gaps[0].energy_highest_occupied is not None
    assert band_gaps[0].energy_lowest_unoccupied is not None
    assert band_gaps[0].value == pytest.approx((gap1 * ureg.electron_volt).to(ureg.joule).magnitude)
    assert band_gaps[0].type == gap_type
    assert band_gaps[1].energy_highest_occupied is not None
    assert band_gaps[1].energy_lowest_unoccupied is not None
    assert band_gaps[1].value == pytest.approx((gap1 * ureg.electron_volt).to(ureg.joule).magnitude)
    assert band_gaps[1].type == gap_type
    assert bs.segment[0].energies.shape == (2, 100, 2)
    assert bs.segment[0].kpoints.shape == (100, 3)


def test_dos_phonon():
    # DOS with all correct metainfo
    archive = get_template_dos(type="vibrational")
    dos = archive.results.properties.vibrational.dos_phonon
    assert dos.total[0].value.shape == (101,)
    assert dos.energies.shape == (101,)

    # Electronic instead of vibrational
    archive = get_template_dos(type="electronic")
    vibrational = archive.results.properties.vibrational
    assert vibrational is None

    # Empty values
    archive = get_template_dos(type="vibrational", normalize=False)
    archive.run[0].calculation[0].dos_phonon[0].total = None
    archive = run_normalize(archive)
    vibrational = archive.results.properties.vibrational
    assert vibrational is None

    # Empty energies
    archive = get_template_dos(type="vibrational", normalize=False)
    archive.run[0].calculation[0].dos_phonon[0].energies = []
    archive = run_normalize(archive)
    vibrational = archive.results.properties.vibrational
    assert vibrational is None


def test_band_structure_phonon():
    # Valid phonon band structure
    archive = get_template_band_structure(type="vibrational")
    bs = archive.results.properties.vibrational.band_structure_phonon
    assert bs.segment[0].energies.shape == (1, 100, 2)
    assert bs.segment[0].kpoints.shape == (100, 3)


def test_energy_free_helmholtz(phonon):
    energy_free = phonon.results.properties.vibrational.energy_free_helmholtz
    assert energy_free.temperatures.shape == (11, )
    assert energy_free.energies.shape == (11, )


def test_heat_capacity_constant_volume(phonon):
    heat_cap = phonon.results.properties.vibrational.heat_capacity_constant_volume
    assert heat_cap.temperatures.shape == (11, )
    assert heat_cap.heat_capacities.shape == (11, )


def test_geometry_optimization(geometry_optimization):
    geo_opt_prop = geometry_optimization.results.properties.geometry_optimization
    assert_structure(geo_opt_prop.structure_optimized)
    n_frames = len(geo_opt_prop.trajectory)
    n_energies = len(geo_opt_prop.energies)
    assert n_frames > 0
    assert n_frames == n_energies
    assert geo_opt_prop.final_energy_difference > 0
    assert geo_opt_prop.type == "ionic"


def test_n_calculations(geometry_optimization):
    assert geometry_optimization.results.properties.n_calculations == 2
