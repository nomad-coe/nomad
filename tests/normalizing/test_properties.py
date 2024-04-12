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
    get_template_dft,
    add_template_dos,
    get_template_dos,
    add_template_band_structure,
    get_template_band_structure,
    add_template_magnetic_shielding,
    add_template_spin_spin_coupling,
    add_template_magnetic_susceptibility,
    run_normalize,
)
from nomad.datamodel.metainfo import simulationworkflowschema, SCHEMA_IMPORT_ERROR


def test_eels(eels):
    assert eels.results.method.method_name == 'EELS'
    assert eels.results.properties.spectroscopic is not None
    spectra = eels.results.properties.spectroscopic.spectra
    assert len(spectra) == 1
    assert spectra[0].type == 'EELS'
    assert spectra[0].label == 'experiment'
    assert spectra[0].n_energies == 101
    assert spectra[0].energies[22].to('eV').magnitude == pytest.approx(122.0)
    assert spectra[0].intensities[22] == pytest.approx(22.0)
    assert spectra[0].intensities_units == 'counts'
    assert spectra[0].m_xpath('provenance.eels')
    provenance = spectra[0].provenance
    assert provenance.label == 'EELSDB'
    assert provenance.eels.detector_type == 'Quantum GIF'
    assert provenance.eels.min_energy.to('eV').magnitude == pytest.approx(100.0)
    assert provenance.eels.max_energy.to('eV').magnitude == pytest.approx(200.0)
    assert provenance.eels.resolution.to('eV').magnitude == pytest.approx(1.0)


@pytest.mark.skipif(simulationworkflowschema is None, reason=SCHEMA_IMPORT_ERROR)
def test_bulk_modulus(mechanical_eos):
    bulk_modulus = mechanical_eos.results.properties.mechanical.bulk_modulus
    assert len(bulk_modulus) == 1
    modulus = bulk_modulus[0]
    assert modulus.type == 'murnaghan'
    assert modulus.value.magnitude == pytest.approx(10000)


@pytest.mark.skipif(simulationworkflowschema is None, reason=SCHEMA_IMPORT_ERROR)
def test_shear_modulus(mechanical_elastic):
    shear_modulus = mechanical_elastic.results.properties.mechanical.shear_modulus
    assert len(shear_modulus) == 3
    modulus_voigt_reuss_hill = shear_modulus[0]
    assert modulus_voigt_reuss_hill.type == 'voigt_reuss_hill_average'
    assert modulus_voigt_reuss_hill.value.magnitude == 10000
    modulus_voigt = shear_modulus[1]
    assert modulus_voigt.type == 'voigt_average'
    assert modulus_voigt.value.magnitude == 10000
    modulus_reuss_hill = shear_modulus[2]
    assert modulus_reuss_hill.type == 'reuss_average'
    assert modulus_reuss_hill.value.magnitude == 10000


@pytest.mark.skipif(simulationworkflowschema is None, reason=SCHEMA_IMPORT_ERROR)
def test_energy_volume_curve(mechanical_eos):
    ev = mechanical_eos.results.properties.mechanical.energy_volume_curve
    assert len(ev) == 2
    ev_raw = ev[0]
    assert ev_raw.type == 'raw'
    assert ev_raw.energies_raw.shape == (10,)
    assert ev_raw.volumes.shape == (10,)
    ev_murhagnan = ev[1]
    assert ev_murhagnan.type == 'murnaghan'
    assert ev_murhagnan.energies_fit.shape == (10,)
    assert ev_murhagnan.volumes.shape == (10,)


def test_band_gap():
    archive = get_template_dft()
    archive = add_template_dos(
        archive, fill=[[[0.0, 1.0], [2.0, 3.0]], [[0.0, 1.0], [1.5, 3.0]]]
    )  # Keep the Fermi level at 1.
    archive = add_template_band_structure(archive, band_gaps=[(1.7, 'direct')])
    archive = run_normalize(archive)

    band_gaps = archive.results.properties.electronic.band_gap
    assert band_gaps[0].value.to('eV').magnitude == pytest.approx(1.0, abs=0.1)
    assert band_gaps[1].value.to('eV').magnitude == pytest.approx(0.5, abs=0.1)
    assert band_gaps[2].value.to('eV').magnitude == pytest.approx(1.7, abs=0.1)


def test_dos_electronic():
    gap_fill = [[0, 1], [2, 3]]

    # DOS without energy references (hence band gap information cannot be extracted)
    archive = get_template_dos()
    dos = archive.results.properties.electronic.dos_electronic_new[0]
    assert dos.spin_polarized is False and len(dos.data) == 1
    assert dos.data[0].total.value.shape == dos.data[0].energies.shape == (101,)

    # Unpolarized DOS with gap:
    efermi = 1.5
    archive = get_template_dos(energy_reference_fermi=efermi)
    dos = archive.results.properties.electronic.dos_electronic_new[0]
    assert dos.spin_polarized is False and len(dos.data) == 1
    assert len(dos.data[0].band_gap) == 1
    assert dos.data[0].band_gap[0].energy_highest_occupied is not None
    assert dos.data[0].band_gap[0].energy_highest_occupied.to(
        'eV'
    ).magnitude == pytest.approx(1.0)
    assert dos.data[0].band_gap[0].energy_lowest_unoccupied is not None
    assert dos.data[0].band_gap[0].energy_lowest_unoccupied.to(
        'eV'
    ).magnitude == pytest.approx(1.9)
    assert dos.data[0].band_gap[0].value.to('eV').magnitude == pytest.approx(0.9)
    assert dos.data[0].total.value.shape == dos.data[0].energies.shape == (101,)

    # Polarized DOS
    efermi = 1.5
    archive = get_template_dos(fill=[gap_fill, gap_fill], energy_reference_fermi=efermi)
    dos = archive.results.properties.electronic.dos_electronic_new[0]
    assert dos.spin_polarized is True and len(dos.data) == 2
    assert dos.data[0].spin_channel == 0 and dos.data[1].spin_channel == 1
    for nspin in range(len(dos.data)):
        assert (
            dos.data[nspin].total.value.shape
            == dos.data[nspin].energies.shape
            == (101,)
        )
    assert dos.data[0].band_gap[0].value == dos.data[1].band_gap[0].value
    assert dos.data[0].band_gap[0].value.to('eV').magnitude == pytest.approx(0.9)

    # Vibrational instead of electronic
    archive = get_template_dos(type='vibrational')
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
    archive = get_template_band_structure([(1, 'direct')], has_references=False)
    bs = archive.results.properties.electronic.band_structure_electronic[0]
    band_gaps = bs.band_gap
    assert bs.reciprocal_cell.shape == (3, 3)
    assert bs.spin_polarized is False
    assert len(band_gaps) == 0
    assert bs.segment[0].energies.shape == (1, 100, 2)
    assert bs.segment[0].kpoints.shape == (100, 3)

    # Unpolarized band structure with no gap
    archive = get_template_band_structure([None])
    bs = archive.results.properties.electronic.band_structure_electronic[0]
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
    bs = archive.results.properties.electronic.band_structure_electronic[0]
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
    gap_type = 'direct'
    archive = get_template_band_structure([(gap, gap_type)])
    bs = archive.results.properties.electronic.band_structure_electronic[0]
    band_gaps = bs.band_gap
    assert bs.reciprocal_cell.shape == (3, 3)
    assert bs.spin_polarized is False
    assert bs.energy_fermi is not None
    assert len(band_gaps) == 1
    assert band_gaps[0].energy_highest_occupied is not None
    assert band_gaps[0].energy_lowest_unoccupied is not None
    assert band_gaps[0].value == pytest.approx(
        (gap * ureg.electron_volt).to(ureg.joule).magnitude
    )
    assert band_gaps[0].type == gap_type
    assert bs.segment[0].energies.shape == (1, 100, 2)
    assert bs.segment[0].kpoints.shape == (100, 3)

    # Unpolarized band structure with indirect gap
    gap = 1  # eV
    gap_type = 'indirect'
    archive = get_template_band_structure([(gap, gap_type)])
    bs = archive.results.properties.electronic.band_structure_electronic[0]
    band_gaps = bs.band_gap
    assert bs.reciprocal_cell.shape == (3, 3)
    assert bs.spin_polarized is False
    assert bs.energy_fermi is not None
    assert len(band_gaps) == 1
    assert band_gaps[0].energy_highest_occupied is not None
    assert band_gaps[0].energy_lowest_unoccupied is not None
    assert band_gaps[0].value == pytest.approx(
        (gap * ureg.electron_volt).to(ureg.joule).magnitude
    )
    assert band_gaps[0].type == gap_type
    assert bs.segment[0].energies.shape == (1, 100, 2)
    assert bs.segment[0].kpoints.shape == (100, 3)

    # Polarized band structure with direct gap
    gap1 = 1  # eV
    gap2 = 2  # eV
    gap_type = 'direct'
    archive = get_template_band_structure([(gap1, gap_type), (gap2, gap_type)])
    bs = archive.results.properties.electronic.band_structure_electronic[0]
    band_gaps = bs.band_gap
    assert bs.reciprocal_cell.shape == (3, 3)
    assert bs.spin_polarized is True
    assert bs.energy_fermi is not None
    assert len(band_gaps) == 2
    assert band_gaps[0].energy_highest_occupied is not None
    assert band_gaps[0].energy_lowest_unoccupied is not None
    assert band_gaps[0].value == pytest.approx(
        (gap1 * ureg.electron_volt).to(ureg.joule).magnitude
    )
    assert band_gaps[0].type == gap_type
    assert band_gaps[1].energy_highest_occupied is not None
    assert band_gaps[1].energy_lowest_unoccupied is not None
    assert band_gaps[1].value == pytest.approx(
        (gap2 * ureg.electron_volt).to(ureg.joule).magnitude
    )
    assert band_gaps[1].type == gap_type
    assert bs.segment[0].energies.shape == (2, 100, 2)
    assert bs.segment[0].kpoints.shape == (100, 3)

    # Polarized band structure with indirect gap
    gap1 = 1  # eV
    gap2 = 2  # eV
    gap_type = 'indirect'
    archive = get_template_band_structure([(gap1, gap_type), (gap2, gap_type)])
    bs = archive.results.properties.electronic.band_structure_electronic[0]
    band_gaps = bs.band_gap
    assert bs.reciprocal_cell.shape == (3, 3)
    assert bs.spin_polarized is True
    assert bs.energy_fermi is not None
    assert len(band_gaps) == 2
    assert band_gaps[0].energy_highest_occupied is not None
    assert band_gaps[0].energy_lowest_unoccupied is not None
    assert band_gaps[0].value == pytest.approx(
        (gap1 * ureg.electron_volt).to(ureg.joule).magnitude
    )
    assert band_gaps[0].type == gap_type
    assert band_gaps[1].energy_highest_occupied is not None
    assert band_gaps[1].energy_lowest_unoccupied is not None
    assert band_gaps[1].value == pytest.approx(
        (gap1 * ureg.electron_volt).to(ureg.joule).magnitude
    )
    assert band_gaps[1].type == gap_type
    assert bs.segment[0].energies.shape == (2, 100, 2)
    assert bs.segment[0].kpoints.shape == (100, 3)


def test_magnetic_properties():
    archive = get_template_dft()
    archive = add_template_magnetic_shielding(archive, n_atoms=2)
    archive = add_template_spin_spin_coupling(archive, n_atoms=2)
    archive = add_template_magnetic_susceptibility(archive)
    archive = run_normalize(archive)
    assert archive.results.properties.magnetic
    magnetic_properties = archive.results.properties.magnetic
    # Magnetic shielding testing
    assert magnetic_properties.magnetic_shielding
    magnetic_shielding = magnetic_properties.magnetic_shielding
    assert len(magnetic_shielding) == 1
    assert magnetic_shielding[0].value.shape == (2, 3, 3)
    assert (
        magnetic_shielding[0].value[0].magnitude
        == [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    ).all()
    assert magnetic_shielding[0].value[1].magnitude == pytest.approx(
        2.0 * magnetic_shielding[0].value[0].magnitude
    )
    # Spin spin coupling
    assert magnetic_properties.spin_spin_coupling
    spin_spin_coupling = magnetic_properties.spin_spin_coupling
    assert len(spin_spin_coupling) == 1
    assert spin_spin_coupling[0].source == 'simulation'
    assert spin_spin_coupling[0].contribution == 'total'
    assert spin_spin_coupling[0].value.shape == (2, 2, 3, 3)
    # Magnetic susceptibility
    assert magnetic_properties.magnetic_susceptibility
    magnetic_susceptibility = magnetic_properties.magnetic_susceptibility
    assert len(magnetic_susceptibility) == 1
    assert magnetic_susceptibility[0].source == 'simulation'
    assert magnetic_susceptibility[0].scale_dimension == 'macroscopic'
    assert (
        magnetic_susceptibility[0].value
        == [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    ).all()


@pytest.mark.skipif(simulationworkflowschema is None, reason=SCHEMA_IMPORT_ERROR)
def test_dos_phonon():
    # DOS with all correct metainfo
    archive = get_template_dos(type='vibrational')
    dos = archive.results.properties.vibrational.dos_phonon
    assert dos.total[0].value.shape == (101,)
    assert dos.energies.shape == (101,)

    # Electronic instead of vibrational
    archive = get_template_dos(type='electronic')
    vibrational = archive.results.properties.vibrational
    assert vibrational is None

    # Empty values
    archive = get_template_dos(type='vibrational', normalize=False)
    archive.run[0].calculation[0].dos_phonon[0].total = None
    archive = run_normalize(archive)
    vibrational = archive.results.properties.vibrational
    assert vibrational is None

    # Empty energies
    archive = get_template_dos(type='vibrational', normalize=False)
    archive.run[0].calculation[0].dos_phonon[0].energies = []
    archive = run_normalize(archive)
    vibrational = archive.results.properties.vibrational
    assert vibrational is None


@pytest.mark.skipif(simulationworkflowschema is None, reason=SCHEMA_IMPORT_ERROR)
def test_band_structure_phonon():
    # Valid phonon band structure
    archive = get_template_band_structure(type='vibrational')
    bs = archive.results.properties.vibrational.band_structure_phonon
    assert bs.segment[0].energies.shape == (1, 100, 2)
    assert bs.segment[0].kpoints.shape == (100, 3)


@pytest.mark.skipif(simulationworkflowschema is None, reason=SCHEMA_IMPORT_ERROR)
def test_energy_free_helmholtz(phonon):
    energy_free = phonon.results.properties.vibrational.energy_free_helmholtz
    assert energy_free.temperatures.shape == (11,)
    assert energy_free.energies.shape == (11,)


@pytest.mark.skipif(simulationworkflowschema is None, reason=SCHEMA_IMPORT_ERROR)
def test_heat_capacity_constant_volume(phonon):
    heat_cap = phonon.results.properties.vibrational.heat_capacity_constant_volume
    assert heat_cap.temperatures.shape == (11,)
    assert heat_cap.heat_capacities.shape == (11,)


@pytest.mark.skipif(simulationworkflowschema is None, reason=SCHEMA_IMPORT_ERROR)
def test_geometry_optimization(geometry_optimization):
    geo_opt_prop = geometry_optimization.results.properties.geometry_optimization
    n_frames = len(geo_opt_prop.trajectory)
    n_energies = len(geo_opt_prop.energies)
    assert geo_opt_prop.system_optimized is not None
    assert n_frames > 0
    assert n_frames == n_energies
    assert geo_opt_prop.final_energy_difference > 0
    assert geo_opt_prop.type == 'atomic'


@pytest.mark.skipif(simulationworkflowschema is None, reason=SCHEMA_IMPORT_ERROR)
def test_trajectory(molecular_dynamics):
    trajectories = molecular_dynamics.results.properties.thermodynamic.trajectory
    n_trajectories = len(trajectories)
    assert n_trajectories == 1
    n_steps = 10
    trajectory = trajectories[0]

    assert 'trajectory' in molecular_dynamics.results.properties.available_properties
    assert trajectory.pressure.value.size == trajectory.pressure.time.size == n_steps
    assert trajectory.volume.value.size == trajectory.volume.time.size == n_steps
    assert (
        trajectory.temperature.value.size == trajectory.temperature.time.size == n_steps
    )
    assert (
        trajectory.energy_potential.value.size
        == trajectory.energy_potential.time.size
        == n_steps
    )
    assert trajectory.provenance.molecular_dynamics.time_step == 0.5 * ureg('fs')
    assert trajectory.provenance.molecular_dynamics.ensemble_type == 'NVT'
    assert set(trajectory.available_properties) == set(
        ['pressure', 'volume', 'temperature', 'energy_potential']
    )


@pytest.mark.skipif(simulationworkflowschema is None, reason=SCHEMA_IMPORT_ERROR)
def test_rgs(molecular_dynamics):
    rgs = molecular_dynamics.results.properties.structural.radius_of_gyration
    n_rgs = len(rgs)
    assert n_rgs == 1
    n_steps = 10
    rg = rgs[0]

    assert (
        'radius_of_gyration'
        in molecular_dynamics.results.properties.available_properties
    )
    assert rg.kind == 'molecular'
    assert rg.value.size == rg.time.size == n_steps
    assert rg.label == 'MOL'


@pytest.mark.skipif(simulationworkflowschema is None, reason=SCHEMA_IMPORT_ERROR)
def test_rdfs(molecular_dynamics):
    rdfs = molecular_dynamics.results.properties.structural.radial_distribution_function
    n_rdfs = len(rdfs)
    assert n_rdfs == 1
    rdf = rdfs[0]

    assert (
        'radial_distribution_function'
        in molecular_dynamics.results.properties.available_properties
    )
    assert rdf.type == 'molecular'
    assert rdf.label == 'MOL-MOL'
    assert np.array_equal(rdf.bins.to(ureg.meter).magnitude, [0, 1, 2])
    assert np.array_equal(rdf.value, [0, 1, 2])
    assert rdf.n_bins == len(rdf.bins)
    assert rdf.frame_start == 0
    assert rdf.frame_end == 100


@pytest.mark.skipif(simulationworkflowschema is None, reason=SCHEMA_IMPORT_ERROR)
def test_msds(molecular_dynamics):
    msds = molecular_dynamics.results.properties.dynamical.mean_squared_displacement
    n_msds = len(msds)
    assert n_msds == 1
    msd = msds[0]

    assert (
        'mean_squared_displacement'
        in molecular_dynamics.results.properties.available_properties
    )
    assert msd.type == 'molecular'
    assert msd.direction == 'xyz'
    assert msd.error_type == 'bootstrapping'
    assert msd.label == 'MOL'
    assert np.array_equal(msd.times.to(ureg.second).magnitude, [0, 1, 2])
    assert np.array_equal(msd.value.to(ureg.meter * ureg.meter).magnitude, [0, 1, 2])
    assert msd.n_times == len(msd.times)
    assert np.array_equal(msd.errors, [0, 1, 2])
    assert np.array_equal(
        msd.diffusion_constant_value.to(
            ureg.meter * ureg.meter / ureg.second
        ).magnitude,
        2.1,
    )
    assert msd.diffusion_constant_error_type == 'Pearson correlation coefficient'
    assert np.array_equal(msd.diffusion_constant_errors, [0.98])


@pytest.mark.skipif(simulationworkflowschema is None, reason=SCHEMA_IMPORT_ERROR)
def test_n_calculations(geometry_optimization):
    assert geometry_optimization.results.properties.n_calculations == 2
