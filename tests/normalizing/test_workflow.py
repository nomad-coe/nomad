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

import tests
from .conftest import run_normalize


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


@pytest.fixture(scope='function')
def workflow_archive():
    def _archive(parser_name, filepath):
        archive = tests.parsing.test_parsing.parse_file((parser_name, filepath))
        return run_normalize(archive)

    return _archive


def test_no_workflow(workflow_archive):
    vasp_archive = workflow_archive(
        'parsers/vasp', 'tests/data/parsers/vasp_outcar/OUTCAR_broken')
    assert not vasp_archive.workflow[0].calculations_ref


def test_single_point_workflow(workflow_archive):
    vasp_archive = workflow_archive(
        'parsers/vasp', 'tests/data/normalizers/workflow/vasp/vasprun.xml.static')
    sec_workflow = vasp_archive.workflow[0]
    assert sec_workflow.type == 'single_point'
    assert sec_workflow.calculations_ref is not None
    assert sec_workflow.calculation_result_ref.m_def.name == 'Calculation'
    assert sec_workflow.single_point.method == 'DFT'
    assert sec_workflow.single_point.n_scf_steps == 9
    assert sec_workflow.single_point.final_scf_energy_difference > 0
    assert sec_workflow.single_point.with_density_of_states
    assert not sec_workflow.single_point.with_bandstructure
    assert sec_workflow.single_point.with_eigenvalues
    assert not sec_workflow.single_point.with_volumetric_data
    assert not sec_workflow.single_point.with_spectra
    assert sec_workflow.single_point.is_converged
    sec_workflow2 = vasp_archive.workflow2
    assert sec_workflow2.method.method == 'DFT'
    assert sec_workflow2.results.n_scf_steps == 9
    assert sec_workflow2.results.final_scf_energy_difference > 0
    assert sec_workflow2.results.dos is not None
    assert sec_workflow2.results.band_structure is None
    assert sec_workflow2.results.eigenvalues is not None
    assert sec_workflow2.results.density_charge is None
    assert sec_workflow2.results.spectra is None
    assert sec_workflow2.results.is_converged


def test_geometry_optimization_workflow(workflow_archive):
    vasp_archive = workflow_archive(
        'parsers/vasp', 'tests/data/normalizers/workflow/vasp/vasprun.xml')
    sec_workflow = vasp_archive.workflow[0]

    assert sec_workflow.type == 'geometry_optimization'
    assert sec_workflow.calculations_ref is not None
    assert sec_workflow.calculation_result_ref.m_def.name == 'Calculation'
    assert sec_workflow.geometry_optimization.type == 'cell_shape'
    assert sec_workflow.geometry_optimization.final_energy_difference > 0.0
    assert sec_workflow.geometry_optimization.optimization_steps == 3
    assert sec_workflow.geometry_optimization.final_force_maximum > 0.0
    assert sec_workflow.geometry_optimization.is_converged_geometry

    task = sec_workflow.task
    assert len(task) == len(sec_workflow.calculations_ref) + 1
    assert task[1].input_calculation == sec_workflow.calculations_ref[0]
    assert task[-1].output_workflow == sec_workflow

    sec_workflow2 = vasp_archive.workflow2
    assert sec_workflow2.method.type == 'cell_shape'
    assert sec_workflow2.results.calculation_result_ref.m_def.name == 'Calculation'
    assert sec_workflow2.results.final_energy_difference > 0.0
    assert sec_workflow2.results.optimization_steps == 3
    assert sec_workflow2.results.final_force_maximum > 0.0
    assert sec_workflow2.results.is_converged_geometry


def test_elastic_workflow(workflow_archive):
    elastic_archive = workflow_archive(
        'parsers/elastic', "tests/data/normalizers/workflow/elastic/INFO_ElaStic")
    sec_workflow = elastic_archive.workflow[0]

    assert sec_workflow.type == 'elastic'
    assert sec_workflow.calculations_ref is not None
    assert sec_workflow.calculation_result_ref.m_def.name == 'Calculation'
    assert sec_workflow.elastic.calculation_method == 'energy'
    assert sec_workflow.elastic.elastic_constants_order == 2
    assert sec_workflow.elastic.is_mechanically_stable
    assert sec_workflow.elastic.fitting_error_maximum > 0.0
    assert sec_workflow.elastic.strain_maximum > 0.0

    sec_workflow2 = elastic_archive.workflow2
    sec_workflow2.results.calculation_result_ref.m_def.name == 'Calculation'
    sec_workflow2.method.calculation_method == 'energy'
    sec_workflow2.method.elastic_constants_order == 2
    sec_workflow2.results.is_mechanically_stable
    sec_workflow2.method.fitting_error_maximum > 0.0
    sec_workflow2.method.strain_maximum > 0.0


def test_phonon_workflow(workflow_archive):
    phonopy_archive = workflow_archive(
        'parsers/phonopy',
        'tests/data/normalizers/workflow/phonopy/phonopy-FHI-aims-displacement-01/control.in')

    sec_workflow = phonopy_archive.workflow[0]
    assert sec_workflow.type == 'phonon'
    assert sec_workflow.calculations_ref is not None
    assert sec_workflow.calculation_result_ref.m_def.name == 'Calculation'
    assert sec_workflow.phonon.force_calculator == 'fhi-aims'
    assert sec_workflow.phonon.mesh_density > 0.0
    assert sec_workflow.phonon.n_imaginary_frequencies > 0
    assert not sec_workflow.phonon.random_displacements
    assert not sec_workflow.phonon.with_non_analytic_correction
    assert not sec_workflow.phonon.with_grueneisen_parameters

    sec_workflow2 = phonopy_archive.workflow2
    assert sec_workflow2.results.calculation_result_ref.m_def.name == 'Calculation'
    assert sec_workflow2.method.force_calculator == 'fhi-aims'
    assert sec_workflow2.method.mesh_density > 0.0
    assert sec_workflow2.results.n_imaginary_frequencies > 0
    assert not sec_workflow2.method.random_displacements
    assert not sec_workflow2.method.with_non_analytic_correction
    assert not sec_workflow2.method.with_grueneisen_parameters


def test_molecular_dynamics_workflow(workflow_archive):
    lammmps_archive = workflow_archive(
        'parsers/lammps', 'tests/data/normalizers/workflow/lammps/log.lammps')

    sec_workflow = lammmps_archive.workflow[0]
    assert sec_workflow.type == 'molecular_dynamics'
    assert sec_workflow.calculations_ref is not None
    assert sec_workflow.calculation_result_ref.m_def.name == 'Calculation'
    assert sec_workflow.molecular_dynamics.finished_normally
    assert sec_workflow.molecular_dynamics.with_trajectory

    sec_workflow2 = lammmps_archive.workflow2
    sec_workflow2.results.calculation_result_ref.m_def.name == 'Calculation'
    assert sec_workflow2.results.finished_normally
    assert sec_workflow2.results.trajectory


def test_rdf_and_msd(workflow_archive):
    archive = workflow_archive(
        'parsers/lammps', 'tests/data/parsers/lammps/hexane_cyclohexane/log.hexane_cyclohexane_nvt')

    sec_workflow2 = archive.workflow2
    section_md = sec_workflow2.results

    assert section_md.radial_distribution_functions[0].type == 'molecular'
    assert section_md.radial_distribution_functions[0].n_smooth == 2
    assert section_md.radial_distribution_functions[0].variables_name[0] == 'distance'

    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[0].label == '0-0'
    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[0].n_bins == 198
    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[0].bins[122].magnitude == approx(6.923255643844605 * 10**(-10))
    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[0].bins[122].units == 'meter'
    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[0].value[96] == approx(0.0)
    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[0].frame_start == 0
    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[0].frame_end == 40

    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[3].label == '0-0'
    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[3].n_bins == 198
    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[3].bins[65].magnitude == approx(3.727906885147095 * 10**(-10))
    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[3].bins[65].units == 'meter'
    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[3].value[52] == approx(0.0)
    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[3].frame_start == 120
    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[3].frame_end == 201

    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[5].label == '1-0'
    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[5].n_bins == 198
    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[5].bins[102].magnitude == approx(5.802080640792847 * 10**(-10))
    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[5].bins[102].units == 'meter'
    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[5].value[55] == approx(0.0)
    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[5].frame_start == 40
    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[5].frame_end == 201

    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[10].label == '1-1'
    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[10].n_bins == 198
    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[10].bins[44].magnitude == approx(2.550673131942749 * 10**(-10))
    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[10].bins[44].units == 'meter'
    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[10].value[101] == approx(1.4750986777470825)
    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[10].frame_start == 80
    assert section_md.radial_distribution_functions[0].radial_distribution_function_values[10].frame_end == 201

    assert section_md.mean_squared_displacements[0].type == 'molecular'
    assert section_md.mean_squared_displacements[0].direction == 'xyz'

    assert section_md.mean_squared_displacements[0].mean_squared_displacement_values[0].label == '0'
    assert section_md.mean_squared_displacements[0].mean_squared_displacement_values[0].n_times == 54
    assert section_md.mean_squared_displacements[0].mean_squared_displacement_values[0].times[13].magnitude == approx(1.3 * 10**(-11))
    assert section_md.mean_squared_displacements[0].mean_squared_displacement_values[0].times[13].units == 'second'
    assert section_md.mean_squared_displacements[0].mean_squared_displacement_values[0].value[32].magnitude == approx(8.98473539965496 * 10**(-19))
    assert section_md.mean_squared_displacements[0].mean_squared_displacement_values[0].value[32].units == 'meter^2'
    assert section_md.mean_squared_displacements[0].mean_squared_displacement_values[0].diffusion_constant.value.magnitude == approx(6.09812270414572 * 10**(-9))
    assert section_md.mean_squared_displacements[0].mean_squared_displacement_values[0].diffusion_constant.value.units == 'meter^2/second'
    assert section_md.mean_squared_displacements[0].mean_squared_displacement_values[0].diffusion_constant.error_type == 'Pearson correlation coefficient'
    assert section_md.mean_squared_displacements[0].mean_squared_displacement_values[0].diffusion_constant.errors == approx(0.9924847048341159)

    assert section_md.mean_squared_displacements[0].mean_squared_displacement_values[1].label == '1'
    assert section_md.mean_squared_displacements[0].mean_squared_displacement_values[1].n_times == 54
    assert section_md.mean_squared_displacements[0].mean_squared_displacement_values[1].times[13].magnitude == approx(1.3 * 10**(-11))
    assert section_md.mean_squared_displacements[0].mean_squared_displacement_values[1].times[13].units == 'second'
    assert section_md.mean_squared_displacements[0].mean_squared_displacement_values[1].value[32].magnitude == approx(8.448369705677565 * 10**(-19))
    assert section_md.mean_squared_displacements[0].mean_squared_displacement_values[1].value[32].units == 'meter^2'
    assert section_md.mean_squared_displacements[0].mean_squared_displacement_values[1].diffusion_constant.value.magnitude == approx(5.094072039759048 * 10**(-9))
    assert section_md.mean_squared_displacements[0].mean_squared_displacement_values[1].diffusion_constant.value.units == 'meter^2/second'
    assert section_md.mean_squared_displacements[0].mean_squared_displacement_values[1].diffusion_constant.error_type == 'Pearson correlation coefficient'
    assert section_md.mean_squared_displacements[0].mean_squared_displacement_values[1].diffusion_constant.errors == approx(0.9965870174917716)


def test_radius_of_gyration(workflow_archive):
    archive = workflow_archive(
        'parsers/gromacs', 'tests/data/parsers/gromacs/protein_fsfg/nvt.log')

    sec_calc = archive.run[0].calculation[4]
    sec_rg = sec_calc.radius_of_gyration[0]
    sec_rgvals = sec_rg.radius_of_gyration_values[0]

    assert sec_rg.kind == 'molecular'

    assert sec_rgvals.label == 'Protein-index_0'
    assert sec_rgvals.value.magnitude == approx(5.464423436523278e-10)
    assert sec_rgvals.value.units == 'meter'

    sec_calc = archive.run[0].calculation[7]
    sec_rg = sec_calc.radius_of_gyration[0]
    sec_rgvals = sec_rg.radius_of_gyration_values[1]

    assert sec_rg.kind == 'molecular'
    assert sec_rgvals.label == 'Protein-index_1'
    assert sec_rgvals.value.magnitude == approx(7.326346215313874e-10)
    assert sec_rgvals.value.units == 'meter'

    sec_workflow2 = archive.workflow2
    sec_rg = sec_workflow2.results.radius_of_gyration[0]
    frame = 4

    assert sec_rg.type == 'molecular'

    assert sec_rg.label == 'Protein-index_0'
    assert sec_rg.value[frame].magnitude == approx(5.464423436523278e-10)
    assert sec_rg.value[frame].units == 'meter'

    frame = 7
    sec_rg = sec_workflow2.results.radius_of_gyration[1]
    sec_calc = archive.run[0].calculation[7]

    assert sec_rg.type == 'molecular'
    assert sec_rg.label == 'Protein-index_1'
    assert sec_rg.value[frame].magnitude == approx(7.326346215313874e-10)
    assert sec_rg.value[frame].units == 'meter'
