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
    assert not vasp_archive.workflow2.results.calculations_ref


def test_single_point_workflow(workflow_archive):
    vasp_archive = workflow_archive(
        'parsers/vasp', 'tests/data/normalizers/workflow/vasp/vasprun.xml.static')
    sec_workflow = vasp_archive.workflow2
    assert sec_workflow.method.method == 'DFT'
    assert sec_workflow.results.n_scf_steps == 9
    assert sec_workflow.results.final_scf_energy_difference > 0
    assert sec_workflow.results.dos is not None
    assert sec_workflow.results.band_structure is None
    assert sec_workflow.results.eigenvalues is not None
    assert sec_workflow.results.density_charge is None
    assert sec_workflow.results.spectra is None
    assert sec_workflow.results.is_converged


def test_gw_workflow(gw_workflow):
    """Testing GW workflow entry"""
    workflow = gw_workflow.workflow2
    assert workflow.name == 'GW'
    assert workflow.method.gw_method_ref.type == 'G0W0'
    assert workflow.method.electrons_representation.type == 'plane waves'
    assert workflow.method.starting_point.name == 'GGA_X_PBE'
    results = gw_workflow.results
    assert results.method.method_name == 'GW'
    assert results.method.workflow_name == 'GW'
    assert results.method.simulation.program_name == 'VASP'
    assert results.method.simulation.program_version == '4.6.35'
    assert results.method.simulation.gw.type == 'G0W0'
    assert results.method.simulation.gw.starting_point_type == 'GGA'
    assert results.method.simulation.gw.starting_point_names == ['GGA_X_PBE']
    assert results.method.simulation.gw.basis_set_type == 'plane waves'
    assert not results.properties.electronic.band_gap
    assert not results.properties.electronic.greens_functions_electronic
    assert len(results.properties.electronic.dos_electronic) == 2
    assert len(results.properties.electronic.band_structure_electronic) == 2
    assert results.properties.electronic.dos_electronic[0].label == 'DFT'
    assert results.properties.electronic.dos_electronic[1].label == 'GW'


def test_geometry_optimization_workflow(workflow_archive):
    vasp_archive = workflow_archive(
        'parsers/vasp', 'tests/data/normalizers/workflow/vasp/vasprun.xml')
    sec_workflow = vasp_archive.workflow2
    assert sec_workflow.method.type == 'cell_shape'
    assert sec_workflow.results.calculation_result_ref.m_def.name == 'Calculation'
    assert sec_workflow.results.final_energy_difference > 0.0
    assert sec_workflow.results.optimization_steps == 3
    assert sec_workflow.results.final_force_maximum > 0.0
    assert sec_workflow.results.is_converged_geometry

    tasks = sec_workflow.tasks
    assert len(tasks) == len(vasp_archive.run[0].calculation)
    assert tasks[0].inputs[0].section == vasp_archive.run[0].system[0]
    assert tasks[-1].outputs[0].section == vasp_archive.run[0].calculation[-1]


def test_elastic_workflow(workflow_archive):
    elastic_archive = workflow_archive(
        'parsers/elastic', "tests/data/normalizers/workflow/elastic/INFO_ElaStic")
    sec_workflow = elastic_archive.workflow2
    sec_workflow.results.calculation_result_ref.m_def.name == 'Calculation'
    sec_workflow.method.calculation_method == 'energy'
    sec_workflow.method.elastic_constants_order == 2
    sec_workflow.results.is_mechanically_stable
    sec_workflow.method.fitting_error_maximum > 0.0
    sec_workflow.method.strain_maximum > 0.0


def test_phonon_workflow(workflow_archive):
    phonopy_archive = workflow_archive(
        'parsers/phonopy',
        'tests/data/normalizers/workflow/phonopy/phonopy-FHI-aims-displacement-01/control.in')

    sec_workflow = phonopy_archive.workflow2
    assert sec_workflow.results.calculation_result_ref.m_def.name == 'Calculation'
    assert sec_workflow.method.force_calculator == 'fhi-aims'
    assert sec_workflow.method.mesh_density > 0.0
    assert sec_workflow.results.n_imaginary_frequencies > 0
    assert not sec_workflow.method.random_displacements
    assert not sec_workflow.method.with_non_analytic_correction
    assert not sec_workflow.method.with_grueneisen_parameters


def test_molecular_dynamics_workflow(workflow_archive):
    lammmps_archive = workflow_archive(
        'parsers/lammps', 'tests/data/normalizers/workflow/lammps/log.lammps')

    sec_workflow = lammmps_archive.workflow2
    sec_workflow.results.calculation_result_ref.m_def.name == 'Calculation'
    assert sec_workflow.results.finished_normally
    assert sec_workflow.results.trajectory


def test_rdf_and_msd(workflow_archive):
    archive = workflow_archive(
        'parsers/lammps', 'tests/data/parsers/lammps/hexane_cyclohexane/log.hexane_cyclohexane_nvt')

    sec_workflow = archive.workflow2
    section_md = sec_workflow.results

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
        'parsers/gromacs', 'tests/data/parsers/gromacs/protein_small_nowater/md.log')

    sec_calc = archive.run[0].calculation[4]
    sec_rg = sec_calc.radius_of_gyration[0]
    sec_rgvals = sec_rg.radius_of_gyration_values[0]

    assert sec_rg.kind == 'molecular'

    assert sec_rgvals.label == 'Protein_chain_X-index_0'
    assert sec_rgvals.value.magnitude == approx(5.081165959952965e-10)
    assert sec_rgvals.value.units == 'meter'

    sec_calc = archive.run[0].calculation[1]
    sec_rg = sec_calc.radius_of_gyration[0]
    sec_rgvals = sec_rg.radius_of_gyration_values[0]

    assert sec_rg.kind == 'molecular'
    assert sec_rgvals.label == 'Protein_chain_X-index_0'
    assert sec_rgvals.value.magnitude == approx(5.036762961380965e-10)
    assert sec_rgvals.value.units == 'meter'

    sec_workflow = archive.workflow2
    sec_rg = sec_workflow.results.radius_of_gyration[0]
    frame = 4

    assert sec_rg.type == 'molecular'

    assert sec_rg.label == 'Protein_chain_X-index_0'
    assert sec_rg.value[frame].magnitude == approx(5.081165959952965e-10)
    assert sec_rg.value[frame].units == 'meter'

    frame = 1
    sec_rg = sec_workflow.results.radius_of_gyration[0]
    sec_calc = archive.run[0].calculation[1]

    assert sec_rg.type == 'molecular'
    assert sec_rg.label == 'Protein_chain_X-index_0'
    assert sec_rg.value[frame].magnitude == approx(5.036762961380965e-10)
    assert sec_rg.value[frame].units == 'meter'
