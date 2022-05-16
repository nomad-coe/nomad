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


@pytest.fixture(scope='function')
def workflow_archive():
    def _archive(parser_name, filepath):
        archive = tests.parsing.test_parsing.parse_file((parser_name, filepath))
        return run_normalize(archive)

    return _archive


def test_no_workflow(workflow_archive):
    vasp_archive = workflow_archive(
        'parsers/vasp', 'tests/data/parsers/vasp_outcar/OUTCAR_broken')
    assert not vasp_archive.workflow


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
    assert not sec_workflow.single_point.with_excited_states
    assert sec_workflow.single_point.is_converged


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


def test_molecular_dynamics_workflow(workflow_archive):
    lammmps_archive = workflow_archive(
        'parsers/lammps', 'tests/data/normalizers/workflow/lammps/log.lammps')

    sec_workflow = lammmps_archive.workflow[0]
    assert sec_workflow.type == 'molecular_dynamics'
    assert sec_workflow.calculations_ref is not None
    assert sec_workflow.calculation_result_ref.m_def.name == 'Calculation'
    assert sec_workflow.molecular_dynamics.finished_normally
    assert sec_workflow.molecular_dynamics.with_trajectory
