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

from tests.parsing.test_parsing import parse_file
from .conftest import run_normalize


@pytest.fixture(scope='function')
def workflow_archive():
    def _archive(parser_name, filepath):
        archive = parse_file((parser_name, filepath))
        return run_normalize(archive)

    return _archive


def test_no_workflow(workflow_archive):
    vasp_archive = workflow_archive(
        'parsers/vasp', 'tests/data/parsers/vasp_outcar/OUTCAR_broken')
    assert vasp_archive.section_workflow is None


def test_geometry_optimization_workflow(workflow_archive):
    vasp_archive = workflow_archive(
        'parsers/vasp', 'tests/data/normalizers/workflow/vasp/vasprun.xml')
    sec_workflow = vasp_archive.section_workflow

    assert sec_workflow.workflow_type == 'geometry_optimization'
    assert sec_workflow.calculations_ref is not None
    assert sec_workflow.calculation_result_ref.m_def.name == 'SingleConfigurationCalculation'
    assert sec_workflow.section_geometry_optimization.geometry_optimization_type == 'cell_shape'
    assert sec_workflow.section_geometry_optimization.final_energy_difference > 0.0
    assert sec_workflow.section_geometry_optimization.optimization_steps == 3
    assert sec_workflow.section_geometry_optimization.final_force_maximum > 0.0


def test_elastic_workflow(workflow_archive):
    elastic_archive = workflow_archive(
        'parsers/elastic', "tests/data/normalizers/workflow/elastic/INFO_ElaStic")
    sec_workflow = elastic_archive.section_workflow

    assert sec_workflow.workflow_type == 'elastic'
    assert sec_workflow.calculations_ref is not None
    assert sec_workflow.calculation_result_ref.m_def.name == 'SingleConfigurationCalculation'
    assert sec_workflow.section_elastic.elastic_calculation_method == 'energy'
    assert sec_workflow.section_elastic.elastic_constants_order == 2
    assert sec_workflow.section_elastic.is_mechanically_stable
    assert sec_workflow.section_elastic.fitting_error_maximum > 0.0
    assert sec_workflow.section_elastic.strain_maximum > 0.0


def test_phonon_workflow(workflow_archive):
    phonopy_archive = workflow_archive(
        'parsers/phonopy',
        'tests/data/normalizers/workflow/phonopy/phonopy-FHI-aims-displacement-01/control.in')

    sec_workflow = phonopy_archive.section_workflow
    assert sec_workflow.workflow_type == 'phonon'
    assert sec_workflow.calculations_ref is not None
    assert sec_workflow.calculation_result_ref.m_def.name == 'SingleConfigurationCalculation'
    assert sec_workflow.section_phonon.force_calculator == 'fhi-aims'
    assert sec_workflow.section_phonon.mesh_density > 0.0
    assert sec_workflow.section_phonon.n_imaginary_frequencies > 0
    assert not sec_workflow.section_phonon.random_displacements
    assert not sec_workflow.section_phonon.with_non_analytic_correction
    assert not sec_workflow.section_phonon.with_grueneisen_parameters


def test_molecular_dynamics_workflow(workflow_archive):
    lammmps_archive = workflow_archive(
        'parsers/lammps', 'tests/data/normalizers/workflow/lammps/log.lammps')

    sec_workflow = lammmps_archive.section_workflow
    assert sec_workflow.workflow_type == 'molecular_dynamics'
    assert sec_workflow.calculations_ref is not None
    assert sec_workflow.calculation_result_ref.m_def.name == 'SingleConfigurationCalculation'
    assert sec_workflow.section_molecular_dynamics.finished_normally
    assert sec_workflow.section_molecular_dynamics.with_trajectory
    assert sec_workflow.section_molecular_dynamics.with_thermodynamics
