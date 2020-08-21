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


def test_vasp_wokflow(dos_si_vasp):
    sec_workflow = dos_si_vasp.section_workflow
    assert sec_workflow.workflow_type == 'geometry_optimization'
    assert sec_workflow.section_relaxation.relaxation_type == 'cell_shape'
    assert sec_workflow.section_relaxation.final_calculation_ref.m_def.name == 'section_single_configuration_calculation'
    assert sec_workflow.section_relaxation.final_energy_difference > 0.0
