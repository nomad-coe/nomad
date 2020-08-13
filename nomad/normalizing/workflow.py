# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the 'License');
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an'AS IS' BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import numpy as np

from nomad.normalizing.normalizer import Normalizer
from nomad.datamodel.metainfo.public import Workflow, Relaxation


class WorkflowNormalizer(Normalizer):
    '''
    This normalizer performs all produces a section all data necessary for the Optimade API.
    It assumes that the :class:`SystemNormalizer` was run before.
    '''
    def __init__(self, backend):
        super().__init__(backend)

    def _get_relaxation_type(self):
        sec_system = self.section_run.section_system
        if not sec_system:
            return

        def compare_cell(cell1, cell2):
            if (cell1 == cell2).all():
                return None
            else:
                cell1_normed = cell1 / np.linalg.norm(cell1)
                cell2_normed = cell2 / np.linalg.norm(cell2)
                if (cell1_normed == cell2_normed).all():
                    return 'cell_volume'
                else:
                    return 'cell_shape'

        if len(sec_system) < 2:
            return 'static'

        else:
            cell_init = np.array(sec_system[0].lattice_vectors)
            cell_final = np.array(sec_system[-1].lattice_vectors)

            cell_relaxation = compare_cell(cell_init, cell_final)

            if cell_relaxation is not None:
                return cell_relaxation

            atom_pos_init = np.array(sec_system[0].atom_positions)
            atom_pos_final = np.array(sec_system[-1].atom_positions)

            if (atom_pos_init == atom_pos_final).all():
                return 'static'

            return 'ionic'

    def normalize_relaxation(self):
        sec_relaxation = self.entry_archive.section_workflow.section_relaxation
        if not sec_relaxation:
            sec_relaxation = self.entry_archive.section_workflow.m_create(Relaxation)

        if not sec_relaxation.relaxation_type:
            sec_relaxation.relaxation_type = self._get_relaxation_type()

        if not sec_relaxation.final_calculation_ref:
            scc = self.section_run.section_single_configuration_calculation
            if scc:
                sec_relaxation.final_calculation_ref = scc[-1]

        if not sec_relaxation.final_energy_difference:
            energies = []
            for scc in self.section_run.section_single_configuration_calculation:
                if scc.energy_total_T0:
                    energies.append(scc.energy_total_T0)

            delta_energy = None
            if len(energies) > 1:
                delta_energy = abs(energies[-1] - energies[-2])

            if delta_energy == 0.0:
                try:
                    delta_energy = abs(energies[-1] - energies[-3])
                except Exception:
                    pass

            if delta_energy:
                sec_relaxation.final_energy_difference = delta_energy

        if not sec_relaxation.final_force_maximum:
            scc = self.section_run.section_single_configuration_calculation
            max_force = None
            if scc:
                if scc[-1].atom_forces is not None:
                    forces = np.array(scc[-1].atom_forces)
                    max_force = np.max(np.linalg.norm(forces, axis=1))
            if max_force is not None:
                sec_relaxation.final_force_maximum = max_force

    def normalize(self, logger=None) -> None:
        # Setup logger
        if logger is not None:
            self.logger = logger.bind(normalizer=self.__class__.__name__)

        # Do nothing if section_run is not present
        if self.section_run is None:
            return

        sec_workflow = self.entry_archive.section_workflow
        if not sec_workflow:
            sec_workflow = self.entry_archive.m_create(Workflow)

        if not sec_workflow.workflow_type:
            sec_sampling_method = self.section_run.section_sampling_method
            if sec_sampling_method:
                workflow_type = sec_sampling_method[-1].sampling_method
                # TODO imho geometry_optimization is not an appropriate name
                # if workflow_type == 'geometry_optimization':
                #     workflow_type = 'relaxation'
                sec_workflow.workflow_type = workflow_type

        if sec_workflow.workflow_type in ['geometry_optimization', 'relaxation']:
            self.normalize_relaxation()
