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
import pint

from nomad.normalizing.normalizer import Normalizer
from nomad.datamodel.metainfo.public import Workflow, Relaxation, Phonon, Elastic


class RelaxationNormalizer(Normalizer):
    def __init__(self, entry_archive):
        super().__init__(entry_archive)

    def _to_numpy_array(self, quantity):
        try:
            return np.array(quantity)
        except pint.UnitStrippedWarning as e:
            self.logger.warn(e)

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
            cell_init = self._to_numpy_array(sec_system[0].lattice_vectors)
            cell_final = self._to_numpy_array(sec_system[-1].lattice_vectors)

            cell_relaxation = compare_cell(cell_init, cell_final)

            if cell_relaxation is not None:
                return cell_relaxation

            atom_pos_init = self._to_numpy_array(sec_system[0].atom_positions)
            atom_pos_final = self._to_numpy_array(sec_system[-1].atom_positions)

            if (atom_pos_init == atom_pos_final).all():
                return 'static'

            return 'ionic'

    def normalize(self):
        self.section = self.entry_archive.section_workflow.section_relaxation
        if not self.section:
            self.section = self.entry_archive.section_workflow.m_create(Relaxation)

        if not self.section.relaxation_type:
            self.section.relaxation_type = self._get_relaxation_type()

        scc = self.section_run.section_single_configuration_calculation
        if not self.section.final_calculation_ref:
            if scc:
                self.section.final_calculation_ref = scc[-1]

        if not self.section.n_relaxation_steps:
            self.section.n_relaxation_steps = len(scc)

        if not self.section.calculations_ref:
            self.section.calculations_ref = scc

        if not self.section.final_energy_difference:
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
                self.section.final_energy_difference = delta_energy

        if not self.section.final_force_maximum:
            scc = self.section_run.section_single_configuration_calculation
            max_force = None
            if scc:
                if scc[-1].atom_forces is not None:
                    forces = self._to_numpy_array(scc[-1].atom_forces)
                    max_force = np.max(np.linalg.norm(forces, axis=1))
            if max_force is not None:
                self.section.final_force_maximum = max_force


class PhononNormalizer(Normalizer):
    def __init__(self, entry_archive):
        super().__init__(entry_archive)

    def _get_n_imaginary_frequencies(self):
        scc = self.entry_archive.section_run[0].section_single_configuration_calculation
        sec_band = scc[0].section_k_band
        result = 0
        for band_segment in sec_band[0].section_k_band_segment:
            freq = band_segment.band_energies
            result += np.count_nonzero(np.array(freq) < 0)
        return result

    def normalize(self):
        self.section = self.entry_archive.section_workflow.section_phonon
        if not self.section:
            self.section = self.entry_archive.section_workflow.m_create(Phonon)

        if not self.section.n_imaginary_frequencies:
            # get number from bands (not complete as this is not the whole mesh)
            self.section.n_imaginary_frequencies = self._get_n_imaginary_frequencies()


class ElasticNormalizer(Normalizer):
    def __init__(self, entry_archive):
        super().__init__(entry_archive)

    def _resolve_mechanical_stability(self):
        spacegroup = self.entry_archive.section_run[-1].section_system[-1].x_elastic_space_group_number
        order = self.entry_archive.section_run[-1].section_method[-1].x_elastic_elastic_constant_order

        if spacegroup is None or order != 2:
            return

        scc = self.entry_archive.section_run[-1].section_single_configuration_calculation[-1]
        c = scc.x_elastic_2nd_order_constants_matrix

        # see Phys. Rev B 90, 224104 (2014)
        res = False
        if spacegroup <= 2:  # Triclinic
            res = np.count_nonzero(c < 0)
        elif spacegroup <= 15:  # Monoclinic
            res = np.count_nonzero(c < 0)
        elif spacegroup <= 74:  # Orthorhombic
            res =\
                c[0][0] > 0 and c[0][0] * c[1][1] > c[0][1] ** 2 and\
                c[0][0] * c[1][1] * c[2][2] + 2 * c[0][1] * c[0][2] * c[1][2] -\
                c[0][0] * c[1][2] ** 2 - c[1][1] * c[0][2] ** 2 - c[2][2] * c[0][1] ** 2 > 0 and\
                c[3][3] > 0 and c[4][4] > 0 and c[5][5] > 0
        elif spacegroup <= 88:  # Tetragonal II
            res =\
                c[0][0] > abs(c[0][1]) and\
                2 * c[0][2] ** 2 < c[2][2] * (c[0][0] + c[0][1])
        elif spacegroup <= 142:  # Tetragonal I
            res =\
                c[0][0] > abs(c[0][1]) and\
                2 * c[0][2] ** 2 < c[2][2] * (c[0][0] + c[0][1]) and\
                c[3][3] > 0 and c[5][5] > 0
        elif spacegroup <= 148:  # rhombohedral II
            res =\
                c[0][0] > abs(c[0][1]) and c[3][3] > 0 and\
                c[0][2] ** 2 < (0.5 * c[2][2] * (c[0][0] + c[0][1])) and\
                c[0][3] ** 2 + c[0][4] ** 2 < 0.5 * c[3][3] * (c[0][0] - c[0][1])
        elif spacegroup <= 167:  # rhombohedral I
            res =\
                c[0][0] > abs(c[0][1]) and c[3][3] > 0 and\
                c[0][2] ** 2 < 0.5 * c[2][2] * (c[0][0] + c[0][1]) and\
                c[0][3] ** 2 < 0.5 * c[3][3] * (c[0][0] - c[0][1])
        elif spacegroup <= 194:  # hexagonal I
            res =\
                c[0][0] > abs(c[0][1]) and\
                2 * c[0][2] ** 2 < c[2][2] * (c[0][0] + c[0][1]) and\
                c[3][3] > 0 and c[5][5] > 0
        else:  # cubic
            res = c[0][0] - c[0][1] > 0 and c[0][0] + 2 * c[0][1] > 0 and c[3][3] > 0

        return res

    def _get_maximum_fit_error(self):
        scc = self.entry_archive.section_run[-1].section_single_configuration_calculation[-1]

        max_error = 0.0
        for diagram in scc.x_elastic_section_strain_diagrams:
            if diagram.x_elastic_strain_diagram_type == 'cross-validation':
                error = np.amax(diagram.x_elastic_strain_diagram_values)
                max_error = error if error > max_error else max_error

        return max_error

    def normalize(self):
        self.section = self.entry_archive.section_workflow.section_elastic
        if not self.section:
            self.section = self.entry_archive.section_workflow.m_create(Elastic)

        self.section.is_mechanically_stable = bool(self._resolve_mechanical_stability())

        self.section.fitting_error_maximum = self._get_maximum_fit_error()


class WorkflowNormalizer(Normalizer):
    '''
    This normalizer performs all produces a section all data necessary for the Optimade API.
    It assumes that the :class:`SystemNormalizer` was run before.
    '''
    def __init__(self, entry_archive):
        super().__init__(entry_archive)

    def normalize(self, logger=None) -> None:
        # Setup logger
        if logger is not None:
            self.logger = logger.bind(normalizer=self.__class__.__name__)

        # Do nothing if section_run is not present
        if self.section_run is None:
            return

        workflow_type = None
        if self.entry_archive.section_workflow:
            workflow_type = self.entry_archive.section_workflow.workflow_type

        if not workflow_type:
            sec_sampling_method = self.section_run.section_sampling_method
            if sec_sampling_method:
                # TODO imho geometry_optimization is not an appropriate name
                # if workflow_type == 'geometry_optimization':
                #     workflow_type = 'relaxation'
                workflow_type = sec_sampling_method[-1].sampling_method

        if not workflow_type:
            return

        workflow = self.entry_archive.section_workflow
        if not workflow:
            workflow = self.entry_archive.m_create(Workflow)
        workflow.workflow_type = workflow_type

        if workflow.workflow_type in ['geometry_optimization', 'relaxation']:
            RelaxationNormalizer(self.entry_archive).normalize()

        elif workflow.workflow_type == 'phonon':
            PhononNormalizer(self.entry_archive).normalize()

        elif workflow.workflow_type == 'elastic':
            ElasticNormalizer(self.entry_archive).normalize()
