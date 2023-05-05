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

from nomad.normalizing.normalizer import Normalizer
from nomad.datamodel.metainfo.simulation.workflow import (
    SinglePoint, GeometryOptimization, MolecularDynamics, Phonon,
    Elastic, GW, GWResults
)
from nomad.datamodel import EntryArchive


class WorkflowNormalizer(Normalizer):
    '''
    This normalizer produces information specific to a workflow.
    '''

    def __init__(self, entry_archive: EntryArchive):
        super().__init__(entry_archive)
        self._elastic_programs = ['elastic']
        self._phonon_programs = ['phonopy']
        self._molecular_dynamics_programs = ['lammps']

    def _resolve_workflow(self):
        if not self.entry_archive.run:
            return

        # resolve it from parser
        workflow = None
        try:
            program_name = self.entry_archive.run[-1].program.name
        except Exception:
            return

        if program_name:
            program_name = program_name.lower()

        if program_name in self._elastic_programs:
            workflow = Elastic()

        elif program_name in self._molecular_dynamics_programs:
            workflow = MolecularDynamics()

        elif program_name in self._phonon_programs:
            workflow = Phonon()

        # resolve if from scc
        if workflow is None:
            # workflow references always to the last run
            # TODO decide if workflow should map to each run
            if len(self.entry_archive.run[-1].calculation) == 1:
                workflow = SinglePoint()
            else:
                workflow = GeometryOptimization()

        return workflow

    def gw_workflow_results(self, workflow: GW):
        """Populates the self.entry_archive.workflow2.results section in the GW workflow
        EntryArchive with references to the normalized calculations in each of these
        tasks: DFT SinglePoint and GW SinglePoint.
        """
        dft_task = workflow.tasks[0]
        gw_task = workflow.tasks[1]

        workflow_results = workflow.results
        if workflow_results is None:  # TODO ask @ladinesa if this is necessary
            workflow_results = workflow.m_create(GWResults)

        for name, section in workflow_results.m_def.all_quantities.items():
            calc_name = '_'.join(name.split('_')[:-1])
            if calc_name in ['dos', 'band_structure']:
                calc_name = f'{calc_name}_electronic'
            calc_section = []
            if 'dft' in name:
                calc_section = getattr(dft_task.outputs[-1].section, calc_name)
            elif 'gw' in name:
                calc_section = getattr(gw_task.outputs[-1].section, calc_name)
            if calc_section:
                workflow_results.m_set(section, calc_section)

    def normalize(self, logger=None) -> None:
        super().normalize(logger)

        # Do nothing if section_run is not present
        if not self.entry_archive.run:
            return

        if not self.entry_archive.workflow2:
            self.entry_archive.workflow2 = self._resolve_workflow()

        # Update workflow2.results with the normalized results for DFT and beyondDFT methods
        # (their normalizers are applied first always)
        workflow = self.entry_archive.workflow2
        workflow_name = workflow.m_def.name
        if workflow_name == 'GW':
            self.gw_workflow_results(workflow)
