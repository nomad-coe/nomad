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
import random

from nomad.utils import get_logger
from nomad.datamodel import EntryArchive, EntryMetadata
from nomad.datamodel.metainfo.workflow import Task
from nomad.datamodel.metainfo.simulation import Run, Calculation
from nomad.datamodel.metainfo.simulation.workflow import (
    SimulationWorkflow, SimulationWorkflowMethod, SimulationWorkflowResults)


class TestSimulationWorkflow:
    '''
    Tests for the base simulation workflow class.
    '''

    n_calc = 10
    logger = get_logger(__name__)

    @pytest.fixture(autouse=True)
    def serial_simulation(self) -> EntryArchive:
        '''
        Simulation with calculations done in serial.
        '''
        archive = EntryArchive()
        archive.metadata = EntryMetadata(entry_type='Workflow')
        archive.workflow2 = SimulationWorkflow(
            method=SimulationWorkflowMethod(), results=SimulationWorkflowResults())

        archive.run.append(Run(calculation=[
            Calculation(time_physical=t, time_calculation=1) for t in range(1, self.n_calc + 1)]))
        return archive

    def test_tasks_serial(self, serial_simulation):
        '''
        Test tasks creation of purely serial calculation.
        '''
        workflow = serial_simulation.workflow2
        workflow.normalize(serial_simulation, self.logger)
        assert len(workflow.tasks) == self.n_calc

        assert workflow.inputs[0].section == workflow.tasks[0].inputs[0].section
        for n, task in enumerate(workflow.tasks[:-1]):
            assert task.name == f'Step {n + 1}'
            assert len(task.inputs) == 1
            assert len(task.outputs) == 1
            assert task.outputs[0].section == workflow.tasks[n + 1].inputs[0].section
        assert workflow.outputs[0].section == workflow.tasks[-1].outputs[0].section

    def test_tasks_defined(self, serial_simulation):
        '''
        Test tasks creation skipped if tasks are predefined
        '''
        workflow = SimulationWorkflow(tasks=[Task(name='1')])
        serial_simulation.workflow2 = workflow
        workflow.normalize(serial_simulation, self.logger)

        assert len(serial_simulation.workflow2.tasks) == 1
        assert serial_simulation.workflow2.tasks[0].name == '1'

    def test_tasks_no_time(self, serial_simulation):
        '''
        Test tasks creation skipped if at least one calculation has no time info.
        '''

        for key in ['time_physical', 'time_calculation']:
            calc = serial_simulation.run[0].calculation[random.randint(0, self.n_calc - 1)]
            calc.m_set(calc.m_get_quantity_definition(key), None)
            serial_simulation.workflow2.normalize(serial_simulation, self.logger)
            assert not serial_simulation.workflow2.tasks

    @pytest.mark.parametrize('calculation_indices', [
        # parallel (0 to 3), 4, 5, parallel (6 to 9)
        [[0, 1, 2, 3], [4], [5], [6, 7, 8, 9]],
        # 0, parallel (1 to 2), 4, 5, 6, 7, 8, 8
        [[0], [1, 2], [3], [4], [5], [6], [7], [8], [9]],
        # parallel (0 to 8), 9
        [[0, 1, 2, 3, 4, 5, 6, 7, 8], [9]]
    ])
    def test_task_not_serial(self, serial_simulation, calculation_indices):
        '''
        Test creation for mixed serial and parallel tasks.
        '''
        def _create_times(indices, start_time=0):
            times = []
            for n in indices:
                if not isinstance(n, int):
                    times.extend(sorted(_create_times(n, start_time=times[-1][1] if times else 0), key=lambda x: x[1]))
                else:
                    calc_time = random.random()
                    dt = random.random() * 0.1  # small perturbation
                    times.append([n, calc_time + start_time + dt, calc_time])
            return times

        for n, time_physical, time_calculation in _create_times(calculation_indices):
            serial_simulation.run[-1].calculation[n].time_physical = time_physical
            serial_simulation.run[-1].calculation[n].time_calculation = time_calculation

        workflow = serial_simulation.workflow2
        workflow.normalize(serial_simulation, self.logger)
        assert len(workflow.tasks) == 10

        # workflow inputs as inputs to first parallel tasks
        for n in calculation_indices[0]:
            assert workflow.tasks[n].name == 'Step 1'
            assert workflow.tasks[n].inputs[0].section == workflow.inputs[0].section

        # outputs of previous tasks are inputs of succeeding tasks in series
        for i in range(1, len(calculation_indices)):
            for n1 in calculation_indices[i]:
                assert workflow.tasks[n1].name == f'Step {i + 1}'
                inputs = [input.section for input in workflow.tasks[n1].inputs]
                for n0 in calculation_indices[i - 1]:
                    assert workflow.tasks[n0].outputs[-1].section in inputs

        # last parallel tasks oututs as workflow outputs
        for n in calculation_indices[-1]:
            assert workflow.tasks[n].outputs[0].section in [output.section for output in workflow.outputs]
