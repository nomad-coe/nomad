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

from nomad.metainfo import Quantity, SubSection, Section

from nomad.datamodel.data import ArchiveSection, EntryData, WorkflowsElnCategory


class Link(ArchiveSection):
    '''
    Instances of Link are used to represent either a single input or single
    output of a Task. Using a separate section for links allows to put
    additional information (e.g. a name) on an input or output.
    '''

    name = Quantity(type=str, description=(
        'Name of the link. Will be used as a label for the input or output in workflow '
        'representations.'),
        a_eln=dict(component='StringEditQuantity'))
    section = Quantity(type=ArchiveSection, description=(
        'A reference to the section that contains the actual input or output data.'),
        a_eln=dict(component='ReferenceEditQuantity'))

    def normalize(self, archive, logger):
        super(Link, self).normalize(archive, logger)
        if not self.name and self.section:
            self.name = getattr(self.section, 'name', None)


class Task(ArchiveSection):
    '''
    Instances of Task are used to represent an activity that happened during workflow
    execution and that was acting on inputs to produce outputs.
    '''

    name = Quantity(type=str, description=(
        'A name of the task. Will be used as a label for the input or output in workflow '
        'representations.'),
        a_eln=dict(component='StringEditQuantity'))
    inputs = SubSection(sub_section=Link, repeats=True, description=(
        'All the links to sections that represent the inputs for this task.'))
    outputs = SubSection(sub_section=Link, repeats=True, description=(
        'All the links to sections that represent the outputs for this task.'))


class TaskReference(Task):
    '''
    A proxy section that can be used to compose a workflow of tasks that are contained
    in a different entry or workflow.
    '''

    task = Quantity(type=Task, description=(
        'A reference to the task that this section is a proxy for.'),
        a_eln=dict(component='ReferenceEditQuantity'))

    def normalize(self, archive, logger):
        super(TaskReference, self).normalize(archive, logger)
        if not self.name and self.task:
            self.name = self.task.name


class Workflow(Task, EntryData):
    '''
    Instances of Workflow are used to represent a set of Tasks that connect input and
    output data objects to produce a provenance graph for those data.

    Workflows themselves can be tasks. This allows to build nested workflows where some
    of the workflow tasks are workflows themselves.
    '''
    m_def = Section(categories=[WorkflowsElnCategory])

    tasks = SubSection(sub_section=Task, repeats=True, description=(
        'The tasks of this workflow as a repeating sub section. Use TaskReference if '
        'tasks cannot be contained.'))

    def normalize(self, archive, logger):
        super(Workflow, self).normalize(archive, logger)

        from nomad.datamodel import EntryArchive
        if isinstance(self.m_parent, EntryArchive):
            archive.workflow2 = self
