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

from nomad.metainfo.metainfo import MObject, Section, Quantity


class TestM3:
    def test_section(self):
        assert Section.m_section == Section.m_section.m_section
        assert Section.m_section.name == 'Section'

    def test_quantity(self):
        assert Quantity.m_section.m_section == Section.m_section
        assert Quantity.m_section.name == 'Quantity'
        assert Quantity.m_section.parent == Section.m_section


class TestReflection:
    def test_m_get_success(self):
        assert Section.m_section.m_get('name') == 'Section'

    def test_m_get_fail(self):
        try:
            Section.m_section.m_get('doesnotexist')
            assert False
        except KeyError:
            pass

    def test_m_set_success(self):
        Section.m_section.m_set('name', 'newname')
        assert Section.m_section.m_get('name') == 'newname'


class TestPureReflection:

    def test_instantiation(self):
        test_section_def = Section(name='TestSection')
        test_section_def.m_create(Quantity, name='test_quantity')

        obj = MObject(m_section=test_section_def)
        assert obj.m_section.name == 'TestSection'
        # FIXME assert obj.m_get('test_quantity') is None
        obj.m_set('test_quantity', 'test_value')
        assert obj.m_get('test_quantity') == 'test_value'


class Run(MObject):
    m_section = Section()
    code_name = Quantity(type=str)


class System(MObject):
    m_section = Section(repeats=True, parent=Run.m_section)
    n_atoms = Quantity(type=int)


class TestM2:

    def test_quantities(self):
        assert len(Run.m_section.quantities) == 1
        assert Run.m_section.quantities['code_name'] == Run.__dict__['code_name']

    def test_sub_sections(self):
        assert len(Run.m_section.sub_sections) == 1
        assert Run.m_section.sub_sections['System'] == System.m_section

    def test_attributes(self):
        assert len(Run.m_section.attributes) == 2
        assert Run.m_section.attributes['System'] == System.m_section
        assert Run.m_section.attributes['code_name'] == Run.__dict__['code_name']


class TestM1:

    def test_run(self):
        class Run(MObject):
            m_section = Section()

        run = Run()

        assert run.m_section == Run.m_section
        assert run.m_section.name == 'Run'
        assert len(run.m_data) == 0

    def test_system(self):
        class System(MObject):
            m_section = Section()
            atom_labels = Quantity(type=str, shape=['1..*'])

        system = System()
        system.atom_labels = ['H']
        assert len(system.atom_labels) == 1
        assert len(system.m_data) == 1
