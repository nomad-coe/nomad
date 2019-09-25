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
import numpy as np

from nomad.metainfo.metainfo import MObject, Section, Quantity, Definition, sub_section


def assert_section_def(section_def: Section):
    assert isinstance(section_def, Section)
    assert section_def.m_section is not None
    assert isinstance(section_def.m_section, Section)
    assert section_def.m_section.name is not None
    assert section_def.m_section.m_section == Section.m_section

    assert section_def.name is not None

    if section_def.parent is not None:
        if section_def.parent != section_def:
            assert_section_def(section_def.parent)

    if section_def.repeats:
        assert section_def.parent is not None


def assert_section_instance(section: MObject):
    assert_section_def(section.m_section)

    if section.m_parent is not None:
        assert section.m_parent.m_sub_section(section.m_section, section.m_parent_index) == section


class TestM3:
    """ Test for meta-info definition that are used to define other definitions. """

    def test_section(self):
        assert Section.m_section == Section.m_section.m_section
        assert Section.m_section.name == 'Section'
        assert Section.name is not None
        assert Section.name == Definition.name
        assert Section.name.m_section == Quantity.m_section
        assert Section.description.description is not None

        assert Section.m_section.m_sub_section(Quantity, 0).name in Section.m_section.attributes

        assert_section_instance(Section.m_section)

    def test_quantity(self):
        assert Quantity.m_section.m_section == Section.m_section
        assert Quantity.m_section.name == 'Quantity'
        assert Quantity.m_section.parent == Section.m_section

        assert_section_instance(Quantity.m_section)


class TestPureReflection:
    """ Test for using meta-info instances without knowing/using the respective definitions. """

    def test_instantiation(self):
        test_section_def = Section(name='TestSection')
        test_section_def.m_create(Quantity, name='test_quantity')

        obj = MObject(m_section=test_section_def)
        assert obj.m_section.name == 'TestSection'
        # FIXME assert obj.m_get('test_quantity') is None
        setattr(obj, 'test_quantity', 'test_value')
        assert getattr(obj, 'test_quantity') == 'test_value'


class Run(MObject):
    """ This is the description.

    And some more description.
    """

    code_name = Quantity(
        type=str, description='''
        The code_name description.
        ''')


class System(MObject):
    m_section = Section(repeats=True, parent=Run.m_section)
    n_atoms = Quantity(type=int, default=0)
    atom_label = Quantity(type=str, shape=['n_atoms'])
    atom_positions = Quantity(type=np.dtype('f8'), shape=['n_atoms', 3])


class Parsing(MObject):
    m_section = Section(parent=Run.m_section)


class TestM2:
    """ Test for meta-info definitions. """

    def test_basics(self):
        assert_section_def(Run.m_section)
        assert_section_def(System.m_section)

    def test_default_section_def(self):
        """ A section class without an explicit section def must set a default section def. """
        assert Run.m_section is not None
        assert Run.m_section.name == 'Run'
        assert not Run.m_section.repeats
        assert Run.m_section.parent is None

    def test_quantities(self):
        assert len(Run.m_section.quantities) == 1
        assert Run.m_section.quantities['code_name'] == Run.__dict__['code_name']

    def test_sub_sections(self):
        assert len(Run.m_section.sub_sections) == 2
        assert Run.m_section.sub_sections['System'] == System.m_section

    def test_attributes(self):
        assert len(Run.m_section.attributes) == 3
        assert Run.m_section.attributes['System'] == System.m_section
        assert Run.m_section.attributes['code_name'] == Run.__dict__['code_name']

    def test_get_quantity_def(self):
        assert System.n_atoms == System.m_section.attributes['n_atoms']

    def test_add_quantity(self):
        System.m_section.add_quantity(Quantity(name='test', type=str))

        system = System()
        system.test = 'test_value'

        assert 'test' in system.m_data
        assert system.test == 'test_value'
        assert getattr(System, 'test') == System.m_section.quantities['test']

    def test_section_name(self):
        assert Run.m_section.name == 'Run'

    def test_quantity_name(self):
        assert Run.code_name.name == 'code_name'

    def test_section_description(self):
        assert Run.m_section.description is not None
        assert Run.m_section.description.strip() == Run.m_section.description.strip()

    def test_quantity_description(self):
        assert Run.code_name.description is not None
        assert Run.code_name.description == 'The code_name description.'
        assert Run.code_name.description.strip() == Run.code_name.description.strip()


class TestM1:
    """ Test for meta-info instances. """

    def test_run(self):
        class Run(MObject):
            pass

        run = Run()

        assert run.m_section == Run.m_section
        assert run.m_section.name == 'Run'
        assert len(run.m_data) == 0

        assert_section_instance(run)

    def test_system(self):
        class System(MObject):
            m_section = Section()
            atom_labels = Quantity(type=str, shape=['1..*'])

        system = System()
        system.atom_labels = ['H']
        assert len(system.atom_labels) == 1
        assert len(system.m_data) == 1

        assert_section_instance(system)

    def test_defaults(self):
        assert System().n_atoms == 0
        assert System().atom_label is None
        try:
            System().does_not_exist
            assert False, 'Supposed unreachable'
        except AttributeError:
            pass
        else:
            assert False, 'Expected AttributeError'

    def test_m_section(self):
        assert Run().m_section == Run.m_section

    def test_children_parent(self):
        run = Run()
        system = run.m_create(System)

        assert run.system[0] == system  # pylint: disable=E1101
        assert run.m_sub_section(System, 0) == system

    def test_children_sub_section(self):
        setattr(Run, 'a_system_sub_section', sub_section(System))
        run = Run()
        system = run.m_create(System)

        assert run.a_system_sub_section[0] == system  # pylint: disable=E1101
        assert run.system[0] == system  # pylint: disable=E1101
        assert run.m_sub_section(System, 0) == system

    def test_parent_repeats(self):
        run = Run()
        system = run.m_create(System)

        assert system.m_parent == run
        assert system.m_parent_index == 0

    def test_parent_not_repeats(self):
        run = Run()
        parsing = run.m_create(Parsing)

        assert parsing.m_parent == run
        assert parsing.m_parent_index == -1

    def test_wrong_type(self):
        try:
            Run().code_name = 1
            assert False, 'Supposed unreachable'
        except TypeError:
            pass
        else:
            assert False, 'Expected TypeError'

    def test_wrong_shape_1(self):
        try:
            Run().code_name = ['name']
            assert False, 'Supposed unreachable'
        except TypeError:
            pass
        else:
            assert False, 'Expected TypeError'

    def test_wrong_shape_2(self):
        try:
            System().atom_label = 'label'
            assert False, 'Supposed unreachable'
        except TypeError:
            pass
        else:
            assert False, 'Expected TypeError'

    def test_np(self):
        system = System()
        system.atom_positions = [[1, 2, 3]]
        assert type(system.atom_positions) == np.ndarray
