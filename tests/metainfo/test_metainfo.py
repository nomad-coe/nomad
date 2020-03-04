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

import pytest
import numpy as np
import pint.quantity
import datetime

from nomadcore.local_meta_info import InfoKindEl, InfoKindEnv

from nomad.metainfo.metainfo import (
    MSection, MCategory, Section, Quantity, SubSection, Definition, Package, DeriveError,
    MetainfoError, Environment, MResource, Datetime, units, Annotation, SectionAnnotation,
    DefinitionAnnotation)
from nomad.metainfo.example import Run, VaspRun, System, SystemHash, Parsing, m_package as example_package
from nomad.metainfo.legacy import LegacyMetainfoEnvironment
from nomad.parsing.metainfo import MetainfoBackend

from tests.utils import assert_exception


def assert_section_def(section_def: Section):
    assert isinstance(section_def, Section)
    assert section_def.m_def is not None
    assert isinstance(section_def.m_def, Section)
    assert section_def.m_def.name is not None
    assert section_def.m_def.m_def == Section.m_def

    assert section_def.name is not None


def assert_section_instance(section: MSection):
    assert_section_def(section.m_def)

    if section.m_parent is not None:
        assert section.m_parent.m_get_sub_section(section.m_parent_sub_section, section.m_parent_index) == section


class TestM3:
    ''' Test for meta-info definition that are used to define other definitions. '''

    def test_section(self):
        assert Section.m_def == Section.m_def.m_def
        assert Section.m_def.name == 'Section'
        assert Section.name is not None
        assert Section.name == Definition.name
        assert Section.name.m_def == Quantity.m_def
        assert Section.description.description is not None

        for quantity in iter(Section.m_def.quantities):
            assert quantity.name in Section.m_def.all_properties
            assert quantity.name in Section.m_def.all_quantities
            assert quantity.m_parent == Section.m_def

        for sub_section in iter(Section.m_def.sub_sections):
            assert sub_section.name in Section.m_def.all_properties
            assert sub_section.name in Section.m_def.all_sub_sections
            assert sub_section.sub_section in Section.m_def.all_sub_sections_by_section
            assert sub_section.m_parent == Section.m_def

        assert 'quantities' in Section.m_def.all_sub_sections
        assert 'sub_sections' in Section.m_def.all_sub_sections

        assert_section_instance(Section.m_def)

    def test_quantity(self):
        assert Quantity.m_def.m_def == Section.m_def
        assert Quantity.m_def.name == 'Quantity'

        assert_section_instance(Quantity.m_def)

    def test_definition(self):
        assert len(Section.m_def.base_sections) == 1
        assert len(Section.m_def.all_base_sections) == 1
        assert Section.m_def.m_follows(Definition.m_def)


class TestPureReflection:
    ''' Test for using meta-info instances without knowing/using the respective definitions. '''

    def test_instantiation(self):
        test_section_def = Section(name='TestSection')
        test_section_def.m_create(Quantity, name='test_quantity')

        obj = MSection(m_def=test_section_def)
        assert obj.m_def.name == 'TestSection'
        # FIXME assert obj.m_get('test_quantity') is None
        setattr(obj, 'test_quantity', 'test_value')
        assert getattr(obj, 'test_quantity') == 'test_value'


class MaterialDefining(MCategory):
    '''Quantities that add to what constitutes a different material.'''
    pass


class TestM2:
    ''' Test for meta-info definitions. '''

    def test_basics(self):
        assert_section_def(Run.m_def)
        assert_section_def(System.m_def)

    def test_default_section_def(self):
        ''' A section class without an explicit section def must set a default section def. '''
        assert Run.m_def is not None
        assert Run.m_def.name == 'Run'

    def test_quantities(self):
        assert len(Run.m_def.quantities) == 3
        assert Run.m_def.all_quantities['code_name'] in Run.m_def.quantities
        assert Run.m_def.all_quantities['code_name'] == Run.__dict__['code_name']

    def test_sub_sections(self):
        assert len(Run.m_def.sub_sections) == 3
        assert Run.m_def.all_sub_sections['systems'] in Run.m_def.sub_sections
        assert Run.m_def.all_sub_sections['systems'].sub_section == System.m_def
        assert len(Run.m_def.all_sub_sections_by_section[System.m_def]) == 1
        assert Run.m_def.all_sub_sections_by_section[System.m_def][0].sub_section == System.m_def

    def test_properties(self):
        assert len(Run.m_def.all_properties) == 6

    def test_get_quantity_def(self):
        assert System.n_atoms == System.m_def.all_properties['n_atoms']

    def test_section_name(self):
        assert Run.m_def.name == 'Run'

    def test_quantity_name(self):
        assert Run.code_name.name == 'code_name'

    def test_section_description(self):
        assert Run.m_def.description is not None
        assert Run.m_def.description.strip() == Run.m_def.description.strip()

    def test_quantity_description(self):
        assert Run.code_name.description is not None
        assert Run.code_name.description == 'The name of the code that was run.'
        assert Run.code_name.description.strip() == Run.code_name.description.strip()

    def test_direct_category(self):
        assert len(System.atom_labels.categories) == 1
        assert SystemHash.m_def in System.atom_labels.categories
        assert System.atom_labels in SystemHash.m_def.definitions

    def test_package(self):
        assert example_package.name == 'nomad.metainfo.example'
        assert example_package.description == 'An example metainfo package.'
        assert example_package.m_sub_section_count(Package.section_definitions) == 5
        assert example_package.m_sub_section_count(Package.category_definitions) == 1
        assert len(example_package.all_definitions) == 6

    def test_base_sections(self):
        assert Definition.m_def in iter(Section.m_def.base_sections)
        assert 'name' in Section.m_def.all_quantities
        assert 'name' in Quantity.m_def.all_quantities

    def test_unit(self):
        assert System.lattice_vectors.unit is not None

    def test_extension(self):
        assert getattr(Run, 'x_vasp_raw_format', None) is not None
        assert 'x_vasp_raw_format' in Run.m_def.all_quantities

    def test_constraints(self):
        assert len(Run.m_def.constraints) > 0

    def test_unique_names(self):
        class TestBase(MSection):
            name = Quantity(type=str)

        with assert_exception(MetainfoError):
            class TestSection(TestBase):  # pylint: disable=unused-variable
                name = Quantity(type=int)

    def test_unique_names_extends(self):
        class TestBase(MSection):
            name = Quantity(type=str)

        with assert_exception(MetainfoError):
            class TestSection(TestBase):  # pylint: disable=unused-variable
                m_def = Section(extends_base_section=True)
                name = Quantity(type=int)

    def test_multiple_sub_sections(self):
        class TestSection(MSection):  # pylint: disable=unused-variable
            one = SubSection(sub_section=System)
            two = SubSection(sub_section=System)

        assert len(TestSection.m_def.all_sub_sections_by_section[System.m_def]) == 2

    def test_dimension_exists(self):
        with assert_exception(MetainfoError):
            class TestSection(MSection):  # pylint: disable=unused-variable
                test = Quantity(type=str, shape=['does_not_exist'])

    def test_dimension_is_int(self):
        with assert_exception(MetainfoError):
            class TestSection(MSection):  # pylint: disable=unused-variable
                dim = Quantity(type=str)
                test = Quantity(type=str, shape=['dim'])

    def test_dimension_is_shapeless(self):
        with assert_exception(MetainfoError):
            class TestSection(MSection):  # pylint: disable=unused-variable
                dim = Quantity(type=int, shape=[1])
                test = Quantity(type=str, shape=['dim'])

    def test_higher_shapes_require_dtype(self):
        with assert_exception(MetainfoError):
            class TestSection(MSection):  # pylint: disable=unused-variable
                test = Quantity(type=int, shape=[3, 3])

    def test_only_extends_one_base(self):
        with assert_exception(MetainfoError):
            class TestSection(Run, System):  # pylint: disable=unused-variable
                m_def = Section(extends_base_section=True)

    def test_parent_section_sub_section_defs(self):
        len(System.m_def.parent_section_sub_section_defs) > 0

    def test_qualified_name(self):
        assert System.m_def.qualified_name() == 'nomad.metainfo.example.System'

    def test_derived_virtual(self):
        assert System.n_atoms.virtual

    def test_annotations(self):
        class TestSectionAnnotation(SectionAnnotation):
            def init_annotation(self, definition):
                super().init_annotation(definition)
                section_cls = definition.section_cls
                assert definition.name == 'TestSection'
                assert 'test_quantity' in definition.all_quantities
                assert section_cls.test_quantity.m_get_annotations('test').initialized
                assert section_cls.test_quantity.a_test.initialized
                assert section_cls.test_quantity.m_get_annotations('test', as_list=True)[0].initialized
                assert section_cls.test_quantity.m_get_annotations(Annotation).initialized
                assert all(a.initialized for a in section_cls.list_test_quantity.a_test)
                assert all(a.initialized for a in section_cls.list_test_quantity.m_get_annotations(Annotation))
                self.initialized = True

            def new(self, section):
                return dict(test='test annotation')

        class TestQuantityAnnotation(DefinitionAnnotation):
            def init_annotation(self, definition):
                super().init_annotation(definition)
                assert definition.name in ['test_quantity', 'list_test_quantity']
                assert definition.m_parent is not None
                self.initialized = True

        class TestSection(MSection):
            m_def = Section(a_test=TestSectionAnnotation())

            test_quantity = Quantity(type=str, a_test=TestQuantityAnnotation())
            list_test_quantity = Quantity(
                type=str,
                a_test=[TestQuantityAnnotation(), TestQuantityAnnotation()])

        assert TestSection.m_def.a_test.initialized
        assert TestSection.m_def.m_get_annotations(TestSectionAnnotation).initialized

        assert TestSection().a_test == 'test annotation'


class TestM1:
    ''' Test for meta-info instances. '''

    def test_run(self):
        class Run(MSection):
            pass

        run = Run()

        assert run.m_def == Run.m_def
        assert run.m_def.name == 'Run'

        assert_section_instance(run)

    def test_system(self):
        class System(MSection):
            m_def = Section()
            atom_labels = Quantity(type=str, shape=['1..*'])

        system = System()
        system.atom_labels = ['H']
        assert len(system.atom_labels) == 1

        assert_section_instance(system)

    def test_set_none(self):
        run = Run()
        run.code_name = 'test'
        assert run.code_name is not None

        run.code_name = None
        assert run.code_name is None

    def test_set_subsection(self):
        run = Run()
        first = Parsing()
        run.parsing = first
        assert first.m_parent == run
        assert run.parsing == first

        second = Parsing()
        run.parsing = second
        assert first.m_parent is None
        assert second.m_parent == run
        assert run.parsing == second

        run.parsing = None
        assert run.parsing is None

    def test_defaults(self):
        assert len(System().periodic_dimensions) == 3
        assert System().atom_labels is None

        with assert_exception(AttributeError):
            getattr(System(), 'does_not_exist')

    def test_m_section(self):
        assert Run().m_def == Run.m_def

    def test_children_parent(self):
        run = Run()
        system1 = run.m_create(System)
        run.m_create(System)

        assert run.systems[0] == system1  # pylint: disable=E1101
        assert run.m_get_sub_section(Run.systems, 0) == system1
        assert run.m_sub_section_count(Run.systems) == 2

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
        with assert_exception(TypeError):
            Run().code_name = 1

    def test_wrong_shape_1(self):
        with assert_exception(TypeError):
            Run().code_name = ['name']

    def test_wrong_shape_2(self):
        with assert_exception(TypeError):
            System().atom_labels = 'label'

    def test_np_array(self):
        system = System()
        system.atom_positions = [[1, 2, 3]]
        assert isinstance(system.atom_positions, pint.quantity._Quantity)

    def test_np_scalar(self):
        class TestSection(MSection):
            test_quantity = Quantity(type=np.dtype('int16'))

        test_section = TestSection()
        test_section.test_quantity = 12
        assert test_section.test_quantity == 12
        assert type(test_section.test_quantity) == np.int16

    def test_unit_conversion(self):
        system = System()
        system.atom_positions = [[1, 2, 3]] * units.angstrom
        assert system.atom_positions.units == units.meter
        assert system.atom_positions[0][0] < 0.1 * units.meter

    def test_synonym(self):
        system = System()
        system.lattice_vectors = [[1.2e-10, 0, 0], [0, 1.2e-10, 0], [0, 0, 1.2e-10]]
        assert isinstance(system.lattice_vectors, pint.quantity._Quantity)
        assert isinstance(system.unit_cell, pint.quantity._Quantity)
        assert np.array_equal(system.unit_cell, system.lattice_vectors)

    @pytest.fixture(scope='function')
    def example_data(self):
        run = Run()
        run.code_name = 'test code name'
        run.m_create(Parsing)
        system: System = run.m_create(System)
        system.atom_labels = ['H', 'H', 'O']
        system.atom_positions = np.array([[1.2e-10, 0, 0], [0, 1.2e-10, 0], [0, 0, 1.2e-10]])

        return run

    def assert_example_data(self, data: Run):
        assert_section_instance(data)
        assert data.m_def == Run.m_def
        assert data.code_name == 'test code name'
        system: System = data.systems[0]
        assert_section_instance(system)
        assert system.m_def == System.m_def
        assert system.n_atoms == 3
        assert system.atom_labels == ['H', 'H', 'O']
        assert isinstance(system.atom_positions, pint.quantity._Quantity)

    def test_to_dict(self, example_data):
        dct = example_data.m_to_dict()
        new_example_data = Run.m_from_dict(dct)

        self.assert_example_data(new_example_data)

    def test_to_dict_defaults(self, example_data):
        dct = example_data.m_to_dict()
        assert 'nomad_version' not in dct['parsing']
        assert 'n_atoms' not in dct['systems'][0]

        dct = example_data.m_to_dict(include_defaults=True)
        assert 'nomad_version' in dct['parsing']
        assert 'n_atoms' not in dct['systems'][0]

    def test_derived(self):
        system = System()

        with assert_exception(DeriveError):
            assert system.n_atoms == 3

        system.atom_labels = ['H', 'H', 'O']
        assert system.n_atoms == 3
        pass

    def test_extension(self):
        run = Run()
        run.m_as(VaspRun).x_vasp_raw_format = 'outcar'
        assert run.m_as(VaspRun).x_vasp_raw_format == 'outcar'

    def test_resolve(self):
        run = Run()
        system = run.m_create(System)

        assert run.m_resolve('/systems/0') == system
        assert system.m_resolve('/systems/0') == system

    def test_validate(self):
        run = Run()
        run.m_create(System)

        errors = run.m_validate()
        assert len(errors) == 1

    def test_validate_dimension(self):
        system = System()
        system.atom_labels = ['H']
        system.atom_positions = []
        assert len(system.m_validate()) > 0

    def test_multiple_sub_sections(self):
        class TestSection(MSection):  # pylint: disable=unused-variable
            one = SubSection(sub_section=System)
            two = SubSection(sub_section=System)

        test_section = TestSection()
        with assert_exception():
            test_section.m_create(System)

        test_section.m_create(System, TestSection.one)

        assert test_section.one is not None

    def test_resource(self):
        resource = MResource()
        run = resource.create(Run)
        run.m_create(System)
        run.m_create(System)

        assert len(resource.all(System)) == 2

    def test_mapping(self):
        run = Run()
        run.m_create(Parsing).parser_name = 'test'
        system = run.m_create(System)
        system.atom_labels = ['H', 'O']

        assert run.systems[0].atom_labels == ['H', 'O']
        assert run['systems.0.atom_labels'] == ['H', 'O']
        assert run['systems/0/atom_labels'] == ['H', 'O']
        assert run['parsing.parser_name'] == 'test'


class TestDatatypes:

    def test_datetime(self):
        class TestSection(MSection):
            datetime = Quantity(type=Datetime)

        obj = TestSection()
        assert obj.datetime is None
        assert 'datetime' not in obj.m_to_dict()

        obj.datetime = datetime.datetime.now()
        assert obj.datetime is not None
        assert isinstance(obj.m_to_dict()['datetime'], str)

        obj.datetime = obj.datetime.isoformat()
        assert obj.datetime is not None
        assert isinstance(obj.m_to_dict()['datetime'], str)

        obj.datetime = None
        assert obj.datetime is None
        assert obj.m_to_dict()['datetime'] is None


class TestEnvironment:

    @pytest.fixture
    def env(self) -> Environment:
        env = Environment()
        env.m_add_sub_section(Environment.packages, example_package)
        return env

    def test_create(self, env):
        assert env is not None

    def test_resolve(self, env: Environment):
        sub_section_system = env.resolve_definition('systems', SubSection)
        assert sub_section_system.m_def == SubSection.m_def
        assert sub_section_system.name == 'systems'


class TestLegacy:

    @pytest.fixture(scope='session')
    def legacy_example(self):
        return {
            'type': 'nomad_meta_info_1_0',
            'metaInfos': [
                dict(name='section_run', kindStr='type_section'),
                dict(
                    name='program_name', dtypeStr='C', superNames=['program_info']),
                dict(name='bool_test', dtypeStr='b', superNames=['section_run']),

                dict(name='section_system', kindStr='type_section', superNames=['section_run']),
                dict(
                    name='atom_labels', dtypeStr='C', shape=['number_of_atoms'],
                    superNames=['section_system']),
                dict(
                    name='atom_positions', dtypeStr='f', shape=['number_of_atoms', 3],
                    superNames=['section_system']),
                dict(
                    name='number_of_atoms', dtypeStr='i', kindStr='type_dimension',
                    superNames=['section_system', 'system_info']),

                dict(name='section_method', kindStr='type_section', superNames=['section_run']),
                dict(
                    name='method_system_ref', dtypeStr='i',
                    referencedSections=['section_system'],
                    superNames=['section_method']),

                dict(
                    name='program_info', kindStr='type_abstract_document_content',
                    superNames=['section_run']),
                dict(name='system_info', kindStr='type_abstract_document_content')
            ]
        }

    @pytest.fixture(scope='session')
    def legacy_env(self, legacy_example):
        env = InfoKindEnv()
        for definition in legacy_example.get('metaInfos'):
            env.addInfoKindEl(InfoKindEl(
                description='test_description', package='test_package', **definition))
        return env

    @pytest.fixture(scope='session')
    def env(self, legacy_env):
        return LegacyMetainfoEnvironment(legacy_env)

    def test_environment(self, env, no_warn):
        assert env.all_defs.get('section_system') is not None
        assert env.all_defs.get('section_run') is not None
        assert env.all_defs.get('section_method') is not None
        assert env.all_defs.get('method_system_ref').type.target_section_def == \
            env.all_defs.get('section_system')
        assert env.all_defs.get('atom_labels').type == str
        assert env.all_defs.get('atom_positions').type == np.dtype('f')

        assert env.env.resolve_definition('section_system', Section) == \
            env.all_defs.get('section_system', '<does not exist>')
        assert env.env.resolve_definition('section_system', SubSection) is not None

        assert env.env.resolve_definition('bool_test', Quantity).type == bool
        assert env.env.resolve_definition('program_name', Quantity).m_parent.name == 'section_run'

    def test_backend(self, env, no_warn):
        backend = MetainfoBackend(env)
        run = backend.openSection('section_run')
        backend.addValue('program_name', 'vasp')

        system_0 = backend.openSection('section_system')
        assert system_0 == 0

        backend.addArrayValues('atom_labels', np.array(['H', 'H', 'O']))
        backend.addValue('number_of_atoms', 3)
        backend.addArrayValues('atom_positions', [[0, 0, 0], [0, 0, 0], [0, 0, 0]])
        backend.closeSection('section_system', system_0)

        system_1 = backend.openSection('section_system')
        assert system_1 == 1
        backend.closeSection('section_system', system_1)

        method = backend.openSection('section_method')
        backend.addValue('method_system_ref', system_0)
        backend.closeSection('section_method', method)

        backend.closeSection('section_run', run)

        assert backend.get_sections('section_system') == [0, 1]
        assert len(backend.get_value('atom_labels', 0)) == 3
        assert backend.get_value('method_system_ref', 0) == 0
        assert backend.get_value('program_name', 0) == 'vasp'

        backend.openContext('section_run/0/section_system/1')
        backend.addValue('number_of_atoms', 10)
        backend.closeContext('section_run/0/section_system/1')
        assert backend.get_value('number_of_atoms', 1) == 10
