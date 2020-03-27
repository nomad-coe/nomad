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

from nomad.metainfo.metainfo import (
    MSection, MCategory, Section, Quantity, SubSection, Definition, Package, DeriveError,
    MetainfoError, Environment, MResource, Datetime, units, Annotation, SectionAnnotation,
    DefinitionAnnotation, Reference, MProxy, derived)
from nomad.metainfo.example import Run, VaspRun, System, SystemHash, Parsing, SCC, m_package as example_package


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
        assert len(Run.m_def.extending_sections) == 1
        assert len(Run.m_def.all_quantities) == 3
        assert Run.m_def.all_quantities['code_name'] in Run.m_def.quantities
        assert Run.m_def.all_quantities['code_name'] == Run.__dict__['code_name']

    def test_sub_sections(self):
        assert len(Run.m_def.sub_sections) == 3
        assert Run.m_def.all_sub_sections['systems'] in Run.m_def.sub_sections
        assert Run.m_def.all_sub_sections['systems'].sub_section == System.m_def
        assert len(Run.m_def.all_sub_sections_by_section[System.m_def]) == 1
        assert Run.m_def.all_sub_sections_by_section[System.m_def][0].sub_section == System.m_def

    def test_unset_sub_section(self):
        run = Run()
        assert run.systems == []
        assert run.parsing is None

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

        with pytest.raises(MetainfoError):
            class TestSection(TestBase):  # pylint: disable=unused-variable
                name = Quantity(type=int)

    def test_unique_names_extends(self):
        class TestBase(MSection):
            name = Quantity(type=str)

        with pytest.raises(MetainfoError):
            class TestSection(TestBase):  # pylint: disable=unused-variable
                m_def = Section(extends_base_section=True)
                name = Quantity(type=int)

    def test_multiple_sub_sections(self):
        class TestSection(MSection):  # pylint: disable=unused-variable
            one = SubSection(sub_section=System)
            two = SubSection(sub_section=System)

        assert len(TestSection.m_def.all_sub_sections_by_section[System.m_def]) == 2

    def test_dimension_exists(self):
        class TestSection(MSection):  # pylint: disable=unused-variable
            test = Quantity(type=str, shape=['does_not_exist'])

        assert len(TestSection.m_def.warnings) > 0

    def test_dimension_is_int(self):
        class TestSection(MSection):  # pylint: disable=unused-variable
            dim = Quantity(type=str)
            test = Quantity(type=str, shape=['dim'])

        assert len(TestSection.m_def.warnings) > 0

    def test_dimension_is_shapeless(self):
        class TestSection(MSection):  # pylint: disable=unused-variable
            dim = Quantity(type=int, shape=[1])
            test = Quantity(type=str, shape=['dim'])

        assert len(TestSection.m_def.warnings) > 0

    def test_higher_shapes_require_dtype(self):
        with pytest.raises(MetainfoError):
            class TestSection(MSection):  # pylint: disable=unused-variable
                test = Quantity(type=int, shape=[3, 3])

    def test_only_extends_one_base(self):
        with pytest.raises(MetainfoError):
            class TestSection(Run, System):  # pylint: disable=unused-variable
                m_def = Section(extends_base_section=True)

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

        with pytest.raises(AttributeError):
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
        with pytest.raises(TypeError):
            Run().code_name = 1

    def test_wrong_shape_1(self):
        with pytest.raises(TypeError):
            Run().code_name = ['name']

    def test_wrong_shape_2(self):
        with pytest.raises(TypeError):
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
        assert np.array_equal(system.unit_cell.magnitude, system.lattice_vectors.magnitude)

    @pytest.fixture(scope='function')
    def example_data(self):
        run = Run()
        run.code_name = 'test code name'
        run.m_create(Parsing)
        system: System = run.m_create(System)
        system.atom_labels = ['H', 'H', 'O']
        system.atom_positions = np.array([[1.2e-10, 0, 0], [0, 1.2e-10, 0], [0, 0, 1.2e-10]])
        system.atom_labels = np.array(['H', 'H', 'O'])
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

    def test_to_dict_category_filter(self, example_data: Run):
        system = example_data.systems[0]
        system.system_type = 'bulk'
        dct = system.m_to_dict(categories=[SystemHash])
        assert 'atom_labels' in dct
        assert 'n_atoms' not in dct  # derived
        assert 'system_type' not in dct  # not system hash

    def test_to_dict_defaults(self, example_data):
        dct = example_data.m_to_dict()
        assert 'nomad_version' not in dct['parsing']
        assert 'n_atoms' not in dct['systems'][0]

        dct = example_data.m_to_dict(include_defaults=True)
        assert 'nomad_version' in dct['parsing']
        assert 'n_atoms' not in dct['systems'][0]

    def test_derived(self):
        system = System()

        with pytest.raises(DeriveError):
            assert system.n_atoms == 3

        system.atom_labels = ['H', 'H', 'O']
        assert system.n_atoms == 3
        pass

    def test_derived_cached(self):
        class TestSection(MSection):
            value = Quantity(type=str)
            list = Quantity(type=str, shape=[1])

            @derived(cached=True)
            def derived(self):
                return self.value + self.list[0]

        assert TestSection.derived.cached
        test_section = TestSection(value='test', list=['1'])
        assert test_section.derived == 'test1'
        test_section.value = '2'
        assert test_section.derived == '21'
        test_section.list[0] = '2'
        assert test_section.derived == '21'

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

        errors, _ = run.m_validate()
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
        with pytest.raises(Exception):
            test_section.m_create(System)

        test_section.m_create(System, TestSection.one)

        assert test_section.one is not None

    def test_resource(self):
        resource = MResource()
        run = resource.create(Run)
        run.m_create(System)
        run.m_create(System)

        assert len(resource.all(System)) == 2

    def test_resource_move(self):
        resource = MResource()
        run = resource.create(Run)
        system = run.m_create(System)

        run = Run()
        run.m_add_sub_section(Run.systems, system)

    def test_mapping(self):
        run = Run()
        run.m_create(Parsing).parser_name = 'test'
        system = run.m_create(System)
        system.atom_labels = ['H', 'O']

        assert run.systems[0].atom_labels == ['H', 'O']
        assert run['systems.0.atom_labels'] == ['H', 'O']
        assert run['systems/0/atom_labels'] == ['H', 'O']
        assert run['parsing.parser_name'] == 'test'

    def test_np_dtype(self):
        scc = SCC()
        scc.energy_total_0 = 1.0
        scc.an_int = 1
        assert scc.energy_total_0.m == 1.0  # pylint: disable=no-member
        assert scc.energy_total_0 == 1.0 * units.J
        assert scc.m_to_dict()['energy_total_0'] == 1.0
        assert scc.an_int == 1
        assert scc.an_int.__class__ == np.int32
        assert scc.an_int.item() == 1  # pylint: disable=no-member

    def test_proxy(self):
        class OtherSection(MSection):
            name = Quantity(type=str)

        class ReferencingSection(MSection):
            proxy = Quantity(type=Reference(OtherSection.m_def))
            sub = SubSection(sub_section=OtherSection.m_def)

        obj = ReferencingSection()
        referenced = obj.m_create(OtherSection)
        referenced.name = 'test_value'
        obj.proxy = referenced

        assert obj.proxy == referenced
        assert obj.m_to_dict()['proxy'] == '/sub'
        assert obj.m_resolve('sub') == referenced
        assert obj.m_resolve('/sub') == referenced

        obj.proxy = MProxy('doesnotexist', m_proxy_section=obj, m_proxy_quantity=ReferencingSection.proxy)
        with pytest.raises(ReferenceError):
            obj.proxy.name

        obj.proxy = MProxy('sub', m_proxy_section=obj, m_proxy_quantity=ReferencingSection.proxy)
        assert obj.proxy.name == 'test_value'
        assert not isinstance(obj.proxy, MProxy)

        obj = ReferencingSection.m_from_dict(obj.m_to_dict(with_meta=True))
        assert obj.proxy.name == 'test_value'

    def test_copy(self):
        run = Run()
        run.m_create(Parsing).parser_name = 'test'
        system = run.m_create(System)
        system.atom_labels = ['H', 'O']

        copy = run.m_copy()
        assert copy is not run
        assert copy.m_def is run.m_def
        assert copy.systems is run.systems

        copy = run.m_copy(deep=True)
        assert copy is not run
        assert copy.systems is not run.systems
        assert copy.systems[0] is not run.systems[0]
        assert copy.systems[0].m_parent_index == 0
        assert copy.systems[0].m_parent_sub_section is run.systems[0].m_parent_sub_section

    def test_default_defaults(self):
        class TestSection(MSection):
            int_quantity = Quantity(type=int)
            float_quantity = Quantity(type=float)
            bool_quantity = Quantity(type=bool)

        section = TestSection()
        assert section.int_quantity == 0
        assert section.float_quantity == 0.0
        assert section.bool_quantity is False


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
