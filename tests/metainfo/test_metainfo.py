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

# Contains more general test cases that are replaced continiously by more specialized
# in-depth tests in test_* files of the same module.

import pytest
import numpy as np
import pandas as pd
import pint.quantity

from nomad.metainfo.metainfo import (
    MSection, MCategory, Section, Quantity, SubSection, Definition, Package, DeriveError,
    MetainfoError, Environment, Annotation, AnnotationModel, SectionAnnotation, Context,
    DefinitionAnnotation, derived, MTypes)
from nomad.metainfo.example import Run, VaspRun, System, SystemHash, Parsing, SCC, m_package as example_package
from nomad import utils
from nomad.units import ureg

from tests import utils as test_utils


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
        assert run.systems == []  # pylint: disable=use-implicit-booleaness-not-comparison
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
        assert Quantity.m_def in Definition.m_def.all_inheriting_sections

    def test_unit(self):
        assert System.lattice_vectors.unit is not None

    @pytest.mark.parametrize('unit', [
        pytest.param('delta_degC / hr'),
        pytest.param('ΔdegC / hr'),
        pytest.param(ureg.delta_degC / ureg.hour),
    ])
    def test_unit_explicit_delta(self, unit):
        '''Explicit delta values are not allowed when setting or de-serializing.
        '''
        with pytest.raises(TypeError):
            Quantity(type=np.dtype(np.float64), unit=unit)
        with pytest.raises(TypeError):
            Quantity.m_from_dict({'m_def': 'nomad.metainfo.metainfo.Quantity', 'unit': str(unit)})

    @pytest.mark.parametrize('unit', [
        pytest.param('degC / hr'),
        pytest.param(ureg.degC / ureg.hour),
    ])
    def test_unit_implicit_delta(self, unit):
        '''Implicit delta values are allowed in setting and deserializing, delta
        prefixes are not serialized.
        '''
        quantity = Quantity(type=np.dtype(np.float64), unit=unit)
        serialized = quantity.m_to_dict()
        assert serialized['unit'] == 'degree_Celsius / hour'
        Quantity.m_from_dict(serialized)

    @pytest.mark.parametrize('dtype', [
        pytest.param(np.longlong),
        pytest.param(np.ulonglong)
    ])
    def test_unsupported_type(self, dtype):
        with pytest.raises(MetainfoError):
            Quantity(type=dtype)

    def test_extension(self):
        assert getattr(Run, 'x_vasp_raw_format', None) is not None
        assert 'x_vasp_raw_format' in Run.m_def.all_quantities

    def test_constraints(self):
        assert len(Run.m_def.constraints) > 0

    def test_unique_names(self):
        class TestBase(MSection):
            name = Quantity(type=str)

        # this is possible, can overwrite existing quantity
        class TestSection(TestBase):  # pylint: disable=unused-variable
            name = Quantity(type=int)

        pkg = Package(section_definitions=[TestBase.m_def, TestSection.m_def])
        pkg.init_metainfo()

    def test_unique_names_extends(self):
        class TestBase(MSection):
            name = Quantity(type=str)

        class TestSection(TestBase):  # pylint: disable=unused-variable
            m_def = Section(extends_base_section=True)
            name = Quantity(type=int)

        pkg = Package(section_definitions=[TestBase.m_def, TestSection.m_def])

        # this is not possible, cant replace existing quantity
        with pytest.raises(MetainfoError):
            pkg.init_metainfo()

    def test_alias(self):
        class SubTest(MSection):
            pass

        class Test(MSection):
            one = Quantity(type=str, aliases=['two', 'three'])
            section_one = SubSection(sub_section=SubTest, aliases=['section_two'])

        t = Test()
        t.one = 'value 1'
        assert t.one == 'value 1'
        assert t.two == 'value 1'
        assert t.three == 'value 1'

        t.two = 'value 2'
        assert t.one == 'value 2'
        assert t.two == 'value 2'
        assert t.three == 'value 2'

        sub_section = t.m_create(SubTest)
        assert t.section_one == sub_section
        assert t.section_two == sub_section

    def test_multiple_sub_sections(self):
        class TestSection(MSection):  # pylint: disable=unused-variable
            one = SubSection(sub_section=System)
            two = SubSection(sub_section=System)

        assert len(TestSection.m_def.all_sub_sections_by_section[System.m_def]) == 2

    def test_dimension_exists(self):
        class TestSection(MSection):  # pylint: disable=unused-variable
            test = Quantity(type=str, shape=['does_not_exist'])

        pkg = Package(section_definitions=[TestSection.m_def])
        pkg.init_metainfo()
        assert len(pkg.warnings) > 0

    def test_dimension_is_int(self):
        class TestSection(MSection):  # pylint: disable=unused-variable
            dim = Quantity(type=str)
            test = Quantity(type=str, shape=['dim'])

        pkg = Package(section_definitions=[TestSection.m_def])
        pkg.init_metainfo()
        assert len(pkg.warnings) > 0

    def test_dimension_is_shapeless(self):
        class TestSection(MSection):  # pylint: disable=unused-variable
            dim = Quantity(type=int, shape=[1])
            test = Quantity(type=str, shape=['dim'])

        pkg = Package(section_definitions=[TestSection.m_def])
        pkg.init_metainfo()
        assert len(pkg.warnings) > 0

    # TODO
    @pytest.mark.skip(reason=(
        'We disabled the constraint that is tested here, because some Nexus definitions '
        'are violating it.'))
    def test_higher_shapes_require_dtype(self):
        class TestSection(MSection):  # pylint: disable=unused-variable
            test = Quantity(type=int, shape=[3, 3])

        with pytest.raises(MetainfoError):
            pkg = Package(section_definitions=[TestSection.m_def])
            pkg.init_metainfo()

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

        class TestDefinitionAnnotation(DefinitionAnnotation):
            def init_annotation(self, definition):
                super().init_annotation(definition)
                assert definition.name in ['test_quantity', 'list_test_quantity', 'test_sub_section']
                assert definition.m_parent is not None
                self.initialized = True

        class TestSection(MSection):
            m_def = Section(a_test=TestSectionAnnotation())

            test_quantity = Quantity(type=str, a_test=TestDefinitionAnnotation())
            list_test_quantity = Quantity(
                type=str,
                a_test=[TestDefinitionAnnotation(), TestDefinitionAnnotation()])

            test_sub_section = SubSection(sub_section=System, a_test=TestDefinitionAnnotation())

        assert TestSection.m_def.a_test.initialized
        assert TestSection.m_def.m_get_annotations(TestSectionAnnotation).initialized
        assert TestSection().a_test == 'test annotation'

        assert TestSection.test_quantity.a_test is not None
        assert len(TestSection.list_test_quantity.m_get_annotations(TestDefinitionAnnotation)) == 2
        assert TestSection.test_sub_section.a_test is not None

    @pytest.mark.parametrize('annotation, passes', [
        pytest.param(dict(string='test_value'), True, id='passes-string'),
        pytest.param(dict(integer=1), True, id='passes-int'),
        pytest.param(dict(integer='string'), False, id='fails')
    ])
    def test_annotation_models(self, annotation, passes):
        class TestAnnotation(AnnotationModel):
            string: str = 'default'
            integer: int = 0
            no_default: str = None

        AnnotationModel.m_registry['test'] = TestAnnotation

        if passes:
            section_def = Section(name='test', a_test=annotation)
            assert isinstance(section_def.a_test, TestAnnotation)
            assert section_def.a_test.m_definition is not None
            as_dict = section_def.m_to_dict(with_meta=True)
            assert 'm_definition' not in as_dict
            section_def = Section.m_from_dict(as_dict)
            assert isinstance(section_def.a_test[0], TestAnnotation)
            assert section_def.a_test[0].m_definition is not None
        else:
            section_def = Section(name='test', a_test=annotation)
            assert isinstance(section_def.a_test, AnnotationModel)
            assert section_def.a_test.m_error is not None
            errors, _ = section_def.m_all_validate()
            assert len(errors) == 1

    def test_more_property(self):
        class TestSection(MSection):
            m_def = Section(this_does_not_exist_in_metainfo='value')
            test_quantity = Quantity(type=str, also_no_metainfo_quantity=1, one_more=False)
            test_delayed_more_quantity = Quantity(type=str)
            another_test_quantity = Quantity(type=str)

        TestSection.test_delayed_more_quantity.more = dict(one_more='test')

        assert TestSection.m_def.more['this_does_not_exist_in_metainfo'] == 'value'
        assert TestSection.test_quantity.more['also_no_metainfo_quantity'] == 1
        assert not TestSection.test_quantity.more['one_more']
        assert TestSection.test_delayed_more_quantity.more['one_more'] == 'test'
        assert len(TestSection.another_test_quantity.more) == 0

        assert TestSection.m_def.this_does_not_exist_in_metainfo == 'value'
        assert TestSection.test_quantity.also_no_metainfo_quantity == 1
        assert not TestSection.test_quantity.one_more
        with pytest.raises(AttributeError):
            assert TestSection.not_even_in_more is None

        serialized = TestSection.m_def.m_to_dict()
        assert 'more' in serialized
        assert 'this_does_not_exist_in_metainfo' in serialized['more']
        assert 'this_does_not_exist_in_metainfo' not in serialized


existing_repeating = Run()
existing_repeating.systems.append(System())
existing_nonrepeating = Run()
existing_nonrepeating.parsing = Parsing()
existing_multiple = Run()
existing_multiple.systems = [System(), System()]


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
        def run_test(section, quantity):
            assert getattr(section, quantity.name) is not None
            assert section.m_is_set(quantity)

            for _ in range(0, 2):
                setattr(section, quantity.name, None)
                assert not section.m_is_set(quantity)
                assert getattr(section, quantity.name) is None
                assert quantity.name not in section.m_to_dict()

        run = Run()
        run.code_name = 'test'
        run_test(run, Run.code_name)

        system = System()
        system.atom_positions = [[0, 0, 0]]
        run_test(system, System.atom_positions)

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

    def test_sub_section_lst(self):
        run = Run()
        assert run.systems == []  # pylint: disable=use-implicit-booleaness-not-comparison
        run.systems.append(System())

        assert len(run.systems) == 1
        assert run.systems[0].m_parent == run
        assert run.systems[0].m_parent_index == 0
        assert run.systems[0].m_parent_sub_section == Run.systems

        with pytest.raises(NotImplementedError):
            run.systems[0] = System()

        run.systems.append(System())
        first = run.systems[0]
        del(run.systems[0])
        assert first.m_parent is None
        assert run.systems[0].m_parent_index == 0

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

    def test_pd_dataframe(self):
        system = System()
        system.atom_positions = pd.DataFrame([[1, 2], [3, 4]])
        assert isinstance(system.atom_positions, pint.quantity._Quantity)
        assert np.all(system.atom_positions.m == [[1, 2], [3, 4]])

    def test_np_scalar(self):
        class TestSection(MSection):
            test_quantity = Quantity(type=np.dtype('int32'))

        test_section = TestSection()
        test_section.test_quantity = 12
        assert test_section.test_quantity == 12
        assert type(test_section.test_quantity) == np.int32

    def test_pd_dataframe_quantity(self):
        class TestSection(MSection):
            test_quantity = Quantity(type=np.dtype('int32'))

        test_section = TestSection()
        test_section.test_quantity = pd.DataFrame([[1, 2]])
        assert np.all(test_section.test_quantity == [1, 2])

    def test_unit_conversion(self):
        system = System()
        system.atom_positions = [[1, 2, 3]] * ureg.angstrom
        assert system.atom_positions.units == ureg.meter
        assert system.atom_positions[0][0] < 0.1 * ureg.meter

    @pytest.mark.parametrize('dtype', MTypes.num)
    @pytest.mark.parametrize('shape', [None, [1, 2]])
    def test_setting_with_dimensionless_unit(self, dtype, shape):
        if dtype not in MTypes.numpy:
            shape = None
        class TestSection(MSection):
            test_quantity = Quantity(type=dtype, shape=shape)

        test_section = TestSection()
        if dtype in MTypes.int:
            value = 42
        elif dtype in MTypes.float:
            value = 3.14
        elif dtype in MTypes.complex:
            value = 1+2j
        else:
            raise Exception('Unsupported type')
        if shape:
            value = np.full(shape, value, dtype=dtype)
        test_section.test_quantity = value * ureg.dimensionless
        if shape:
            assert np.all(test_section.test_quantity == value)
        elif dtype == np.float16:
            assert test_section.test_quantity == pytest.approx(value, 1e-2)
        elif dtype == np.float32:
            assert test_section.test_quantity == pytest.approx(value, 1e-7)
        else:
            assert test_section.test_quantity == value

    def test_synonym(self):
        system = System()
        system.lattice_vectors = [[1.2e-10, 0, 0], [0, 1.2e-10, 0], [0, 0, 1.2e-10]]
        assert isinstance(system.lattice_vectors, pint.quantity._Quantity)
        assert isinstance(system.unit_cell, pint.quantity._Quantity)
        assert np.array_equal(system.unit_cell.magnitude, system.lattice_vectors.magnitude)  # pylint: disable=no-member

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

    def test_derived(self):
        system = System()

        with pytest.raises(DeriveError):
            assert system.n_atoms == 3

        system.atom_labels = ['H', 'H', 'O']
        assert system.n_atoms == 3

        assert 'n_atoms' not in system.m_to_dict()
        assert 'n_atoms' in system.m_to_dict(include_derived=True)
        pass

    def test_derived_cached(self):
        class TestSection(MSection):
            value = Quantity(type=str)
            list = Quantity(type=str, shape=[1])

            @derived(cached=True)
            def derived(self):
                return self.value + self.list[0]

        assert TestSection.derived.cached  # pylint: disable=no-member
        test_section = TestSection(value='test', list=['1'])
        assert test_section.derived == 'test1'
        test_section.value = '2'
        assert test_section.derived == '21'
        test_section.list[0] = '2'
        assert test_section.derived == '21'

    def test_derived_deserialize(self):
        class TestSection(MSection):
            value = Quantity(type=str, derived=lambda _: 'test_value')

        test_section = TestSection()
        assert test_section.value == 'test_value'

        with pytest.raises(MetainfoError):
            test_section = TestSection(value='wrong')

        test_section = TestSection.m_from_dict({'test_value': 'wrong'})
        assert test_section.value == 'test_value'

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
        assert scc.energy_total_0 == 1.0 * ureg.J
        assert scc.m_to_dict()['energy_total_0'] == 1.0
        assert scc.an_int == 1
        assert scc.an_int.__class__ == np.int32
        assert scc.an_int.item() == 1  # pylint: disable=no-member

        scc.energy_total_0 = 1
        assert scc.energy_total_0.m == 1.0  # pylint: disable=no-member

        system = System()
        value = [[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]
        system.lattice_vectors = value
        assert isinstance(system.lattice_vectors.m, np.ndarray)  # pylint: disable=no-member
        assert system.lattice_vectors.m[1][2] == 1.0  # pylint: disable=no-member
        assert system.m_to_dict()['lattice_vectors'] == value

        class TestSection(MSection):
            f32_default = Quantity(
                type=np.dtype(np.float32),
                shape=[], default=1.0)
            f32 = Quantity(
                type=np.dtype(np.float32),
                shape=[])
            f64 = Quantity(
                type=np.dtype(np.float64),
                shape=[])

        section = TestSection()
        section.f32_default = -200
        section.f32 = -200
        section.f64 = -200

    def test_np_allow_wrong_shape(self, caplog):
        class MyContext(Context):
            def warning(self, event, **kwargs):
                utils.get_logger(__name__).warn(event, **kwargs)

        scc = SCC(m_context=MyContext())
        scc.energy_total_0 = np.array([1.0, 1.0, 1.0])
        scc.m_to_dict()
        test_utils.assert_log(caplog, 'WARNING', 'numpy quantity has wrong shape')

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

    def test_copy_keeps_m_sub_section_list(self):
        run = Run()
        run.m_create(Parsing).parser_name = 'test'
        system = run.m_create(System)
        system.atom_labels = ['H', 'O']

        copy = run.m_copy(deep=True)
        copy.systems.append(System())

        copy2 = copy.m_copy(deep=True)
        copy2.m_to_dict()

    def test_not_default_defaults(self):
        class TestSection(MSection):
            int_quantity = Quantity(type=int)
            float_quantity = Quantity(type=float)
            bool_quantity = Quantity(type=bool)

        section = TestSection()
        assert section.int_quantity is None
        assert section.float_quantity is None
        assert section.bool_quantity is None

    @pytest.mark.filterwarnings("ignore")
    def test_xpath(self):
        run = Run()
        run.code_name = 'amazingX'
        system = run.m_create(System)
        system.atom_labels = ['H', 'O']
        system.system_type = 'molecule'
        calc = run.m_create(SCC)
        calc.energy_total = -1.20E-23
        calc.system = system

        assert run.m_xpath('code_name') == 'amazingX'
        assert run.m_xpath('systems[-1].system_type') == 'molecule'
        assert run.m_xpath('sccs[0].system.atom_labels') == ['H', 'O']
        assert run.m_xpath('systems[?system_type == `molecule`].atom_labels') == [['H', 'O']]
        assert run.m_xpath('sccs[?energy_total < `1.0E-23`].system') == [{'atom_labels': ['H', 'O'], 'system_type': 'molecule'}]

    def test_m_update(self):
        class Child(MSection):
            pass

        class Parent(MSection):
            quantity = Quantity(type=str)
            single_sub_section = SubSection(sub_section=Child)
            many_sub_section = SubSection(sub_section=Child, repeats=True)

        parent = Parent()
        parent.m_update(
            quantity='Hello',
            single_sub_section=Child(),
            many_sub_section=[Child(), Child()])

        assert parent.quantity == 'Hello'
        assert parent.single_sub_section is not None
        assert len(parent.many_sub_section) == 2

    @pytest.mark.parametrize('root,path,exception', [
        pytest.param(Run(), 'parsing', None, id="non-existing non-repeating section"),
        pytest.param(Run(), 'systems', None, id="non-existing repeating section"),
        pytest.param(existing_nonrepeating, 'parsing', None, id="existing non-repeating section"),
        pytest.param(existing_repeating, 'systems', None, id="existing repeating section"),
        pytest.param(Run(), 'code_name', 'Could not find section definition for path "code_name"', id="cannot target quantity"),
        pytest.param(Run(), 'missing', 'Could not find section definition for path "missing"', id="invalid path"),
        pytest.param(existing_multiple, 'systems', 'Cannot resolve "systems" as several instances were found', id="ambiguous path"),
    ])
    def test_m_setdefault(self, root, path, exception):
        if not exception:
            system = root.m_setdefault(path)
            assert system
        else:
            with pytest.raises(Exception) as exc_info:
                system = root.m_setdefault(path)
            assert exception in str(exc_info.value)

    def test_m_traverse(self):
        expected = [
            ['quantity'],
            ['child', 'quantity'],
            ['child'],
            ['child_repeated', 0, 'quantity'],
            ['child_repeated'],
            ['child_repeated', 1, 'quantity'],
            ['child_repeated'],
        ]

        class Child(MSection):
            quantity = Quantity(type=int)

        class ChildRepeated(MSection):
            quantity = Quantity(type=int)

        class Parent(MSection):
            quantity = Quantity(type=int)
            child = SubSection(sub_section=Child)
            child_repeated = SubSection(sub_section=ChildRepeated, repeats=True)

        section = Parent(
            quantity=1,
            child=Child(quantity=2),
            child_repeated=[ChildRepeated(quantity=3), ChildRepeated(quantity=4)]
        )
        for i, [_, _, _, path] in enumerate(section.m_traverse()):
            assert path == expected[i]


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


@pytest.mark.parametrize('as_dict', [True, False])
@pytest.mark.parametrize('add_key', [True, False])
@pytest.mark.parametrize('str_type', [True, False])
@pytest.mark.parametrize('dup_key', [True, False])
def test_serialise_as_dict(as_dict, add_key, str_type, dup_key):
    class TestSection(MSection):
        q = Quantity(type=str if str_type else int)

    class TestContainer(MSection):
        s = SubSection(sub_section=TestSection, repeats=True, key_quantity='q' if add_key else None)

    container = TestContainer()

    def __key(_i):
        if str_type:
            return 'abc' if dup_key else f'abc{_i}'

        return 0 if dup_key else _i

    for i in range(3):
        section = TestSection(q=__key(i))
        container.m_add_sub_section(TestContainer.s, section)

    kwarg = {'subsection_as_dict': as_dict}
    if not str_type and add_key:
        with pytest.raises(TypeError):
            for i in range(3):
                _ = container.s[i].m_key
    else:
        for i in range(3):
            assert container.s[i].q == __key(i)
            if add_key and not dup_key:
                assert container.s[f'abc{i}'].q == f'abc{i}'
            else:
                with pytest.raises(KeyError):
                    _ = container.s[f'abc{i}'].q

        json_dict = container.m_to_dict(**kwarg)

        if not as_dict or add_key and dup_key:
            assert isinstance(json_dict['s'], list)
        else:
            assert isinstance(json_dict['s'], dict)
        assert json_dict == TestContainer.m_from_dict(json_dict).m_to_dict(**kwarg)
