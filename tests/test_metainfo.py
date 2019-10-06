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

from nomad.metainfo.metainfo import MSection, MCategory, Section, Quantity, SubSection, Definition, Package, DeriveError, MetainfoError, units
from nomad.metainfo.example import Run, VaspRun, System, SystemHash, Parsing, m_package as example_package

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
    """ Test for meta-info definition that are used to define other definitions. """

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
    """ Test for using meta-info instances without knowing/using the respective definitions. """

    def test_instantiation(self):
        test_section_def = Section(name='TestSection')
        test_section_def.m_create(Quantity, name='test_quantity')

        obj = MSection(m_def=test_section_def)
        assert obj.m_def.name == 'TestSection'
        # FIXME assert obj.m_get('test_quantity') is None
        setattr(obj, 'test_quantity', 'test_value')
        assert getattr(obj, 'test_quantity') == 'test_value'


class MaterialDefining(MCategory):
    """Quantities that add to what constitutes a different material."""
    pass


class TestM2:
    """ Test for meta-info definitions. """

    def test_basics(self):
        assert_section_def(Run.m_def)
        assert_section_def(System.m_def)

    def test_default_section_def(self):
        """ A section class without an explicit section def must set a default section def. """
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
        assert Run.m_def.all_sub_sections_by_section[System.m_def].sub_section == System.m_def

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

    def test_unique_sub_sections(self):
        with assert_exception(MetainfoError):
            class TestSection(MSection):  # pylint: disable=unused-variable
                one = SubSection(sub_section=System)
                two = SubSection(sub_section=System)

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


class TestM1:
    """ Test for meta-info instances. """

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

    def test_defaults(self):
        assert len(System().periodic_dimensions) == 3
        assert System().atom_labels is None
        try:
            System().does_not_exist
            assert False, 'Supposed unreachable'
        except AttributeError:
            pass
        else:
            assert False, 'Expected AttributeError'

    def test_m_section(self):
        assert Run().m_def == Run.m_def

    def test_children_parent(self):
        run = Run()
        system = run.m_create(System)

        assert run.systems[0] == system  # pylint: disable=E1101
        assert run.m_get_sub_section(Run.systems, 0) == system

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
            System().atom_labels = 'label'
            assert False, 'Supposed unreachable'
        except TypeError:
            pass
        else:
            assert False, 'Expected TypeError'

    def test_np(self):
        system = System()
        system.atom_positions = [[1, 2, 3]]
        assert isinstance(system.atom_positions, pint.quantity._Quantity)

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

    def test_derived(self):
        system = System()

        try:
            assert system.n_atoms == 3
            assert False, 'supposed unreachable'
        except DeriveError:
            pass
        else:
            assert False, 'supposed unreachable'

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
