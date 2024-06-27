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
import numpy as np
import yaml

from nomad.app.v1.routers.metainfo import (
    get_package_by_section_definition_id,
    store_package_definition,
)
from nomad.metainfo import MSection, MCategory, Quantity, SubSection
from nomad.metainfo.metainfo import Datetime, Package, MEnum, Reference, Definition

# resolve_references are tested in .test_references
# type specific serialization is tested in .test_quantities


class Category(MCategory):
    pass


class Abstract(MSection):
    scalar = Quantity(type=str, categories=[Category])
    many = Quantity(type=str, shape=['*'])
    matrix = Quantity(type=np.dtype(np.float64), shape=['3', '3'])


class Child(Abstract):
    pass


class Root(Abstract):
    quantity = Quantity()
    default = Quantity(type=str, default='test_value')
    derived = Quantity(type=str, derived=lambda *args, **kwargs: 'test_value')

    child = SubSection(sub_section=Child.m_def, categories=[Category])
    children = SubSection(sub_section=Child.m_def, repeats=True, categories=[Category])

    abstract = SubSection(sub_section=Abstract.m_def, categories=[Category])


values = dict(
    scalar='test_value',
    many=['test_value_1', 'test_value_2'],
    matrix=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
)

expected_child = dict(**values)
expected_root = dict(
    child=expected_child,
    children=[expected_child, expected_child],
    abstract=dict(m_def='tests.metainfo.test_to_dict.Child', **expected_child),
    **values,
)


@pytest.fixture
def example():
    root = Root(**values)
    root.m_create(Child, Root.child, **values)
    for _ in range(0, 2):
        root.m_create(Child, Root.children, **values)
    root.abstract = Child(**values)

    return root


def test_plain(example):
    assert example.m_to_dict() == expected_root


def test_m_from_dict(example):
    assert Root.m_from_dict(example.m_to_dict()).m_to_dict() == expected_root


@pytest.mark.parametrize(
    'metainfo_data',
    [
        pytest.param(
            {
                'm_def': 'nomad.metainfo.metainfo.Package',
                'name': 'test.Package',
                'section_definitions': [{'name': 'MySection'}],
            },
            id='python',
        )
    ],
)
def test_from_dict(metainfo_data, monkeypatch, mongo_module):
    assert (
        MSection.from_dict(metainfo_data).m_to_dict(
            with_root_def=True, with_out_meta=True
        )
        == metainfo_data
    )

    monkeypatch.setattr('nomad.config.process.add_definition_id_to_reference', True)

    metainfo_data['m_def'] += f'@{Package.m_def.definition_id}'

    package = MSection.from_dict(metainfo_data)

    assert package.m_to_dict(with_root_def=True, with_out_meta=True) == metainfo_data

    store_package_definition(package, with_root_def=True, with_out_meta=True)

    section_id = package.section_definitions[0].definition_id

    pkg_definition = get_package_by_section_definition_id(section_id)
    del pkg_definition['entry_id_based_name']
    assert pkg_definition == metainfo_data


def test_with_meta(example):
    assert example.m_to_dict(with_meta=True) == dict(
        m_def='tests.metainfo.test_to_dict.Root',
        child=dict(
            m_def='tests.metainfo.test_to_dict.Child',
            m_parent_sub_section='child',
            **expected_child,
        ),
        children=[
            dict(
                m_def='tests.metainfo.test_to_dict.Child',
                m_parent_sub_section='children',
                m_parent_index=0,
                **expected_child,
            ),
            dict(
                m_def='tests.metainfo.test_to_dict.Child',
                m_parent_sub_section='children',
                m_parent_index=1,
                **expected_child,
            ),
        ],
        abstract=dict(
            m_def='tests.metainfo.test_to_dict.Child',
            m_parent_sub_section='abstract',
            **expected_child,
        ),
        **values,
    )


def test_include_defaults(example):
    assert example.m_to_dict(include_defaults=True) == dict(
        default='test_value', **expected_root
    )


def test_derived(example):
    assert example.m_to_dict(include_derived=True) == dict(
        derived='test_value', **expected_root
    )


@pytest.mark.parametrize('include', [True, False])
def test_exclude_include(example, include: bool):
    def filter_function(prop, section):
        if isinstance(prop, Quantity) and section.m_def == Root.m_def:
            return not include

        if prop == Root.children or prop == Root.abstract:
            return not include

        return include

    if include:
        kwargs = dict(include=filter_function)
    else:
        kwargs = dict(exclude=filter_function)

    assert example.m_to_dict(**kwargs) == dict(child=expected_child)


def test_categories(example):
    root = dict(**expected_root)
    del root['many']
    del root['matrix']

    assert example.m_to_dict(categories=[Category]) == root


@pytest.mark.parametrize(
    'type, serialized_type',
    [
        pytest.param(str, dict(type_kind='python', type_data='str'), id='primitive'),
        pytest.param(
            Reference(Definition),
            dict(type_kind='reference', type_data='/section_definitions/0'),
            id='reference',
        ),
        pytest.param(
            MEnum('A', 'B', m_descriptions=dict(A='a', B='b')),
            dict(
                type_kind='enum',
                type_data=['A', 'B'],
                type_descriptions=dict(A='a', B='b'),
            ),
            id='enum',
        ),
        pytest.param(
            np.float64, dict(type_kind='numpy', type_data='float64'), id='numpy'
        ),
    ],
)
def test_quantity_type(type, serialized_type):
    class MySection(MSection):
        my_quantity = Quantity(type=type)

    assert MySection.m_def.m_to_dict()['quantities'][0]['type'] == serialized_type


def test_transform(example):
    def transform(quantity, section, value, path):
        if quantity == Abstract.scalar and section.m_def == Root.m_def:
            return 'other_value'

        return value

    root = dict(**expected_root)
    root.update(scalar='other_value')
    assert example.m_to_dict(transform=transform) == root


@pytest.fixture(scope='function')
def schema_yaml():
    return """
name: advanced_metainfo_example
section_definitions:
  - name: BaseSection
    quantities:
      - name: notes
        type: str
        description: |
          Additional notes about this section in Markdown
        links:
          - http://markdown.org
        format: Markdown

    inner_section_definitions:
      - name: 'User'
        quantities:
          - name: first_name
            type: str
          - name: last_name
            type: str
          - name: email
            type: str
            format: email
            optional: true

    sub_sections:
      - name: authors
        sub_section: User
        description: The user that authored this section
        repeats: true

  - name: ApplicationSection
    inner_section_definitions:
      - name: User
        base_sections:
          - BaseSection.User
        quantities:
          - name: user_id
            type: int
          - name: email
            type: str
            deprecated: 'Use user_id as replacement'
    sub_sections:
      - name: authors
        sub_section: User
      - name: data
        sub_section: ApplicationData
        repeats: true

  - name: ApplicationData
    quantities:
      - name: name
        type: str
      - name: value
        type: str
      - name: related
        type: ApplicationData
      - name: time
        type: Datetime
      - name: floaty
        type: np.float64
"""


def test_schema_deserialization(schema_yaml):
    schema_dict = yaml.load(schema_yaml, yaml.SafeLoader)
    pkg = Package.m_from_dict(schema_dict)
    pkg.init_metainfo()

    inner_section = pkg.section_definitions[1].inner_section_definitions[0]  # pylint: disable=all
    inner_section_base = inner_section.base_sections[0].m_resolved()
    base_inner_section = pkg.section_definitions[0].inner_section_definitions[0]
    application_data = pkg.all_definitions['ApplicationData']

    assert inner_section_base == base_inner_section
    assert (
        pkg.all_definitions['ApplicationSection']
        .all_properties['authors']
        .section_def.m_resolved()
        == inner_section
    )
    assert (
        pkg.all_definitions['ApplicationSection']
        .all_properties['data']
        .section_def.m_resolved()
        == pkg.all_definitions['ApplicationData']
    )
    assert application_data.all_properties['name'].type.standard_type() == 'str'
    assert (
        application_data.all_properties['related'].type.target_section_def.m_resolved()
        == application_data
    )
    assert application_data.all_properties['time'].type.standard_type() == 'datetime'
    assert application_data.all_properties['floaty'].type.standard_type() == 'float64'


def test_schema_definition_id(schema_yaml):
    """
    Test if the definition id is correctly generated.
    """
    schema_dict = yaml.load(schema_yaml, yaml.SafeLoader)
    pkg = Package.m_from_dict(schema_dict)
    pkg.init_metainfo()

    def check_dict(value):
        if not isinstance(value, (dict, list)):
            return value
        if isinstance(value, list):
            return [check_dict(v) for v in value]

        if 'm_def' in value:
            assert value['m_def_id'] is not None
            assert value['definition_id'] is not None

        return {
            k: check_dict(v)
            for k, v in value.items()
            if k not in ('m_def_id', 'definition_id')
        }

    check_dict(pkg.m_to_dict(with_def_id=True))
