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

from nomad.metainfo import (
    MSection, MCategory, Quantity, SubSection)

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


values = dict(
    scalar='test_value',
    many=['test_value_1', 'test_value_2'],
    matrix=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])

expected_child = dict(**values)
expected_root = dict(
    child=expected_child,
    children=[expected_child, expected_child],
    **values)


@pytest.fixture
def example():
    root = Root(**values)
    root.m_create(Child, Root.child, **values)
    for _ in range(0, 2):
        root.m_create(Child, Root.children, **values)

    return root


def test_plain(example):
    assert example.m_to_dict() == expected_root


def test_with_meta(example):
    assert example.m_to_dict(with_meta=True) == dict(
        m_def='Root',
        child=dict(m_def='Child', m_parent_sub_section='child', **expected_child),
        children=[
            dict(m_def='Child', m_parent_sub_section='children', m_parent_index=0, **expected_child),
            dict(m_def='Child', m_parent_sub_section='children', m_parent_index=1, **expected_child)],
        **values)


def test_include_defaults(example):
    assert example.m_to_dict(include_defaults=True) == dict(
        default='test_value', **expected_root)


def test_derived(example):
    assert example.m_to_dict(include_derived=True) == dict(
        derived='test_value', **expected_root)


@pytest.mark.parametrize('include', [True, False])
def test_exclude_include(example, include: bool):
    def filter_function(prop, section):
        if isinstance(prop, Quantity) and section.m_def == Root.m_def:
            return not include

        if prop == Root.children:
            return not include

        return include

    if include:
        kwargs = dict(include=filter_function)
    else:
        kwargs = dict(exclude=filter_function)

    assert example.m_to_dict(**kwargs) == dict(
        child=expected_child)


def test_categories(example):
    root = dict(**expected_root)
    del(root['many'])
    del(root['matrix'])

    assert example.m_to_dict(categories=[Category]) == root


def test_transform(example):
    def transform(quantity, section, value, path):
        if quantity == Abstract.scalar and section.m_def == Root.m_def:
            return 'other_value'

        return value

    root = dict(**expected_root)
    root.update(scalar='other_value')
    assert example.m_to_dict(transform=transform) == root
