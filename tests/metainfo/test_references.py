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

from typing import cast
import pytest
import os.path

from nomad.datamodel import UserReference, AuthorReference
from nomad.metainfo import (
    MSection,
    Quantity,
    Section,
    SubSection,
    MProxy,
    Reference,
    QuantityReference,
    File,
    MetainfoReferenceError,
    Package as MetainfoPackage,
    Context,
    Package,
)


class Referenced(MSection):
    str_quantity = Quantity(type=str)


class Referencing(MSection):
    section_reference = Quantity(type=Referenced)
    section_reference_list = Quantity(type=Referenced, shape=['*'])

    quantity_reference = Quantity(type=Referenced.str_quantity)
    file_reference = Quantity(type=File)


class Root(MSection):
    referenced = SubSection(sub_section=Referenced)
    referenceds = SubSection(sub_section=Referenced, repeats=True)
    referencing = SubSection(sub_section=Referencing)


@pytest.fixture(scope='function', params=['reference', 'proxy'])
def definitions(request):
    """Alters the type of the references used in the reference properties definitions."""
    reference_type = request.param
    if reference_type == 'reference':
        Referencing.section_reference.type = Reference(Referenced.m_def)
        Referencing.section_reference_list.type = Reference(Referenced.m_def)
        Referencing.quantity_reference.type = QuantityReference(Referenced.str_quantity)

    elif reference_type == 'proxy':
        Referencing.section_reference.type = Reference(Referenced.m_def)
        Referencing.section_reference_list.type = Reference(Referenced.m_def)
        Referencing.quantity_reference.type = QuantityReference(Referenced.str_quantity)

    else:
        raise NotImplementedError()


@pytest.fixture(scope='function')
def example_data(definitions):
    def create_referenced():
        referenced = Referenced()
        referenced.str_quantity = 'test_value'
        return referenced

    referenced = create_referenced()
    referenced_1 = create_referenced()
    referenced_2 = create_referenced()

    root = Root()
    root.referenced = referenced
    root.referenceds.append(referenced_1)
    root.referenceds.append(referenced_2)

    referencing = Referencing()
    referencing.section_reference = referenced
    referencing.section_reference_list = [referenced_1, referenced_2]
    referencing.quantity_reference = referenced
    root.referencing = referencing

    return root


def assert_data(example_data):
    def assert_properties(example_data):
        assert (
            example_data.referencing.section_reference.m_resolved()
            == example_data.referenced
        )
        assert (
            example_data.referencing.m_to_dict()['section_reference'] == '/referenced'
        )
        assert (
            example_data.referencing.section_reference_list[1].m_resolved()
            == example_data.referenceds[1]
        )
        assert example_data.referencing.m_to_dict()['section_reference_list'] == [
            '/referenceds/0',
            '/referenceds/1',
        ]
        assert example_data.referencing.quantity_reference == 'test_value'
        assert (
            example_data.referencing.m_to_dict()['quantity_reference']
            == '/referenced/str_quantity'
        )

        assert example_data.referencing.m_is_set(Referencing.section_reference)
        assert example_data.referencing.m_is_set(Referencing.section_reference_list)
        assert example_data.referencing.m_is_set(Referencing.quantity_reference)

    assert_properties(example_data)

    example_data_serialized = example_data.m_to_dict(with_meta=True)
    example_data = Root.m_from_dict(example_data_serialized)
    assert_properties(example_data)


def test_references(example_data):
    assert_data(example_data)


def test_section_proxy(example_data):
    example_data.referencing.section_reference = MProxy(
        'doesnotexist',
        m_proxy_section=example_data.referencing,
        m_proxy_type=Referencing.section_reference.type,
    )
    with pytest.raises(MetainfoReferenceError):
        example_data.referencing.section_reference.str_quantity

    example_data.referencing.section_reference = MProxy(
        '/referenced',
        m_proxy_section=example_data.referencing,
        m_proxy_type=Referencing.section_reference.type,
    )

    assert_data(example_data)


def test_quantity_proxy(example_data):
    with pytest.raises(MetainfoReferenceError):
        example_data.referencing.quantity_reference = MProxy(
            'doesnotexist',
            m_proxy_section=example_data.referencing,
            m_proxy_type=Referencing.section_reference.type,
        )

    example_data.referencing.quantity_reference = MProxy(
        '/referenced',
        m_proxy_section=example_data.referencing,
        m_proxy_type=Referencing.section_reference.type,
    )
    assert example_data.referencing.quantity_reference == 'test_value'

    assert_data(example_data)


def test_resolve_references(example_data):
    assert example_data.m_to_dict(resolve_references=True) == {
        'referenced': {'str_quantity': 'test_value'},
        'referenceds': [{'str_quantity': 'test_value'}, {'str_quantity': 'test_value'}],
        'referencing': {
            'section_reference': {'str_quantity': 'test_value'},
            'section_reference_list': [
                {'str_quantity': 'test_value'},
                {'str_quantity': 'test_value'},
            ],
            'quantity_reference': 'test_value',
        },
    }


def test_quantity_references_serialize():
    source = {
        'referenced': {'str_quantity': 'test_value'},
        'referencing': {'quantity_reference': '/referenced/str_quantity'},
    }
    root = Root.m_from_dict(source)
    assert source == root.m_to_dict()


def test_quantity_references_unset_serialize():
    referencing = Referencing()
    referencing.quantity_reference = Referenced()

    assert not referencing.m_is_set(Referencing.quantity_reference)
    assert 'quantity_reference' not in referencing.m_to_dict()


def test_section_reference_serialize():
    class TargetSection(MSection):
        test_quantity = Quantity(type=str)

    class SourceSection(MSection):
        section_ref = Quantity(type=Reference(TargetSection.m_def))
        quantity_ref = Quantity(type=QuantityReference(TargetSection.test_quantity))

    pkg = MetainfoPackage(
        section_definitions=[TargetSection.m_def, SourceSection.m_def]
    )

    json_data = pkg.m_to_dict()
    assert 'base_sections' not in json_data['section_definitions'][0]
    assert 'constraints' not in json_data['section_definitions'][0]
    assert MetainfoPackage.m_from_dict(json_data).m_to_dict() == json_data


@pytest.mark.parametrize(
    'ref',
    [
        pytest.param('/section_definitions/0/inner_section_definitions/0', id='base'),
        pytest.param('tests.metainfo.test_references.Referenced', id='python'),
        pytest.param('/TestSection/Referenced', id='metainfo'),
    ],
)
def test_section_reference_deserialize(ref):
    json_data = {
        'm_def': 'nomad.metainfo.metainfo.Package',
        'section_definitions': [
            {
                'base_sections': [ref],
                'name': 'TestSection',
                'inner_section_definitions': [{'name': 'Referenced'}],
                'quantities': [
                    {
                        'name': 'test_quantity',
                        'type': {'type_kind': 'reference', 'type_data': ref},
                    }
                ],
                'sub_sections': [{'name': 'test_sub_section', 'section_def': ref}],
            }
        ],
    }

    pkg = cast(Package, MSection.from_dict(json_data))
    section_def = pkg.section_definitions[0]

    def assert_section(section):
        assert section is not None
        section = section.m_resolved()
        assert isinstance(section, Section)
        assert not isinstance(section, MProxy)

    assert_section(section_def.base_sections[0])
    assert_section(section_def.quantities[0].type.target_section_def)
    assert_section(section_def.sub_sections[0].section_def)


@pytest.mark.parametrize(
    'url,value',
    [
        pytest.param('/referenced', '/referenced', id='archive-plain'),
        pytest.param('#/referenced', '/referenced', id='archive-anchor'),
        pytest.param('/api#/referenced', None, id='api-no-support'),
        pytest.param(
            '../upload/archive/my_entry_id#/referenced',
            '/referenced',
            id='upload-entry-id',
        ),
        pytest.param(
            '../upload/archive/mainfile/my/main/file#/referenced',
            '/referenced',
            id='upload-mainfile',
        ),
        pytest.param(
            '../uploads/my_upload_id/archive/my_entry_id#/referenced',
            '/referenced',
            id='api',
        ),
    ],
)
def test_reference_urls(example_data, url, value):
    class MyContext(Context):
        def resolve_archive_url(self, url):
            if url == '../upload/archive/my_entry_id':
                return example_data
            if url == '../upload/archive/mainfile/my/main/file':
                return example_data
            if url == '../uploads/my_upload_id/archive/my_entry_id':
                return example_data

            raise MetainfoReferenceError()

    if value:
        context = MyContext()
        example_data.m_context = context

    example_data.referencing.section_reference = url

    if value:
        value = example_data.m_resolve(value).m_resolved()
        assert example_data.referencing.section_reference.m_resolved() == value

    else:
        with pytest.raises(MetainfoReferenceError):
            example_data.referencing.section_reference.m_resolved()


def test_file_references(example_data):
    example_data.referencing.file_reference = '../upload/raw/a_file.txt'

    assert example_data.referencing.file_reference == '../upload/raw/a_file.txt'


def test_def_reference():
    definitions = MetainfoPackage.m_from_dict(
        {
            'section_definitions': [
                {
                    'name': 'TestSection',
                    'quantities': [{'name': 'test_quantity', 'type': 'str'}],
                }
            ]
        }
    )

    multiple_definitions = MetainfoPackage.m_from_dict(
        {
            'section_definitions': [
                {
                    'name': 'TestSection',
                    'quantities': [{'name': 'test_quantity', 'type': 'str'}],
                },
                {
                    'name': 'TestSection',
                    'quantities': [{'name': 'test_quantity', 'type': 'str'}],
                },
            ]
        }
    )

    class TestContext(Context):
        def resolve_archive_url(self, url):
            assert url == 'definitions'
            return definitions

    class TestContextMulti(Context):
        def resolve_archive_url(self, url):
            assert url == 'definitions'
            return multiple_definitions

    data_with_section_index = {
        'm_def': 'definitions#section_definitions/0',
        'test_quantity': 'TestValue',
    }

    data_with_section_name = {
        'm_def': 'definitions#section_definitions/TestSection',
        'test_quantity': 'TestValue',
    }

    result_with_section_index = MSection.from_dict(
        data_with_section_index, m_context=TestContext()
    )
    result_with_section_name = MSection.from_dict(
        data_with_section_name, m_context=TestContext()
    )

    assert result_with_section_index.m_to_dict() == {'test_quantity': 'TestValue'}

    assert result_with_section_name.m_to_dict() == {'test_quantity': 'TestValue'}

    with pytest.raises(MetainfoReferenceError):
        MSection.from_dict(data_with_section_name, m_context=TestContextMulti())


@pytest.mark.parametrize('mainfile', ['intra-entry', 'inter-entry'])
def test_parse_with_references(mainfile):
    from nomad.client import parse, normalize_all

    entry_archive = parse(
        os.path.join(
            os.path.dirname(__file__), f'../data/metainfo/{mainfile}.archive.json'
        )
    )[0]
    normalize_all(entry_archive)

    m_def = entry_archive.m_to_dict()['data']['m_def']
    assert '#/definitions/' in m_def


@pytest.mark.parametrize(
    'def_type, value, expected_name',
    [
        pytest.param(
            UserReference, '00000000-0000-0000-0000-000000000001', 'Sheldon Cooper'
        ),
        pytest.param(
            AuthorReference,
            {'first_name': 'Mohammad', 'last_name': 'Nakhaee'},
            'Mohammad Nakhaee',
        ),
        pytest.param(
            AuthorReference, '00000000-0000-0000-0000-000000000001', 'Sheldon Cooper'
        ),
    ],
)
def test_user_author(def_type, value, expected_name):
    class UserAuthorSection(MSection):
        quantity = Quantity(type=def_type)

    section = UserAuthorSection()

    # test assignment
    section.quantity = value
    quantity = section.quantity
    resolved_quantity = quantity.m_resolved()

    assert quantity.m_proxy_value == value
    assert (
        quantity.m_proxy_type.target_section_def.name
        == def_type().target_section_def.name
    )
    assert quantity.m_proxy_section == section
    assert resolved_quantity.name == expected_name

    # test serialization
    serialized_section = section.m_to_dict()
    assert serialized_section['quantity'] == value

    # test deserialization
    deserialized_section = UserAuthorSection().m_from_dict(serialized_section)
    deserialized_quantity = deserialized_section.quantity
    resolved_deserialized_quantity = deserialized_quantity.m_resolved()

    assert deserialized_quantity.m_proxy_value == value
    assert (
        deserialized_quantity.m_proxy_type.target_section_def.name
        == def_type().target_section_def.name
    )
    assert resolved_deserialized_quantity.name == expected_name
