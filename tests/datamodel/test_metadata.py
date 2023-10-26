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
from datetime import datetime
import pytz

from nomad.metainfo import Quantity, MSection, SubSection, Datetime, MEnum
from nomad.datamodel.datamodel import EntryMetadata, SearchableQuantity, EntryArchive
from nomad.metainfo.elasticsearch_extension import schema_separator
from nomad.datamodel import EntryData
# from tests.data.schemas.nomadschemaexample.schema import MySchema
from tests.config import python_schema_name


@pytest.mark.parametrize('source_quantity, source_value, target_quantity, target_value', [
    pytest.param(Quantity(type=str), 'test', SearchableQuantity.str_value, 'test', id='text'),
    pytest.param(Quantity(type=MEnum('one', 'two')), 'two', SearchableQuantity.str_value, 'two', id='keyword'),
    pytest.param(Quantity(type=np.float64), 1.2, SearchableQuantity.float_value, 1.2, id='np-double'),
    pytest.param(Quantity(type=float), 1.2, SearchableQuantity.float_value, 1.2, id='python-double'),
    pytest.param(Quantity(type=int), 1, SearchableQuantity.int_value, 1, id='long'),
    pytest.param(Quantity(type=Datetime), datetime.fromtimestamp(0, tz=pytz.UTC), SearchableQuantity.datetime_value, datetime(1970, 1, 1, 0, 0, tzinfo=pytz.UTC), id='date'),
    pytest.param(Quantity(type=str, shape=['*']), ['test'], None, None, id='shape'),
    pytest.param(Quantity(type=str), None, None, None, id='None'),
    pytest.param(Quantity(type=str, default='test'), None, None, None, id='default'),
])
def test_search_quantities(source_quantity, source_value, target_quantity, target_value):

    class Base(MSection):
        test_quantity = source_quantity

    class TestSchema(EntryData, Base):
        pass

    schema = TestSchema()
    if source_value:
        schema.test_quantity = source_value

    archive = EntryArchive(data=schema, metadata=EntryMetadata())
    archive.metadata.apply_archive_metadata(archive)

    if target_quantity:
        assert len(archive.metadata.search_quantities) == 1
    else:
        assert len(archive.metadata.search_quantities) == 0
        return

    searchable_quantity = archive.metadata.search_quantities[0]
    assert searchable_quantity.m_get(target_quantity) == target_value
    assert searchable_quantity.id == f'data.test_quantity{schema_separator}test_metadata.TestSchema'
    assert searchable_quantity.definition == 'test_metadata.Base.test_quantity'


def test_search_quantities_nested():

    class MySubSection(MSection):
        value = Quantity(type=str)

    class TestSchema(EntryData, MySubSection):
        children = SubSection(section=MySubSection, repeats=True)

    archive = EntryArchive(metadata=EntryMetadata(entry_name='test'))
    archive.data = TestSchema(value='root')
    archive.data.children.append(MySubSection(value='child1'))
    archive.data.children.append(MySubSection(value='child2'))

    archive.metadata.apply_archive_metadata(archive)
    assert len(archive.metadata.search_quantities) == 3
    data = [
        (item.id, item.str_value)
        for item in archive.metadata.search_quantities]
    assert data == [
        (f'data.value{schema_separator}test_metadata.TestSchema', 'root'),
        (f'data.children.value{schema_separator}test_metadata.TestSchema', 'child1'),
        (f'data.children.value{schema_separator}test_metadata.TestSchema', 'child2')
    ]


def populate_child(data):
    from nomadschemaexample.schema import MySection
    data.child = MySection(name='test')


def populate_reference(data):
    from nomadschemaexample.schema import MySection
    data.child = MySection()
    data.reference_section = data.child


@pytest.mark.parametrize('source_quantity, source_value, target_quantity, target_value, definition', [
    pytest.param('name', 'test', SearchableQuantity.str_value, 'test', f'{python_schema_name}.name', id='str'),
    pytest.param('count', 1, SearchableQuantity.int_value, 1, f'{python_schema_name}.count', id='int'),
    pytest.param('frequency', 1.0, SearchableQuantity.float_value, 1.0, f'{python_schema_name}.frequency', id='float'),
    pytest.param('timestamp', datetime.fromtimestamp(0, tz=pytz.UTC), SearchableQuantity.datetime_value, datetime(1970, 1, 1, 0, 0, tzinfo=pytz.UTC), f'{python_schema_name}.timestamp', id='datetime'),
    pytest.param('non_scalar', np.eye(3), None, np.eye(3), None, id='non-scalar'),
    pytest.param('reference_section', populate_reference, None, datetime(1970, 1, 1, 0, 0, tzinfo=pytz.UTC), None, id='references are skipped'),
    pytest.param('child.name', populate_child, SearchableQuantity.str_value, 'test', 'nomadschemaexample.schema.MySection.name', id='child quantity'),
])
def test_search_quantities_plugin(plugin_schema, source_quantity, source_value, target_quantity, target_value, definition):
    '''Tests that different types of search quantities are loaded correctly from
    plugin and saved into the search_quantities field.'''
    from nomadschemaexample.schema import MySchema

    data = MySchema()

    if callable(source_value):
        source_value(data)
    else:
        setattr(data, source_quantity, source_value)
    archive = EntryArchive(data=data, metadata=EntryMetadata())
    archive.metadata.apply_archive_metadata(archive)

    if target_quantity:
        assert len(archive.metadata.search_quantities) == 1
    else:
        assert len(archive.metadata.search_quantities) == 0
        return

    searchable_quantity = archive.metadata.search_quantities[0]
    assert searchable_quantity.m_get(target_quantity) == target_value
    assert searchable_quantity.id == f'data.{source_quantity}{schema_separator}{python_schema_name}'
    assert searchable_quantity.definition == definition
