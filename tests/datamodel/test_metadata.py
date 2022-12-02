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

from nomad.metainfo import MSection, Quantity, Section, Datetime, MEnum, SubSection
from nomad.datamodel.datamodel import EntryMetadata, SearchableQuantity, EntryArchive
from nomad.datamodel import EntryData


@pytest.mark.parametrize('source_quantity,source_value,target_quantity,target_value', [
    pytest.param(Quantity(type=str), 'test', SearchableQuantity.text_value, 'test', id='text'),
    pytest.param(Quantity(type=MEnum('one', 'two')), 'two', SearchableQuantity.keyword_value, 'two', id='keyword'),
    pytest.param(Quantity(type=np.float64), 1.2, SearchableQuantity.double_value, 1.2, id='np-double'),
    pytest.param(Quantity(type=float), 1.2, SearchableQuantity.double_value, 1.2, id='python-double'),
    pytest.param(Quantity(type=int), 1, SearchableQuantity.long_value, 1, id='long'),
    pytest.param(Quantity(type=Datetime), datetime.fromtimestamp(0, tz=pytz.UTC), SearchableQuantity.date_value, "1970-01-01T00:00:00+00:00", id='date'),
    pytest.param(Quantity(type=str, shape=['*']), ['test'], None, None, id='shape'),
    pytest.param(Quantity(type=str), None, None, None, id='None'),
    pytest.param(Quantity(type=str, default='test'), None, None, None, id='default'),
])
def test_searchable_quantities(source_quantity, source_value, target_quantity, target_value):
    source_quantity.name = 'test_quantity'
    section_def = Section(name='MySection', quantities=[source_quantity], base_sections=[EntryData])
    section = section_def.section_cls()  # pylint: disable=not-callable
    if source_value:
        section.test_quantity = source_value

    archive = EntryArchive(data=section, metadata=EntryMetadata())
    archive.metadata.apply_archive_metadata(archive)

    if target_quantity:
        assert len(archive.metadata.searchable_quantities) == 1
    else:
        assert len(archive.metadata.searchable_quantities) == 0
        return

    searchable_quantity = archive.metadata.searchable_quantities[0]
    assert searchable_quantity.m_get(target_quantity) == target_value
    assert searchable_quantity.path == 'data.test_quantity'
    assert searchable_quantity.section_definition == 'MySection'
    assert searchable_quantity.quantity_name == 'test_quantity'


def test_searchable_quantities_nested():

    class MySubSection(MSection):
        value = Quantity(type=str)

    class MySection(EntryData, MySubSection):
        children = SubSection(section=MySubSection, repeats=True)

    archive = EntryArchive(metadata=EntryMetadata(entry_name='test'))
    archive.data = MySection(value='root')
    archive.data.children.append(MySubSection(value='child1'))
    archive.data.children.append(MySubSection(value='child2'))

    archive.metadata.apply_archive_metadata(archive)
    assert len(archive.metadata.searchable_quantities) == 3
    data = [
        (item.path, item.text_value)
        for item in archive.metadata.searchable_quantities]
    assert data == [('data.value', 'root'), ('data.children.value', 'child1'), ('data.children.value', 'child2')]
