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
import json
import datetime
import pytz
from nomad.metainfo import MSection, Quantity, Unit, units, JSON, Dimension, Datetime, Capitalized, Bytes


@pytest.mark.parametrize('def_type, value', [
    pytest.param(str, 'hello', id='str'),
    pytest.param(int, 23, id='int'),
    pytest.param(float, 3.14e23, id='float'),
    pytest.param(bool, True, id='bool'),
    pytest.param(JSON, dict(key='value'), id='JSON'),
    pytest.param(Unit, units.parse_units('m*m/s'), id='Unit'),
    pytest.param(Dimension, '*', id='Dimension-*'),
    pytest.param(Dimension, 1, id='Dimension-1'),
    pytest.param(Dimension, 'quantity', id='Dimension-quantity'),
    pytest.param(Datetime, datetime.datetime.now(datetime.timezone.utc), id='Datetime'),
    pytest.param(Datetime, datetime.datetime.now(pytz.timezone('America/Los_Angeles')), id='Datetime'),
    pytest.param(Datetime, datetime.date.today(), id='Date'),
    pytest.param(Capitalized, 'Hello', id='Capitalize'),
    pytest.param(Bytes, b'hello', id='Bytes')
])
def test_basic_types(def_type, value):
    class TestSectionA(MSection):
        quantity = Quantity(type=def_type)

    section = TestSectionA()
    assert section.quantity is None
    section.quantity = value
    if not isinstance(value, datetime.datetime) and isinstance(value, datetime.date):
        assert section.quantity == datetime.datetime.combine(value, datetime.datetime.min.time()).replace(
            tzinfo=pytz.utc)
    else:
        assert section.quantity == value

    section_serialized = section.m_to_dict()
    json.dumps(section_serialized)
    section = TestSectionA.m_from_dict(section_serialized)
    if not isinstance(value, datetime.datetime) and isinstance(value, datetime.date):
        assert section.quantity == datetime.datetime.combine(value, datetime.datetime.min.time()).replace(
            tzinfo=pytz.utc)
    else:
        assert section.quantity == value

    class TestSectionB(MSection):
        quantity = Quantity(type=def_type, default=value)

    section = TestSectionB()
    assert section.quantity == value
    assert 'quantity' not in section.m_to_dict()


@pytest.mark.parametrize('def_type, orig_value, normalized_value', [
    pytest.param(Unit, 'm*m/s', units.parse_units('m*m/s'), id='Unit'),
    pytest.param(Datetime, '1970-01-01 01:00:00', None, id='Datetime-str'),
    pytest.param(Datetime, '1970-01-01 01:00+01', None, id='Datetime-str-tz'),
    pytest.param(Datetime, '1970-01-01 01:00:00.0000', None, id='Datetime-str-ms'),
    pytest.param(Datetime, 'Wed, 01 Jan 1970 00:00:00 -0100', None, id='Datetime-rfc822'),
    pytest.param(Datetime, '1970-01-01T00:00:00Z', None, id='Datetime-aniso861-time'),
    pytest.param(Datetime, '1970-01-01', None, id='Datetime-aniso861-date'),
    pytest.param(Datetime, '2022-05-19T05:16:32.237914-07:00', None, id='Datetime-conversion-from-localtime-to-UTC'),
    pytest.param(Capitalized, 'hello', 'Hello', id='Capitalize')
])
def test_value_normalization(def_type, orig_value, normalized_value):
    class TestSection(MSection):
        quantity = Quantity(type=def_type)

    section = TestSection()
    assert section.quantity is None
    section.quantity = orig_value
    assert normalized_value is None or section.quantity == normalized_value
