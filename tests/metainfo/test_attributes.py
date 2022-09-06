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

import datetime

import pytest
import numpy as np
import pytz

from nomad.metainfo import MSection, Quantity, Attribute, MEnum, Reference, Datetime, Property
from nomad.metainfo.metainfo import SubSection, MetainfoError
from nomad.metainfo.nx_unit import validate_allowable_list
from nomad.units import ureg


@pytest.mark.parametrize('type,errors,value', [
    pytest.param(str, 0, 'test_value', id='str'),
    pytest.param(np.float64, 0, 1.1, id='numpy'),
    pytest.param(MEnum('value1'), 0, 'value1', id='enum'),
    pytest.param(Datetime, 0, datetime.datetime.now(tz=pytz.utc), id='datetime'),
    pytest.param(Reference(Quantity.m_def), 1, None, id='reference')
])
def test_attributes(type, errors, value):
    class MySection(MSection):
        my_quantity = Quantity(
            type=str,
            attributes=[
                Attribute(name='my_quantity_attribute', type=type)
            ]
        )
        my_section = SubSection(
            section=Quantity.m_def,
            attributes=[
                Attribute(name='my_section_attribute', type=type)
            ]
        )

    assert Property.all_attributes.derived is not None
    assert len(MySection.m_def.m_all_validate()[0]) == errors * 2

    assert MySection.my_quantity.attributes[0].name == 'my_quantity_attribute'
    assert MySection.my_quantity.attributes[0].type == type
    assert MySection.my_section.attributes[0].name == 'my_section_attribute'
    assert MySection.my_section.attributes[0].type == type

    if errors > 0:
        return

    section = MySection()
    attributes = [
        (MySection.my_quantity, 'my_quantity_attribute'),
        (MySection.my_section, 'my_section_attribute')
    ]
    for property, attribute in attributes:
        assert section.m_get_attribute(property, attribute) is None
        section.m_set_attribute(property, attribute, value)
        assert section.m_get_attribute(property, attribute) == value

    json_data = section.m_to_dict()
    assert json_data == {}

    section.my_quantity = 'test'
    section.my_section = MySection.my_quantity

    json_data = section.m_to_dict()
    section = MySection.m_from_dict(json_data)
    for property, attribute in attributes:
        assert section.m_get_attribute(property, attribute) == value


@pytest.mark.parametrize('name,value', [
    pytest.param(None, None, id='none'),
    pytest.param('m_source_unit', ureg('m').units, id='source-unit')
])
def test_m_attributes(name, value):
    class MySection(MSection):
        my_quantity = Quantity(type=float)

    section = MySection(my_quantity=1)
    if name:
        section.m_set_attribute('my_quantity', name, value)

    json_data = section.m_to_dict()
    section = MySection.m_from_dict(json_data)

    if name:
        assert section.m_get_attribute('my_quantity', name) == value
    else:
        for key in json_data.keys():
            assert '@' not in key


def test_variable_name():
    class MySection(MSection):
        MY_quantity = Quantity(
            type=str, variable=True,
            attributes=[
                Attribute(name='MY_attribute', type=str, variable=True)
            ]
        )

    section = MySection()

    section.MY_quantity = 'v1'
    section.m_set_attribute('MY_quantity', 'MY_attribute', 'v1')
    assert section.MY_quantity == 'v1'
    assert section.m_get_attribute('MY_quantity', 'MY_attribute') == 'v1'

    section.test_quantity = 'v2'
    section.m_set_attribute('test_quantity', 'test_attribute', 'v2')
    assert section.MY_quantity == 'v2'
    assert section.m_get_attribute('MY_quantity', 'm_source_name') == 'test_quantity'
    assert section.m_get_attribute('MY_quantity', 'MY_attribute') == 'v2'
    assert section.test_quantity == 'v2'
    assert section.m_get_attribute('test_quantity', 'test_attribute') == 'v2'

    with pytest.raises(MetainfoError):
        section.completely_off = 'v1'


@pytest.mark.parametrize('token,units,result', [
    pytest.param('[length]', ['m', 'cm', 'm^2/m'], True, id='length_true'),
    pytest.param('[length]', ['m/m/m', '1/m'], False, id='length_false'),
    pytest.param('dimensionless', ['1', 'm/m', 'kg*m/s/s/m^2/MPa'], True, id='dimensionless_true')
])
def test_nx_unit_compatibility(token, units, result):
    assert validate_allowable_list(token, units) == result

    if result:
        class MySection(MSection):
            numerical = Quantity(
                type=np.dtype(np.float64), dimensionality=token, unit=units[0]
            )

        section = MySection()
        for u in units:
            section.numerical = 1 * ureg.parse_units(u)
    else:
        with pytest.raises(TypeError):
            class MySection(MSection):
                numerical = Quantity(
                    type=np.dtype(np.float64), dimensionality=token, unit=units[0]
                )

            section = MySection()
            for u in units:
                section.numerical = 1 * ureg.parse_units(u)
