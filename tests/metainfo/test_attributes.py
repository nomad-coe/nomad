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

from nomad.metainfo import MSection, Quantity, Attribute, MEnum, Reference, Datetime
from nomad.metainfo.metainfo import MQuantity, Definition
from nomad.metainfo.util import validate_allowable_unit
from nomad.units import ureg


@pytest.mark.parametrize(
    'type,errors,value',
    [
        pytest.param(str, 0, 'test_value', id='str'),
        pytest.param(np.float64, 0, 1.1, id='numpy'),
        pytest.param(MEnum('value1'), 0, 'value1', id='enum'),
        pytest.param(Datetime, 0, datetime.datetime.now(tz=pytz.utc), id='datetime'),
        pytest.param(Reference(Quantity.m_def), 1, None, id='reference'),
    ],
)
def test_attributes(type, errors, value):
    class MySection(MSection):
        my_quantity = Quantity(
            type=str, attributes=[Attribute(name='my_quantity_attribute', type=type)]
        )

    assert Definition.all_attributes.derived is not None
    assert len(MySection.m_def.m_all_validate()[0]) == errors

    assert MySection.my_quantity.attributes[0].name == 'my_quantity_attribute'
    assert MySection.my_quantity.attributes[0].type == type

    if errors > 0:
        return

    section = MySection()

    assert (
        section.m_get_quantity_attribute('my_quantity', 'my_quantity_attribute') is None
    )
    section.m_set_quantity_attribute('my_quantity', 'my_quantity_attribute', value)
    assert (
        section.m_get_quantity_attribute('my_quantity', 'my_quantity_attribute')
        == value
    )

    section.my_quantity = 'test'

    json_data = section.m_to_dict()
    section = MySection.m_from_dict(json_data)

    assert (
        section.m_get_quantity_attribute('my_quantity', 'my_quantity_attribute')
        == value
    )


@pytest.mark.parametrize(
    'name,value',
    [
        pytest.param(None, None, id='none'),
        pytest.param('m_source_unit', ureg('m').units, id='source-unit'),
    ],
)
def test_m_attributes(name, value):
    class MySection(MSection):
        my_quantity = Quantity(
            type=float, attributes=[Attribute(name='m_source_unit', type=str)]
        )

    section = MySection(my_quantity=1)
    if name:
        section.m_set_quantity_attribute('my_quantity', name, value)

    json_data = section.m_to_dict()
    section = MySection.m_from_dict(json_data)

    if name:
        assert section.m_get_quantity_attribute('my_quantity', name) == value


def test_variable_name():
    class MySection(MSection):
        MY_quantity = Quantity(
            type=str,
            variable=True,
            attributes=[Attribute(name='MY_attribute', type=str, variable=True)],
        )

    section = MySection()

    section.MY_quantity = MQuantity.wrap('v1', 'MY_quantity')
    section.m_set_quantity_attribute('MY_quantity', 'MY_attribute', 'v1')
    assert section.MY_quantity == 'v1'
    assert section.m_get_quantity_attribute('MY_quantity', 'MY_attribute') == 'v1'

    section.test_quantity = MQuantity.wrap('v2', 'test_quantity')
    section.m_set_quantity_attribute('test_quantity', 'test_attribute', 'v2')
    assert section.test_quantity == 'v2'
    assert section.m_get_quantity_attribute('test_quantity', 'test_attribute') == 'v2'

    with pytest.raises(ValueError):
        section.completely_off = 'v1'


@pytest.mark.parametrize(
    'token,units,result',
    [
        pytest.param('[length]', ['m', 'cm', 'm^2/m'], True, id='length_true'),
        pytest.param('[length]', ['m/m/m', '1/m'], False, id='length_false'),
        pytest.param(
            'dimensionless',
            ['1', 'm/m', 'kg*m/s/s/m^2/MPa'],
            True,
            id='dimensionless_true',
        ),
    ],
)
def test_unit_compatibility(token, units, result):
    assert validate_allowable_unit(token, units) == result

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


def test_repeating_quantity():
    class MySection(MSection):
        TEST_repeat = Quantity(repeats=True, type=float, unit='m')
        TEST_nonrepeat = Quantity(variable=True, type=float)

    my_section = MySection()

    my_section.TEST_repeat = MQuantity.wrap(ureg.Quantity(1.0, 'cm'))

    assert my_section.TEST_repeat.m == 0.01  # pylint: disable=E1101

    my_section.TEST_repeat = None

    assert my_section.TEST_repeat is None

    with pytest.raises(ValueError):
        _ = my_section.instance_repeat

    my_section.instance_nonrepeat = MQuantity('instance_nonrepeat', 1.23)

    assert my_section.instance_nonrepeat == 1.23

    my_section.instance_nonrepeat = None

    _ = my_section.instance_nonrepeat
