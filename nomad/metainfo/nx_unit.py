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

from typing import Union

import pint

from nomad.units import ureg


def validate_allowable_list(
        dimensionality: Union[str, None], allowable_list: Union[str, list, pint.Unit, pint.Quantity]) -> bool:
    '''
    For a given list of units, e.g., ['m', 'cm', 'mm'], and a target NX unit token such as 'NX_LENGTH',
    this function check the compatibility of the target unit with the list of units.

    Returns:
        True if ALL units are compatible with the unit token (dimensionality).
        False if at least one unit cannot be represented by the unit token (dimensionality).
    '''
    if not dimensionality:
        return True

    if isinstance(allowable_list, str):
        if dimensionality in ('1', 'dimensionless'):
            return ureg.Quantity(1, allowable_list).dimensionless

        try:
            return ureg.Quantity(1, allowable_list).check(dimensionality)
        except KeyError:
            return False

    if isinstance(allowable_list, (pint.Unit, pint.Quantity)):
        if dimensionality == ('1', 'dimensionless'):
            return allowable_list.dimensionless

        return allowable_list.dimensionality == dimensionality

    for unit in allowable_list:
        if not validate_allowable_list(dimensionality, unit):
            return False

    return True
