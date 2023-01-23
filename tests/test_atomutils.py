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
from nomad.atomutils import Formula


@pytest.mark.parametrize('format, formula, expected_formula', [
    pytest.param('hill', 'BrBaHCoCIH', 'CH2BaBrCoI', id='large number of elements'),
    pytest.param('hill', 'HCO', 'CHO', id='carbon, hydrogen and others'),
    pytest.param('hill', 'HC', 'CH', id='carbon and hydrogen'),
    pytest.param('hill', 'OC', 'CO', id='carbon and others'),
    pytest.param('hill', 'HAl', 'AlH', id='hydrogen and others'),
    pytest.param('hill', 'TiO', 'OTi', id='no carbon or hydrogen'),
    pytest.param('hill', 'Cu2ZnSn0.5Ge0.5Se4', 'Cu4GeSe8SnZn2', id='non-integer coefficients'),
    pytest.param('hill', 'C2HOZa', 'C2HOX', id='unknown species'),
    pytest.param('iupac', 'AsGa', 'GaAs', id='ordering across periodic table groups'),
    pytest.param('iupac', 'CSi', 'SiC', id='ordering within same periodic table group'),
    pytest.param('iupac', 'ClHe', 'HeCl', id='halogens in the beginning'),
    pytest.param('iupac', 'PoH', 'HPo', id='hydrogen exception'),
    pytest.param('iupac', 'HN', 'NH', id='hydrogen exception'),
    pytest.param('iupac', 'CAc', 'AcC', id='actinide'),
    pytest.param('iupac', 'CLa', 'LaC', id='lanthanide'),
    pytest.param('iupac', 'H2O', 'H2O', id='iupac gcd 1'),
    pytest.param('iupac', 'H2O2', 'HO', id='iupac gcd 2'),
    pytest.param('iupac', 'Mg6O3', 'Mg2O', id='iupac gcd 3'),
    pytest.param('iupac', 'Cu2ZnSn0.5Ge0.5Se4', 'Cu4Zn2SnGeSe8', id='non-integer coefficients'),
    pytest.param('iupac', 'C2HOZa', 'C2HOX', id='unknown species'),
    pytest.param('reduced', 'H2O', 'H2O', id='reduced gcd 1'),
    pytest.param('reduced', 'H2O2', 'HO', id='reduced gcd 2'),
    pytest.param('reduced', 'Mg6O3', 'Mg2O', id='reduced gcd 3'),
    pytest.param('reduced', 'H2O', 'H2O', id='original order'),
    pytest.param('reduced', 'OH2', 'H2O', id='reversed order'),
    pytest.param('reduced', 'Cu2ZnSn0.5Ge0.5Se4', 'Cu4GeSe8SnZn2', id='non-integer coefficients'),
    pytest.param('reduced', 'C2HOZa', 'C2HOX', id='unknown species'),
    pytest.param('anonymous', 'HO', 'AB', id='anonymous gcd 1'),
    pytest.param('anonymous', 'O2H2', 'AB', id='anonymous gcd 2'),
    pytest.param('anonymous', 'Mg6O3', 'A2B', id='anonymous gcd 3'),
    pytest.param('anonymous', 'H2O', 'A2B', id='original order'),
    pytest.param('anonymous', 'OH2', 'A2B', id='reversed order'),
    pytest.param('anonymous', 'Cu2ZnSn0.5Ge0.5Se4', 'A8B4C2DE', id='non-integer coefficients'),
    pytest.param('anonymous', 'C2HOZa', 'A2BCD', id='unknown species'),
])
def test_formula(formula, format, expected_formula):
    assert Formula(formula).format(format) == expected_formula
