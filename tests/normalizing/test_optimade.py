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
from ase import Atoms

from tests.normalizing.conftest import get_template_for_structure


@pytest.mark.parametrize('formula, expected', [
    ('NaClHC', 'CHClNa'),
    ('NaClH2', 'ClH2Na')
])
def test_chemical_formula_hill(formula, expected):
    archive = get_template_for_structure(Atoms(formula))
    assert archive.metadata.optimade.chemical_formula_hill == expected


@pytest.mark.parametrize('formula, expected', [
    ('Na3Cl2H', 'A3B2C'),
    ('NaNaNaClClHH', 'A3B2C2')
])
def test_chemical_formula_anonymous(formula, expected):
    archive = get_template_for_structure(Atoms(formula))
    assert archive.metadata.optimade.chemical_formula_anonymous == expected


@pytest.mark.parametrize('formula, expected', [
    ('Na3Cl2H', 'Cl2HNa3'),
    ('NaNaNaClClHH', 'Cl2H2Na3')
])
def test_chemical_formula_reduced(formula, expected):
    archive = get_template_for_structure(Atoms(formula))
    assert archive.metadata.optimade.chemical_formula_reduced == expected
