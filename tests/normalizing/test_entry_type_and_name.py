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
import tests
from .conftest import run_normalize, run_processing


@pytest.fixture(scope='session')
def unknown_material_archive():
    archive = tests.parsing.test_parsing.parse_file(
        ('parsers/lammps', 'tests/data/normalizers/workflow/lammps/log.lammps')
    )
    return run_normalize(archive)


@pytest.mark.parametrize(
    'fixture, entry_type, entry_name',
    [
        pytest.param(
            'dmft',
            'w2dynamics DMFT SinglePoint',
            'KSi2Br w2dynamics DMFT SinglePoint simulation',
            id='DMFT SinglePoint for inorganic material',
        ),
        pytest.param(
            'geometry_optimization',
            'VASP GeometryOptimization',
            'Si VASP GeometryOptimization simulation',
            id='Workflow for inorganic material',
        ),
        pytest.param(
            'unknown_material_archive',
            'LAMMPS MolecularDynamics',
            'X LAMMPS MolecularDynamics simulation',
            id='Unknown material',
        ),
        pytest.param(
            'organic_formula',
            'VASP GeometryOptimization',
            'CHCl3 VASP GeometryOptimization simulation',
            id='Organic material',
        ),
        pytest.param(
            'organic_carbonyl_formula',
            'VASP GeometryOptimization',
            'CAgO VASP GeometryOptimization simulation',
            id='Organic carbonyl material',
        ),
        pytest.param(
            'inorganic_carbonyl_formula',
            'VASP GeometryOptimization',
            'FeC5O5 VASP GeometryOptimization simulation',
            id='Inorganic carbonyl material',
        ),
        pytest.param(
            'inorganic_special_formula',
            'VASP SinglePoint',
            'KHCO3 VASP SinglePoint simulation',
            id='Inorganic material with special formula',
        ),
        pytest.param(
            'unknown_program',
            'not processed SinglePoint',
            'Si not processed SinglePoint simulation',
            id='Unknown program name',
        ),
    ],
)
def test_entry_type_and_name(fixture, entry_type, entry_name, request):
    """Tests if entry_type and entry_name for simulations are properly defined."""
    archive = request.getfixturevalue(fixture)
    assert archive.metadata.entry_type == entry_type
    assert archive.metadata.entry_name == entry_name
