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

import re
from io import BytesIO, StringIO

import pytest
import numpy as np
import ase.io
from ase import Atoms as ASEAtoms

from nomad.units import ureg
from nomad.normalizing.common import ase_atoms_from_nomad_atoms
from nomad.datamodel.datamodel import EntryArchive
from nomad.datamodel.results import Results, Material, System
from nomad.datamodel.metainfo.simulation.run import Run
from nomad.datamodel.metainfo.simulation.system import System as SystemRun, Atoms
from nomad.utils.exampledata import ExampleData
from nomad.app.v1.routers.systems import format_map, FormatFeature, WrapModeEnum

from .common import assert_response, assert_browser_download_headers


def ase_atoms(content, format):
    '''Creates an ase.Atoms object given given file contents and format.'''
    format = {
        'pdb': 'proteindatabank'
    }.get(format, format)
    atoms = ase.io.read(StringIO(content), format=format)
    return atoms


def assert_contents(a: str, b: str):
    '''Compares two file contents to each other. Differences in trailing
    whitespace are ignored.
    '''
    def normalize(a):
        return re.sub(' +\n', '\n', a)

    assert normalize(a) == normalize(b)


def assert_atoms(a: ASEAtoms, b: ASEAtoms, compare_cell: bool = True, compare_pbc: bool = False, atol: float = 0, rtol: float = 1e-3):
    '''Compares two different ase.Atoms objects to see if they have the same
    structure.
    '''
    assert np.allclose(a.get_positions(), b.get_positions(), atol=atol, rtol=rtol)
    if compare_cell:
        assert np.allclose(a.get_cell(), b.get_cell(), atol=atol, rtol=rtol)
    if compare_pbc:
        assert np.array_equal(a.get_pbc(), b.get_pbc())
    assert np.array_equal(a.get_atomic_numbers(), b.get_atomic_numbers())
    assert np.array_equal(a.get_chemical_symbols(), b.get_chemical_symbols())


atoms_with_cell = Atoms(
    n_atoms=2,
    labels=["C", "H"],
    species=[6, 1],
    positions=np.array([[0, 0, 0], [1, 1, 1]]) * ureg.angstrom,
    lattice_vectors=np.array([[5, 0, 0], [0, 5, 0], [0, 0, 5]]) * ureg.angstrom,
    periodic=[True, True, True],
)

atoms_without_cell = Atoms(
    n_atoms=2,
    labels=["N", "O"],
    species=[7, 8],
    positions=np.array([[0, 0, 0], [1, 1, 1]]) * ureg.angstrom,
)

atoms_missing_positions = Atoms(
    n_atoms=2,
    labels=["N", "O"],
    species=[7, 8],
)

atoms_wrap_mode = Atoms(
    n_atoms=2,
    labels=["C", "H"],
    species=[6, 1],
    positions=np.array([[-15, -15, -15], [17, 17, 17]]) * ureg.angstrom,
    lattice_vectors=np.array([[5, 0, 0], [0, 5, 0], [0, 0, 5]]) * ureg.angstrom,
    periodic=[True, True, True],
)

atoms_wrap_mode_no_pbc = atoms_wrap_mode.m_copy()
atoms_wrap_mode_no_pbc.periodic = [False, False, False]


@pytest.fixture(scope="module")
def example_data_systems(elastic_module, mongo_module, test_user):
    data = ExampleData(main_author=test_user)
    upload_id = "systems_upload"

    data.create_upload(
        upload_id=upload_id,
        published=True
    )
    archive = EntryArchive(run=[Run(system=[
        SystemRun(atoms=atoms_with_cell),
        SystemRun(atoms=atoms_missing_positions),
        SystemRun(atoms=atoms_without_cell),
        SystemRun(atoms=atoms_wrap_mode),
        SystemRun(atoms=atoms_wrap_mode_no_pbc),
    ])])
    archive.results = Results(
        material=Material(
            topology=[
                System(atoms=atoms_without_cell),
                System(atoms_ref=archive.run[0].system[0].atoms),
                System(atoms_ref=archive.run[0].system[0].atoms, indices=[0]),
                System(atoms_ref=archive.run[0].system[0].atoms, indices=[[0]]),
                System(atoms_ref=archive.run[0].system[2].atoms),
            ]
        )
    )

    data.create_entry(
        upload_id=upload_id,
        entry_id="systems_entry_1",
        mainfile="test_content/test_entry/main-file.json",
        entry_archive=archive
    )

    data.save()

    yield

    data.delete()
    from nomad.search import search
    assert search(query=dict(upload_id=upload_id)).pagination.total == 0


def run_query(entry_id, path, format, client, wrap_mode=None):
    response = client.get(f'systems/{entry_id}/?path={path}&format={format}{f"&wrap_mode={wrap_mode}" if wrap_mode else ""}', headers={})
    return response


@pytest.mark.parametrize("path, status_code", [
    pytest.param('run/0/system/0', 200, id='explicit indexing'),
    pytest.param('run/0/system/-1', 200, id='negative indexing'),
    pytest.param('/run/0/system/0', 200, id='start with slash'),
    pytest.param('/run/0/system/0', 200, id='end with slash'),
    pytest.param('results/material/topology/0', 200, id='saved in topology'),
    pytest.param('results/material/topology/1', 200, id='referenced in topology'),
    pytest.param('run/0/system/1', 500, id='cannot serialize'),
    pytest.param('results/does_not_exist', 404, id='invalid path'),
    pytest.param('results/material/topology/100', 404, id='not found'),
    pytest.param('run/100/system/0', 404, id='not found'),
])
def test_paths(path, status_code, client, example_data_systems, indices=None):
    response = run_query('systems_entry_1', path, 'pdb', client)
    assert_response(response, status_code)


@pytest.mark.parametrize('format, content_expected, filename', [
    pytest.param(
        'pdb',
        '''TITLE     NOMAD ENTRY ID: systems_entry_1
REMARK 285 LATTICE VECTORS
REMARK 285  A: 5.000, 0.000, 0.000
REMARK 285  B: 0.000, 5.000, 0.000
REMARK 285  C: 0.000, 0.000, 5.000
REMARK 285 PBC (A, B, C): TRUE, TRUE, TRUE
CRYST1    5.000    5.000    5.000  90.00  90.00  90.00 P 1           1
ATOM      1  C   UNK X   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  H   UNK X   1       1.000   1.000   1.000  1.00  0.00           H
END
''',
        'CH.pdb',
        id='pdb'
    ),
    pytest.param(
        'cif',
        '''# NOMAD ENTRY ID: systems_entry_1
data_CH
_cell_length_a       5
_cell_length_b       5
_cell_length_c       5
_cell_angle_alpha    90
_cell_angle_beta     90
_cell_angle_gamma    90

_symmetry_space_group_name_H-M    "P 1"
_symmetry_int_tables_number       1

loop_
  _symmetry_equiv_pos_as_xyz
  'x, y, z'

loop_
  _atom_site_label
  _atom_site_occupancy
  _atom_site_fract_x
  _atom_site_fract_y
  _atom_site_fract_z
  _atom_site_thermal_displace_type
  _atom_site_B_iso_or_equiv
  _atom_site_type_symbol
  C1       1.0000 0.00000  0.00000  0.00000  Biso   1.000  C
  H1       1.0000 0.20000  0.20000  0.20000  Biso   1.000  H
''',
        'CH.cif',
        id='cif'
    ),
    pytest.param(
        'xyz',
        '''2
Lattice="5.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 5.0" Properties=species:S:1:pos:R:3 pbc="T T T" nomad_entry_id="systems_entry_1"
C        0.00000000       0.00000000       0.00000000
H        1.00000000       1.00000000       1.00000000
''',
        'CH.xyz',
        id='xyz'
    ),
])
def test_formats_with_cell(format, content_expected, filename, client, example_data_systems):
    '''Test that writing a structure with valid unit cell information produces
    the expected output.
    '''
    format_info = format_map[format]
    response = run_query(
        'systems_entry_1',
        'run/0/system/0',
        format_info['label'],
        client
    )
    assert_response(response, 200)
    assert_browser_download_headers(
        response,
        format_info['mime_type'],
        filename
    )
    content = response.content.decode('utf-8')
    assert_contents(content, content_expected)
    atoms = ase_atoms(content, format)
    assert_atoms(
        atoms,
        ase_atoms_from_nomad_atoms(atoms_with_cell),
        compare_cell=True,
        compare_pbc=format_info['features'][FormatFeature.PBC],
    )


@pytest.mark.parametrize('format, content_expected, filename', [
    pytest.param(
        'pdb',
        '''TITLE     NOMAD ENTRY ID: systems_entry_1
REMARK 285 UNITARY VALUES FOR THE UNIT CELL SET BECAUSE UNIT CELL INFORMATION
REMARK 285 WAS MISSING. PROTEIN DATA BANK CONVENTIONS REQUIRE THAT CRYST1
REMARK 285 RECORD IS INCLUDED, BUT THE VALUES ON THIS RECORD ARE MEANINGLESS.
CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1
ATOM      1  N   UNK X   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  O   UNK X   1       1.000   1.000   1.000  1.00  0.00           O
END
''',
        'NO.pdb',
        id='pdb'
    ),
    pytest.param(
        'cif',
        '''# NOMAD ENTRY ID: systems_entry_1
data_NO
loop_
  _atom_site_label
  _atom_site_occupancy
  _atom_site_Cartn_x
  _atom_site_Cartn_y
  _atom_site_Cartn_z
  _atom_site_thermal_displace_type
  _atom_site_B_iso_or_equiv
  _atom_site_type_symbol
  N1       1.0000 0.00000  0.00000  0.00000  Biso   1.000  N
  O1       1.0000 1.00000  1.00000  1.00000  Biso   1.000  O
''',
        'NO.cif',
        id='cif'
    ),
    pytest.param(
        'xyz',
        '''2
Properties=species:S:1:pos:R:3 pbc="F F F" nomad_entry_id="systems_entry_1"
N        0.00000000       0.00000000       0.00000000
O        1.00000000       1.00000000       1.00000000
''',

        'NO.xyz',
        id='xyz'
    ),
])
def test_formats_without_cell(format, content_expected, filename, client, example_data_systems):
    '''Test that writing a structure without unit cell information produces the
    expected output. Note that certains formats cannot be serialized without a
    dummy placeholder cell.
    '''
    format_info = format_map[format]
    response = run_query('systems_entry_1', '/run/0/system/2', format, client)
    assert_response(response, 200)
    assert_browser_download_headers(
        response,
        format_info['mime_type'],
        filename
    )
    content = response.content.decode('utf-8')
    assert_contents(content, content_expected)
    atoms = ase_atoms(content, format)
    assert_atoms(
        atoms,
        ase_atoms_from_nomad_atoms(atoms_without_cell),
        compare_cell=format_info['features'][FormatFeature.NO_UNIT_CELL],
        compare_pbc=format_info['features'][FormatFeature.PBC],
    )


@pytest.mark.parametrize("path, filename, n_atoms", [
    pytest.param('results/material/topology/2', 'C.cif', 1, id='1D indices'),
    pytest.param('results/material/topology/3', 'C.cif', 1, id='2D indices'),
])
def test_indices(path, filename, n_atoms, client, example_data_systems):
    '''Test that systems where indices have been specified are returned
    correctly by including only a subset of atoms.
    '''
    format_info = format_map['cif']
    response = run_query('systems_entry_1', path, 'cif', client)
    assert_response(response, 200)
    assert_browser_download_headers(
        response,
        format_info['mime_type'],
        filename
    )
    atoms = ase.io.read(BytesIO(response.content), format='cif')
    assert len(atoms) == n_atoms


@pytest.mark.parametrize("path, wrap_mode, expected_positions", [
    pytest.param('/run/0/system/3', None, [[-15, -15, -15], [17, 17, 17]], id='default'),
    pytest.param('/run/0/system/3', WrapModeEnum.original, [[-15, -15, -15], [17, 17, 17]], id='original'),
    pytest.param('/run/0/system/3', WrapModeEnum.wrap, [[0, 0, 0], [2, 2, 2]], id='wrap, pbc=[1, 1, 1]'),
    pytest.param('/run/0/system/4', WrapModeEnum.wrap, [[-15, -15, -15], [17, 17, 17]], id='wrap, pbc=[0, 0, 0]'),
    pytest.param('/run/0/system/3', WrapModeEnum.unwrap, [[1.5, 1.5, 1.5], [3.5, 3.5, 3.5]], id='unwrap, pbc=[1, 1, 1]'),
    pytest.param('/run/0/system/4', WrapModeEnum.unwrap, [[-15, -15, -15], [17, 17, 17]], id='unwrap, pbc=[0, 0, 0]'),
])
def test_wrap_mode(path, wrap_mode, expected_positions, client, example_data_systems):
    '''Test that the wrap_mode parameter is handled correctly.
    '''
    response = run_query('systems_entry_1', path, 'xyz', client, wrap_mode)
    assert_response(response, 200)
    atoms = ase.io.read(StringIO(response.text), format='xyz')
    positions = atoms.get_positions()
    assert np.allclose(positions, expected_positions)
