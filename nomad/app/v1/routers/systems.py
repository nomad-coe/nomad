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
from typing import Union, Dict, List
from io import StringIO, BytesIO
from collections import OrderedDict
from enum import Enum

import numpy as np
from fastapi import APIRouter, Depends, Path, Query, status, HTTPException
from fastapi.responses import Response
import ase.io
import ase.build
from MDAnalysis.lib.util import NamedStream
from MDAnalysis.coordinates.PDB import PDBWriter

from nomad.units import ureg
from nomad.utils import strip, deep_get, query_list_to_dict
from nomad.atomutils import Formula
from nomad.normalizing.common import ase_atoms_from_nomad_atoms, mda_universe_from_nomad_atoms
from nomad.datamodel.metainfo.simulation.system import Atoms as NOMADAtoms
from .entries import answer_entry_archive_request

from .auth import create_user_dependency
from ..utils import create_responses
from ..models import User, HTTPExceptionModel

router = APIRouter()
default_tag = 'systems'


def write_pdb(atoms: NOMADAtoms, entry_id: str = None, formula: str = None) -> str:
    '''For writing a PDB file.'''

    # Add custom title that contains the entry id.
    lines = []
    if entry_id is not None:
        lines.append(f'TITLE     NOMAD ENTRY ID: {entry_id}\n')

    # PDB files do not contain a field for the full lattice vectors and the PBC.
    # To work around this, they are stored as REMARKs. The info on full lattice
    # vectors and pbc are used by the the GUI visualizer.
    lattice_vectors = atoms.lattice_vectors
    pbc = atoms.periodic
    if lattice_vectors is not None:
        cell = lattice_vectors.to(ureg.angstrom).magnitude
        lines.append('REMARK 285 LATTICE VECTORS\n')
        lines.append(f'REMARK 285  A: {cell[0, 0]:.3f}, {cell[0, 1]:.3f}, {cell[0, 2]:.3f}\n')
        lines.append(f'REMARK 285  B: {cell[1, 0]:.3f}, {cell[1, 1]:.3f}, {cell[1, 2]:.3f}\n')
        lines.append(f'REMARK 285  C: {cell[2, 0]:.3f}, {cell[2, 1]:.3f}, {cell[2, 2]:.3f}\n')
    else:
        lines.append('REMARK 285 UNITARY VALUES FOR THE UNIT CELL SET BECAUSE UNIT CELL INFORMATION\n')
        lines.append('REMARK 285 WAS MISSING. PROTEIN DATA BANK CONVENTIONS REQUIRE THAT CRYST1\n')
        lines.append('REMARK 285 RECORD IS INCLUDED, BUT THE VALUES ON THIS RECORD ARE MEANINGLESS.\n')
    if pbc is not None:
        pbc = ['TRUE' if x else 'FALSE' for x in pbc]
        lines.append(f'REMARK 285 PBC (A, B, C): {pbc[0]}, {pbc[1]}, {pbc[2]}\n')

    mda_string_stream = StringIO()
    mda_named_stream = NamedStream(mda_string_stream, f'temp.{format}', close=False, reset=False)
    writer = PDBWriter(mda_named_stream, remarks='')
    universe = mda_universe_from_nomad_atoms(atoms)
    writer.write(universe)
    writer.close()

    # We skip the title line that is written by MDA (cannot be disabled otherwise)
    mda_string_stream.seek(0)
    for line in mda_string_stream.readlines():
        if not line.startswith(('REMARK', 'TITLE')):
            lines.append(line)

    content = "".join(lines)
    return content


def write_cif(atoms: NOMADAtoms, entry_id: str = None, formula: str = None) -> str:
    '''For writing a CIF file.'''
    # The ASE CIF writer expects a BytesIO, unlike other formats supported by ASE.
    byte_stream = BytesIO()
    atoms = ase_atoms_from_nomad_atoms(atoms)
    ase.io.write(byte_stream, atoms, format='cif')
    byte_stream.seek(0)
    content = ''
    if entry_id is not None:
        content = f'# NOMAD ENTRY ID: {entry_id}\n'
    content += byte_stream.read().decode('utf-8')
    content = content.replace('data_image0', f'data_{formula}')
    return content


def write_xyz(atoms: NOMADAtoms, entry_id: str, formula: str = None) -> str:
    '''For writing an XYZ file.'''
    stream = StringIO()
    atoms = ase_atoms_from_nomad_atoms(atoms)
    ase.io.write(stream, atoms, format='xyz')
    stream.seek(0)
    content = stream.read()
    if entry_id is not None:
        content = content.replace('"\n', f'" nomad_entry_id=\"{entry_id}\"\n', 1)
    return content


class FormatFeature(str, Enum):
    '''Contains features that are relevant for the different file formats used to
    serialize atomic configurations.
    '''
    NO_UNIT_CELL = 'Cartesian positions without unit cell'
    LATTICE_VECTORS = 'Full lattice vectors'
    PBC = 'Periodic boundary conditions (PBC)'


format_map: Dict[str, dict] = OrderedDict({
    'cif': {
        'label': 'cif',
        'description': 'Crystallographic Information File',
        'extension': 'cif',
        'features': {
            FormatFeature.PBC: False,
            FormatFeature.LATTICE_VECTORS: False,
            FormatFeature.NO_UNIT_CELL: True
        },
        'mime_type': 'chemical/x-cif',
        'writer': write_cif
    },
    'xyz': {
        'label': 'xyz',
        'description': '''XYZ file. The comment line contains information that
        complies with the extended XYZ specification.''',
        'extension': 'xyz',
        'features': {
            FormatFeature.PBC: True,
            FormatFeature.LATTICE_VECTORS: True,
            FormatFeature.NO_UNIT_CELL: True
        },
        'mime_type': 'chemical/x-xyz',
        'writer': write_xyz
    },
    'pdb': {
        'label': 'pdb',
        'description': '''Protein Data Bank file. Note that valid PDB files
        require a CRYST1 record, while certains systems in NOMAD may not have a
        unit cell associated with them. In this case the returned structure file
        will contain a dummy CRYST1 record in order to load the atomic
        positions.''',
        'features': {
            FormatFeature.PBC: False,
            FormatFeature.LATTICE_VECTORS: False,
            FormatFeature.NO_UNIT_CELL: False
        },
        'extension': 'pdb',
        'mime_type': 'chemical/x-pdb',
        'writer': write_pdb
    },
})

format_description = "\n".join([f'- `{format["label"]}`: {format["description"]}' for format in format_map.values()])
format_features_list = []
format_features_list.append('|'.join(['Format'] + [feature.value for feature in FormatFeature]))
format_features_list.append('|'.join([':---'] + [':---:'] * (len(FormatFeature))))
for format in format_map.values():
    format_features_list.append('|'.join([format['label']] + ['&#9745;' if format['features'][feature.value] else '&#9744;' for feature in FormatFeature]))
format_features = "\n".join(format_features_list)


class TempFormatEnum(str, Enum):
    pass


FormatEnum = TempFormatEnum("FormatEnum", {format: format for format in format_map.keys()})  # type: ignore


_file_response = status.HTTP_200_OK, {
    'content': {'application/octet-stream': {}},
    'description': strip('''
        A byte stream with file contents. The content length is not known in advance.
        The final mime-type may be more specific depending on the format.
    ''')}

_not_found_response = status.HTTP_404_NOT_FOUND, {
    'model': HTTPExceptionModel,
    'description': strip('''
        Could not find data for given entry id and path. Check that the
        arguments are correct and that you are authorized to access the data.
    ''')}

_serialization_error_response = status.HTTP_500_INTERNAL_SERVER_ERROR, {
    'model': HTTPExceptionModel,
    'description': strip('''
        Could not serialize the system information in the given format. The
        information may be invalid or incomplete for this format.
    ''')}


@router.get(
    '/{entry_id}',
    tags=[default_tag],
    summary=strip('''
    Build and retrieve an atomistic structure file from data within an entry.
    '''),
    response_class=Response,
    responses=create_responses(_file_response, _not_found_response, _serialization_error_response))
async def get_entry_raw_file(
        entry_id: str = Path(..., description='The unique entry id of the entry to retrieve archive data from.'),
        path: str = Query(
            ...,
            example="run/0/system/0",
            description=strip('''
            Path to a NOMAD System inside the archive. The targeted path should
            point to a system in `run.system` or `results.material.topology`.
            The following path types are supported:

            - `run/0/system/0`: Path to system in `run.system`
            - `results/material/topology/0`: Path to system in `results.material.topology`
            - `run/0/system/-1`: Negative indices are supported.''')
        ),
        format: FormatEnum = Query(  # type: ignore
            default='cif',
            description=f'''The file format for the system. The following formats are supported:

{format_description}

Here is a brief rundown of the different features each format supports:

{format_features}'''
        ),
        user: User = Depends(create_user_dependency(signature_token_auth_allowed=True))):
    '''
    Build and retrieve a structure file containing an atomistic system stored
    within an entry. Note that some formats are more restricted and cannot fully
    describe certains kinds of systems. For examples some entries within NOMAD
    do not contain a unit cell (e.g. molecules), whereas some formats require it
    to be present.
    '''
    # Remove prefix
    for prefix in ['#/']:
        if path.startswith(prefix):
            path = path[len(prefix):]

    # Add indexing
    query_list: List[Union[str, int]] = []
    paths = [x for x in path.split('/') if x != '']
    i = 0
    while i < len(paths):
        a = paths[i]
        query_list.append(a)
        try:
            b = int(paths[i + 1])
            i += 1
            query_list.append(b)
        except Exception:
            pass
        i += 1

    # We extract both atoms and atoms_ref
    value = {'atoms': '*'}
    if 'topology' in path:
        value['atoms_ref'] = 'include-resolved'
        value['indices'] = '*'

    # Fetch the specific part of the archive. If path not found, raise exception
    required = query_list_to_dict(query_list, value)
    required['resolve-inplace'] = True
    query = {'entry_id': entry_id}
    try:
        archive = answer_entry_archive_request(query, required=required, user=user)['data']['archive']
    except Exception as e:
        raise HTTPException(
            status_code=_not_found_response[0],
            detail=_not_found_response[1]['description']
        ) from e

    # The returned archive contains a single section when an index has been
    # specified, and thus we have to set all original indices to zero when
    # extracting the data.
    try:
        result_dict = deep_get(archive, *[0 if isinstance(x, int) else x for x in query_list])
    except Exception as e:
        raise HTTPException(
            status_code=_not_found_response[0],
            detail='The given path does not exist in the archive.'
        ) from e

    # Extract the atoms: they can be given under 'atoms' or 'atoms_ref'
    atoms = result_dict.get('atoms', result_dict.get('atoms_ref'))
    formula = None

    # Write file into stream in memory
    format_info = format_map[format]
    try:
        atoms = NOMADAtoms.m_from_dict(atoms)

        # Use indices to strip down the atoms to include. When the system has
        # several representative indices, the first one is returned.
        if indices := result_dict.get('indices'):
            indices = np.array(indices)
            if len(indices.shape) == 2:
                indices = indices[0]
            if atoms.atomic_numbers is not None:
                atoms.atomic_numbers = atoms.atomic_numbers[indices]
            if atoms.positions is not None:
                atoms.positions = atoms.positions[indices]
            if atoms.species is not None:
                atoms.species = atoms.species[indices]
            if atoms.labels is not None:
                atoms.labels = np.array(atoms.labels)[indices].tolist()

        try:
            formula = Formula(''.join(atoms.labels)).format('iupac')
        except Exception:
            pass

        content = format_info['writer'](atoms, entry_id, formula)
    except Exception as e:
        raise HTTPException(
            status_code=_serialization_error_response[0],
            detail=_serialization_error_response[1]['description']
        ) from e

    # Return the contents of the stream. A simple Response is used instead of
    # StreamingResponse since the content is already in memory and it is not
    # very big
    return Response(
        content=content,
        media_type=format_info['mime_type'],
        headers={'Content-Disposition': f'attachment; filename="{formula or "system"}.{format_info["extension"]}"'}
    )
