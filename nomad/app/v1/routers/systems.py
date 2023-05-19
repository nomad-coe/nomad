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
from io import StringIO
from collections import OrderedDict
from enum import Enum

from fastapi import APIRouter, Depends, Path, Query, status, HTTPException
from fastapi.responses import StreamingResponse
import ase.io
import ase.build
from MDAnalysis.lib.util import NamedStream
from MDAnalysis.coordinates.PDB import PDBWriter

from nomad.units import ureg
from nomad.utils import strip, deep_get, query_list_to_dict
from nomad.normalizing.common import ase_atoms_from_nomad_atoms, mda_universe_from_nomad_atoms
from nomad.datamodel.metainfo.simulation.system import Atoms as NOMADAtoms
from .entries import answer_entry_archive_request, _bad_id_response

from .auth import create_user_dependency
from ..utils import create_responses
from ..models import User, HTTPExceptionModel

router = APIRouter()
default_tag = 'systems'
format_list: List[dict] = [
    {
        'label': 'pdb',
        'description': 'Protein Data Bank file',
        'extension': 'pdb',
        'mime_type': 'chemical/x-pdb',
        'reader': 'mdanalysis',
    }
]
format_map: Dict[str, dict] = OrderedDict()
for fmt in format_list:
    format_map[fmt['label']] = fmt

format_description = "\n".join([f' - `{format["label"]}`: {format["description"]}' for format in format_map.values()])


class TempFormatEnum(str, Enum):
    pass


FormatEnum = TempFormatEnum("FormatEnum", {format: format for format in format_map.keys()})  # type: ignore

_bad_path_response = status.HTTP_404_NOT_FOUND, {
    'model': HTTPExceptionModel,
    'description': strip('The given path in the archive was not found.')}

_system_file_response = 200, {
    'content': {'application/octet-stream': {}},
    'description': strip('''
        A byte stream with system file contents. The content length is not known in advance.
    ''')}


@router.get(
    '/{entry_id}',
    tags=[default_tag],
    summary=strip('''
    Build and retrieve an atomistic structure file from data within an entry.
    '''),
    response_class=StreamingResponse,
    responses=create_responses(_bad_id_response, _bad_path_response, _system_file_response))
async def get_entry_raw_file(
        entry_id: str = Path(..., description='The unique entry id of the entry to retrieve archive data from.'),
        path: str = Query(
            ...,
            example="run/0/system/0",
            description=strip('''
            Path to a NOMAD System inside the archive. The targeted path should
            point to a system in `nomad.datamodel.simulation.system` or
            `nomad.results.material.topology`. The following path types are
            supported:

            - `run/0/system/0`: Path to system in `run`
            - `results/material/topology/0`: Path to system in `results`
            - `run/0/system/-1`: Negative indices are supported.''')
        ),
        format: FormatEnum = Query(  # type: ignore
            default='pdb',
            description=f'The file format for the system. The following formats are supported:\n{format_description}'
        ),
        user: User = Depends(create_user_dependency(signature_token_auth_allowed=True))):
    '''
    Build and retrieve a structure file containing an atomistic system stored
    within an entry. The file is streamed and will contain atomic positions,
    labels and the simulation cell together with periodicity if the file format
    supports them.
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
    query_list.append('atoms')

    # Fetch the specific part of the archive. If path not found, raise exception
    required = query_list_to_dict(query_list, '*')
    required['resolve-inplace'] = False
    query = {'entry_id': entry_id}
    archive = answer_entry_archive_request(query, required=required, user=user)['data']['archive']
    # The returned archive contains a single section when an index has been
    # specified, and thus we have to set all original indices to zero when
    # extracting the data
    result_dict = deep_get(archive, *[0 if isinstance(x, int) else x for x in query_list])

    # Transform the system into ase atoms, catch any errors
    format_info = format_map[format]
    stringio = StringIO()

    try:
        atoms = NOMADAtoms.m_from_dict(result_dict)
        if format_info['reader'] == 'ase':
            atoms = ase_atoms_from_nomad_atoms(atoms)
            ase_format_map = {
                'pdb': 'proteindatabank'
            }
            ase_format = ase_format_map[format]
            ase.io.write(stringio, atoms, ase_format)
        elif format_info['reader'] == 'mdanalysis':
            # For some reason NamedStream tries to close the stream: here we disable
            # closing. The memory will be freed when stringio goes out of scope.
            stringio.close = lambda: None  # type: ignore[assignment]
            namedstream = NamedStream(stringio, f'temp.{format}', close=False, reset=False)
            globals()[f'write_{format}'](atoms, namedstream)
        else:
            raise ValueError("No reader for filetype.")
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail='Could not serialize the system system information in the given format. The information may be invalid or incomplete for this format.'
        ) from e

    # Stream the system file from the in-memory object. A single system
    # should fit easily into memory.
    stringio.seek(0)
    return StreamingResponse(
        stringio,
        media_type=format_info['mime_type'],
        headers={'Content-Disposition': f'filename=system.{format_info["extension"]}'}
    )


def write_pdb(atoms: NOMADAtoms, stream: Union[StringIO, NamedStream]) -> None:
    '''For writing a PDB file.'''
    writer = PDBWriter(stream)
    writer.TITLE('')

    # PDB files do not contain a field for the full lattice vectors and
    # the PBC. To work around this, they are stored as REMARKs.
    remarks = []
    lattice_vectors = atoms.lattice_vectors
    pbc = atoms.periodic
    if lattice_vectors is not None:
        cell = lattice_vectors.to(ureg.angstrom).magnitude
        remarks.append('LATTICE VECTORS')
        remarks.append(f' A: {cell[0, 0]:.3f}, {cell[0, 1]:.3f}, {cell[0, 2]:.3f}')
        remarks.append(f' B: {cell[1, 0]:.3f}, {cell[1, 1]:.3f}, {cell[1, 2]:.3f}')
        remarks.append(f' C: {cell[2, 0]:.3f}, {cell[2, 1]:.3f}, {cell[2, 2]:.3f}')
    if pbc is not None:
        pbc = ['TRUE' if x else 'FALSE' for x in pbc]
        remarks.append(f'PBC (A, B, C): {pbc[0]}, {pbc[1]}, {pbc[2]}')
    if remarks:
        writer.REMARK(*remarks)

    universe = mda_universe_from_nomad_atoms(atoms)
    writer.write(universe)
