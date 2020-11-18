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

import io
import re
import uuid
import json

import numpy as np
import requests
import ase
import bs4
import matid

from nomad import atomutils, config


def write_prototype_data_file(aflow_prototypes: dict, filepath) -> None:
    '''Writes the prototype data file in a compressed format to a python
    module.

    Args:
        aflow_prototypes
    '''
    class NoIndent(object):
        def __init__(self, value):
            self.value = value

    class NoIndentEncoder(json.JSONEncoder):
        '''A custom JSON encoder that can pretty-print objects wrapped in the
        NoIndent class.
        '''
        def __init__(self, *args, **kwargs):
            super(NoIndentEncoder, self).__init__(*args, **kwargs)
            self.kwargs = dict(kwargs)
            del self.kwargs['indent']
            self._replacement_map = {}

        def default(self, o):  # pylint: disable=E0202
            if isinstance(o, NoIndent):
                key = uuid.uuid4().hex
                self._replacement_map[key] = json.dumps(o.value, **self.kwargs)
                return "@@%s@@" % (key,)
            else:
                return super(NoIndentEncoder, self).default(o)

        def encode(self, o):
            result = super(NoIndentEncoder, self).encode(o)
            for k, v in self._replacement_map.items():
                result = result.replace('"@@%s@@"' % (k,), v)
            return result

    prototype_dict = aflow_prototypes["prototypes_by_spacegroup"]
    for prototypes in prototype_dict.values():
        for prototype in prototypes:
            # Save the information back in a prettified form
            prototype["atom_positions"] = NoIndent(prototype["atom_positions"])
            prototype["atom_labels"] = NoIndent(prototype["atom_labels"])
            prototype["lattice_vectors"] = NoIndent(prototype["lattice_vectors"])
            try:
                prototype["normalized_wyckoff_matid"] = NoIndent(prototype["normalized_wyckoff_matid"])
            except KeyError:
                pass

    # Save the updated data
    with io.open(filepath, "w", encoding="utf8") as f:
        json_dump = json.dumps(aflow_prototypes, ensure_ascii=False, indent=4, sort_keys=True, cls=NoIndentEncoder)
        json_dump = re.sub(r"\"(-?\d+(?:[\.,]\d+)?)\"", r'\1', json_dump)  # Removes quotes around numbers
        f.write("# -*- coding: utf-8 -*-\naflow_prototypes = {}\n".format(json_dump))


def update_prototypes(ctx, filepath, matches_only):

    if matches_only:
        from nomad.aflow_prototypes import aflow_prototypes
    else:
        # The basic AFLOW prototype data is available in a Javascript file. Here we
        # retrieve it and read only the prototype list from it.
        prototypes_file_url = 'http://aflowlib.org/CrystalDatabase/js/table_sort.js'
        r = requests.get(prototypes_file_url, allow_redirects=True)
        datastring = r.content.decode("utf-8")
        datastring = datastring.split('];')[0]
        datastring = datastring.split('= [')[1]
        data = json.loads('[' + datastring + ']')

        newdictarray = []
        n_prototypes = 0
        n_missing = 0
        for protodict in data:
            n_prototypes += 1
            newdict = {}

            # Make prototype plaintext
            prototype = bs4.BeautifulSoup(protodict["Prototype"], "html5lib").getText()

            # Add to new dictionary
            newdict['Notes'] = protodict['Notes']
            newdict['Prototype'] = prototype
            newdict['Space Group Symbol'] = protodict['Space Group Symbol']
            newdict['Space Group Number'] = protodict['Space Group Number']
            newdict['Pearsons Symbol'] = protodict['Pearson Symbol']
            newdict['Strukturbericht Designation'] = protodict['Strukturbericht Designation']
            newdict['aflow_prototype_id'] = protodict['AFLOW Prototype']
            newdict['aflow_prototype_url'] = 'http://www.aflowlib.org/CrystalDatabase/' + protodict['href'][2:]

            # Download cif or poscar if possible make ASE ase.Atoms object if possible
            # to obtain labels, positions, cell
            cifurl = 'http://www.aflowlib.org/CrystalDatabase/CIF/' + protodict['href'][2:-5] + '.cif'
            r = requests.get(cifurl, allow_redirects=True)
            cif_str = r.content.decode("utf-8")
            cif_file = io.StringIO()
            cif_file.write(cif_str)
            cif_file.seek(0)
            try:
                atoms = ase.io.read(cif_file, format='cif')
            except Exception:
                print("Error in getting prototype structure from CIF: {}", format(cifurl))
                # Then try to get structure from POSCAR
                try:
                    poscarurl = 'http://www.aflowlib.org/CrystalDatabase/POSCAR/' + protodict['href'][2:-5] + '.poscar'
                    r = requests.get(poscarurl, allow_redirects=True)
                    poscar_str = r.content.decode("utf-8")
                    poscar_file = io.StringIO()
                    poscar_file.write(poscar_str)
                    poscar_file.seek(0)
                    atoms = ase.io.read(poscar_file, format='vasp')
                except Exception:
                    print("Error in getting prototype structure from POSCAR: {}".format(poscarurl))
                    print("Could not read prototype structure from CIF or POSCAR file for prototype: {}, {}, ".format(prototype, newdict['aflow_prototype_url']))
                    n_missing += 1
                    continue

            atom_positions = atoms.get_positions()
            atom_labels = atoms.get_chemical_symbols()
            cell = atoms.get_cell()

            newdict['lattice_vectors'] = cell.tolist()
            newdict['atom_positions'] = atom_positions.tolist()
            newdict['atom_labels'] = atom_labels
            newdictarray.append(newdict)

            print("Processed: {}".format(len(newdictarray)))

        # Sort prototype dictionaries by spacegroup and make dictionary
        structure_types_by_spacegroup = {}
        for i_sg in range(1, 231):
            protos_sg = []
            for newdict in newdictarray:
                if newdict['Space Group Number'] == i_sg:
                    protos_sg.append(newdict)
            structure_types_by_spacegroup[i_sg] = protos_sg

        # Wrap in a dictionary that can hold other data, e.g. the symmemtry tolerance parameter.
        aflow_prototypes = {
            "prototypes_by_spacegroup": structure_types_by_spacegroup
        }
        print(
            "Extracted latest AFLOW prototypes online. Total number of "
            "successfully fetched prototypes: {}, missing: {}"
            .format(n_prototypes, n_missing)
        )

    # Update matches
    n_prototypes = 0
    n_failed = 0
    n_unmatched = 0
    prototype_dict = aflow_prototypes["prototypes_by_spacegroup"]

    for aflow_spg_number, prototypes in prototype_dict.items():
        n_prototypes += len(prototypes)
        for prototype in prototypes:

            # Read prototype structure
            pos = np.array(prototype["atom_positions"])
            labels = prototype["atom_labels"]
            cell = np.array(prototype["lattice_vectors"])
            atoms = ase.Atoms(
                symbols=labels,
                positions=pos,
                cell=cell,
                pbc=True
            )

            # Try to first see if the space group can be matched with the one in AFLOW
            try:
                symm = matid.SymmetryAnalyzer(atoms, config.normalize.prototype_symmetry_tolerance)
                spg_number = symm.get_space_group_number()
                wyckoff_matid = symm.get_wyckoff_letters_conventional()
                norm_system = symm.get_conventional_system()
            except Exception:
                n_failed += 1
            else:
                # If the space group is matched, add the MatID normalized Wyckoff
                # letters to the data.
                if spg_number == aflow_spg_number:
                    atomic_numbers = norm_system.get_atomic_numbers()
                    normalized_wyckoff_matid = atomutils.get_normalized_wyckoff(atomic_numbers, wyckoff_matid)
                    prototype["normalized_wyckoff_matid"] = normalized_wyckoff_matid
                else:
                    n_unmatched += 1
    print(
        "Updated matches in AFLOW prototype library. Total number of "
        "prototypes: {}, unmatched: {}, failed: {}"
        .format(n_prototypes, n_unmatched, n_failed)
    )

    # Write data file to the specified path
    write_prototype_data_file(aflow_prototypes, filepath)
