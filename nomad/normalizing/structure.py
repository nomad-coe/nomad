# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import io
import functools
import fractions
import json
import uuid

from matid import SymmetryAnalyzer
import numpy as np
from ase import Atoms
from nomad.normalizing.data.aflow_prototypes import aflow_prototypes
from nomad import config

# The AFLOW symmetry information is checked once on import
old_symmetry_tolerance = aflow_prototypes["matid_symmetry_tolerance"]
symmetry_tolerance = config.normalize.symmetry_tolerance
if old_symmetry_tolerance != symmetry_tolerance:
    raise AssertionError(
        "The AFLOW prototype information is outdated due to changed "
        "tolerance for symmetry detection. Please update the AFLOW "
        "prototype information by running once the function "
        "'update_aflow_prototype_information'."
    )


def get_normalized_wyckoff(atomic_number: np.array, wyckoff: np.array) -> dict:
    """Returns a normalized Wyckoff sequence for the given atomic numbers and
    corresponding wyckoff letters. In a normalized sequence the chemical
    species are "anonymized" by replacing the with upper case alphabets.

    Args:
        atomic_number: Array of atomic numbers.
        wyckoff: Array of Wyckoff letters as strings.

    Returns:
    """
    atomCount = {}
    for nr in atomic_number:
        atomCount[nr] = atomCount.get(nr, 0) + 1
    wycDict = {}

    for i, wk in enumerate(wyckoff):
        oldVal = wycDict.get(wk, {})
        nr = atomic_number[i]
        oldVal[nr] = oldVal.get(nr, 0) + 1
        wycDict[wk] = oldVal
    sortedWyc = list(wycDict.keys())
    sortedWyc.sort()

    def cmpp(a, b):
        return ((a < b) - (a > b))

    def compareAtNr(at1, at2):
        c = cmpp(atomCount[at1], atomCount[at2])
        if (c != 0):
            return c
        for wk in sortedWyc:
            p = wycDict[wk]
            c = cmpp(p.get(at1, 0), p.get(at2, 0))
            if c != 0:
                return c
        return 0

    sortedAt = list(atomCount.keys())
    sortedAt.sort(key=functools.cmp_to_key(compareAtNr))
    standardAtomNames = {}
    for i, at in enumerate(sortedAt):
        standardAtomNames[at] = ("X_%d" % i)
    standardWyc = {}
    for wk, ats in wycDict.items():
        stdAts = {}
        for at, count in ats.items():
            stdAts[standardAtomNames[at]] = count
        standardWyc[wk] = stdAts
    if standardWyc:
        counts = [c for x in standardWyc.values() for c in x.values()]
        gcd = counts[0]
        for c in counts[1:]:
            gcd = fractions.gcd(gcd, c)
        if gcd != 1:
            for wk, d in standardWyc.items():
                for at, c in d.items():
                    d[at] = c // gcd
    return standardWyc


def search_aflow_prototype(space_group: int, norm_wyckoff: dict) -> dict:
    """Searches the AFLOW prototype library for a match for the given space
    group and normalized Wyckoff sequence. The normalized Wyckoff sequence is
    assumed to come from the MatID symmetry routine.

    Currently only contains Part I of the prototype library (M. J. Mehl, D.
    Hicks, C. Toher, O. Levy, R. M. Hanson, G. L. W. Hart, and S. Curtarolo,
    The AFLOW Library of Crystallographic Prototypes: Part 1, Comp. Mat. Sci.
    136, S1-S828 (2017), 10.1016/j.commatsci.2017.01.017)

    Args:
        space_group_number: Space group number
        norm_wyckoff: Normalized Wyckoff occupations

    Returns:
        Dictionary containing the AFLOW prototype information.
    """
    structure_type_info = None
    type_descriptions = aflow_prototypes["prototypes_by_spacegroup"].get(space_group, [])
    for type_description in type_descriptions:
        current_norm_wyckoffs = type_description.get("normalized_wyckoff_matid")
        if current_norm_wyckoffs and current_norm_wyckoffs == norm_wyckoff:
            structure_type_info = type_description
            break
    return structure_type_info


def update_aflow_prototype_information(filepath: str) -> None:
    """Used to update the AFLOW prototype information. Creates a new python
    module with updated symmetry tolerance parameter and the wyckoff positions
    as detected by MatID.

    This function is relative heavy and should only be run if the symmetry
    tolerance has been changed or the symmetry detection routine has been
    updated.

    Args:
        filepath: Path to the python file in which the new symmetry information
        will be written.
    """
    class NoIndent(object):
        def __init__(self, value):
            self.value = value

    class NoIndentEncoder(json.JSONEncoder):
        """A custom JSON encoder that can pretty-print objects wrapped in the
        NoIndent class.
        """
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

    n_prototypes = 0
    n_failed = 0
    n_unmatched = 0
    prototype_dict = aflow_prototypes["prototypes_by_spacegroup"]
    for aflow_spg_number, prototypes in prototype_dict.items():
        n_prototypes += len(prototypes)
        for prototype in prototypes:

            # Read prototype structure
            pos = np.array(prototype["atom_positions"]) * 1E10
            labels = prototype["atom_labels"]
            cell = np.array(prototype["lattice_vectors"]) * 1E10
            atoms = Atoms(
                symbols=labels,
                positions=pos,
                cell=cell,
                pbc=True
            )

            # Try to first see if the space group can be matched with the one in AFLOW
            tolerance = config.normalize.symmetry_tolerance
            try:
                symm = SymmetryAnalyzer(atoms, tolerance)
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
                    normalized_wyckoff_matid = get_normalized_wyckoff(atomic_numbers, wyckoff_matid)
                    prototype["normalized_wyckoff_matid"] = NoIndent(normalized_wyckoff_matid)
                else:
                    n_unmatched += 1

            # Save the information back in a prettified form
            prototype["atom_positions"] = NoIndent(prototype["atom_positions"])
            prototype["atom_labels"] = NoIndent(prototype["atom_labels"])
            prototype["lattice_vectors"] = NoIndent(prototype["lattice_vectors"])
            try:
                prototype["normalized_wyckoff"] = NoIndent(prototype["normalized_wyckoff"])
            except KeyError:
                pass

    print(f"Updated AFLOW prototype library. Total number of prototypes: {n_prototypes}, unmatched: {n_unmatched}, failed: {n_failed}")

    # Save the updated data
    with io.open(filepath, "w", encoding="utf8") as f:
        json_dump = json.dumps(aflow_prototypes, ensure_ascii=False, indent=4, cls=NoIndentEncoder)
        f.write("# -*- coding: utf-8 -*-\naflow_prototypes = {}\n".format(json_dump))
