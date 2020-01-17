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

from math import gcd as gcd
from typing import List, Dict, Tuple
from functools import reduce
import numpy as np
from ase import Atoms

from nomad.constants import NUMBER_TO_MASS_MAP_KG


def get_summed_atomic_mass(atomic_numbers: np.ndarray) -> float:
    """Calculates the summed atomic mass for the given atomic numbers.

    Args:
        atomic_numbers: Array of valid atomic numbers

    Returns:
        The atomic mass in kilograms.
    """
    # It is assumed that the atomic numbers are valid at this point.
    mass = np.sum(NUMBER_TO_MASS_MAP_KG[atomic_numbers])
    return mass


def get_symmetry_string(space_group: int, wyckoff_sets: Dict) -> str:
    """Used to serialize symmetry information into a string. The Wyckoff
    positions are assumed to be normalized and ordered as is the case if using
    the matid-library.

    Args:
        space_group: 3D space group number
        wyckoff_sets: Wyckoff sets that map a Wyckoff letter to related
            information

    Returns:
        A string that encodes the symmetry properties of an atomistic
        structure.
    """
    wyckoff_strings = []
    for group in wyckoff_sets:
        element = group.element
        wyckoff_letter = group.wyckoff_letter
        n_atoms = len(group.indices)
        i_string = "{} {} {}".format(element, wyckoff_letter, n_atoms)
        wyckoff_strings.append(i_string)
    wyckoff_string = ", ".join(sorted(wyckoff_strings))
    string = "{} {}".format(space_group, wyckoff_string)

    return string


def get_lattice_parameters(normalized_cell: np.ndarray) -> np.ndarray:
    """Calculate the lattice parameters for the normalized cell.

    Args:
        normalized_cell: The normalized cell as a 3x3 array. Each row is a
            basis vector.

    Returns:
        Six parameters a, b, c, alpha, beta, gamma (in this order) as a numpy
        array. Here is an explanation of each parameter:

        a = length of first basis vector
        b = length of second basis vector
        c = length of third basis vector
        alpha = angle between b and c
        beta  = angle between a and c
        gamma = angle between a and b
    """
    if normalized_cell is None:
        return None

    # Lengths
    lengths = np.linalg.norm(normalized_cell, axis=1)
    a, b, c = lengths

    # Angles
    angles = np.zeros(3)
    for i in range(3):
        j = (i + 1) % 3
        k = (i + 2) % 3
        angles[i] = np.dot(
            normalized_cell[j],
            normalized_cell[k]) / (lengths[j] * lengths[k])
    angles = np.clip(angles, -1.0, 1.0)
    alpha, beta, gamma = np.arccos(angles)

    return [a, b, c, alpha, beta, gamma]


def get_hill_decomposition(atom_labels: np.ndarray, reduced: bool = False) -> Tuple[List[str], List[int]]:
    """Given a list of atomic labels, returns the chemical formula using the
    Hill system (https://en.wikipedia.org/wiki/Hill_system).

    Args:
        atom_labels: Atom labels.
        reduced: Whether to divide the number of atoms by the greatest common
            divisor

    Returns:
        An ordered list of chemical symbols and the corresponding counts.
    """
    # Count occurancy of elements
    names = []
    counts = []
    unordered_names, unordered_counts = np.unique(atom_labels, return_counts=True)
    element_count_map = dict(zip(unordered_names, unordered_counts))

    # Apply basic Hill system:
    # 1. Is Carbon part of the system?
    if "C" in element_count_map:
        names.append("C")
        counts.append(element_count_map["C"])
        del element_count_map['C']

        # 1a. add hydrogren
        if "H" in element_count_map:
            names.append("H")
            counts.append(element_count_map["H"])
            del element_count_map["H"]

    # 2. all remaining elements in alphabetic order
    for element in sorted(element_count_map):
        names.append(element)
        counts.append(element_count_map[element])

    # 3. Binary ionic compounds: cation first, anion second
    # If any of the most electronegative elements is first
    # by alphabetic order, we move it to second
    if len(counts) == 2:
        order = {
            "F": 1,
            "O": 2,
            "N": 3,
            "Cl": 4,
            "Br": 5,
            "C": 6,
            "Se": 7,
            "S": 8,
            "I": 9,
            "As": 10,
            "H": 11,
            "P": 12,
            "Ge": 13,
            "Te": 14,
            "B": 15,
            "Sb": 16,
            "Po": 17,
            "Si": 18,
            "Bi": 19
        }
        if (names[0] in order):
            if (names[1] in order):
                if(order[names[0]] < order[names[1]]):
                    # For non-metals:
                    # Swap symbols and counts if first element
                    # is more electronegative than the second one,
                    # because the more electronegative element is the anion
                    names[0], names[1] = names[1], names[0]
                    counts[0], counts[1] = counts[1], counts[0]
            else:
                # Swap symbols and counts always if second element
                # is any other element,i.e.,
                # put non-metal last because it is the anion
                names[0], names[1] = names[1], names[0]
                counts[0], counts[1] = counts[1], counts[0]

    # TODO: implement all further exceptions regarding ordering
    #       in chemical formulas:
    #         - ionic compounds (ordering wrt to ionization)
    #         - oxides, acids, hydroxides...

    # Reduce if requested
    if reduced:
        greatest_common_divisor = reduce(gcd, counts)
        counts = np.array(counts) / greatest_common_divisor

    return names, counts


def get_formula_string(symbols: List[str], counts: List[int]) -> str:
    """Used to form a single formula string from a list of chemical speices and
    their counts.

    Args:
        symbols: List of chemical species
        counts: List of chemical species occurences

    Returns:
        The formula as a string.
    """
    formula = ""
    for symbol, count in zip(symbols, counts):
        if count > 1:
            formula += "%s%d" % (symbol, count)
        else:
            formula += symbol
    return formula


def find_vacuum_directions(system: Atoms, threshold: float) -> np.array:
    """Searches for vacuum gaps that are separating the periodic copies.

    Args:
        system: The structure to analyze
        threshold: Vacuum threshold in angstroms

    Returns:
        np.ndarray: An array with a boolean for each lattice basis
        direction indicating if there is enough vacuum to separate the
        copies in that direction.
    """
    rel_pos = system.get_scaled_positions()
    pbc = system.get_pbc()

    # Find the maximum vacuum gap for all basis vectors
    gaps = np.empty(3, dtype=bool)
    for axis in range(3):
        if not pbc[axis]:
            gaps[axis] = True
            continue
        comp = rel_pos[:, axis]
        ind = np.sort(comp)
        ind_rolled = np.roll(ind, 1, axis=0)
        distances = ind - ind_rolled

        # The first distance is from first to last, so it needs to be
        # wrapped around
        distances[0] += 1

        # Find maximum gap in cartesian coordinates
        max_gap = np.max(distances)
        basis = system.get_cell()[axis, :]
        max_gap_cartesian = np.linalg.norm(max_gap * basis)
        has_vacuum_gap = max_gap_cartesian >= threshold
        gaps[axis] = has_vacuum_gap

    return gaps
