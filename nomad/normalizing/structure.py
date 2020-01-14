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

from typing import Dict
import numpy as np


def create_symmetry_string(space_group: int, wyckoff_sets: Dict) -> str:
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


def get_lattice_parameters(normalized_cell):
    """Calculate the lattice parameters for the normalized cell.

    :param normalized_cell: The normalized cell as a 2D array. Each row is a
        basis vector.
    :type normalized_cell: numpy.array

    :return: Six parameters a, b, c, alpha, beta, gamma (in this order) as a
        list. Here is an explanation of each parameter:

        #. a = length of first basis vector
        #. b = length of second basis vector
        #. c = length of third basis vector
        #. alpha = angle between b and c
        #. beta  = angle between a and c
        #. gamma = angle between a and b

    :rtype: numpy.array
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
