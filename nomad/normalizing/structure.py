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


def get_summed_atomic_mass(atomic_numbers: np.ndarray) -> float:
    """Calculates the summed atomic mass for the given atomic numbers.

    Args:
        atomic_numbers: Array of valid atomic numbers

    Returns:
        The atomic mass in kilograms.
    """
    # It is assumed that the atomic numbers are valid at this point.
    masses = np.array(number_to_mass_map)
    mass = np.sum(masses[atomic_numbers])
    return mass


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


def get_lattice_parameters(normalized_cell: np.ndarray) -> np.ndarray:
    """Calculate the lattice parameters for the normalized cell.

    Args:
        normalized_cell: The normalized cell as a 2D array. Each row is a basis
            vector.

    Returns:
        Six parameters a, b, c, alpha, beta, gamma (in this order) as a numpy array.
        Here is an explanation of each parameter:

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


"""List of atomic masses (natural isotope dist.) in order, atomic mass units.
origin: phonopy-1.11.2.25.tar.gz:phonopy/structure/atoms.py:atom_data
"""
number_to_mass_map_amu = np.array([
    np.nan,       # 0
    1.00794,      # 1
    4.002602,     # 2
    6.941,        # 3
    9.012182,     # 4
    10.811,       # 5
    12.0107,      # 6
    14.0067,      # 7
    15.9994,      # 8
    18.9984032,   # 9
    20.1797,      # 10
    22.98976928,  # 11
    24.3050,      # 12
    26.9815386,   # 13
    28.0855,      # 14
    30.973762,    # 15
    32.065,       # 16
    35.453,       # 17
    39.948,       # 18
    39.0983,      # 19
    40.078,       # 20
    44.955912,    # 21
    47.867,       # 22
    50.9415,      # 23
    51.9961,      # 24
    54.938045,    # 25
    55.845,       # 26
    58.933195,    # 27
    58.6934,      # 28
    63.546,       # 29
    65.38,        # 30
    69.723,       # 31
    72.64,        # 32
    74.92160,     # 33
    78.96,        # 34
    79.904,       # 35
    83.798,       # 36
    85.4678,      # 37
    87.62,        # 38
    88.90585,     # 39
    91.224,       # 40
    92.90638,     # 41
    95.96,        # 42
    98.9062,      # 43 - NIST
    101.07,       # 44
    102.90550,    # 45
    106.42,       # 46
    107.8682,     # 47
    112.411,      # 48
    114.818,      # 49
    118.710,      # 50
    121.760,      # 51
    127.60,       # 52
    126.90447,    # 53
    131.293,      # 54
    132.9054519,  # 55
    137.327,      # 56
    138.90547,    # 57
    140.116,      # 58
    140.90765,    # 59
    144.242,      # 60
    145,          # 61 most stable isotope
    150.36,       # 62
    151.964,      # 63
    157.25,       # 64
    158.92535,    # 65
    162.500,      # 66
    164.93032,    # 67
    167.259,      # 68
    168.93421,    # 69
    173.054,      # 70
    174.9668,     # 71
    178.49,       # 72
    180.94788,    # 73
    183.84,       # 74
    186.207,      # 75
    190.23,       # 76
    192.217,      # 77
    195.084,      # 78
    196.966569,   # 79
    200.59,       # 80
    204.3833,     # 81
    207.2,        # 82
    208.98040,    # 83
    209,          # 84 - NIST
    210,          # 85 - NIST
    222,          # 86 - NIST
    223,          # 87 - NIST
    226,          # 88 - NIST
    227,          # 89 - NIST
    232.03806,    # 90
    231.03588,    # 91
    238.02891,    # 92
    237,          # 93 - NIST
    244,          # 94 - NIST
    243,          # 95 - most stable isotope
    247,          # 96 - most stable isotope
    247,          # 97 - most stable isotope
    251,          # 98 - most stable isotope
    252,          # 99 - most stable isotope
    257,          # 100 - most stable isotope
    258,          # 101 - most stable isotope
    259,          # 102 - most stable isotope
    262,          # 103 - most stable isotope
    261,          # 104 - most stable isotope
    262,          # 105 - most stable isotope
    266,          # 106 - most stable isotope
    264,          # 107 - most stable isotope
    277,          # 108 - most stable isotope
    268,          # 109 - most stable isotope
    281,          # 110
    282,          # 111
    285,          # 112
    286,          # 113
    289,          # 114
    290,          # 115
    293,          # 116
    294,          # 117
    294,          # 118
])


"""List of atomic masses (natural isotope dist.) in order, SI units (kg)"""
number_to_mass_map = number_to_mass_map_amu * 1.6605389e-27
