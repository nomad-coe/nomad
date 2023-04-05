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

from __future__ import annotations

import warnings
import functools
import itertools
from collections import namedtuple
from array import array
from itertools import chain
import math
import re
from string import ascii_uppercase
from functools import reduce
from typing import List, Dict, Tuple, Any, Union, Iterable, cast, Callable, TYPE_CHECKING
import logging
from nptyping import NDArray

import numpy as np
from scipy import sparse
from scipy.spatial import Voronoi  # pylint: disable=no-name-in-module
from scipy.stats import linregress
from ase.utils import pbc2pbc
import ase.geometry
import ase.data
from ase import Atoms
from ase.formula import Formula as ASEFormula
from pymatgen.core.periodic_table import get_el_sp
from pymatgen.core import Composition
import MDAnalysis
from MDAnalysis.core.topology import Topology
from MDAnalysis.core._get_readers import get_reader_for
from MDAnalysis.core.universe import Universe
import MDAnalysis.analysis.rdf as MDA_RDF

from nomad.aflow_prototypes import aflow_prototypes
from nomad.constants import atomic_masses
from nomad.units import ureg

if TYPE_CHECKING:
    from nomad.datamodel.results import Material


def get_summed_atomic_mass(atomic_numbers: NDArray[Any]) -> float:
    '''
    Calculates the summed atomic mass for the given atomic numbers.

    Args:
        atomic_numbers: Array of valid atomic numbers

    Returns:
        The atomic mass in kilograms.
    '''
    # It is assumed that the atomic numbers are valid at this point.
    mass = np.sum(atomic_masses[atomic_numbers])
    return mass


def get_volume(basis: NDArray[Any]) -> float:
    '''
    Calculates the volume of the given parallelepiped.

    Args:
        basis: 3x3 matrix with basis vectors of a parallellepiped as rows.

    Returns:
        Volume of the parallelepiped defined by the basis.
    '''
    return np.abs(np.linalg.det(basis))


def is_valid_basis(basis: NDArray[Any]) -> bool:
    '''
    Checks if the given set of basis vectors are valid. Currently does not
    check for linear independence, only for empty rows.

    Args:
        basis: 3x3 matrix with basis vectors as rows.

    Returns:
        True if the basis is valid, False otherwise.
    '''
    if basis is None:
        return False
    for row in np.asarray(basis):
        if not np.any(row):
            return False
    return True


def wrap_positions(
        positions: NDArray[Any],
        cell: NDArray[Any] = None,
        pbc: Union[bool, NDArray[Any]] = True,
        center: NDArray[Any] = [0.5, 0.5, 0.5],
        eps: float = 1e-12) -> NDArray[Any]:
    '''
    Wraps the given position so that they are within the unit cell. If no
    cell is given, scaled positions are assumed. For wrapping cartesian
    positions you also need to provide the cell.

    Args:
        positions: Positions of the atoms. Accepts both scaled and
            cartesian positions.
        cell: Lattice vectors for wrapping cartesian positions.
        pbc: For each axis in the unit cell decides whether the positions will
            be wrapped along this axis.
        center: The position in fractional coordinates that the wrapped
            positions will be nearest possible to.
        eps: Small number to prevent slightly negative coordinates from being
            wrapped.
    '''
    if not hasattr(center, '__len__'):
        center = (center,) * 3

    pbc = pbc2pbc(pbc)
    shift = np.asarray(center) - 0.5 - eps

    # Don't change coordinates when pbc is False
    shift[np.logical_not(pbc)] = 0.0

    if cell is None:
        fractional = positions
    else:
        fractional = to_scaled(positions, cell)
    fractional -= shift

    for i, periodic in enumerate(pbc):
        if periodic:
            fractional[:, i] %= 1.0
        fractional[:, i] += shift[i]
    if cell is not None:
        return np.dot(fractional, cell)
    else:
        return fractional


def chemical_symbols(atomic_numbers: Iterable[int]) -> List[str]:
    '''
    Converts atomic numbers to chemical_symbols.

    Args:
        atomic_numbers: The atomic numbers to convert.

    Returns:
        Array of chemical symbols.
    '''
    return [ase.data.chemical_symbols[x] for x in atomic_numbers]


def to_scaled(
        positions: NDArray[Any],
        cell: NDArray[Any] = None) -> NDArray[Any]:
    '''
    Converts cartesian positions into scaled position one using the given
    cell lattice vectors as a basis.

    Args:
        positions: Scaled positions.
        cell: Lattice vectors.

    Returns:
        The given positions in scaled coordinates.
    '''
    return np.linalg.solve(complete_cell(cell).T, positions.T).T


def to_cartesian(
        positions: NDArray[Any],
        cell: NDArray[Any] = None) -> NDArray[Any]:
    '''
    Converts scaled positions into cartesian one using the given cell
    lattice vectors as a basis.

    Args:
        positions: Scaled positions.
        cell: Lattice vectors.

    Returns:
        The given positions in cartesian coordinates.
    '''
    cartesian_positions = np.dot(positions, complete_cell(cell))
    return cartesian_positions


def complete_cell(cell: NDArray[Any]) -> NDArray[Any]:
    '''
    Creates placeholder axes for cells with zero-dimensional lattice vectors
    in order to do linear algebra.

    Args:
        cell: Lattice vectors.

    Returns:
        The given cell with zero-dimensional lattice vectors filled with
        placeholder axes.
    '''
    return ase.geometry.complete_cell(cell)


def reciprocal_cell(cell: NDArray[Any]) -> NDArray[Any]:
    '''
    Returns the reciprocal cell without the factor or 2*Pi.

    Args:
        cell: Lattice vectors.

    Returns:
        Reciprocal cell as a 3x3 array.
    '''
    return np.linalg.pinv(cell).transpose()


def find_match(pos: NDArray[Any], positions: NDArray[Any], eps: float) -> Union[int, None]:
    '''
    Attempts to find a position within a larger list of positions.

    Args:
        pos: The point to search for
        positions: The points within which the search is performed.
        eps: Match tolerance.

    Returns:
        Index of the matched position or None if match not found.
    '''
    displacements = positions - pos
    distances = np.linalg.norm(displacements, axis=1)
    min_arg = np.argmin(distances)
    min_value = distances[min_arg]
    if min_value <= eps:
        return cast(int, min_arg)
    else:
        return None


def cellpar_to_cell(cellpar: NDArray[Any], ab_normal: NDArray[Any] = [0, 0, 1], a_direction: NDArray[Any] = None, degrees=False) -> NDArray[Any]:
    '''
    Creates a 3x3 cell from the given lattice_parameters.

    The returned cell is orientated such that a and b are normal to `ab_normal`
    and a is parallel to the projection of `a_direction` in the a-b plane.

    Default `a_direction` is (1,0,0), unless this is parallel to `ab_normal`,
    in which case default `a_direction` is (0,0,1).

    The returned cell has the vectors va, vb and vc along the rows. The cell
    will be oriented such that va and vb are normal to `ab_normal` and va will
    be along the projection of `a_direction` onto the a-b plane.

    Args:
        cellpar: Six lattice parameters: [a, b, c, alpha, beta, gamma].
        The following typical convention is used:

            a = length of first basis vector
            b = length of second basis vector
            c = length of third basis vector
            alpha = angle between b and c in radians
            beta  = angle between a and c in radians
            gamma = angle between a and b in radians
        degrees: Use degrees in place of radians.

    Returns:
        Six parameters (in this order) as a numpy
        array. Here is an explanation of each parameter:
    '''
    if not degrees:
        cellpar[3:6] *= 180.0 / np.pi

    return ase.geometry.cell.cellpar_to_cell(cellpar, ab_normal, a_direction)


def cell_to_cellpar(cell: NDArray[Any], degrees=False) -> NDArray[Any]:
    '''
    Returns lattice parameters for the given cell.

    Args:
        normalized_cell: The normalized cell as a 3x3 array. Each row is a
            basis vector.
        degrees: Use degrees in place of radians.

    Returns:
        Six parameters [a, b, c, alpha, beta, gamma]. The following typical
        convention ik used:

            a = length of first basis vector
            b = length of second basis vector
            c = length of third basis vector
            alpha = angle between b and c in radians
            beta  = angle between a and c in radians
            gamma = angle between a and b in radians
    '''
    # Lengths
    lengths = np.linalg.norm(cell, axis=1)

    # Angles
    angles = np.zeros(3)
    for i in range(3):
        j = (i + 1) % 3
        k = (i + 2) % 3
        angles[i] = np.dot(
            cell[j],
            cell[k]) / (lengths[j] * lengths[k])
    angles = np.arccos(np.clip(angles, -1.0, 1.0))
    if degrees:
        angles *= 180.0 / np.pi

    return np.concatenate((lengths, angles), axis=0)


def get_symmetry_string(space_group: int, wyckoff_sets: List[Any], is_2d: bool = False) -> str:
    '''
    Used to serialize symmetry information into a string. The Wyckoff
    positions are assumed to be normalized and ordered as is the case if using
    the matid-library.

    The symmetry analysis of 2D structures is run by artificially extending the
    system in the non-periodic direction. Due to this the symmetry information
    may overlap with real 3D structures unless we include a distinguishing
    label for 2D structures.

    Args:
        space_group: 3D space group number
        wyckoff_sets: Wyckoff sets that map a Wyckoff letter to related
            information
        is_2d: Whether the symmetry information is analyzed from a 2D
            structure. If true, a prefix is added to the string to distinguish
            2D from 3D.

    Returns:
        A string that encodes the symmetry properties of an atomistic
        structure.
    '''
    wyckoff_strings = []
    for group in wyckoff_sets:
        element = group.element
        wyckoff_letter = group.wyckoff_letter
        n_atoms = len(group.indices)
        i_string = '{} {} {}'.format(element, wyckoff_letter, n_atoms)
        wyckoff_strings.append(i_string)
    wyckoff_string = ', '.join(sorted(wyckoff_strings))
    if is_2d:
        string = '2D {} {}'.format(space_group, wyckoff_string)
    else:
        string = '{} {}'.format(space_group, wyckoff_string)

    return string


def get_hill_decomposition(atom_labels: NDArray[Any], reduced: bool = False) -> Tuple[List[str], List[int]]:
    '''
    Given a list of atomic labels, returns the chemical formula using the
    Hill system (https://en.wikipedia.org/wiki/Hill_system) with an exception
    for binary ionic compounds where the cation is always given first.

    Args:
        atom_labels: Atom labels.
        reduced: Whether to divide the number of atoms by the greatest common
            divisor

    Returns:
        An ordered list of chemical symbols and the corresponding counts.
    '''
    # Count occurancy of elements
    names = []
    counts = []
    unordered_names, unordered_counts = np.unique(atom_labels, return_counts=True)
    element_count_map = dict(zip(unordered_names, unordered_counts))

    # Apply basic Hill system:
    # 1. Is Carbon part of the system?
    if 'C' in element_count_map:
        names.append('C')
        counts.append(element_count_map['C'])
        del element_count_map['C']

        # 1a. add hydrogren
        if 'H' in element_count_map:
            names.append('H')
            counts.append(element_count_map['H'])
            del element_count_map['H']

    # 2. all remaining elements in alphabetic order
    for element in sorted(element_count_map):
        names.append(element)
        counts.append(element_count_map[element])

    # 3. Binary ionic compounds: cation first, anion second
    # If any of the most electronegative elements is first
    # by alphabetic order, we move it to second
    if len(counts) == 2 and names != ['C', 'H']:
        order = {
            'F': 1,
            'O': 2,
            'N': 3,
            'Cl': 4,
            'Br': 5,
            'C': 6,
            'Se': 7,
            'S': 8,
            'I': 9,
            'As': 10,
            'H': 11,
            'P': 12,
            'Ge': 13,
            'Te': 14,
            'B': 15,
            'Sb': 16,
            'Po': 17,
            'Si': 18,
            'Bi': 19
        }
        if (names[0] in order):
            if (names[1] in order):
                if (order[names[0]] < order[names[1]]):
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
        greatest_common_divisor = reduce(math.gcd, counts)
        counts = np.array(counts) / greatest_common_divisor

    return names, counts


def get_formula_string(symbols: Iterable[str], counts: Iterable[int]) -> str:
    '''
    Used to form a single formula string from a list of chemical species and
    their counts.

    Args:
        symbols: List of chemical species
        counts: List of chemical species occurences

    Returns:
        The formula as a string.
    '''
    formula = ''
    for symbol, count in zip(symbols, counts):
        if count > 1:
            formula += '%s%d' % (symbol, count)
        else:
            formula += symbol
    return formula


def get_normalized_wyckoff(atomic_numbers: NDArray[Any], wyckoff_letters: NDArray[Any]) -> Dict[str, Dict[str, int]]:
    '''
    Returns a normalized Wyckoff sequence for the given atomic numbers and
    corresponding wyckoff letters. In a normalized sequence the chemical
    species are 'anonymized' by replacing them with upper case alphabets.

    Args:
        atomic_numbers: Array of atomic numbers.
        wyckoff_letters: Array of Wyckoff letters as strings.

    Returns:
        Returns a dictionary that maps each present Wyckoff letter to a
        dictionary. The dictionary contains the number of atoms for each
        species, where the species names have been anomymized in the form
        'X_<index>'.
    '''
    # Count the occurrence of each chemical species
    atom_count: Dict[int, int] = {}
    for atomic_number in atomic_numbers:
        atom_count[atomic_number] = atom_count.get(atomic_number, 0) + 1

    # Form a dictionary that maps Wyckoff letters to a dictionary with the
    # number of atoms of that Wyckoff letter for each atomic number
    wyc_dict: dict = {}
    for i, wyckoff_letter in enumerate(wyckoff_letters):
        old_val = wyc_dict.get(wyckoff_letter, {})
        atomic_number = atomic_numbers[i]
        old_val[atomic_number] = old_val.get(atomic_number, 0) + 1
        wyc_dict[wyckoff_letter] = old_val
    sorted_wyckoff_letters = list(wyc_dict.keys())
    sorted_wyckoff_letters.sort()

    # Anonymize the atomic species to X_<index>, where the index is calculated
    # by ordering the species.
    def compare_atomic_number(at1, at2):
        def cmpp(a, b):
            return ((a < b) - (a > b))

        c = cmpp(atom_count[at1], atom_count[at2])
        if (c != 0):
            return c
        for wyckoff_letter in sorted_wyckoff_letters:
            p = wyc_dict[wyckoff_letter]
            c = cmpp(p.get(at1, 0), p.get(at2, 0))
            if c != 0:
                return c
        return 0
    sorted_species = list(atom_count.keys())
    sorted_species.sort(key=functools.cmp_to_key(compare_atomic_number))
    standard_atom_names = {}
    for i, at in enumerate(sorted_species):
        standard_atom_names[at] = ('X_%d' % i)

    # Rename with anonymized species labels
    standard_wyc: dict = {}
    for wk, ats in wyc_dict.items():
        std_ats = {}
        for at, count in ats.items():
            std_ats[standard_atom_names[at]] = count
        standard_wyc[wk] = std_ats

    # Divide atom counts with greatest common divisor
    if standard_wyc:
        counts = [c for x in standard_wyc.values() for c in x.values()]
        gcd = counts[0]
        for c in counts[1:]:
            gcd = math.gcd(gcd, c)
        if gcd != 1:
            for d in standard_wyc.values():
                for at, c in d.items():
                    d[at] = c // gcd

    return standard_wyc


def search_aflow_prototype(space_group: int, norm_wyckoff: dict) -> dict:
    '''
    Searches the AFLOW prototype library for a match for the given space
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
    '''
    structure_type_info = None
    type_descriptions: Any = aflow_prototypes['prototypes_by_spacegroup'].get(space_group, [])
    for type_description in type_descriptions:
        current_norm_wyckoffs = type_description.get('normalized_wyckoff_matid')
        if current_norm_wyckoffs and current_norm_wyckoffs == norm_wyckoff:
            structure_type_info = type_description
            break
    return structure_type_info


def get_brillouin_zone(reciprocal_lattice: NDArray[Any]) -> dict:
    '''
    Calculates the Brillouin Zone information from the given reciprocal
    lattice.

    This function uses the crystallographic definition, so there is no factor
    of 2*Pi.

    Args:
        primitive_lattice: The primitive cell as a matrix where rows are the
            cell basis vectors.

    Returns:
        A dictionary containing:
        'vertices': The vertices of the first Brillouin zone
        'faces': The indices of the vertices that make up the faces on the
            first Brillouin zone. The order of these indices matter, because
            only when combined sequentially they form the correct face.
    '''
    # Create the near lattice points that surround the origin
    b1 = reciprocal_lattice[0, :]
    b2 = reciprocal_lattice[1, :]
    b3 = reciprocal_lattice[2, :]

    list_k_points = []
    for i, j, k in itertools.product([-1, 0, 1], [-1, 0, 1], [-1, 0, 1]):
        list_k_points.append(i * b1 + j * b2 + k * b3)

    # Create the first Brillouin zone by calculating a Voronoi cell starting
    # from the reciprocal cell origin.
    voronoi = Voronoi(list_k_points)
    origin_index = 13

    # Get the vertices. The regions attribute will contain a list of
    # different regions that were found during the Voronoi creation. We want
    # the Voronoi region for the point at the origin.
    point_region = voronoi.point_region[13]
    vertice_indices = voronoi.regions[point_region]
    vertices = voronoi.vertices[vertice_indices].tolist()

    # Create a mapping between the original index and an index in the new list
    index_map = {
        old_id: new_id for (new_id, old_id) in enumerate(vertice_indices)
    }

    # The ridges are the faces of a 3D Voronoi cell. Here we search for ridges
    # that are placed between the origin and some other point. These form the
    # BZ faces.
    faces = []
    for key in voronoi.ridge_dict:
        if key[0] == origin_index or key[1] == origin_index:
            ridge_indices = voronoi.ridge_dict[key]
            new_ridge_indices = [index_map[i] for i in ridge_indices]
            faces.append(new_ridge_indices)
    faces = faces

    brillouin_zone = {
        'vertices': vertices,
        'faces': faces,
    }
    return brillouin_zone


def get_minimized_structure(atoms: Atoms):
    '''
    Reduce cell size to just fit the system in the non-periodic dimensions.

    Args:
        atoms: The structure to minimize

    Returns:
        A new structure where the non-periodic dimension have been minimized.
    '''
    min_atoms = atoms.copy()
    pos = atoms.get_scaled_positions(wrap=False)
    cell = atoms.get_cell()
    new_cell = np.array(cell)
    translation = np.zeros(3)
    pbc = atoms.get_pbc()
    for index, periodic in enumerate(pbc):
        if not periodic:
            imin = np.min(pos[:, index])
            imax = np.max(pos[:, index])
            translation -= cell[index, :] * imin
            new_cell[index] = cell[index, :] * (imax - imin)

    min_atoms.translate(translation)
    min_atoms.set_cell(new_cell)

    return min_atoms


def swap_basis(atoms, a, b):
    cell_old = atoms.get_cell()
    pbc_old = atoms.get_pbc()
    cell_new = np.array(cell_old)
    cell_new[a] = cell_old[b]
    cell_new[b] = cell_old[a]
    pbc_new = np.array(pbc_old)
    pbc_new[a] = pbc_old[b]
    pbc_new[b] = pbc_old[a]
    atoms.set_cell(cell_new)
    atoms.set_pbc(pbc_new)


class Formula():
    '''Helper class for extracting formulas used by NOMAD.
    '''
    formula_parentheses = re.compile(r'([A-Z][a-z]?)\((\d+)\)')
    symbols = set(ase.data.chemical_symbols[1:])
    placeholder_symbol = 'X'

    def __init__(self, formula: str, unknown: str = 'replace'):
        '''
        Args:
            formula: Chemical formula to work on
            unknown: The handling of unknown species. Can be one of:
            - 'replace': Replaces the unknown species with X.
            - 'remove': Removes the unknown species.
            - 'keep': Keeps the labels for the unknown species.

        Raises:
            ValueError if no meaningful formula can be extracted.
        '''
        # Try if ASE can directly make sense of the formula
        try:
            self.ase_formula = ASEFormula(formula)
        # Try if can be interpreted as fractional formula
        except Exception:
            try:
                self.ase_formula = ASEFormula(Composition(formula).get_integer_formula_and_factor()[0])
            # Try if formula contains parentheses that can be removed
            except Exception:
                self.ase_formula = ASEFormula(self._remove_parentheses(formula))

        count = self.ase_formula.count()
        if unknown == 'remove':
            for key in list(count.keys()):
                if key not in self.symbols:
                    del count[key]
        elif unknown == 'replace':
            for key, value in list(count.items()):
                if key not in self.symbols and key != self.placeholder_symbol:
                    if self.placeholder_symbol not in count:
                        count[self.placeholder_symbol] = value
                    else:
                        count[self.placeholder_symbol] += value
                    del count[key]
        elif unknown == 'keep':
            pass
        else:
            raise ValueError('Invalid option for the argument "unknown"')
        self._count = count
        if len(count) == 0:
            raise ValueError(f'Could not extract any species from the formula "{formula}"')
        self._original_formula = formula

    def count(self) -> Dict[str, int]:
        '''Return dictionary that maps chemical symbol to number of atoms.
        '''
        return self._count.copy()

    def format(self, fmt: str) -> str:
        '''
        Args:
            fmt: The used format. Available options:
            - hill: Formula in Hill notation (see chemical_formula_hill)
            - iupac: The IUPAC formula (see chemical_formula_iupac)
            - reduced: Reduced formula (see chemical_formula_reduced)
            - anonymous: Anonymized formula (see chemical_formula_anonymous)
            - original: The originally supplied formula format
        '''
        if fmt == 'hill':
            return self._formula_hill()
        if fmt == 'iupac':
            return self._formula_iupac()
        if fmt == 'reduced':
            return self._formula_reduced()
        if fmt == 'anonymous':
            return self._formula_anonymous()
        if fmt == 'original':
            return self._original_formula
        else:
            raise ValueError(f'Invalid format option "{fmt}"')

    def elements(self) -> List[str]:
        '''Returns the list of chemical elements present in the formula.
        '''
        return sorted(self.count().keys())

    def atomic_fractions(self) -> Dict[str, float]:
        '''Returns dictionary that maps chemical symbol to atomic fraction.

        Returns:
            Dict[str, float]: Dictionary with chemical symbol as key and the atomic
            fraction as value.
        '''
        count = self.count()
        total_count = sum(count.values())
        atomic_fractions = {key: value / total_count for key, value in count.items()}
        return atomic_fractions

    def populate_material(self, material: Material,
                          descriptive_format: Union[str, None] = 'original') -> None:
        '''Populates the supplied material object with the list of elements as well as the
        formulae in the formats: hill, reduced, iupac, anonymous and descriptive.
        The descriptive formula defaults to the originally supplied formula but can be
        changed to any other format by supplying the `descriptive_format` argument as:
        'hill', 'reduced', 'iupac' or 'anonymous'. If descriptive_format is None, no
        descriptive formula will be added to the material.

        Args:
            material (Material): The material object to be populated with elements and
            formulae.
            descriptive_format (Union[str, None], optional): The format used for the
            descriptive formula (see `format` method for details). If None, the materials
            descriptive formula is not set. Defaults to 'original'.
        '''
        material.elements = self.elements()
        material.chemical_formula_hill = self.format('hill')
        material.chemical_formula_reduced = self.format('reduced')
        material.chemical_formula_iupac = self.format('iupac')
        material.chemical_formula_anonymous = self.format('anonymous')
        if descriptive_format:
            material.chemical_formula_descriptive = self.format(descriptive_format)

    def _remove_parentheses(self, formula: str) -> str:
        '''Used to remove parentheses from a formula. E.g. C(2)O(1) becomes C2O1
        Args:
            formula: the original formula

        Returns:
            The formula without parentheses.
        '''
        matches = list(re.finditer(self.formula_parentheses, formula))
        if matches:
            n_matched_chars = sum([len(match[0]) for match in matches])
            n_formula = len(formula.strip())
            if n_matched_chars == n_formula:
                formula = ''.join(['{}{}'.format(match[1], match[2]) for match in matches])
        return formula

    def _formula_hill(self) -> str:
        '''Returns the Hill formula.
        '''
        count = self.count()
        count_hill = {}
        if 'C' in count:
            count_hill['C'] = count.pop('C')
            if 'H' in count:
                count_hill['H'] = count.pop('H')
        for symb, n in sorted(count.items()):
            count_hill[symb] = n

        return self._dict2str(count_hill)

    def _formula_iupac(self) -> str:
        '''Returns the IUPAC formula.
        '''
        count = self.count()
        counts = count.values()
        symbols = list(count.keys())
        gcd = reduce(math.gcd, counts)
        symbols_sorted = sorted(
            symbols,
            key=lambda x: float('inf') if x == self.placeholder_symbol else get_el_sp(x).iupac_ordering)
        count_iupac = {symbol: int(count[symbol] / gcd) for symbol in symbols_sorted}

        return self._dict2str(count_iupac)

    def _formula_reduced(self) -> str:
        '''Returns the reduced formula.
        '''
        count = self.count()
        counts = count.values()
        gcd = reduce(math.gcd, counts)
        count_reduced = {key: int(value / gcd) for key, value in sorted(count.items())}

        return self._dict2str(count_reduced)

    def _formula_anonymous(self) -> str:
        '''Returns the anonymous formula.
        '''
        count = self.count()
        counts = count.values()
        gcd = reduce(math.gcd, counts)
        count_anonymous = {
            ascii_uppercase[index]: int(value / gcd)
            for index, value in enumerate(reversed(sorted(count.values())))
        }
        return self._dict2str(count_anonymous)

    def _dict2str(self, dct: Dict[str, int]) -> str:
        '''Convert symbol-to-count dict to a string. Omits the chemical
        proportion number 1.

        Args:
            dct: The dictionary of symbol-to-count items to convert.

        Returns:
            Chemical formula as a string.
        '''
        return ''.join(symb + (str(n) if n > 1 else '') for symb, n in dct.items())


def create_empty_universe(n_atoms: int, n_frames: int = 1, n_residues: int = 1, n_segments: int = 1,
                          atom_resindex: List[int] = None, residue_segindex: List[int] = None,
                          flag_trajectory: bool = False, flag_velocities: bool = False, flag_forces: bool = False):
    '''Create a blank Universe

    This function was adapted from the function empty() within the MDA class Universe().
    The only difference is that the Universe() class is imported directly here, whereas in the
    original function is is passed as a function argument, since the function there is a classfunction.

    Useful for building a Universe without requiring existing files,
    for example for system building.

    If `flag_trajectory` is set to True, a
    :class:`MDAnalysis.coordinates.memory.MemoryReader` will be
    attached to the Universe.

    Parameters
    ----------
    n_atoms: int
      number of Atoms in the Universe
    n_residues: int, default 1
      number of Residues in the Universe, defaults to 1
    n_segments: int, default 1
      number of Segments in the Universe, defaults to 1
    atom_resindex: array like, optional
      mapping of atoms to residues, e.g. with 6 atoms,
      `atom_resindex=[0, 0, 1, 1, 2, 2]` would put 2 atoms
      into each of 3 residues.
    residue_segindex: array like, optional
      mapping of residues to segments
    flag_trajectory: bool, optional
      if True, attaches a :class:`MDAnalysis.coordinates.memory.MemoryReader`
      allowing coordinates to be set and written.  Default is False
    flag_velocities: bool, optional
      include velocities in the :class:`MDAnalysis.coordinates.memory.MemoryReader`
    flag_forces: bool, optional
      include forces in the :class:`MDAnalysis.coordinates.memory.MemoryReader`

    Returns
    -------
    MDAnalysis.Universe object

    Examples
    --------
    For example to create a new Universe with 6 atoms in 2 residues, with
    positions for the atoms and a mass attribute:

    >>> u = mda.Universe.empty(6, 2,
                                atom_resindex=np.array([0, 0, 0, 1, 1, 1]),
                                flag_trajectory=True,
            )
    >>> u.add_TopologyAttr('masses')

    .. versionadded:: 0.17.0
    .. versionchanged:: 0.19.0
        The attached Reader when flag_trajectory=True is now a MemoryReader
    .. versionchanged:: 1.0.0
        Universes can now be created with 0 atoms
    '''

    if not n_atoms:
        n_residues = 0
        n_segments = 0

    if atom_resindex is None:
        warnings.warn(
            'Residues specified but no atom_resindex given.  '
            'All atoms will be placed in first Residue.',
            UserWarning)

    if residue_segindex is None:
        warnings.warn(
            'Segments specified but no segment_resindex given.  '
            'All residues will be placed in first Segment',
            UserWarning)

    topology = Topology(n_atoms, n_residues, n_segments, atom_resindex=atom_resindex,
                        residue_segindex=residue_segindex)

    universe = Universe(topology)

    if flag_trajectory:
        coords = np.zeros((n_frames, n_atoms, 3), dtype=np.float32)
        vels = np.zeros_like(coords) if flag_velocities else None
        forces = np.zeros_like(coords) if flag_forces else None

        # grab and attach a MemoryReader
        universe.trajectory = get_reader_for(coords)(
            coords, order='fac', n_atoms=n_atoms,
            velocities=vels, forces=forces)

    return universe


def archive_to_universe(archive, system_index: int = 0, method_index: int = -1, model_index: int = -1):
    '''Extract the topology from a provided run section of an archive entry

        Input:

            archive_sec_run: section run of an EntryArchive

            system_index: list index of archive.run[].system to be used for topology extraction

            method_index: list index of archive.run[].method to be used for atom parameter (charges and masses) extraction

            model_index: list index of archive.run[].method[].force_field.model for bond list extraction

        Variables:

            n_frames (int):

            n_atoms (int):

            atom_names (str, shape=(n_atoms)):

            atom_types (str, shape=(n_atoms)):

            atom_resindex (str, shape=(n_atoms)):

            atom_segids (str, shape=(n_atoms)):

            n_segments (int): Segments correspond to a group of the same type of molecules.

            n_residues (int): The number of distinct residues (nb - individual molecules are also denoted as a residue).

            resnames (str, shape=(n_residues)): The name of each residue.

            residue_segindex (int, shape=(n_residues)): The segment index that each residue belongs to.

            residue_molnums (int, shape=(n_residues)): The molecule index that each residue belongs to.

            residue_moltypes (int, shape=(n_residues)): The molecule type of each residue.

            n_molecules (int):

            masses (float, shape=(n_atoms)):  atom masses, units = amu

            charges (float, shape=(n_atoms)): atom partial charges, units = e

            positions (float, shape=(n_frames,n_atoms,3)): atom positions

            velocities (float, shape=(n_frames,n_atoms,3)): atom velocities

            dimensions (float, shape=(n_frames,6)): box dimensions (nb - currently assuming a cubic box!)

            bonds (tuple, shape=([])): list of tuples with the atom indices of each bond
    '''

    try:
        sec_run = archive.run[-1]
        sec_system = sec_run.system
        sec_system_top = sec_run.system[system_index]
        sec_atoms = sec_system_top.atoms
        sec_atoms_group = sec_system_top.atoms_group
        sec_method = sec_run.method[method_index] if sec_run.get('method') is not None else None
        sec_force_field = sec_method.force_field if sec_method is not None else None
        sec_model = sec_force_field.model[model_index] if sec_force_field is not None else None
    except IndexError:
        logging.warning('Supplied indices or necessary sections do not exist in archive. Cannot build the MDA universe.')
        return None

    n_atoms = sec_atoms.get('n_atoms')
    if n_atoms is None:
        logging.warning('No atoms found in the archive. Cannot build the MDA universe.')
        return None

    n_frames = len(sec_system) if sec_system is not None else None
    atom_names = sec_atoms.get('labels')
    model_atom_parameters = sec_method.get('atom_parameters')
    atom_types = [atom.label for atom in model_atom_parameters] if model_atom_parameters else atom_names
    atom_resindex = np.arange(n_atoms)
    atoms_segindices = np.empty(n_atoms)
    atom_segids = np.array(range(n_atoms), dtype='object')
    molecule_groups = sec_atoms_group
    n_segments = len(molecule_groups)

    n_residues = 0
    n_molecules = 0
    residue_segindex = []
    resnames = []
    residue_moltypes = []
    residue_min_atom_index = []
    residue_n_atoms = []
    molecule_n_res = []
    for mol_group_ind, mol_group in enumerate(molecule_groups):
        atoms_segindices[mol_group.atom_indices] = mol_group_ind
        atom_segids[mol_group.atom_indices] = mol_group.label
        molecules = mol_group.atoms_group if mol_group.atoms_group is not None else []
        for mol in molecules:
            monomer_groups = mol.atoms_group
            mol_res_counter = 0
            if monomer_groups:
                for mon_group in monomer_groups:
                    monomers = mon_group.atoms_group
                    for mon in monomers:
                        resnames.append(mon.label)
                        residue_segindex.append(mol_group_ind)
                        residue_moltypes.append(mol.label)
                        residue_min_atom_index.append(np.min(mon.atom_indices))
                        residue_n_atoms.append(len(mon.atom_indices))
                        n_residues += 1
                        mol_res_counter += 1
            else:  # no monomers => whole molecule is it's own residue
                resnames.append(mol.label)
                residue_segindex.append(mol_group_ind)
                residue_moltypes.append(mol.label)
                residue_min_atom_index.append(np.min(mol.atom_indices))
                residue_n_atoms.append(len(mol.atom_indices))
                n_residues += 1
                mol_res_counter += 1
            molecule_n_res.append(mol_res_counter)
            n_molecules += 1
    n_residues = n_residues

    # reorder the residues by atom_indices
    residue_data = np.array([
        [residue_min_atom_index[i], residue_n_atoms[i], residue_segindex[i], residue_moltypes[i], resnames[i]]
        for i in range(len(residue_min_atom_index))], dtype=object)
    residue_data = np.array(sorted(residue_data, key=lambda x: x[0], reverse=False)).T
    residue_n_atoms = residue_data[1].astype(int)
    residue_segindex = residue_data[2].astype(int)
    residue_moltypes = residue_data[3]
    resnames = residue_data[4]
    res_index_counter = 0
    for i_residue, res_n_atoms in enumerate(residue_n_atoms):
        atom_resindex[res_index_counter:res_index_counter + res_n_atoms] = i_residue
        res_index_counter += res_n_atoms
    residue_molnums = np.array(range(n_residues))
    mol_index_counter = 0
    for i_molecule, n_res in enumerate(molecule_n_res):
        residue_molnums[mol_index_counter:mol_index_counter + n_res] = i_molecule
        mol_index_counter += n_res

    # get the atom masses and charges
    masses = np.empty(n_atoms)
    charges = np.empty(n_atoms)
    atom_parameters = sec_method.get('atom_parameters') if sec_method is not None else []
    atom_parameters = atom_parameters if atom_parameters is not None else []
    for atom_ind, atom in enumerate(atom_parameters):
        if atom.get('mass'):
            masses[atom_ind] = ureg.convert(atom.mass.magnitude, atom.mass.units, ureg.amu)
        if atom.get('charge'):
            charges[atom_ind] = ureg.convert(atom.charge.magnitude, atom.charge.units, ureg.e)

    # get the atom positions, velocites, and box dimensions
    positions = np.empty(shape=(n_frames, n_atoms, 3))
    velocities = np.empty(shape=(n_frames, n_atoms, 3))
    dimensions = np.empty(shape=(n_frames, 6))
    for frame_ind, frame in enumerate(sec_system):
        sec_atoms_fr = frame.get('atoms')
        if sec_atoms_fr is not None:
            positions_frame = sec_atoms_fr.positions
            positions[frame_ind] = ureg.convert(positions_frame.magnitude, positions_frame.units,
                                                ureg.angstrom) if positions_frame is not None else None
            velocities_frame = sec_atoms_fr.velocities
            velocities[frame_ind] = ureg.convert(velocities_frame.magnitude, velocities_frame.units,
                                                 ureg.angstrom / ureg.picosecond) if velocities_frame is not None else None
            latt_vec_tmp = sec_atoms_fr.get('lattice_vectors')
            if latt_vec_tmp is not None:
                length_conversion = ureg.convert(1.0, sec_atoms_fr.lattice_vectors.units, ureg.angstrom)
                dimensions[frame_ind] = [
                    sec_atoms_fr.lattice_vectors.magnitude[0][0] * length_conversion,
                    sec_atoms_fr.lattice_vectors.magnitude[1][1] * length_conversion,
                    sec_atoms_fr.lattice_vectors.magnitude[2][2] * length_conversion,
                    90, 90, 90]  # TODO: extend to non-cubic boxes

    # get the bonds
    bonds = []
    contributions = sec_model.get('contributions') if sec_model is not None else []
    contributions = contributions if contributions is not None else []
    for contribution in contributions:
        if contribution.type == 'bond':  # and contribution.atom_indices is not None:
            bonds.append(tuple(contribution.atom_indices))

    # create the Universe
    metainfo_universe = create_empty_universe(
        n_atoms, n_frames=n_frames, n_residues=n_residues, n_segments=n_segments, atom_resindex=atom_resindex,
        residue_segindex=residue_segindex, flag_trajectory=True, flag_velocities=True)

    # set the positions and velocities
    for frame_ind, frame in enumerate(metainfo_universe.trajectory):
        metainfo_universe.atoms.positions = positions[frame_ind]
        metainfo_universe.atoms.velocities = velocities[frame_ind]

    # add the atom attributes
    metainfo_universe.add_TopologyAttr('name', atom_names)
    metainfo_universe.add_TopologyAttr('type', atom_types)
    metainfo_universe.add_TopologyAttr('mass', masses)
    metainfo_universe.add_TopologyAttr('charge', charges)
    if n_segments != 0:
        metainfo_universe.add_TopologyAttr('segids', np.unique(atom_segids))
    if n_residues != 0:
        metainfo_universe.add_TopologyAttr('resnames', resnames)
        metainfo_universe.add_TopologyAttr('resids', np.unique(atom_resindex) + 1)
        metainfo_universe.add_TopologyAttr('resnums', np.unique(atom_resindex) + 1)
    if len(residue_molnums) > 0:
        metainfo_universe.add_TopologyAttr('molnums', residue_molnums)
    if len(residue_moltypes) > 0:
        metainfo_universe.add_TopologyAttr('moltypes', residue_moltypes)

    # add the box dimensions
    for frame_ind, frame in enumerate(metainfo_universe.trajectory):
        metainfo_universe.atoms.dimensions = dimensions[frame_ind]

    # add the bonds
    if hasattr(metainfo_universe, 'bonds'):
        logging.warning('archive_to_universe() failed, universe already has bonds.')
        return None
    metainfo_universe.add_TopologyAttr('bonds', bonds)

    return metainfo_universe


class BeadGroup(object):
    # see https://github.com/MDAnalysis/mdanalysis/issues/1891#issuecomment-387138110
    # by @richardjgowers with performance improvements
    def __init__(self, atoms, compound='fragments'):
        '''Initialize with an AtomGroup instance.
        Will split based on keyword 'compounds' (residues or fragments).
        '''
        self._atoms = atoms
        self.compound = compound
        self._nbeads = len(getattr(self._atoms, self.compound))
        # for caching
        self._cache = {}
        self._cache['positions'] = None
        self.__last_frame = None

    def __len__(self):
        return self._nbeads

    @property
    def positions(self):
        # cache positions for current frame
        if self.universe.trajectory.frame != self.__last_frame:
            self._cache['positions'] = self._atoms.center_of_mass(
                unwrap=True, compound=self.compound)
            self.__last_frame = self.universe.trajectory.frame
        return self._cache['positions']

    @property  # type: ignore
    @MDAnalysis.lib.util.cached('universe')
    def universe(self):
        return self._atoms.universe


def __log_indices(first: int, last: int, num: int = 100):
    ls = np.logspace(0, np.log10(last - first + 1), num=num)
    return np.unique(np.int_(ls) - 1 + first)


def __correlation(function, positions: List[float]):
    iterator = iter(positions)
    start_frame = next(iterator)
    return map(lambda f: function(start_frame, f), chain([start_frame], iterator))


def shifted_correlation_average(function: Callable, times: NDArray, positions: NDArray,
                                index_distribution: Callable = __log_indices, correlation: Callable = __correlation,
                                segments: int = 10, window: float = 0.5, skip: int = 0):

    '''
    Code adapted from MDevaluate module: https://github.com/mdevaluate/mdevaluate.git

    Calculate the time series for a correlation function.

    The times at which the correlation is calculated are determined automatically by the
    function given as ``index_distribution``. The default is a logarithmic distribution.

    The function has been edited so that the average is always calculated, i.e., average=True below.

    Args:
        function:   The function that should be correlated
        positions:     The coordinates of the simulation data
        index_distribution (opt.):
                    A function that returns the indices for which the timeseries
                    will be calculated
        correlation (function, opt.):
                    The correlation function
        segments (int, opt.):
                    The number of segments the time window will be shifted
        window (float, opt.):
                    The fraction of the simulation the time series will cover
        skip (float, opt.):
                    The fraction of the trajectory that will be skipped at the beginning,
                    if this is None the start index of the frames slice will be used,
                    which defaults to 0.
        counter (bool, opt.):
                    If True, returns length of frames (in general number of particles specified)
        average (bool, opt.):
                    If True,
    Returns:
        tuple:
            A list of length N that contains the indices of the frames at which
            the time series was calculated and a numpy array of shape (segments, N)
            that holds the (non-avaraged) correlation data

            if has_counter == True: adds number of counts to output tupel.
                                    if average is returned it will be weighted.

    Example:
        Calculating the mean square displacement of a coordinates object named ``coords``:

        >>> indices, data = shifted_correlation(msd, coords)
    '''
    if window + skip >= 1:
        warnings.warn('Invalid parameters for shifted_correlation(), '
                      'resetting to defaults.', UserWarning)
        window = 0.5
        skip = 0

    start_frames = np.unique(np.linspace(
        len(positions) * skip, len(positions) * (1 - window),
        num=segments, endpoint=False, dtype=int
    ))
    num_frames = int(len(positions) * (window))

    idx = index_distribution(0, num_frames)

    def correlate(start_frame):
        shifted_idx = idx + start_frame
        return correlation(function, map(positions.__getitem__, shifted_idx))

    correlation_times = np.array([times[i] for i in idx]) - times[0]

    result: NDArray
    for i_start_frame, start_frame in enumerate(start_frames):
        if i_start_frame == 0:
            result = np.array(list(correlate(start_frame)))
        else:
            result += np.array(list(correlate(start_frame)))
    result = np.array(result)
    result = result / len(start_frames)

    return correlation_times, result


def __calc_diffusion_constant(times: NDArray, values: NDArray, dim: int = 3):
    '''
    Determines the diffusion constant from a fit of the mean squared displacement
    vs. time according to the Einstein relation.
    '''
    linear_model = linregress(times, values)
    slope = linear_model.slope
    error = linear_model.rvalue
    return slope * 1 / (2 * dim), error


def __get_molecular_bead_groups(universe: MDAnalysis.Universe, moltypes: List[str] = None):

    if moltypes is None:
        atoms_moltypes = getattr(universe.atoms, 'moltypes', [])
        moltypes = np.unique(atoms_moltypes)
    bead_groups = {}
    for moltype in moltypes:
        ags_by_moltype = universe.select_atoms('moltype ' + moltype)
        ags_by_moltype = ags_by_moltype[ags_by_moltype.masses > abs(1e-2)]  # remove any virtual/massless sites (needed for, e.g., 4-bead water models)
        bead_groups[moltype] = BeadGroup(ags_by_moltype, compound='fragments')

    return bead_groups


def calc_molecular_rdf(universe: MDAnalysis.Universe, n_traj_split: int = 10, n_prune: int = 1, interval_indices=None):
    '''
    Calculates the radial distribution functions between for each unique pair of
    molecule types as a function of their center of mass distance.

    interval_indices: 2D array specifying the groups of the n_traj_split intervals to be averaged
    '''

    if not universe or not universe.trajectory or universe.trajectory[0].dimensions is None:
        return

    n_frames = universe.trajectory.n_frames
    if n_frames < n_traj_split:
        n_traj_split = 1
        frames_start = np.array([0])
        frames_end = np.array([n_frames])
        n_frames_split = np.array([n_frames])
        interval_indices = [[0]]
    else:
        run_len = int(n_frames / n_traj_split)
        frames_start = np.arange(n_traj_split) * run_len
        frames_end = frames_start + run_len
        frames_end[-1] = n_frames
        n_frames_split = frames_end - frames_start
        if np.sum(n_frames_split) != n_frames:
            logging.error('Something went wrong with input parameters in calc_molecular_rdf().'
                          'Radial distribution functions will not be calculated.')
            return
        if not interval_indices:
            interval_indices = [[i] for i in range(n_traj_split)]

    atoms_moltypes = getattr(universe.atoms, 'moltypes', [])
    moltypes = np.unique(atoms_moltypes)
    bead_groups = __get_molecular_bead_groups(universe, moltypes=moltypes)

    min_box_dimension = np.min(universe.trajectory[0].dimensions[:3])
    max_rdf_dist = min_box_dimension / 2
    n_bins = 200
    n_smooth = 2

    rdf_results: Dict[str, Any] = {}
    rdf_results['n_smooth'] = n_smooth
    rdf_results['n_prune'] = n_prune
    rdf_results['type'] = 'molecular'
    rdf_results['types'] = []
    rdf_results['variables_name'] = []
    rdf_results['bins'] = []
    rdf_results['value'] = []
    rdf_results['frame_start'] = []
    rdf_results['frame_end'] = []
    for i, moltype_i in enumerate(moltypes):
        for j, moltype_j in enumerate(moltypes):
            if j > i:
                continue
            elif i == j and bead_groups[moltype_i].positions.shape[0] == 1:  # skip if only 1 mol in group
                continue

            if i == j:
                exclusion_block = (1, 1)  # remove self-distance
            else:
                exclusion_block = None
            pair_type = f'{moltype_i}-{moltype_j}'
            rdf_results_interval: Dict[str, Any] = {}
            rdf_results_interval['types'] = []
            rdf_results_interval['variables_name'] = []
            rdf_results_interval['bins'] = []
            rdf_results_interval['value'] = []
            rdf_results_interval['frame_start'] = []
            rdf_results_interval['frame_end'] = []
            for i_interval in range(n_traj_split):
                rdf_results_interval['types'].append(pair_type)
                rdf_results_interval['variables_name'].append(['distance'])
                rdf = MDA_RDF.InterRDF(
                    bead_groups[moltype_i], bead_groups[moltype_j], range=(0, max_rdf_dist),
                    exclusion_block=exclusion_block, nbins=n_bins).run(
                    frames_start[i_interval], frames_end[i_interval], n_prune)
                rdf_results_interval['frame_start'].append(frames_start[i_interval])
                rdf_results_interval['frame_end'].append(frames_end[i_interval])

                rdf_results_interval['bins'].append(rdf.results.bins[int(n_smooth / 2):-int(n_smooth / 2)] * ureg.angstrom)
                rdf_results_interval['value'].append(np.convolve(
                    rdf.results.rdf, np.ones((n_smooth,)) / n_smooth,
                    mode='same')[int(n_smooth / 2):-int(n_smooth / 2)])

            flag_logging_error = False
            for interval_group in interval_indices:
                split_weights = n_frames_split[np.array(interval_group)] / np.sum(n_frames_split[np.array(interval_group)])
                if abs(np.sum(split_weights) - 1.0) > 1e-6:
                    flag_logging_error = True
                    continue
                rdf_values_avg = split_weights[0] * rdf_results_interval['value'][interval_group[0]]
                for i_interval, interval in enumerate(interval_group[1:]):
                    if rdf_results_interval['types'][interval] != rdf_results_interval['types'][interval - 1]:
                        flag_logging_error = True
                        continue
                    if rdf_results_interval['variables_name'][interval] != rdf_results_interval['variables_name'][interval - 1]:
                        flag_logging_error = True
                        continue
                    if not (rdf_results_interval['bins'][interval] == rdf_results_interval['bins'][interval - 1]).all():
                        flag_logging_error = True
                        continue
                    rdf_values_avg += split_weights[i_interval + 1] * rdf_results_interval['value'][interval]
                if flag_logging_error:
                    logging.error('Something went wrong in calc_molecular_rdf(). Some interval groups were skipped.')
                rdf_results['types'].append(rdf_results_interval['types'][interval_group[0]])
                rdf_results['variables_name'].append(rdf_results_interval['variables_name'][interval_group[0]])
                rdf_results['bins'].append(rdf_results_interval['bins'][interval_group[0]])
                rdf_results['value'].append(rdf_values_avg)
                rdf_results['frame_start'].append(int(rdf_results_interval['frame_start'][interval_group[0]]))
                rdf_results['frame_end'].append(int(rdf_results_interval['frame_end'][interval_group[-1]]))

    return rdf_results


def calc_molecular_mean_squared_displacements(universe: MDAnalysis.Universe):
    '''
    Calculates the mean squared displacement for the center of mass of each
    molecule type.
    '''

    def mean_squared_displacement(start: NDArray, current: NDArray):
        '''
        Calculates mean square displacement between current and initial (start) coordinates.
        '''
        vec = start - current
        return (vec ** 2).sum(axis=1).mean()

    if not universe or not universe.trajectory or universe.trajectory[0].dimensions is None:
        return

    n_frames = universe.trajectory.n_frames
    if n_frames < 50:
        warnings.warn('Less than 50 frames in trajectory, not calculating molecular'
                      'mean squared displacements.', UserWarning)
        return

    atoms_moltypes = getattr(universe.atoms, 'moltypes', [])
    moltypes = np.unique(atoms_moltypes)
    bead_groups = __get_molecular_bead_groups(universe, moltypes=moltypes)
    times = np.arange(n_frames) * universe.trajectory.dt

    msd_results: Dict[str, Any] = {}
    msd_results['type'] = 'molecular'
    msd_results['direction'] = 'xyz'
    msd_results['value'] = []
    msd_results['times'] = []
    msd_results['diffusion_constant'] = []
    msd_results['error_diffusion_constant'] = []
    for moltype in moltypes:
        positions = __get_nojump_positions(universe, bead_groups[moltype])
        results = shifted_correlation_average(mean_squared_displacement, times, positions)
        if results:
            msd_results['value'].append(results[1])
            msd_results['times'].append(results[0])
            diffusion_constant, error = __calc_diffusion_constant(*results)
            msd_results['diffusion_constant'].append(diffusion_constant)
            msd_results['error_diffusion_constant'].append(error)

    msd_results['types'] = moltypes
    msd_results['times'] = np.array(msd_results['times']) * ureg.picosecond
    msd_results['value'] = np.array(msd_results['value']) * ureg.angstrom**2
    msd_results['diffusion_constant'] = (np.array(
        msd_results['diffusion_constant']) * ureg.angstrom**2 / ureg.picosecond)
    msd_results['error_diffusion_constant'] = np.array(msd_results['error_diffusion_constant'])

    return msd_results


def __parse_jumps(universe: MDAnalysis.Universe, selection: MDAnalysis.AtomGroup):
    __ = universe.trajectory[0]
    prev = np.array(selection.positions)
    box = universe.trajectory[0].dimensions[:3]
    sparse_data = namedtuple('SparseData', ['data', 'row', 'col'])  # type: ignore[name-match]
    jump_data = (
        sparse_data(data=array('b'), row=array('l'), col=array('l')),
        sparse_data(data=array('b'), row=array('l'), col=array('l')),
        sparse_data(data=array('b'), row=array('l'), col=array('l'))
    )

    for i_frame, _ in enumerate(universe.trajectory[1:]):
        curr = np.array(selection.positions)
        delta = ((curr - prev) / box).round().astype(np.int8)
        prev = np.array(curr)
        for d in range(3):
            col, = np.where(delta[:, d] != 0)
            jump_data[d].col.extend(col)
            jump_data[d].row.extend([i_frame] * len(col))
            jump_data[d].data.extend(delta[col, d])

    return jump_data


def __generate_nojump_matrices(universe: MDAnalysis.Universe, selection: MDAnalysis.AtomGroup):
    jump_data = __parse_jumps(universe, selection)
    n_frames = len(universe.trajectory)
    n_atoms = selection.positions.shape[0]

    nojump_matrices = tuple(
        sparse.csr_matrix((np.array(m.data), (m.row, m.col)), shape=(n_frames, n_atoms)) for m in jump_data
    )
    return nojump_matrices


def __get_nojump_positions(universe: MDAnalysis.Universe, selection: MDAnalysis.AtomGroup):
    nojump_matrices = __generate_nojump_matrices(universe, selection)
    box = universe.trajectory[0].dimensions[:3]

    nojump_positions = []
    for i_frame, __ in enumerate(universe.trajectory):
        delta = np.array(np.vstack(
            [m[:i_frame, :].sum(axis=0) for m in nojump_matrices]
        ).T) * box
        nojump_positions.append(selection.positions - delta)

    return np.array(nojump_positions)


def calc_radius_of_gyration(universe: MDAnalysis.Universe, molecule_atom_indices: NDArray):
    '''
    Calculates the radius of gyration as a function of time for the atoms 'molecule_atom_indices'.
    '''

    if not universe or not universe.trajectory or universe.trajectory[0].dimensions is None:
        return

    selection = ' '.join([str(i) for i in molecule_atom_indices])
    selection = f'index {selection}'
    molecule = universe.select_atoms(selection)
    rg_results: Dict[str, Any] = {}
    rg_results['type'] = 'molecular'
    rg_results['times'] = []
    rg_results['value'] = []
    for __ in universe.trajectory:
        rg_results['times'].append(universe.trajectory.time)
        rg_results['value'].append(molecule.radius_of_gyration())
    rg_results['n_frames'] = len(rg_results['times'])
    rg_results['times'] = np.array(rg_results['times']) * ureg.picosecond
    rg_results['value'] = np.array(rg_results['value']) * ureg.angstrom

    return rg_results


def calc_molecular_radius_of_gyration(universe: MDAnalysis.Universe, system_topology):
    '''
    Calculates the radius of gyration as a function of time for each polymer in the system.
    '''

    if not system_topology:
        return []

    rg_results = []
    for molgroup in system_topology:
        for molecule in molgroup.get('atoms_group'):
            sec_monomer_groups = molecule.get('atoms_group')
            group_type = sec_monomer_groups[0].type if sec_monomer_groups else None
            if group_type != 'monomer_group':
                continue
            rg_result = calc_radius_of_gyration(universe, molecule.atom_indices)
            rg_result['label'] = molecule.label + '-index_' + str(molecule.index)
            rg_result['atomsgroup_ref'] = molecule
            rg_results.append(rg_result)

    return rg_results
