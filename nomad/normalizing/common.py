#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the 'License');
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an 'AS IS' BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
import numpy as np
from ase import Atoms
from typing import List, Set, Any, Optional
from nptyping import NDArray
from matid import SymmetryAnalyzer
from matid.symmetry.wyckoffset import WyckoffSet as WyckoffSetMatID
import matid.geometry

from nomad import atomutils
from nomad import config
from nomad.units import ureg
from nomad.datamodel.metainfo.simulation.system import System, Atoms as NOMADAtoms
from nomad.datamodel.optimade import Species
from nomad.datamodel.results import (
    Cell,
    Structure,
    LatticeParameters,
    WyckoffSet,
)


def wyckoff_sets_from_matid(wyckoff_sets: List[WyckoffSetMatID]) -> List[WyckoffSet]:
    '''Given a dictionary of wyckoff sets, returns the metainfo equivalent.

    Args:
        wyckoff_sets: List of Wyckoff sets as returned by MatID.

    Returns:
        List of NOMAD WyckoffSet objects.
    '''
    wsets = []
    for group in wyckoff_sets:
        wset = WyckoffSet()
        if group.x is not None or group.y is not None or group.z is not None:
            if group.x is not None:
                wset.x = float(group.x)
            if group.y is not None:
                wset.y = float(group.y)
            if group.z is not None:
                wset.z = float(group.z)
        wset.indices = group.indices
        wset.element = group.element
        wset.wyckoff_letter = group.wyckoff_letter
        wsets.append(wset)
    return wsets


def species(labels: List[str], atomic_numbers: List[int], logger=None) -> Optional[List[Species]]:
    '''Given a list of atomic labels and atomic numbers, returns the
    corresponding list of Species objects.

    Args:
        labels: List of atomic labels.
        atomic_numbers: List of atomic numbers.  If the atomic number does not
            correspond to an actual element, an error message is logged and the
            species for thatitem is not created.

    Returns:
        List of Species objects.
    '''
    if labels is None or atomic_numbers is None:
        return None
    species: Set[str] = set()
    species_list = []
    for label, atomic_number in zip(labels, atomic_numbers):
        if label not in species:
            species.add(label)
            i_species = Species()
            i_species.name = label
            try:
                symbol = atomutils.chemical_symbols([atomic_number])[0]
            except ValueError:
                if logger:
                    logger.info(f'could not identify chemical symbol for atomic number {atomic_number}')
            else:
                i_species.chemical_symbols = [symbol]
            i_species.concentration = [1.0]
            species_list.append(i_species)

    return species_list


def lattice_parameters_from_array(lattice_vectors: NDArray[Any]) -> LatticeParameters:
    '''Converts the given 3x3 numpy array into metainfo LatticeParameters.
    Undefined angle values are not stored.

    Args:
        lattice_vectors: 3x3 array where the lattice vectors are the rows.
        Should be given in meters.

    Returns:
        LatticeParameters object.
    '''
    param_values = atomutils.cell_to_cellpar(lattice_vectors)
    params = LatticeParameters()
    params.a = float(param_values[0])
    params.b = float(param_values[1])
    params.c = float(param_values[2])
    alpha = float(param_values[3])
    params.alpha = None if np.isnan(alpha) else alpha
    beta = float(param_values[4])
    params.beta = None if np.isnan(beta) else beta
    gamma = float(param_values[5])
    params.gamma = None if np.isnan(gamma) else gamma
    return params


def cell_from_ase_atoms(atoms: Atoms) -> Cell:
    '''Extracts Cell metainfo from the given ASE Atoms.
    Undefined angle values are not stored.

    Args:
        atoms: The system from which the information is extracted from.

    Returns:
        Cell object.
    '''
    param_values = atomutils.cell_to_cellpar(atoms.cell)
    cell = Cell()
    cell.a = float(param_values[0]) * ureg.angstrom
    cell.b = float(param_values[1]) * ureg.angstrom
    cell.c = float(param_values[2]) * ureg.angstrom
    alpha = float(param_values[3])
    cell.alpha = None if np.isnan(alpha) else alpha * ureg.radian
    beta = float(param_values[4])
    cell.beta = None if np.isnan(beta) else beta * ureg.radian
    gamma = float(param_values[5])
    cell.gamma = None if np.isnan(gamma) else gamma * ureg.radian

    volume = atoms.cell.volume * ureg.angstrom ** 3
    mass = atomutils.get_summed_atomic_mass(atoms.get_atomic_numbers()) * ureg.kg
    mass_density = mass / volume
    number_of_atoms = atoms.get_number_of_atoms()
    atomic_density = number_of_atoms / volume
    cell.volume = volume
    cell.atomic_density = atomic_density
    cell.mass_density = mass_density

    return cell


def cell_from_structure(structure: Structure) -> Cell:
    '''Extracts Cell metainfo from the given Structure.
    Undefined angle values are not stored.

    Args:
        structure: The system from which the information is extracted from.

    Returns:
        Cell object.
    '''
    return Cell(
        a=structure.lattice_parameters.a,
        b=structure.lattice_parameters.b,
        c=structure.lattice_parameters.c,
        alpha=structure.lattice_parameters.alpha,
        beta=structure.lattice_parameters.beta,
        gamma=structure.lattice_parameters.gamma,
        volume=structure.cell_volume,
        atomic_density=structure.atomic_density,
        mass_density=structure.mass_density,
    )


def structure_from_ase_atoms(system: Atoms, wyckoff_sets: List[WyckoffSetMatID] = None, logger=None) -> Structure:
    '''Returns a populated NOMAD Structure instance from an ase.Atoms-object.

    Args:
        atoms: The system to transform.

    Returns:
        NOMAD Structure instance.
    '''
    if system is None:
        return None
    struct = Structure()
    labels = system.get_chemical_symbols()
    atomic_numbers = system.get_atomic_numbers()
    struct.species_at_sites = system.get_chemical_symbols()
    struct.cartesian_site_positions = system.get_positions() * ureg.angstrom
    struct.species = species(labels, atomic_numbers, logger)
    lattice_vectors = system.get_cell()
    if lattice_vectors is not None:
        lattice_vectors = (lattice_vectors * ureg.angstrom).to(ureg.meter).magnitude
        struct.dimension_types = [1 if x else 0 for x in system.get_pbc()]
        struct.lattice_vectors = lattice_vectors
        cell_volume = atomutils.get_volume(lattice_vectors)
        struct.cell_volume = cell_volume
        if system.get_pbc().all() and cell_volume:
            mass = atomutils.get_summed_atomic_mass(atomic_numbers)
            struct.mass_density = mass / cell_volume
            struct.atomic_density = len(system) / cell_volume
        if wyckoff_sets:
            struct.wyckoff_sets = wyckoff_sets_from_matid(wyckoff_sets)
        struct.lattice_parameters = lattice_parameters_from_array(lattice_vectors)
    return struct


def structure_from_nomad_system(system: System, logger=None) -> Structure:
    '''Returns a populated NOMAD Structure instance from a NOMAD System.

    Args:
        system: The system to transform.
        logger: Optional logger to use

    Returns:
        NOMAD Structure instance.
    '''
    if system is None:
        return None
    struct = Structure()
    struct.cartesian_site_positions = system.atoms.positions
    struct.species_at_sites = system.atoms.labels
    labels = system.atoms.labels
    atomic_numbers = system.atoms.species
    struct.species = species(labels, atomic_numbers, logger)
    lattice_vectors = system.atoms.lattice_vectors
    if lattice_vectors is not None:
        lattice_vectors = lattice_vectors.magnitude
        struct.dimension_types = np.array(system.atoms.periodic).astype(int)
        struct.lattice_vectors = lattice_vectors
        cell_volume = atomutils.get_volume(lattice_vectors)
        struct.cell_volume = cell_volume
        if all(system.atoms.periodic) and cell_volume:
            if system.atoms.species is not None:
                mass = atomutils.get_summed_atomic_mass(atomic_numbers)
                struct.mass_density = mass / cell_volume
            struct.atomic_density = len(system.atoms.labels) / cell_volume
        struct.lattice_parameters = lattice_parameters_from_array(lattice_vectors)
    return struct


def nomad_atoms_from_ase_atoms(system: Atoms) -> NOMADAtoms:
    '''Returns a populated NOMAD Atoms instance from an ase.Atoms-object.

    Args:
        atoms: The system to transform.

    Returns:
        NOMAD Atoms instance.
    '''
    if system is None:
        return None

    atoms = NOMADAtoms()
    atoms.positions = system.get_positions() * ureg.angstrom
    atoms.labels = system.get_chemical_symbols()
    atoms.atomic_numbers = system.get_atomic_numbers()
    atoms.species = atoms.atomic_numbers
    atoms.lattice_vectors = system.get_cell() * ureg.angstrom
    atoms.periodic = system.get_pbc()

    return atoms


def ase_atoms_from_nomad_atoms(system: NOMADAtoms) -> Atoms:
    '''Returns a populated ase.Atoms instance from a NOMAD Atoms instance.

    Args:
        system: The system to transform.

    Returns:
        ase.Atoms instance.
    '''
    return Atoms(
        positions=system.positions.to(ureg.angstrom),
        numbers=system.atomic_numbers,
        cell=system.lattice_vectors.to(ureg.angstrom),
        pbc=system.periodic
    )


def ase_atoms_from_structure(system: Structure) -> Atoms:
    '''Returns an instance of ASE.Atoms from a NOMAD Structure-section.

    Args:
        system: The system to transform

    Returns:
        A new ASE.Atoms created from the given data.
    '''
    symbol_map = {}
    for species in system.species:
        assert len(species.chemical_symbols) == 1, "Unable to transform system with multi-species sites as ASE.Atoms."
        symbol_map[species.name] = species.chemical_symbols[0]
    symbols = [symbol_map[x] for x in system.species_at_sites]

    return Atoms(
        positions=system.cartesian_site_positions.to(ureg.angstrom),
        symbols=symbols,
        cell=system.lattice_vectors.to(ureg.angstrom),
        pbc=np.array(system.dimension_types, dtype=bool)
    )


def structures_2d(original_atoms, logger=None):
    conv_atoms = None
    prim_atoms = None
    wyckoff_sets = None
    spg_number = None
    try:
        # Get dimension of system by also taking into account the covalent radii
        dimensions = matid.geometry.get_dimensions(original_atoms, [True, True, True])
        basis_dimensions = np.linalg.norm(original_atoms.get_cell(), axis=1)
        gaps = basis_dimensions - dimensions
        periodicity = gaps <= config.normalize.cluster_threshold

        # If two axis are not periodic, return. This only happens if the vacuum
        # gap is not aligned with a cell vector or if the linear gap search is
        # unsufficient (the structure is "wavy" making also the gap highly
        # nonlinear).
        if sum(periodicity) != 2:
            # TODO: @ Lauri: I'm not sure if this is smart/wise, but the periodicity calculation doesn't seem to work for conventional cell systems. The calculated periodicity is always (True, True, True). At least that's my interpretation of what's happening. This is my quick and dirty and at the moment only solution for this problem.
            periodicity = original_atoms.pbc
            # if logger:
            #     logger.error("could not detect the periodic dimensions in a 2D system")
            # return conv_atoms, prim_atoms, wyckoff_sets, spg_number

        # Center the system in the non-periodic direction, also taking
        # periodicity into account. The get_center_of_mass()-function in MatID
        # takes into account periodicity and can produce the correct CM unlike
        # the similar function in ASE.
        pbc_cm = matid.geometry.get_center_of_mass(original_atoms)
        cell_center = 0.5 * np.sum(original_atoms.get_cell(), axis=0)
        translation = cell_center - pbc_cm
        symm_system = original_atoms.copy()
        symm_system.translate(translation)
        symm_system.wrap()

        # Set the periodicity according to detected periodicity in order
        # for SymmetryAnalyzer to use the symmetry analysis designed for 2D
        # systems.
        symm_system.set_pbc(periodicity)
        symmetry_analyzer = SymmetryAnalyzer(
            symm_system,
            config.normalize.symmetry_tolerance,
            config.normalize.flat_dim_threshold
        )

        spg_number = symmetry_analyzer.get_space_group_number()
        wyckoff_sets = symmetry_analyzer.get_wyckoff_sets_conventional(return_parameters=False)
        conv_atoms = symmetry_analyzer.get_conventional_system()
        prim_atoms = symmetry_analyzer.get_primitive_system()

        # Reduce cell size to just fit the system in the non-periodic
        # dimensions.
        conv_atoms = atomutils.get_minimized_structure(conv_atoms)
        prim_atoms = atomutils.get_minimized_structure(prim_atoms)

        # Swap the cell axes so that the non-periodic one is always the
        # last basis (=c)
        swap_dim = 2
        for i, periodic in enumerate(conv_atoms.get_pbc()):
            if not periodic:
                non_periodic_dim = i
                break
        if non_periodic_dim != swap_dim:
            atomutils.swap_basis(conv_atoms, non_periodic_dim, swap_dim)
            atomutils.swap_basis(prim_atoms, non_periodic_dim, swap_dim)
    except Exception as e:
        if logger:
            logger.error(
                'could not construct a conventional system for a 2D material',
                exc_info=e
            )
    return conv_atoms, prim_atoms, wyckoff_sets, spg_number
