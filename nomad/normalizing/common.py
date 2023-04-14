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
import MDAnalysis as mda
from matid import SymmetryAnalyzer  # pylint: disable=import-error
from matid.symmetry.wyckoffset import WyckoffSet as WyckoffSetMatID  # pylint: disable=import-error
import matid.geometry  # pylint: disable=import-error

from nomad import atomutils
from nomad import config
from nomad.utils import hash
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
    cell.volume = volume
    mass = atomutils.get_summed_atomic_mass(atoms.get_atomic_numbers()) * ureg.kg
    cell.mass_density = None if volume == 0 else mass / volume
    number_of_atoms = atoms.get_number_of_atoms()
    cell.atomic_density = None if volume == 0 else number_of_atoms / volume
    cell.pbc = np.zeros(3, bool)[:] = atoms.get_pbc()

    return cell


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
        numbers=system.species,
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


def mda_universe_from_nomad_atoms(system: Atoms, logger=None) -> mda.Universe:
    '''Returns an instance of mda.Universe from a NOMAD Atoms-section.

    Args:
        system: The atoms to transform

    Returns:
        A new mda.Universe created from the given data.
    '''
    n_atoms = len(system.positions)
    n_residues = 1
    atom_resindex = [0] * n_atoms
    residue_segindex = [0]

    universe = mda.Universe.empty(
        n_atoms,
        n_residues=n_residues,
        atom_resindex=atom_resindex,
        residue_segindex=residue_segindex,
        trajectory=True
    )

    # Add positions
    universe.atoms.positions = system.positions.to(ureg.angstrom).magnitude

    # Add atom attributes
    atom_names = system.labels
    universe.add_TopologyAttr('name', atom_names)
    universe.add_TopologyAttr('type', atom_names)
    universe.add_TopologyAttr('element', atom_names)

    # Add the box dimensions
    universe.atoms.dimensions = atomutils.cell_to_cellpar(
        system.lattice_vectors.to(ureg.angstrom).magnitude,
        degrees=True
    )

    return universe


def structures_2d(original_atoms, logger=None):
    conv_atoms = None
    prim_atoms = None
    wyckoff_sets = None
    spg_number = None
    symm_system = original_atoms.copy()
    n_pbc = sum(original_atoms.get_pbc())

    try:
        # If the given system if fully periodic, try to extract a 2D system by
        # checking the presence of vacuum gaps.
        if n_pbc == 3:
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
                if logger:
                    logger.error("could not detect the periodic dimensions in a 2D system")
                return conv_atoms, prim_atoms, wyckoff_sets, spg_number

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
            0.4,  # The value is increased here to better match 2D materials.
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
    except Exception as e:
        if logger:
            logger.error(
                'could not construct a conventional system for a 2D material',
                exc_info=e
            )
    return conv_atoms, prim_atoms, wyckoff_sets, spg_number


def material_id_bulk(spg_number: int, wyckoff_sets) -> str:
    if spg_number is None or wyckoff_sets is None:
        return None
    norm_hash_string = atomutils.get_symmetry_string(spg_number, wyckoff_sets)
    return hash(norm_hash_string)


def material_id_2d(spg_number: int, wyckoff_sets) -> str:
    if spg_number is None or wyckoff_sets is None:
        return None
    norm_hash_string = atomutils.get_symmetry_string(spg_number, wyckoff_sets, is_2d=True)
    return hash(norm_hash_string)


def material_id_1d(conv_atoms: Atoms) -> str:
    '''Hash to be used as identifier for a 1D material. Based on Coulomb
    matrix eigenvalues and the Hill formula.

    The fingerprint is based on calculating a discretized version of a
    sorted Coulomb matrix eigenspectrum (Grégoire Montavon, Katja Hansen,
    Siamac Fazli, Matthias Rupp, Franziska Biegler, Andreas Ziehe,
    Alexandre Tkatchenko, Anatole V. Lilienfeld, and Klaus-Robert Müller.
    Learning invariant representations of molecules for atomization energy
    prediction. In F. Pereira, C. J. C. Burges, L. Bottou, and K. Q.
    Weinberger, editors, Advances in Neural Information Processing Systems
    25, pages 440–448. Curran Associates, Inc., 2012.).

    The fingerprints are discretized in order to perform O(n) matching
    between structures (no need to compare fingerprints against each
    other). As regular discretization is susceptible to the 'edge problem',
    a robust discretization is used instead (Birget, Jean-Camille & Hong,
    Dawei & Memon, Nasir. (2003). Robust discretization, with an
    application to graphical passwords. IACR Cryptology ePrint Archive.
    2003. 168.) Basically for the 1-dimensional domain two grids are
    created and the points are mapped to the first grid in which they are
    robust using a minimum tolerance parameter r, with the maximum
    tolerance being 5r.

    There are other robust discretization methods that can guarantee exact
    r-tolerance (e.g. Sonia Chiasson, Jayakumar Srinivasan, Robert Biddle,
    and P. C. van Oorschot. 2008. Centered discretization with application
    to graphical passwords. In Proceedings of the 1st Conference on
    Usability, Psychology, and Security (UPSEC’08). USENIX Association,
    USA, Article 6, 1–9.). This method however requires that a predefined
    'correct' structure exists against which the search is done.
    '''
    if conv_atoms is None:
        return None

    # Calculate charge part
    q = conv_atoms.get_atomic_numbers()
    qiqj = np.sqrt(q[None, :] * q[:, None])

    # Calculate distance part. Notice that the minimum image convention
    # must be used. Without it, differently oriented atoms in the same cell
    # may be detected as the same material.
    pos = conv_atoms.get_positions()
    cell = conv_atoms.get_cell()
    cmat = 10 - matid.geometry.get_distance_matrix(pos, pos, cell, pbc=True, mic=True)
    cmat = np.clip(cmat, a_min=0, a_max=None)
    np.fill_diagonal(cmat, 0)
    cmat = qiqj * cmat

    # Calculate eigenvalues
    eigval, _ = np.linalg.eigh(cmat)

    # Sort eigenvalues
    eigval = np.array(sorted(eigval))

    # Perform robust discretization (see function docstring for details). r
    # = 0.5 ensures that all grids are integers which can be uniquely
    # mapped to strings. If finer grid is needed adjust the eigenvalue scale
    # instead.
    eigval /= 25  # Go to smaller scale where integer numbers are meaningful
    dimension = 1
    r = 0.5
    spacing = 2 * r * (dimension + 1)
    phi_k = 2 * r * np.array(range(dimension + 1))
    t = np.mod((eigval[None, :] + phi_k[:, None]), (2 * r * (dimension + 1)))
    grid_mask = (r <= t) & (t < r * (2 * dimension + 1))
    safe_grid_k = np.argmax(grid_mask == True, axis=0)   # noqa: E712
    discretization = spacing * np.floor((eigval + (2 * r * safe_grid_k)) / spacing)
    discretization[safe_grid_k == 1] += 2 * r

    # Construct formula
    names, counts = atomutils.get_hill_decomposition(conv_atoms.get_chemical_symbols(), reduced=False)
    formula = atomutils.get_formula_string(names, counts)

    # Form hash
    strings = []
    for number in discretization:
        num_str = str(int(number))
        strings.append(num_str)
    fingerprint = ';'.join(strings)
    id_strings = []
    id_strings.append(formula)
    id_strings.append(fingerprint)
    hash_seed = ', '.join(id_strings)
    hash_val = hash(hash_seed)

    return hash_val
