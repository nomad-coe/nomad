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

from ase import neighborlist, geometry
from ase.data import chemical_symbols, covalent_radii, atomic_numbers
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.local_env import JmolNN
from nomad.datamodel.results import Coordination

"""
Function for deconstructing MOFs into building units
There are three types of building units
1) The organic linkers, which contains the all atoms of the organic ligands
2) The metal ase_atom, which contains the metal cluster found in the MOF
3) The organic ase_atom, which is the fragment of the organic linker cut at the point of extension
"""


def transition_metals():
    """
    Function containing list of symbols of all metals found in the periodic table
    """
    metal = [
        symbol
        for symbol in chemical_symbols
        if symbol
        not in [
            chemical_symbols[main_index]
            for main_index in [
                1,
                2,
                5,
                6,
                7,
                8,
                9,
                10,
                14,
                15,
                16,
                17,
                18,
                33,
                34,
                35,
                36,
                52,
                53,
                54,
                85,
                86,
            ]
        ]
    ]
    return metal


def inter_atomic_distance_check(ase_atom):
    """
    Check that no two atoms are within a distance of 1.0 Amstrong unless it is an X-H bond

    As simple script to convert from ase atom object to pybel
    Parameters
    ----------
    ase_atom : ASE atoms object

    Returns
    -------
    boolean
    """
    valid = True
    distances = ase_atom.get_all_distances(mic=True)
    for i in range(len(distances)):
        if ase_atom[i].symbol != 'H':
            for j in range(len(distances[i])):
                if i != j:
                    if distances[i, j] < 0.90:
                        valid = False
                        break
    return valid


def covalent_radius(element):
    """
    Extract ase covelent radii for an element
    """
    a_n = atomic_numbers[element]
    return covalent_radii[a_n]


def ase2xyz(atoms):
    """
    Create and xyz string from an ase atom object to be compatible with
    pybel in order to perfom some cheminformatics.

    """

    if any(atoms.get_pbc()):
        raise RuntimeError(' Does not support periodic systems!')
    num_atoms = len(atoms)
    elements = atoms.get_chemical_symbols()
    all_atoms = zip(elements, atoms.get_positions())
    a_str = str(num_atoms) + '\n' + '\n'
    for atom in all_atoms:
        a_str += atom[0] + ' ' + ' '.join([str(x) for x in atom[1]]) + '\n'
    return a_str[:-1]


# def compute_inchis(obmol):
#     '''
#     Openbabel function to extract inchi and inchikey for molecules
#     '''
#     conv = ob.OBConversion()
#     conv.SetOutFormat("inchi")
#     inchi = conv.WriteString(obmol).rstrip()
#     inchi = inchi.split('InChI=')[1]
#     conv.SetOptions("K", conv.OUTOPTIONS)
#     inchikey = conv.WriteString(obmol).rstrip()
#     return inchi, inchikey


# def compute_smi(obmol):
#     '''
#     Openbabel function to extract smile strings for molecules
#     '''
#     conv = ob.OBConversion()
#     conv.SetOutFormat("smi")
#     smi = conv.WriteString(obmol).rstrip()
#     return smi


# def ase2pybel(atoms):
#     """
#     As simple script to convert from ase atom object to pybel
#     Parameters
#     ----------
#     atoms : ASE atoms object

#     Returns
#     -------
#     pybel: pybel molecular object.
#     """
#     a_str = ase2xyz(atoms)
#     pybel_mol = pb.readstring("xyz", a_str)

#     return pybel_mol


def max_index(lists):
    for i, elt in enumerate(lists):
        if elt == max(lists):
            return i


# def compute_openbabel_cheminformatic(ase_atom):
#     '''
#     A procedure that uses openbabel to generate
#     computer readable file formats
#     '''

#     new_ase_atom = wrap_systems_in_unit_cell(ase_atom)

#     new_ase_atom.set_pbc(False)
#     pybel_mol = ase2pybel(new_ase_atom)

#     obmol = pybel_mol.OBMol
#     inchi, inchikey = compute_inchis(obmol)
#     smi = compute_smi(obmol)
#     return smi, inchi, inchikey


def compute_ase_neighbour(ase_atom):
    """
    Create a connectivity graph using ASE neigbour list.

    Parameters:
    -----------
    ASE atoms

    Returns
    -------
    Returns a python dictionary, wherein each atom index is key and the value
    are the indices of it neigbours.
    e.g.
    atom_neighbors ={0:[1,2,3,4], 1:[3,4,5]...}
    """
    atom_neighbors = {}
    cutOff = neighborlist.natural_cutoffs(ase_atom)

    neighborList = neighborlist.NeighborList(
        cutOff, self_interaction=False, bothways=True
    )
    neighborList.update(ase_atom)
    matrix = neighborList.get_connectivity_matrix(sparse=False)

    for atoms in ase_atom:
        connectivity, _ = neighborList.get_neighbors(atoms.index)
        atom_neighbors[atoms.index] = connectivity

    return atom_neighbors, matrix


def matrix2dict(bond_matrix):
    """
    A simple procedure to convert an adjacent matrix to a python dictionary
    Parameters:
    -----------
    Bond matrix
    type: nxn ndarray

    Returns
    -------
    python dictionary
    """
    graph = {}
    for idx, row in enumerate(bond_matrix):
        temp = []
        for r in range(len(row)):
            if row[r] != 0:
                temp.append(r)
        graph[idx] = temp
    return graph


def dfsutil_graph_method(graph, temp, node, visited):
    """
    Depth-first search graph algorithm for traversing graph data structures.
    I starts at the root 'node' and explores as far as possible along
    each branch before backtracking.
    It is used here a a util for searching connected components in the MOF graph
    Parameters:
    -----------
    graph: any python dictionary
    temp: a python list to hold nodes that have been visited
    node: a key in the python dictionary (graph), which is used as the starting or root node
    visited: python list containing nodes that have been traversed.
    Returns
    -------
    python dictionary
    """
    visited[node] = True
    temp.append(node)
    for i in graph[node]:
        if visited[i] is False:
            temp = dfsutil_graph_method(graph, temp, i, visited)
    return temp


def longest_list(lst):
    """
    return longest list in list of list
    """
    return max(lst, key=len)


def remove_unbound_guest(ase_atom):
    """
    A simple script to remove guest from a metal organic framework.
    1)It begins by computing a connected graph component of all the fragments in the system
    using ASE neighbour list.
    2) Secondly it selects indicies of connected components which contain a metal
    3)if the there are two or more components, we create a pytmagen graph for each components and filter out all components that are not polymeric
    4) If there are two or more polymeric components, we check wether these systems there are identical or different
    and select only unique polymeric components

    Parameters:
    -----------
    ASE atoms

    Returns
    -------
    mof_indices : indinces of the guest free system. The guest free ase_atom object
    can be obtain as follows;
    E.g.
    guest_free_system = ase_atom[mof_indices]
    """
    atom_neighbors, _ = compute_ase_neighbour(ase_atom)
    fragments = connected_components(atom_neighbors)
    if len(fragments) == 1:
        return [atom.index for atom in ase_atom]
    else:
        polymeric_indices = []
        for i in range(len(fragments)):
            super_cell = ase_atom[fragments[i]] * (2, 1, 1)
            coordination_graph, _ = compute_ase_neighbour(super_cell)
            pymat_graph = connected_components(coordination_graph)
            if len(pymat_graph) == 1:
                polymeric_indices.append(i)
        if len(polymeric_indices) > 0:
            Graphs = [
                StructureGraph.with_local_env_strategy(
                    AseAtomsAdaptor.get_structure(ase_atom[fragments[i]]), JmolNN()
                )
                for i in polymeric_indices
            ]
            temp_indices = [polymeric_indices[0]]
            unique = [Graphs[0]]
            for k in range(len(polymeric_indices[1:])):
                graph = Graphs[k]
                if True in [graph.diff(j) for j in unique]:
                    unique.append(Graphs[k])
                    temp_indices.append(polymeric_indices[i])
            mof_indices = []
            for frag_indices in temp_indices:
                mof_indices.extend(fragments[frag_indices])
            if len(mof_indices) == 0:
                return longest_list(fragments)
            else:
                return mof_indices
        else:
            return sum(fragments, [])


def connected_components(graph):
    """
    Find the connected fragments in a graph. Should work for any graph defined as a dictionary
    Parameters:
    -----------
    A graph in the form of dictionary
    e.g.
    graph = {1:[0,1,3], 2:[2,4,5]}

    Returns
    -------
    Returns a python list of list of connected components
    These correspond to individual molecular fragments.
    list_of_connected_components = [[1,2],[1,3,4]]
    """
    visited = []
    list_of_connected_components = []
    for _ in list(graph.keys()):
        visited.append(False)
    for v in list(graph.keys()):
        if visited[v] is False:
            temp = []
            list_of_connected_components.append(
                dfsutil_graph_method(graph, temp, v, visited)
            )
    return list_of_connected_components


def check_planarity(p1, p2, p3, p4):
    """
    A simple procedure to check whether a point is planar to three other points.
    Important to distinguish porphyrin type metals
    Parameters:
    -----------
    p1, p2, p3, p4 : ndarray containing x,y,z values

    Returns
    -------
    Boolean
    True: planar
    False: noneplanar
    """
    planar = False
    a1 = p2[0] - p1[0]
    b1 = p2[1] - p1[1]
    c1 = p2[2] - p1[2]
    a2 = p3[0] - p1[0]
    b2 = p3[1] - p1[1]
    c2 = p3[2] - p1[2]
    # ------------------------------------------------
    a = b1 * c2 - b2 * c1
    b = a2 * c1 - a1 * c2
    c = a1 * b2 - b1 * a2
    d = round((-a * p1[0] - b * p1[1] - c * p1[2]), 0)
    # ------------------------------------------------
    factor = round((a * p4[0] + b * p4[1] + c * p4[2]), 0)
    verify = factor + d
    # ------------------------------------------------
    if verify == 0:
        planar = True
    return planar


def metal_in_porphyrin(ase_atom, graph):
    """
    https://en.wikipedia.org/wiki/Transition_metal_porphyrin_complexes
    Check whether a metal is found at the centre of a phorphirin

    Parameters:
    -----------
    ase_atom: ASE atom
    graph: python dictionary containing neigbours

    Returns
    -------
    list of indices consisting of index of metal atoms found in the ASE atom
    """
    all_porphyrin = []
    all_metal_symbols = [
        atom.index for atom in ase_atom if atom.symbol in transition_metals()
    ]
    for idx in all_metal_symbols:
        connected = graph[idx]
        all_nitrogens = [i for i in connected if ase_atom[i].symbol == 'N']
        if len(all_nitrogens) == 4:
            p1, p2, p3, p4 = (
                ase_atom[all_nitrogens[0]].position,
                ase_atom[all_nitrogens[1]].position,
                ase_atom[all_nitrogens[2]].position,
                ase_atom[all_nitrogens[3]].position,
            )
            planarity = check_planarity(p1, p2, p3, p4)
            if planarity:
                all_porphyrin.append(idx)
                all_porphyrin.extend(all_nitrogens)

    return all_porphyrin


def metal_in_porphyrin2(ase_atom, graph):
    """
    https://en.wikipedia.org/wiki/Transition_metal_porphyrin_complexes
    Check whether a metal is found at the centre of a phorphirin

    Parameters:
    -----------
    ase_atom: ASE atom
    graph: python dictionary containing neigbours

    Returns
    -------
    list of indices consisting of index of metal atoms found in the ASE atom
    """
    old_list_of_connected_components = connected_components(graph)
    all_porphyrin = []
    all_metal_symbols = [
        atom.index for atom in ase_atom if atom.symbol in transition_metals()
    ]
    metal_tmp = []
    N_tmp = []
    for idx in all_metal_symbols:
        connected = graph[idx]
        all_nitrogens = [i for i in connected if ase_atom[i].symbol == 'N']
        if len(all_nitrogens) == 4:
            metal_tmp.append(idx)
            N_tmp.extend(all_nitrogens)
    new_atom_indices = [i.index for i in ase_atom if i.index not in metal_tmp]
    tmp_atom = ase_atom[new_atom_indices]
    atom_neighbors, _ = compute_ase_neighbour(tmp_atom)
    list_of_connected_components = connected_components(atom_neighbors)
    if len(list_of_connected_components) == len(old_list_of_connected_components):
        all_porphyrin = metal_tmp + N_tmp
    return all_porphyrin


def move2front(index_value, coords):
    """
    Parameters:
    -----------
    index_value: index of item to move to the from
    coords: list of coords

    Returns
    -------
    Move an index from any position in the list to the front
    The function is important to set the cell of a rodmof to point in the
    a-axis. Such that the system can be grow along this axis
    """
    if any(isinstance(el, list) for el in coords):
        for data in coords:
            data.insert(0, data.pop(index_value))
    else:
        coords.insert(0, coords.pop(index_value))
    return coords


def mof_regions(ase_atom, list_of_connected_components, atom_pairs_at_breaking_point):
    """
    A function to map all atom indices to exact position in which the find themselves in the MOF.
    This function is used to partition a MOF into regions that correspond to unique
    unique building units.

    Parameters:
    -----------
    ase_atom: ASE atom
    list_of_connected_components : list of list, wherein each list correspond to atom indices of a specific building unit
    atom_pairs_at_breaking_point: dictionary containing pairs of atoms from which the bonds were broken

    Returns
    -------
    Move an index from any position in the list to the front
    The function is important to set the cell of a rodmof to point in the
    a-axis. Such that the system can be grow along this axis

    """
    all_regions = {}
    all_pm_structures = [
        sorted(ase_atom[i].symbols) for i in list_of_connected_components
    ]
    for i in range(len(all_pm_structures)):
        temp = []
        for j in range(len(all_pm_structures)):
            if all_pm_structures[i] == all_pm_structures[j]:
                temp.append(j)
        if temp not in all_regions.values():
            all_regions[i] = temp
    Xis_regions = {}
    for idx in range(len(all_regions.keys())):
        frag = list(all_regions.keys())[idx]
        components = list_of_connected_components[all_regions[frag][0]]
        Xis = [
            [
                atom_pairs_at_breaking_point[j]
                for j in list(atom_pairs_at_breaking_point.keys())
                if j in comp
            ]
            for comp in components
        ]
        Xis_regions[idx] = Xis
    return all_regions, Xis_regions


def find_carboxylates(ase_atom, graph):
    """
    A simple aglorimth to search for carboxylates found in the system.
    Parameters:
    -----------
    ase_atom: ASE atom

    Returns
    -------
    dictionary of key = carbon index and values = oxygen index
    """
    carboxyl = {}
    for atoms in ase_atom:
        if atoms.symbol == 'C':
            index = atoms.index
            oxygen = [i for i in graph[index] if ase_atom[i].symbol == 'O']
            if len(oxygen) == 2:
                oxy_metal = sum(
                    [
                        [
                            j
                            for j in graph[i]
                            if ase_atom[j].symbol in transition_metals()
                        ]
                        for i in oxygen
                    ],
                    [],
                )
                if len(oxy_metal) > 0:
                    carboxyl[index] = oxygen
    return carboxyl


def find_phosphite(ase_atom, graph):
    """
    A simple aglorimth to search for sulfides.
     P
     |
    -C
     |
     P

    Parameters:
    -----------
    ase_atom: ASE atom

    Returns
    -------
    dictionary of key = carbon index and values = phosphorous index
    """
    phosphorous = {}
    for atoms in ase_atom:
        if atoms.symbol == 'C':
            index = atoms.index
            phosphorous_atoms = [i for i in graph[index] if ase_atom[i].symbol == 'P']
            if len(phosphorous_atoms) == 2:
                phosphorous_to_metal = sum(
                    [
                        [
                            j
                            for j in graph[i]
                            if ase_atom[j].symbol in transition_metals()
                        ]
                        for i in phosphorous_atoms
                    ],
                    [],
                )
                if len(phosphorous_to_metal) > 0:
                    phosphorous[index] = phosphorous_atoms
    return phosphorous


def find_carbonyl_sulphate(ase_atom, graph):
    """
    A simple algorithm to search for Carbonyl sulphate found in the system.
       O
       |
    -C-S
       |
       O
    Parameters:
    -----------
    ase_atom: ASE atom

    Returns
    -------
    dictionary of key = carbon index and values = oxygen index
    """
    sulphate = {}
    for atoms in ase_atom:
        if atoms.symbol == 'S':
            index = atoms.index
            oxygen = [i for i in graph[index] if ase_atom[i].symbol == 'O']
            carbon = [i for i in graph[index] if ase_atom[i].symbol == 'C']
            if len(oxygen) >= 1:
                oxy_metal = sum(
                    [
                        [
                            j
                            for j in graph[i]
                            if ase_atom[j].symbol in transition_metals()
                        ]
                        for i in oxygen
                    ],
                    [],
                )
                if len(oxy_metal) > 0:
                    if len(carbon) > 0:
                        sulphate[index] = oxygen
    return sulphate


def find_sulfides(ase_atom, graph):
    """
    A simple aglorimth to search for sulfides.
     S
     |
    -C
     |
     S

    Parameters:
    -----------
    ase_atom: ASE atom

    Returns
    -------
    dictionary of key = carbon index and values = sulphur index
    """
    sulfides = {}
    for atoms in ase_atom:
        if atoms.symbol == 'C':
            index = atoms.index
            sulphure_atoms = [i for i in graph[index] if ase_atom[i].symbol == 'S']
            if len(sulphure_atoms) == 2:
                sulphure_to_metal = sum(
                    [
                        [
                            j
                            for j in graph[i]
                            if ase_atom[j].symbol in transition_metals()
                        ]
                        for i in sulphure_atoms
                    ],
                    [],
                )
                if len(sulphure_to_metal) > 0:
                    sulfides[index] = sulphure_atoms
    return sulfides


def find_phosphate(ase_atom, graph):
    """
    A simple algorithm to search for Carbonyl sulphate found in the system.
       O
       |
      -P-o
       |
       O
    Parameters:
    -----------
    ase_atom: ASE atom

    Returns
    -------
    dictionary of key = carbon index and values = oxygen index
    """
    phosphate = {}
    for atoms in ase_atom:
        if atoms.symbol == 'P':
            index = atoms.index
            oxygen = [i for i in graph[index] if ase_atom[i].symbol == 'O']
            if len(oxygen) >= 1:
                oxy_metal = sum(
                    [
                        [
                            j
                            for j in graph[i]
                            if ase_atom[j].symbol in transition_metals()
                        ]
                        for i in oxygen
                    ],
                    [],
                )
                if len(oxy_metal) > 0:
                    phosphate[index] = oxygen
    return phosphate


def find_COS(ase_atom, graph):
    """
    A simple aglorimth to search for COS.
     O
     |
    -C
     |
     S

    Parameters:
    -----------
    ase_atom: ASE atom

    Returns
    -------
    dictionary of key = carbon index and values = sulphur index
    """
    sulfides = {}
    for atoms in ase_atom:
        if atoms.symbol == 'C':
            index = atoms.index
            if len(graph[index]) == 3:
                sulphure_atoms = [i for i in graph[index] if ase_atom[i].symbol == 'S']
                oxygen_atoms = [i for i in graph[index] if ase_atom[i].symbol == 'O']
                if len(sulphure_atoms) + len(oxygen_atoms) == 2:
                    sulphure_to_metal = sum(
                        [
                            [
                                j
                                for j in graph[i]
                                if ase_atom[j].symbol in transition_metals()
                            ]
                            for i in sulphure_atoms
                        ],
                        [],
                    )
                    oxygen_to_metal = sum(
                        [
                            [
                                j
                                for j in graph[i]
                                if ase_atom[j].symbol in transition_metals()
                            ]
                            for i in oxygen_atoms
                        ],
                        [],
                    )
                    if len(sulphure_to_metal) > 0 and len(oxygen_to_metal) > 0:
                        sulfides[index] = sulphure_atoms + oxygen_atoms
    return sulfides


def secondary_building_units(ase_atom):
    """
    1) Search for all carboxylate that are connected to a metal.
       Cut at the position between the carboxylate carbon and the adjecent carbon.
    2) Find all Nitrogen connected to metal. Check whether the nitrogen is in the
       centre of a porphirin ring. If no, cut at nitrogen metal bond.
    3) Look for oxygen that is connected to metal and two carbon. cut at metal oxygen bond
    Parameters:
    -----------
    ase_atom: ASE atom

    Returns
    -------
     list_of_connected_components : list of connected components, in which each list contains atom indices
     atom_pairs_at_breaking_point  : Dictionary containing point of disconnection
     Porphyrin_checker : Boolean showing whether the metal is in the centre of a porpherin
     Regions : Dictionary of regions.
    """
    atom_pairs_at_breaking_point = {}
    all_regions = {}
    bonds_to_break = []
    seen_phosphorous = []
    seen_carbon = []
    graph, bond_matrix = compute_ase_neighbour(ase_atom)
    porphyrin_checker = metal_in_porphyrin2(ase_atom, graph)
    carboxylates = find_carboxylates(ase_atom, graph)
    all_sulphates = find_carbonyl_sulphate(ase_atom, graph)
    all_sulfides = find_sulfides(ase_atom, graph)
    all_phosphates = find_phosphate(ase_atom, graph)
    cos_group = find_COS(ase_atom, graph)
    all_phosphites = find_phosphite(ase_atom, graph)
    ferocene_metal = all_ferrocene_metals(ase_atom, graph)
    all_metals = [i.index for i in ase_atom if i.symbol in transition_metals()]
    all_metals = [i for i in all_metals if i not in ferocene_metal + porphyrin_checker]
    for atoms in graph:
        if atoms in list(carboxylates.keys()):
            seen_carbon.append(atoms)
            connected = graph[atoms]
            all_carbon_indices = [i for i in connected if ase_atom[i].symbol == 'C']
            all_nitrogens = [i for i in connected if ase_atom[i].symbol == 'N']
            S_indx = [i for i in connected if ase_atom[i].symbol == 'S']
            if len(all_carbon_indices) == 1:
                bonds_to_break.append([atoms] + all_carbon_indices)
                atom_pairs_at_breaking_point[atoms] = all_carbon_indices[0]
            if len(all_nitrogens) == 1:
                bonds_to_break.append([atoms] + all_nitrogens)
                atom_pairs_at_breaking_point[atoms] = all_nitrogens[0]
            if len(S_indx) == 1:
                bonds_to_break.append([atoms] + S_indx)
                atom_pairs_at_breaking_point[atoms] = S_indx[0]

        if atoms in list(all_sulfides.keys()):
            seen_carbon.append(atoms)
            connected = graph[atoms]
            all_carbon_indices = [i for i in connected if ase_atom[i].symbol == 'C']
            all_nitrogens = [i for i in connected if ase_atom[i].symbol == 'N']
            # S_indx = [i for i in connected if ase_atom[i].symbol == 'S']
            if len(all_carbon_indices) == 1:
                bonds_to_break.append([atoms] + all_carbon_indices)
                atom_pairs_at_breaking_point[atoms] = all_carbon_indices[0]
            if len(all_nitrogens) == 1:
                bonds_to_break.append([atoms] + all_nitrogens)
                atom_pairs_at_breaking_point[atoms] = all_nitrogens[0]
            # if len(S_indx) == 1:
            #     bonds_to_break.append([atoms] + S_indx)
            #     atom_pairs_at_breaking_point[atoms] = S_indx[0]
        if ase_atom[atoms].symbol == 'C':
            if atoms not in seen_carbon:
                connected = graph[atoms]
                oxygens = [i for i in connected if ase_atom[i].symbol == 'O']
                carbon_metal = [
                    i for i in connected if ase_atom[i].symbol in transition_metals()
                ]
                carbon_metal = [i for i in carbon_metal if i not in ferocene_metal]
                if len(oxygens) == 1:
                    oxy_metal = [
                        i for i in oxygens if ase_atom[i].symbol in transition_metals()
                    ]
                    oxy_metal = [i for i in oxy_metal if i not in porphyrin_checker]
                    oxy_metal = [i for i in oxy_metal if i not in ferocene_metal]
                    if len(oxy_metal) == 1:
                        atom_pairs_at_breaking_point[oxygens[0]] = oxy_metal[0]
                        bonds_to_break.append([oxygens] + oxy_metal)
                if len(carbon_metal) > 0:
                    for met in carbon_metal:
                        atom_pairs_at_breaking_point[atoms] = met
                        bonds_to_break.append([atoms] + [met])

        if atoms in list(all_sulphates.keys()):
            connected = graph[atoms]
            all_carbon_indices = [i for i in connected if ase_atom[i].symbol == 'C']
            oxygen = all_sulphates[atoms]
            if len(all_carbon_indices) == 1:
                atom_pairs_at_breaking_point[atoms] = all_carbon_indices[0]
                bonds_to_break.append([atoms] + all_carbon_indices)
            if len(all_carbon_indices) > 1:
                for oxy in oxygen:
                    metal = [
                        i
                        for i in graph[oxy]
                        if ase_atom[i].symbol in transition_metals()
                    ]
                    metal = [i for i in metal if i not in porphyrin_checker]
                    metal = [i for i in metal if i not in ferocene_metal]
                    if len(metal) > 0:
                        for met in metal:
                            atom_pairs_at_breaking_point[oxy] = met
                            bonds_to_break.append([oxy] + [met])

        if atoms in list(all_phosphites.keys()):
            all__n_indices = all_phosphites[atoms]
            connected = [i for i in graph[atoms] if i not in all__n_indices]
            for neigbour in connected:
                atom_pairs_at_breaking_point[atoms] = neigbour
                bonds_to_break.append([atoms] + [neigbour])

        if atoms in list(cos_group.keys()):
            connected = graph[atoms]
            all_carbon_indices = [i for i in connected if ase_atom[i].symbol == 'C']
            non_metals = cos_group[atoms]
            if len(all_carbon_indices) == 1:
                atom_pairs_at_breaking_point[atoms] = all_carbon_indices[0]
                bonds_to_break.append([atoms] + all_carbon_indices)
            if len(all_carbon_indices) > 1:
                for atom_idx in non_metals:
                    metal = [
                        i
                        for i in graph[atom_idx]
                        if ase_atom[i].symbol in transition_metals()
                    ]
                    metal = [i for i in metal if i not in ferocene_metal]
                    if len(metal) > 0:
                        for met in metal:
                            atom_pairs_at_breaking_point[atom_idx] = met
                            bonds_to_break.append([atom_idx] + [met])

        if ase_atom[atoms].symbol == 'O':
            seen = sum(
                list(carboxylates.values())
                + list(all_sulphates.values())
                + list(all_phosphates.values())
                + list(cos_group.values()),
                [],
            )
            if atoms not in seen:
                connected = graph[atoms]
                metal = [
                    i for i in connected if ase_atom[i].symbol in transition_metals()
                ]
                metal = [i for i in metal if i not in porphyrin_checker]
                metal = [i for i in metal if i not in ferocene_metal]
                Nitrogen = [i for i in connected if ase_atom[i].symbol == 'N']
                carbon = [
                    i
                    for i in connected
                    if ase_atom[i].symbol == 'C' and i not in list(carboxylates.keys())
                ]
                if len(metal) >= 1 and len(carbon) == 1:
                    atom_pairs_at_breaking_point[atoms] = carbon[0]
                    bonds_to_break.append([atoms] + carbon)
                if len(metal) == 1 and len(Nitrogen) == 1:
                    n_carbon = [
                        i
                        for i in graph[Nitrogen[0]]
                        if ase_atom[i].symbol == 'C'
                        and i not in list(carboxylates.keys())
                    ]
                    n_nitrogen = [
                        i for i in graph[Nitrogen[0]] if ase_atom[i].symbol == 'N'
                    ]
                    n_sulphur = [
                        i
                        for i in graph[Nitrogen[0]]
                        if ase_atom[i].symbol == 'S'
                        and i not in list(all_sulphates.keys())
                    ]
                    if len(n_carbon) > 1:
                        atom_pairs_at_breaking_point[atoms] = Nitrogen[0]
                        bonds_to_break.append([atoms] + Nitrogen)

                    elif len(n_nitrogen) > 1:
                        atom_pairs_at_breaking_point[atoms] = Nitrogen[0]
                        bonds_to_break.append([atoms] + Nitrogen)
                    elif len(n_sulphur) > 1:
                        atom_pairs_at_breaking_point[atoms] = Nitrogen[0]
                        bonds_to_break.append([atoms] + Nitrogen)

        if ase_atom[atoms].symbol == 'N':
            connected = graph[atoms]
            metal = [i for i in connected if ase_atom[i].symbol in transition_metals()]
            metal = [i for i in metal if i not in ferocene_metal]
            if atoms not in porphyrin_checker:
                if len(metal) >= 1:
                    for met in metal:
                        atom_pairs_at_breaking_point[atoms] = met
                        bonds_to_break.append([atoms] + [met])

        if ase_atom[atoms].symbol == 'S':
            seen = sum(list(all_sulfides.values()) + list(cos_group.values()), [])
            if atoms not in seen:
                connected = graph[atoms]
                metal = [
                    i for i in connected if ase_atom[i].symbol in transition_metals()
                ]
                if len(metal) > 0:
                    for met in metal:
                        atom_pairs_at_breaking_point[atoms] = met
                        bonds_to_break.append([atoms, met])

        if atoms in list(all_phosphates.keys()):
            seen_phosphorous.append(atoms)
            connected = graph[atoms]
            all_carbon_indices = [i for i in connected if ase_atom[i].symbol == 'C']
            all_nitrogens = [i for i in connected if ase_atom[i].symbol == 'N']
            S_indx = [i for i in connected if ase_atom[i].symbol == 'S']
            if len(all_carbon_indices) == 1:
                bonds_to_break.append([atoms] + all_carbon_indices)
                atom_pairs_at_breaking_point[atoms] = all_carbon_indices[0]
            if len(all_nitrogens) == 1:
                bonds_to_break.append([atoms] + all_nitrogens)
                atom_pairs_at_breaking_point[atoms] = all_nitrogens[0]
            if len(S_indx) == 1:
                bonds_to_break.append([atoms] + S_indx)
                atom_pairs_at_breaking_point[atoms] = S_indx[0]

        if ase_atom[atoms].symbol == 'P':
            """
            Find the carbon closest to P, which is not bonded to a metal and cut
            1) Look for phosphorous atoms
            2) Look for all it's neigbours
            3) Look for neigbours that are not connected to metal or hydrogen.
            """
            # seen = list(all_phosphates.keys())
            seen = seen_phosphorous + sum(list(all_phosphites.values()), [])
            if atoms not in seen:
                connected = graph[atoms]
                # not_connected_to_metal_or_hygrogen = [[i for i in graph[j] if ase_atom[i].symbol not in transition_metals() or ase_atom[i].symbol != 'H'] for j in connected]

                metal_oxy = [
                    [i for i in graph[j] if ase_atom[i].symbol in transition_metals()]
                    for j in connected
                ]

                metal = sum(metal_oxy, [])
                metal = [i for i in metal if i not in porphyrin_checker]
                metal = [i for i in metal if i not in ferocene_metal]
                closest_atoms = sum(
                    [
                        [
                            i
                            for i in graph[j]
                            if i != atoms
                            and ase_atom[i].symbol not in transition_metals()
                        ]
                        for j in connected
                    ],
                    [],
                )

                if len(metal) > 0:
                    all_carbon_indices = sum(
                        [
                            [i for i in graph[j] if i in connected]
                            for j in closest_atoms
                        ],
                        [],
                    )
                    for frag in all_carbon_indices:
                        atom_pairs_at_breaking_point[atoms] = frag
                        atom_pairs_at_breaking_point[frag] = atoms
                        bonds_to_break.append([atoms, frag])

        if ase_atom[atoms].symbol == 'B':
            """
            Find the carbon closest to P, which is  bonded to a metal and cut
            """
            connected = [
                i for i in graph[atoms] if ase_atom[i].symbol not in transition_metals()
            ]
            metal_oxy = [
                [[i, j] for i in graph[j] if ase_atom[i].symbol in transition_metals()]
                for j in connected
            ]

            metal_connect = sum(metal_oxy, [])
            metal = [i for i in metal if i not in porphyrin_checker]
            metal = [i for i in metal if i not in ferocene_metal]
            if len(metal_connect) > 0:
                for frag in metal_connect:
                    atom_pairs_at_breaking_point[frag[0]] = frag[1]
                    bonds_to_break.append(frag)

    # In special cases some carbon and hydrogen get very closed to metals
    # So in such case it is important to clear this
    for metal in all_metals:
        connected = graph[metal]
        c_H = [i for i in connected if ase_atom[i].symbol in ['H', 'C']]
        if len(c_H) > 0:
            for c_h in c_H:
                bonds_to_break.append([c_h, metal])
                atom_pairs_at_breaking_point[metal] = c_h

    for bonds in bonds_to_break:
        bond_matrix[bonds[0], bonds[1]] = 0
        bond_matrix[bonds[1], bonds[0]] = 0

    new_ase_graph = matrix2dict(bond_matrix)
    try:
        list_of_connected_components = connected_components(new_ase_graph)
    except Exception:
        import networkx as nx

        N_Graph = nx.from_dict_of_lists(new_ase_graph)
        list_of_connected_components = [
            list(i) for i in list(nx.connected_components(N_Graph))
        ]

    all_pm_structures = [
        sorted(ase_atom[i].symbols) for i in list_of_connected_components
    ]
    for i in range(len(all_pm_structures)):
        temp = []
        for j in range(len(all_pm_structures)):
            if all_pm_structures[i] == all_pm_structures[j]:
                temp.append(j)
        if temp not in all_regions.values():
            all_regions[i] = temp

    return (
        list_of_connected_components,
        atom_pairs_at_breaking_point,
        porphyrin_checker,
        all_regions,
    )


def ligands_and_metal_clusters(ase_atom):
    """
    Start by checking whether there are more than 2 layers
    if yes, select one
    Here we select the largest connected component
    Parameters:
    -----------
    ase_atom: ASE atom

    Returns
    -------
    list_of_connected_components  : list of connected components, in which each list contains atom indices
     atom_pairs_at_breaking_point  : Dictionary containing point of disconnection
     Porpyrin_checker : Boolean showing whether the metal is in the centre of a porpherin
     Regions : Dictionary of regions.
    """
    graph, bond_matrix = compute_ase_neighbour(ase_atom)
    porphyrin_checker = metal_in_porphyrin2(ase_atom, graph)
    all_regions = {}
    atom_pairs_at_breaking_point = {}
    bonds_to_break = []
    carboxylates = find_carboxylates(ase_atom, graph)
    all_sulphates = find_carbonyl_sulphate(ase_atom, graph)
    all_sulfides = find_sulfides(ase_atom, graph)
    all_phosphates = find_phosphate(ase_atom, graph)
    all_phosphites = find_phosphite(ase_atom, graph)
    seen_sulphure = sum(list(all_sulfides.values()), [])
    ferocene_metal = all_ferrocene_metals(ase_atom, graph)
    all_metals = [i.index for i in ase_atom if i.symbol in transition_metals()]
    all_metals = [i for i in all_metals if i not in ferocene_metal + porphyrin_checker]
    for atoms in graph:
        if atoms in list(carboxylates.keys()):
            oxygen = carboxylates[atoms]
            for oxy in oxygen:
                metal = [
                    i for i in graph[oxy] if ase_atom[i].symbol in transition_metals()
                ]
                metal = [i for i in metal if i not in porphyrin_checker]
                metal = [i for i in metal if i not in ferocene_metal]
                if len(metal) > 0:
                    for met in metal:
                        bonds_to_break.append([atoms] + [met])
                        atom_pairs_at_breaking_point[atoms] = met

        if atoms in list(all_sulphates.keys()):
            oxygen = all_sulphates[atoms]
            connected = graph[atoms]
            # all_carbon_indices = [i for i in connected if ase_atom[i].symbol == 'C']
            for oxy in oxygen:
                metal = [
                    i for i in graph[oxy] if ase_atom[i].symbol in transition_metals()
                ]
                metal = [i for i in metal if i not in porphyrin_checker]
                if len(metal) > 0:
                    for met in metal:
                        bonds_to_break.append([oxy] + [met])
                        atom_pairs_at_breaking_point[oxy] = met

        if atoms in list(all_phosphites.keys()):
            all__n_indices = all_phosphites[atoms]
            for phos in all__n_indices:
                metals = [
                    i for i in graph[phos] if ase_atom[i].symbol in transition_metals()
                ]
                for met in metals:
                    atom_pairs_at_breaking_point[phos] = met
                    bonds_to_break.append([phos, met])

        if atoms in list(all_sulfides.keys()):
            sulphure = all_sulfides[atoms]
            connected = graph[atoms]
            for sulf in sulphure:
                metal = [
                    i for i in graph[sulf] if ase_atom[i].symbol in transition_metals()
                ]
                metal = [i for i in metal if i not in porphyrin_checker]
                metal = [i for i in metal if i not in ferocene_metal]
                if len(metal) > 0:
                    for met in metal:
                        bonds_to_break.append([sulf, met])
                        atom_pairs_at_breaking_point[sulf] = met

        if ase_atom[atoms].symbol == 'N':
            if atoms not in porphyrin_checker:
                connected = graph[atoms]
                metal = [
                    i for i in connected if ase_atom[i].symbol in transition_metals()
                ]
                metal = [i for i in metal if i not in ferocene_metal]
                if len(metal) > 0 and atoms not in porphyrin_checker:
                    for met in metal:
                        bonds_to_break.append([atoms, met])
                        atom_pairs_at_breaking_point[atoms] = met
                # atom_pairs_at_breaking_point [metal[0]] = atoms

        if ase_atom[atoms].symbol == 'S':
            if atoms not in seen_sulphure:
                connected = graph[atoms]
                metal = [
                    i for i in connected if ase_atom[i].symbol in transition_metals()
                ]
                metal = [i for i in metal if i not in porphyrin_checker]
                metal = [i for i in metal if i not in ferocene_metal]
                if len(metal) > 0:
                    for met in metal:
                        bonds_to_break.append([atoms] + [met])
                        atom_pairs_at_breaking_point[atoms] = met

        if ase_atom[atoms].symbol == 'O':
            connected = graph[atoms]
            metal = [i for i in connected if ase_atom[i].symbol in transition_metals()]
            metal = [i for i in metal if i not in porphyrin_checker]
            metal = [i for i in metal if i not in ferocene_metal]
            Nitrogen = [i for i in connected if ase_atom[i].symbol == 'N']
            carbon = [i for i in connected if ase_atom[i].symbol == 'C']
            if len(metal) > 0 and len(carbon) == 1:
                for met in metal:
                    atom_pairs_at_breaking_point[atoms] = met
                    bonds_to_break.append([atoms, met])
            if len(metal) > 0 and len(Nitrogen) == 1:
                n_carbon = [
                    i
                    for i in graph[Nitrogen[0]]
                    if ase_atom[i].symbol in ['C', 'S', 'N']
                ]
                if len(n_carbon) > 1:
                    for met in metal:
                        atom_pairs_at_breaking_point[atoms] = met
                    bonds_to_break.append([atoms, met])

            if len(carbon) > 1 and len(metal) > 0:
                for met in metal:
                    bonds_to_break.append([atoms] + [met])
                    atom_pairs_at_breaking_point[atoms] = met

        if ase_atom[atoms].symbol == 'P':
            """
            Find the carbon closest to P, which is not bonded to a metal and cut
            """
            seen_phosphurous = list(all_phosphates.keys()) + sum(
                list(all_phosphites.values()), []
            )
            if atoms not in seen_phosphurous:
                connected = [
                    i
                    for i in graph[atoms]
                    if ase_atom[i].symbol not in transition_metals()
                ]
                metal_oxy = [
                    [i for i in graph[j] if ase_atom[i].symbol in transition_metals()]
                    for j in connected
                ]

                metal = sum(metal_oxy, [])
                metal = [i for i in metal if i not in porphyrin_checker]
                metal = [i for i in metal if i not in ferocene_metal]
                closest_atoms = sum(
                    [
                        [
                            [i, j]
                            for i in graph[j]
                            if i != atoms and ase_atom[i].symbol in transition_metals()
                        ]
                        for j in connected
                    ],
                    [],
                )

                for frag in closest_atoms:
                    bonds_to_break.append(frag)
                    atom_pairs_at_breaking_point[frag[0]] = frag[1]

    # In special cases some carbon and hydrogen get very closed to metals
    # So in such case it is important to clear this
    for metal in all_metals:
        connected = graph[metal]
        c_H = [i for i in connected if ase_atom[i].symbol in ['H', 'C']]
        if len(c_H) > 0:
            for c_h in c_H:
                bonds_to_break.append([c_h, metal])
                atom_pairs_at_breaking_point[metal] = c_h

    for bonds in bonds_to_break:
        bond_matrix[bonds[0], bonds[1]] = 0
        bond_matrix[bonds[1], bonds[0]] = 0

    new_ase_graph = matrix2dict(bond_matrix)
    try:
        list_of_connected_components = connected_components(new_ase_graph)
    except Exception:
        import networkx as nx

        N_Graph = nx.from_dict_of_lists(new_ase_graph)
        list_of_connected_components = [
            list(i) for i in list(nx.connected_components(N_Graph))
        ]

    all_pm_structures = [
        sorted(ase_atom[i].symbols) for i in list_of_connected_components
    ]
    for i in range(len(all_pm_structures)):
        temp = []
        for j in range(len(all_pm_structures)):
            if all_pm_structures[i] == all_pm_structures[j]:
                temp.append(j)
        if temp not in all_regions.values():
            all_regions[i] = temp

    return (
        list_of_connected_components,
        atom_pairs_at_breaking_point,
        porphyrin_checker,
        all_regions,
    )


def is_rodlike(metal_sbu, graph):
    """
    Simple test to check whether a metal sbu is a rodlike MOF
    """
    Rod_check = []
    cells = [(2, 1, 1), (1, 2, 1), (1, 1, 2)]
    for index, ijk in enumerate(cells):
        rod = metal_sbu * ijk
        graph, _ = compute_ase_neighbour(rod)
        list_of_connected_components = connected_components(graph)
        if len(list_of_connected_components) == 1:
            Rod_check.append(index)
    return Rod_check


def all_ferrocene_metals(ase_atom, graph):
    """
    A function to find metals corresponding to ferrocene.
    These metals should not be considered during mof-constructions

    Parameters:
    -----------
    ase_atom: ASE atom
    graph : dictionary containing neigbour lists
    """
    list_of_metals = []
    for atom_index in graph:
        if ase_atom[atom_index].symbol in transition_metals():
            connectivity = graph[atom_index]
            if len(connectivity) >= 10:
                number_of_carbons = 0
                for neigbour in connectivity:
                    if ase_atom[neigbour].symbol == 'C':
                        number_of_carbons += 1
                if number_of_carbons >= 10:
                    list_of_metals.append(atom_index)
    return list_of_metals


def is_ferrocene(metal_sbu, graph):
    """
    A simple script to check whether a metal_sbu is ferrocene
    """
    Check = []
    verdict = False
    all_connectivity = list(graph.values())
    for connectivity in all_connectivity:
        if len(connectivity) >= 10:
            carbons = 0
            for bondedAtomIndex in connectivity:
                if metal_sbu[bondedAtomIndex].symbol == 'C':
                    carbons += 1
            if carbons >= 10:
                verdict = True
            Check.append(verdict)
        Correct = False
        if True in Check:
            Correct = True
    return Correct


def is_paddlewheel(metal_sbu, graph):
    """
    Returns True if the atom is part of a paddlewheel motif
    """
    Check = []
    verdict = False
    All_connectivity = list(graph.values())
    for connectivity in All_connectivity:
        metalNeighbours = 0
        oxygenNeighbours = 0
        for bondedAtomIndex in connectivity:
            if metal_sbu[bondedAtomIndex].symbol in transition_metals():
                metalNeighbours += 1
            if metal_sbu[bondedAtomIndex].symbol == 'O':
                oxygenNeighbours += 1
        if (
            metalNeighbours == 1
            and oxygenNeighbours == 4
            and len(connectivity) >= 5
            and len(connectivity) <= 6
        ):
            verdict = True
        Check.append(verdict)
    Correct = False
    if True in Check:
        Correct = True
    return Correct


def is_paddlewheel_with_water(ase_atom, graph):
    """
    Returns True if the atom is part of a paddle wheel with water motif
    """
    Check = []
    metal = []
    verdict = False
    for atoms in ase_atom:
        if atoms.symbol in transition_metals():
            index = atoms.index
            metal.append(index)
            connectivity = graph[index]
            if (
                len(
                    [
                        ase_atom[i].symbol
                        for i in connectivity
                        if ase_atom[i].symbol == 'O'
                    ]
                )
                == 5
            ):
                verdict = True
                Check.append(verdict)
    Correct = False
    if True in Check and len(metal) == 2:
        Correct = True
    return Correct


def is_uio66(ase_atom, graph):
    """
    Returns True if the atom is part of a UIO66 motif
    """
    Check = []
    verdict = False
    All_connectivity = list(graph.values())
    for connectivity in All_connectivity:
        metalNeighbours = 0
        oxygenNeighbours = 0
        for bondedAtomIndex in connectivity:
            if ase_atom[bondedAtomIndex].symbol in transition_metals():
                metalNeighbours += 1
            if ase_atom[bondedAtomIndex].symbol == 'O':
                oxygenNeighbours += 1
        if (
            metalNeighbours == 4
            and (oxygenNeighbours == 6 or oxygenNeighbours == 8)
            and (len(connectivity) == 10 or len(connectivity) == 12)
        ):
            verdict = True
        Check.append(verdict)
        Correct = False
        if True in Check:
            Correct = True
    return Correct


def is_irmof(ase_atom, graph):
    """
    Returns True if the atom is part of a IRMOF motif
    """
    Check = []
    verdict = False
    for atoms in ase_atom:
        if atoms.symbol == 'O':
            index = atoms.index
            connectivity = graph[index]
            if (
                len(connectivity) == 4
                and len(
                    [
                        ase_atom[i].symbol
                        for i in connectivity
                        if ase_atom[i].symbol in transition_metals()
                    ]
                )
                == 4
            ):
                verdict = True
                Check.append(verdict)
    Correct = False
    if True in Check:
        Correct = True
    return Correct


def is_mof32(ase_atom, graph):
    """
    Returns True if the atom is part of a MOF32 motif
    """
    Check = []
    metal = []
    verdict = False
    for atoms in ase_atom:
        if atoms.symbol in transition_metals():
            index = atoms.index
            metal.append(index)
            connectivity = graph[index]
            if (
                len(
                    [
                        ase_atom[i].symbol
                        for i in connectivity
                        if ase_atom[i].symbol == 'O'
                    ]
                )
                == 8
            ):
                verdict = True
                Check.append(verdict)
    Correct = False
    if True in Check and len(metal) == 1:
        Correct = True
    return Correct


def rod_manipulation(ase_atom, checker):
    """
    Script to adjust Rodlike sbus.
    1) Its collects the axis responsible for expanding the rod
    2) It shifts all coordinates to the axis
    3) It rotates the rod to lie in the directions of expansion,
    """

    cell = ase_atom.get_cell().tolist()
    value_index = checker[0]
    cell_value = cell[value_index]
    new_cell = move2front(checker, cell_value)
    coords = ase_atom.positions - cell_value
    new_position = move2front(checker, coords.tolist())
    ase_atom.positions = new_position
    return ase_atom, new_cell


def find_unique_building_units(
    list_of_connected_components,
    atom_pairs_at_breaking_point,
    ase_atom,
    porphyrin_checker,
    all_regions,
    wrap_system=True,
    cheminfo=False,
):
    """
    Find Unique components
    Returns a list of unique molecules
    """
    mof_metal = []
    mof_linker = []
    building_unit_regions = {}
    for idx, key in enumerate(all_regions.keys()):
        frag = list(all_regions.keys())[idx]
        components = list_of_connected_components[all_regions[frag][0]]
        all_breaking_point = list(atom_pairs_at_breaking_point.keys()) + list(
            atom_pairs_at_breaking_point.values()
        )
        point_of_extension = [i for i in all_breaking_point if i in components]
        mapped_indices = dict(
            [(i, j) for i, j in zip(components, range(len(components)))]
        )
        molecule_to_write = ase_atom[components]
        molecule_to_write.info['point_of_extension'] = [
            mapped_indices[i] for i in point_of_extension
        ]
        if wrap_system:
            molecule_to_write = wrap_systems_in_unit_cell(molecule_to_write)
        # if cheminfo:
        #     smi, inchi, inchikey = compute_openbabel_cheminformatic(
        #         molecule_to_write)
        #     molecule_to_write.info['smi'] = smi
        #     molecule_to_write.info['inchi'] = str(inchi)
        #     molecule_to_write.info['inchikey'] = str(inchikey)
        molecule_to_write.info['atom_indices_mapping'] = [
            list_of_connected_components[i] for i in all_regions[key]
        ]
        metal = [
            i.index
            for i in molecule_to_write
            if i.symbol in transition_metals() and i.index not in porphyrin_checker
        ]
        non_ferocene_metal = []
        if len(metal) > 0:
            graph_sbu, _ = compute_ase_neighbour(molecule_to_write)
            """
            Check whether the metal sbu is a rod mof. If is it is rod mof,
            we rotate and aligne the sbu such that the axis of rotation will be the a-axis.
            """
            if len(is_rodlike(molecule_to_write, graph_sbu)) == 1:
                molecule_to_write.info['sbu_type'] = 'rodlike'
            elif is_paddlewheel(molecule_to_write, graph_sbu):
                molecule_to_write.info['sbu_type'] = 'paddlewheel'
            elif is_paddlewheel_with_water(molecule_to_write, graph_sbu):
                molecule_to_write.info['sbu_type'] = 'paddlewheel_with_water'
            elif is_uio66(molecule_to_write, graph_sbu):
                molecule_to_write.info['sbu_type'] = 'UIO66_sbu'
            elif is_mof32(molecule_to_write, graph_sbu):
                molecule_to_write.info['sbu_type'] = 'MOF32_sbu'
            elif is_irmof(molecule_to_write, graph_sbu):
                molecule_to_write.info['sbu_type'] = 'IRMOF_sbu'
            elif is_ferrocene(molecule_to_write, graph_sbu):
                molecule_to_write.info['sbu_type'] = 'ferrocenelike'
                ferocene_metal = all_ferrocene_metals(molecule_to_write, graph_sbu)
                non_ferocene_metal = [i for i in metal if i not in ferocene_metal]
                if len(non_ferocene_metal) == 0:
                    mof_linker.append(molecule_to_write)
                    continue
            else:
                molecule_to_write.info['sbu_type'] = 'still checking!'
            mof_metal.append(molecule_to_write)

        else:
            mof_linker.append(molecule_to_write)
        building_unit_regions[idx] = molecule_to_write

    return mof_metal, mof_linker, building_unit_regions


def metal_coordination_number(ase_atom):
    """
    Extract coordination number of central metal
    """
    metal_coordination = []
    graph, _ = compute_ase_neighbour(ase_atom)
    metal_indices = [i.index for i in ase_atom if i.symbol in transition_metals()]
    metal_elt = []
    for i in metal_indices:
        if ase_atom[i].symbol not in metal_elt:
            metal_elt.append(ase_atom[i].symbol)
            metal_coordination.append(
                Coordination(
                    element=ase_atom[i].symbol, coordination_number=len(graph[i])
                )
            )
    return metal_elt, metal_coordination


def wrap_systems_in_unit_cell(ase_atom, max_iter=30, skin=0.3):
    """
    A simple aglorithm to reconnnect all atoms wrapped in a periodic boundary condition such that all atoms outside the box will appear reconnected.
    """
    new_position = geometry.wrap_positions(
        ase_atom.positions,
        ase_atom.cell,
        pbc=True,
        center=(0, 0, 0),
        pretty_translation=True,
        eps=1e-07,
    )
    ase_atom.positions = new_position

    graph, bond_matrix = compute_ase_neighbour(ase_atom)

    for atom in graph:
        connected = graph[atom]
        for nl in connected:
            check = (
                covalent_radius(ase_atom[atom].symbol)
                + covalent_radius(ase_atom[nl].symbol)
                + skin
            )
            bond = round(ase_atom.get_distance(atom, nl), 2)
            if bond > check:
                bond_matrix[atom][nl] = 0

    new_ase_graph = matrix2dict(bond_matrix)
    list_of_connected_components = connected_components(new_ase_graph)
    number_of_iterations = 0
    while len(list_of_connected_components) != 1:
        all_len = [len(i) for i in list_of_connected_components]
        max_index = all_len.index(max(all_len))
        Root = list_of_connected_components[max_index]
        list_of_connected_components.pop(max_index)
        All_sum = sum(list_of_connected_components, [])
        for atom in Root:
            connected = graph[atom]
            for nl in All_sum:
                if nl in connected:
                    v = ase_atom[nl].position - ase_atom[atom].position
                    mic_vector = geometry.find_mic(v, ase_atom.get_cell(), pbc=True)
                    ase_atom[nl].position = mic_vector[0] + ase_atom[atom].position

        graph, bond_matrix = compute_ase_neighbour(ase_atom)
        for atom in graph:
            connected = graph[atom]
            for nl in connected:
                check = (
                    covalent_radius(ase_atom[atom].symbol)
                    + covalent_radius(ase_atom[nl].symbol)
                    + skin
                )
                bond = round(ase_atom.get_distance(atom, nl), 2)
                if bond > check:
                    bond_matrix[atom][nl] = 0
        new_ase_graph = matrix2dict(bond_matrix)
        list_of_connected_components = connected_components(new_ase_graph)
        number_of_iterations += 1
        if number_of_iterations == max_iter:
            break
    return ase_atom
