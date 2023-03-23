#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.'
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

from typing import Dict, List, Tuple, Optional
from collections import defaultdict

from ase import Atoms
from ase.data import chemical_symbols
from ase.geometry.cell import complete_cell
import numpy as np
import matid.geometry  # pylint: disable=import-error
from matid.classification.structureclusterer import StructureClusterer  # pylint: disable=import-error
from matid import Classifier  # pylint: disable=import-error
from matid.classifications import Class0D, Atom, Class1D, Material2D, Surface, Class3D, Class2D, Unknown  # pylint: disable=import-error
from matid.symmetry.symmetryanalyzer import SymmetryAnalyzer  # pylint: disable=import-error

from nomad import atomutils
from nomad.units import ureg
from nomad.atomutils import Formula
from nomad.datamodel.results import Symmetry, Material, System, Relation, Prototype
from nomad.datamodel.metainfo.simulation.system import Atoms as NOMADAtoms
from nomad.normalizing.common import (
    cell_from_ase_atoms,
    ase_atoms_from_nomad_atoms,
    nomad_atoms_from_ase_atoms,
    structures_2d,
    material_id_bulk,
    material_id_2d,
)

SYMBOLS = np.array(chemical_symbols)


def get_topology_id(index: int) -> str:
    '''Retuns a valid topology identifier with the given index.
    Args:
        index: The index of the topology. Must be unique.

    Returns:
        An identifier string that can be stored in topology.system_id.
    '''
    return f'results/material/topology/{index}'


def get_topology_original(material: Material, atoms: NOMADAtoms = None) -> System:
    '''
    Creates a new topology item for the original structure.
    '''
    original = System(
        method='parser',
        label='original',
        description='A representative system chosen from the original simulation.',
        material_id=material.material_id,
        material_name=material.material_name,
        structural_type=material.structural_type,
        functional_type=material.functional_type,
        compound_type=material.compound_type,
        chemical_formula_hill=material.chemical_formula_hill,
        chemical_formula_iupac=material.chemical_formula_iupac,
        chemical_formula_anonymous=material.chemical_formula_anonymous,
        chemical_formula_reduced=material.chemical_formula_reduced,
        elements=material.elements,
        atoms_ref=atoms
    )

    return original


def add_system_info(system: System, topologies: Dict[str, System]) -> None:
    '''Given a system with minimal information, attempts to add all values than
    can be derived.
    '''
    def get_atoms(system):
        '''Resolves the atoms that the system is constructed from.'''
        if system.atoms:
            return system.atoms
        if system.atoms_ref:
            return system.atoms_ref
        if system.indices is not None:
            return get_atoms(topologies[system.parent_system])

    atoms = get_atoms(system)
    if atoms:
        if system.cell is None:
            if system.atoms or system.atoms_ref:
                ase_atoms = ase_atoms_from_nomad_atoms(atoms)
                system.cell = cell_from_ase_atoms(ase_atoms)
        symbols = SYMBOLS[atoms.species if atoms.species is not None else atoms.atomic_numbers]
        if system.indices is not None:
            symbols = symbols[system.indices[0]]
        n_atoms = len(symbols)
        if system.n_atoms is None:
            system.n_atoms = n_atoms
        try:
            formula = Formula(''.join(symbols))
        except Exception:
            pass
        else:
            if system.chemical_formula_hill is None:
                system.chemical_formula_hill = formula.format('hill')
            if system.chemical_formula_iupac is None:
                system.chemical_formula_iupac = formula.format('iupac')
            if system.chemical_formula_reduced is None:
                system.chemical_formula_reduced = formula.format('reduced')
            if system.chemical_formula_anonymous is None:
                system.chemical_formula_anonymous = formula.format('anonymous')
            if not system.elements:
                system.elements = formula.elements()


def add_system(system: System, topologies: Dict[str, System], parent: Optional[System] = None) -> None:
    '''Adds the given system to the topology.
    '''
    index = len(topologies)
    system.system_id = get_topology_id(index)
    if parent:
        children = parent.child_systems if parent.child_systems else []
        children.append(system.system_id)
        if parent.child_systems is not children:
            parent.child_systems = children
        system.parent_system = parent.system_id
    topologies[system.system_id] = system


class TopologyNormalizer():
    '''Handles the creation of topology information.
    '''
    def __init__(self, entry_archive, repr_system, logger):
        self.entry_archive = entry_archive
        self.repr_system = repr_system
        self.structural_type = None
        self.logger = logger

    def topology(self, material) -> Optional[List[System]]:
        '''Returns a dictionary that contains all of the topologies mapped by
        topology id.'''
        # Use the calculation topology primarily
        topology_calc = self.topology_calculation(material)
        if topology_calc:
            return topology_calc
        topology_matid = self.topology_matid(material)
        if topology_matid:
            return topology_matid

        return None

    def topology_calculation(self, material: Material) -> Optional[List[System]]:
        '''Extracts the system topology as defined in the original calculation.
        This topology typically comes from e.g. classical force fields that
        define a topology for the system.
        '''
        try:
            groups = self.entry_archive.run[0].system[0].atoms_group
            if len(groups) == 0:
                return None
        except Exception:
            return None
        try:
            atoms = self.repr_system.atoms
        except Exception:
            atoms = None
        if not atoms:
            return None

        topology: Dict[str, System] = {}
        original = get_topology_original(material, atoms)
        add_system(original, topology)
        label_to_indices: Dict[str, list] = defaultdict(list)

        def add_group(groups, parent=None):
            if not groups: return
            for group in groups:
                label = group.label
                # Groups with the same label are mapped to the same system.
                old_labels = label_to_indices[label]
                instance_indices = group.atom_indices
                if not len(old_labels):
                    description_map = {
                        'molecule': 'Molecule extracted from the calculation topology.',
                        'molecule_group': 'Group of molecules extracted from the calculation topology.',
                        'monomer_group': 'Group of monomers extracted from the calculation topology.',
                        'monomer': 'Monomer extracted from the calculation topology.'
                    }
                    structural_type_map = {
                        'molecule': 'molecule',
                        'molecule_group': 'group',
                        'monomer': 'monomer',
                        'monomer_group': 'group',
                    }
                    system = System(
                        method='parser',
                        description=description_map.get(group.type, None),
                        label=group.label,
                        structural_type=structural_type_map.get(group.type, None),
                        system_relation=Relation(type='subsystem'),
                    )
                    add_system(system, topology, parent)
                    add_group(group.atoms_group, system)
                    old_labels.append(instance_indices)
                else:
                    if len(old_labels[0]) == len(instance_indices):
                        old_labels.append(instance_indices)
                    else:
                        self.logger.warn((
                            'the topology contains entries with the same label but with '
                            'different number of atoms'
                        ))

        add_group(groups, original)

        # Add the derived system information once all indices etc. are gathered.
        for top in topology.values():
            top.indices = label_to_indices.get(top.label)
            add_system_info(top, topology)

        return list(topology.values())

    def topology_matid(self, material: Material) -> Optional[List[System]]:
        # Topology is currently created only for systems that are classified as 2D
        # or surface. TODO: Also other entries should get a topology: all of
        # the symmetry analysis and additional system partitioning through MatID
        # should be stored in the topology.
        '''
        Returns a list of the identified systems with topological relations and
        classification of subsystems
        '''
        if material.structural_type not in {'2D', 'surface', 'unavailable'}:
            return None
        try:
            nomad_atoms = self.repr_system.atoms
            atoms = ase_atoms_from_nomad_atoms(nomad_atoms)
        except Exception:
            return None
        # TODO: MatID does not currently support non-periodic structures
        if not any(atoms.pbc):
            return None
        # TODO: Currently the topology creation is skipped completely. Needs to
        # be re-enabled once we work out performance problems.
        n_atoms = len(atoms)
        if n_atoms > 0:
            return None

        topology: Dict[str, System] = {}
        original = get_topology_original(material, nomad_atoms)
        add_system(original, topology)
        add_system_info(original, topology)

        # Add all meaningful clusters to the topology
        clusters = self._perform_matid_clustering(atoms)
        cluster_indices_list, cluster_symmetries = self._filter_clusters(clusters)
        for indices, symm in zip(cluster_indices_list, cluster_symmetries):
            subsystem_atoms = atoms[indices]
            subsystem = self._create_subsystem(subsystem_atoms, indices)
            structural_type = subsystem.structural_type
            # TODO: Currently only 2D and surface subsystems are included
            if structural_type in {'surface', '2D'}:
                add_system(subsystem, topology, original)
                add_system_info(subsystem, topology)
                conventional_cell = self._create_conv_cell_system(symm, structural_type)
                add_system(conventional_cell, topology, subsystem)
                add_system_info(conventional_cell, topology)

        return list(topology.values())

    def _create_subsystem(self, atoms: Atoms, indices: List[int]) -> System:
        '''
        Creates a new subsystem as detected by MatID.
        '''
        subsystem = System(
            method='matid',
            label='subsystem',
            description='Automatically detected subsystem.',
            system_relation=Relation(type='subsystem'),
            indices=[indices]
        )
        subsystem.structural_type = self._system_type_analysis(atoms)

        return subsystem

    def _create_conv_cell_system(self, symm, structural_type: str):
        '''
        Creates a new topology item for a conventional cell.
        '''
        symmsystem = System(
            method='matid',
            label='conventional cell',
            system_relation=Relation(type='subsystem'),
        )
        conv_system = symm.get_conventional_system()
        if structural_type == 'surface':
            symmsystem.description = 'The conventional cell of the bulk material from which the surface is constructed from.'
            symmsystem.structural_type = 'bulk'
            symmsystem.atoms = nomad_atoms_from_ase_atoms(conv_system)
        elif structural_type == '2D':
            symmsystem.description = 'The conventional cell of the 2D material.'
            symmsystem.structural_type = '2D'
            symmsystem.atoms = nomad_atoms_from_ase_atoms(conv_system)

        if structural_type == 'surface':
            symmsystem = self._create_symmsystem_surface(symm, symmsystem)
        elif structural_type == '2D':
            symmsystem = self._create_symmsystem_2D(symm, symmsystem)
        return symmsystem

    def _perform_matid_clustering(self, atoms: Atoms) -> List:
        '''
        Performs the clustering with MatID
        '''
        clusterer = StructureClusterer()
        clusters = clusterer.get_clusters(
            atoms,
            max_cell_size=6,
            pos_tol=0.70,
            angle_tol=20,
            merge_threshold=0.3,
            merge_radius=5
        )
        return clusters

    def _filter_clusters(self, clusters: StructureClusterer) -> Tuple[List[List[int]], List[Symmetry]]:
        '''
        Filters all clusters < 2 atoms and creates a cluster indices list of the remaining
        clusters
        '''
        # TODO: Add proper decision mechanism for identifying when the detected
        # cluster are OK. For now, we filter out cluster that have very small
        # number of atoms: we do not want to show e.g. random outlier atoms in
        # the topology. This will in general not work e.g. for periodic systems
        # where the cluster size is not indicative of the actual physical size
        # due to periodic boundaries. TODO: Improve this, it is currently based
        # only on the number of atoms in the cluster. Maybe the cluster given
        # out by MatID could already be grouped a bit better, e.g. all outliers
        # would be grouped into one cluster?
        filtered_cluster = filter(lambda x: len(x.indices) > 1, clusters)
        cluster_indices_list: List[List[int]] = []
        cluster_symmetries: List[Symmetry] = []
        for cluster in filtered_cluster:
            indices = list(cluster.indices)
            cluster_indices_list += [indices]
            regions = cluster.regions
            number_of_atoms: List[int] = []
            for region in regions:
                if region:
                    number_of_atoms.append(region.cell.get_number_of_atoms())
                else:
                    number_of_atoms.append(-1)

            # If there are 2 regions with the same size, the one with the smaller index is selected
            largest_region_index = number_of_atoms.index(max(number_of_atoms))
            largest_region_system = regions[largest_region_index].cell

            # TODO: only SymmetryAnalyzer for 2D and surface
            symm = SymmetryAnalyzer(largest_region_system)
            cluster_symmetries += [symm]
        return cluster_indices_list, cluster_symmetries

    def _system_type_analysis(self, atoms: Atoms) -> matid.classifications:
        '''
        Classifies ase.atoms and returns the MatID system type as a string.
        '''
        system_types = {Class3D: 'bulk',
                        Atom: 'atom',
                        Class0D: 'molecule / cluster',
                        Class1D: '1D',
                        Class2D: 'unavailable',
                        Surface: 'surface',
                        Material2D: '2D',
                        Unknown: 'unavailable'}
        try:
            classifier = Classifier(radii='covalent', cluster_threshold=1)
            cls = classifier.classify(atoms)
        except Exception as e:
            self.logger.error(
                'matid project system classification failed', exc_info=e, error=str(e))
            return 'unavailable'
        else:
            classification = type(cls)
            try:
                system_type = system_types[classification]
            except Exception as e:
                self.logger.error(
                    'matid project system classification unavailable', exc_info=e, error=str(e))
                system_type = 'unavailable'
        return system_type

    def _create_symmsystem_surface(self, symm: SymmetryAnalyzer, subsystem: System) -> System:
        '''
        Creates the subsystem with the symmetry information of the conventional cell
        '''
        conv_system = symm.get_conventional_system()
        prototype = self._create_prototype(symm, conv_system)
        spg_number = symm.get_space_group_number()
        subsystem.prototype = prototype
        subsystem.cell = cell_from_ase_atoms(conv_system)
        symmetry = self._create_symmetry(symm)
        wyckoff_sets = symm.get_wyckoff_sets_conventional()
        material_id = material_id_bulk(spg_number, wyckoff_sets)
        subsystem.material_id = material_id
        subsystem.symmetry = symmetry
        return subsystem

    def _create_symmsystem_2D(self, symm: SymmetryAnalyzer, subsystem: System) -> System:
        '''
        Creates the subsystem with the symmetry information of the conventional cell
        '''
        subsystem_atoms = Atoms(
            symbols=subsystem.atoms.labels,
            positions=subsystem.atoms.positions.to(ureg.angstrom),
            cell=complete_cell(subsystem.atoms.lattice_vectors.to(ureg.angstrom)),
            pbc=np.array(subsystem.atoms.periodic, dtype=bool)
        )
        conv_atoms, __, wyckoff_sets, spg_number = structures_2d(subsystem_atoms)
        subsystem.cell = cell_from_ase_atoms(conv_atoms)

        # Here we zero out the irrelevant lattice parameters to correctly handle
        # 2D systems with nonzero thickness (e.g. MoS2).
        if subsystem.cell.c:
            subsystem.cell.c = 0.0
        if subsystem.cell.alpha:
            subsystem.cell.alpha = None
        if subsystem.cell.beta:
            subsystem.cell.beta = None

        prototype = self._create_prototype(symm, conv_atoms)
        subsystem.prototype = prototype
        subsystem.material_id = material_id_2d(spg_number, wyckoff_sets)
        symmetry = None
        subsystem.symmetry = symmetry
        return subsystem

    def _create_symmetry(self, symm: SymmetryAnalyzer) -> Symmetry:
        international_short = symm.get_space_group_international_short()

        sec_symmetry = Symmetry()
        sec_symmetry.symmetry_method = 'MatID'
        sec_symmetry.space_group_number = symm.get_space_group_number()
        sec_symmetry.space_group_symbol = international_short
        sec_symmetry.hall_number = symm.get_hall_number()
        sec_symmetry.hall_symbol = symm.get_hall_symbol()
        sec_symmetry.international_short_symbol = international_short
        sec_symmetry.point_group = symm.get_point_group()
        sec_symmetry.crystal_system = symm.get_crystal_system()
        sec_symmetry.bravais_lattice = symm.get_bravais_lattice()
        sec_symmetry.origin_shift = symm._get_spglib_origin_shift()
        sec_symmetry.transformation_matrix = symm._get_spglib_transformation_matrix()
        return sec_symmetry

    def _create_prototype(self, symm: SymmetryAnalyzer, conv_system: System) -> Prototype:
        spg_number = symm.get_space_group_number()
        atom_species = conv_system.get_atomic_numbers()
        if type(conv_system) == Atoms or conv_system.wyckoff_letters is None:
            wyckoffs = symm.get_wyckoff_letters_conventional()
        else:
            wyckoffs = conv_system.wyckoff_letters
        norm_wyckoff = atomutils.get_normalized_wyckoff(atom_species, wyckoffs)
        protoDict = atomutils.search_aflow_prototype(spg_number, norm_wyckoff)

        if protoDict is not None:
            aflow_prototype_name = protoDict['Prototype']
            aflow_strukturbericht_designation = protoDict['Strukturbericht Designation']
            prototype_label = '%d-%s-%s' % (
                spg_number,
                aflow_prototype_name,
                protoDict.get('Pearsons Symbol', '-')
            )
            prototype = Prototype()
            prototype.label = prototype_label

            prototype.formula = Formula(''.join(protoDict['atom_labels'])).format('hill')
            prototype.aflow_id = protoDict['aflow_prototype_id']
            prototype.aflow_url = protoDict['aflow_prototype_url']
            prototype.assignment_method = 'normalized-wyckoff'
            prototype.m_cache['prototype_notes'] = protoDict['Notes']
            prototype.m_cache['prototype_name'] = aflow_prototype_name
            if aflow_strukturbericht_designation != 'None':
                prototype.m_cache['strukturbericht_designation'] = aflow_strukturbericht_designation
        else:
            prototype = None
        return prototype
