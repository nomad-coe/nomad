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

from typing import Dict, List, Optional
from collections import defaultdict
import pathlib
import json

from ase import Atoms
from ase.data import chemical_symbols
import numpy as np
from matid.clustering import Clusterer, Cluster, Classification  # pylint: disable=import-error
from matid.symmetry.symmetryanalyzer import SymmetryAnalyzer  # pylint: disable=import-error
from matid.classifications import Class0D, Atom, Class1D, Class2D, Material2D, Surface, Class3D  # pylint: disable=import-error

from nomad import utils
from nomad import config
from nomad import atomutils
from nomad.atomutils import Formula
from nomad.datamodel.results import SymmetryNew as Symmetry, Material, System, Relation, structure_name_map
from nomad.datamodel.metainfo.simulation.system import Atoms as NOMADAtoms
from nomad.datamodel.datamodel import EntryArchive
from nomad.normalizing.common import (
    cell_from_ase_atoms,
    ase_atoms_from_nomad_atoms,
    nomad_atoms_from_ase_atoms,
    wyckoff_sets_from_matid,
    structures_2d,
    material_id_bulk,
    material_id_2d,
    material_id_1d,
)

conventional_description = 'The conventional cell of the material from which the subsystem is constructed from.'
subsystem_description = 'Automatically detected subsystem.'
chemical_symbols = np.array(chemical_symbols)
with open(pathlib.Path(__file__).parent / 'data/top_50k_material_ids.json', "r") as fin:
    top_50k_material_ids = json.load(fin)


def get_topology_id(index: int) -> str:
    '''Retuns a valid topology identifier with the given index.
    Args:
        index: The index of the topology. Must be unique.

    Returns:
        An identifier string that can be stored in topology.system_id.
    '''
    return f'results/material/topology/{index}'


def get_topology_original(atoms: NOMADAtoms = None, archive: EntryArchive = None) -> System:
    '''
    Creates a new topology item for the original structure.
    '''
    dimensionality = None
    try:
        classification = archive.run[0].m_cache["classification"]
    except Exception:
        pass
    else:
        dimensionality_map = {
            Class0D: 0,
            Atom: 0,
            Class1D: 1,
            Class2D: 2,
            Surface: 2,
            Material2D: 2,
            Class3D: 3,
        }
        dimension = dimensionality_map.get(classification)
        if dimension is not None:
            dimensionality = f'{dimension}D'

    original = System(
        method='parser',
        label='original',
        description='A representative system chosen from the original simulation.',
        dimensionality=dimensionality,
        system_relation=Relation(type='root'),
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
        atomic_numbers = atoms.species if atoms.species is not None else atoms.atomic_numbers
        symbols = chemical_symbols[atomic_numbers]

        if system.indices is not None:
            system.atoms_ref = atoms
            total_mass = atomutils.get_summed_atomic_mass(atomic_numbers)
            total_atoms = len(atomic_numbers)
            symbols = symbols[system.indices[0]]
            system.atomic_fraction = len(symbols) / total_atoms
            sub_mass = atomutils.get_summed_atomic_mass(atomic_numbers[system.indices[0]])
            system.mass_fraction = sub_mass / total_mass
        n_atoms = len(symbols)
        if system.n_atoms is None:
            system.n_atoms = n_atoms
        try:
            formula = Formula(''.join(symbols))
        except Exception:
            pass
        else:
            formula.populate(system, descriptive_format='descriptive')


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
    def __init__(self, entry_archive, repr_system, repr_symmetry, conv_atoms, logger):
        self.entry_archive = entry_archive
        self.repr_system = repr_system
        self.repr_symmetry = repr_symmetry
        self.structural_type = None
        self.conv_atoms = conv_atoms
        self.logger = logger

    def topology(self, material) -> Optional[List[System]]:
        '''Returns a dictionary that contains all of the topologies mapped by id.
        '''
        # If topology already exists (e.g. written by another normalizer), do
        # not overwrite it.
        topology = self.entry_archive.m_xpath('results.material.topology')
        if topology:
            return None
        # Next use the topology from the calculation
        topology_calc = self.topology_calculation()
        if topology_calc:
            return topology_calc
        # Finally if no other topology exists, try creating one with MatID
        with utils.timer(self.logger, 'calculating topology with matid'):
            topology_matid = self.topology_matid(material)
        if topology_matid:
            return topology_matid

        return None

    def topology_calculation(self) -> Optional[List[System]]:
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
        if atoms is None:
            return None
        elif atoms.positions is None or len(atoms.positions) == 0:
            return None
        elif (atoms.species is None or len(atoms.species) == 0) and (atoms.atomic_numbers is None or len(atoms.atomic_numbers) == 0):
            return None

        topology: Dict[str, System] = {}
        original = get_topology_original(atoms, self.entry_archive)
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
                        'monomer': 'Monomer extracted from the calculation topology.',
                        'projection': 'Atom(s) considered to obtain the projected tight-binding model.',
                        'core_hole': 'Atom with the core-hole state.'
                    }
                    structural_type_map = {
                        'molecule': 'molecule',
                        'molecule_group': 'group',
                        'monomer': 'monomer',
                        'monomer_group': 'group',
                        'projection': 'group',
                        'core_hole': 'group'
                    }
                    building_block_map = {
                        'molecule': 'molecule',
                        'monomer': 'monomer',
                    }
                    relation_map = {
                        'molecule': 'subsystem',
                        'molecule_group': 'group',
                        'monomer': 'subsystem',
                        'monomer_group': 'group',
                        'projection': 'group',
                        'core_hole': 'group'
                    }
                    system = System(
                        method='parser',
                        description=description_map.get(group.type),
                        label=group.label,
                        structural_type=structural_type_map.get(group.type),
                        building_block=building_block_map.get(group.type),
                        system_relation=Relation(type=relation_map.get(group.type)),
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
        '''
        Returns a list of systems that have been identified with MatID.
        '''
        # See if a system is available
        try:
            nomad_atoms = self.repr_system.atoms
            atoms = ase_atoms_from_nomad_atoms(nomad_atoms)
        except Exception:
            return None
        if not atoms or len(atoms) == 0:
            return None

        # Create topology for the original system
        topology: Dict[str, System] = {}
        original = get_topology_original(nomad_atoms, self.entry_archive)
        add_system(original, topology)
        add_system_info(original, topology)

        # Since we still need to run the old classification code
        # (matid.classification.classify), we use it's results to populate the
        # topology for bulk and 1D systems. Also the new clustering cannot
        # currently be run for systems without a cell. In other cases we run the
        # new classification code (matid.clustering.clusterer).
        n_atoms = len(atoms)
        cell = atoms.get_cell()
        if material.structural_type == 'bulk':
            self._topology_bulk(original, topology)
        elif material.structural_type == '1D':
            self._topology_1d(original, topology)
        elif cell is None or cell.volume == 0:
            pass
        # Continue creating topology if system size is not too large
        elif n_atoms <= config.normalize.clustering_size_limit:
            # Add all meaningful clusters to the topology
            clusterer = Clusterer()
            clusters = clusterer.get_clusters(atoms, pos_tol=0.8)
            for cluster in clusters:
                subsystem = self._create_subsystem(cluster)
                if not subsystem:
                    continue
                structural_type = subsystem.structural_type
                # If the found cell has many basis atoms, it is more likely that
                # some of the symmetries were not correctly found than the cell
                # actually being very complicated. Thus we ignore these clusters to
                # minimize false-positive and to limit the time spent on symmetry
                # calculation.
                cell = cluster.cell()
                if len(cell) > 6:
                    self.logger.info(f"cell with many atoms ({len(cell)}) was ignored")
                    continue
                try:
                    conventional_cell = self._create_conv_cell_system(cluster, structural_type)
                except Exception as e:
                    self.logger.error(
                        "conventional cell information could not be created",
                        exc_info=e,
                        error=str(e)
                    )
                    continue
                # We only accept the subsystem if the material id exists in the top
                # 50k materials with most entries attached to them. This ensures
                # that the material_id link points to valid materials and that we
                # don't report anything too weird. The top 50k materials are
                # pre-stored in a pickle file that has been created by using ES
                # terms aggregation.
                if conventional_cell.material_id in top_50k_material_ids:
                    add_system(subsystem, topology, original)
                    add_system_info(subsystem, topology)
                    add_system(conventional_cell, topology, subsystem)
                    add_system_info(conventional_cell, topology)
                else:
                    self.logger.info(f"material_id {conventional_cell.material_id} could not be verified")

        return list(topology.values())

    def _topology_bulk(self, original, topology) -> None:
        '''Creates a topology for bulk structures as detected by the old matid
        classification.'''
        if self.conv_atoms is None:
            return None

        # Subsystem
        subsystem = System(
            method='matid',
            label='subsystem',
            dimensionality='3D',
            structural_type='bulk',
            description=subsystem_description,
            system_relation=Relation(type='subsystem'),
            indices=[list(range(original.n_atoms))]
        )
        add_system(subsystem, topology, original)
        add_system_info(subsystem, topology)

        # Conventional system
        conv_system = System(
            method='matid',
            label='conventional cell',
            system_relation=Relation(type='conventional_cell'),
            dimensionality='3D',
            structural_type='bulk',
            description=conventional_description
        )
        conv_system.atoms = nomad_atoms_from_ase_atoms(self.conv_atoms)
        symmetry_analyzer = self.repr_symmetry.m_cache.get("symmetry_analyzer")
        conv_system.symmetry = self._create_symmetry(symmetry_analyzer)
        conv_system.cell = cell_from_ase_atoms(self.conv_atoms)
        conv_system.material_id = material_id_bulk(
            symmetry_analyzer.get_space_group_number(),
            symmetry_analyzer.get_wyckoff_sets_conventional()
        )
        add_system(conv_system, topology, subsystem)
        add_system_info(conv_system, topology)

    def _topology_1d(self, original, topology):
        '''Creates a topology for 1D structures as detected by the old matid
        classification.'''
        if self.conv_atoms is None:
            return None

        # Subsystem
        subsystem = System(
            method='matid',
            label='subsystem',
            dimensionality='1D',
            structural_type='1D',
            description=subsystem_description,
            system_relation=Relation(type='subsystem'),
            indices=[list(range(original.n_atoms))]
        )
        add_system(subsystem, topology, original)
        add_system_info(subsystem, topology)

        # Conventional system
        conv_system = System(
            method='matid',
            label='conventional cell',
            system_relation=Relation(type='conventional_cell'),
            dimensionality='1D',
            structural_type='1D',
        )
        conv_system.atoms = nomad_atoms_from_ase_atoms(self.conv_atoms)
        conv_system.cell = cell_from_ase_atoms(self.conv_atoms)

        # The lattice parameters that are not well defined for 1D structures are unset
        conv_system.cell.b = None
        conv_system.cell.c = None
        conv_system.cell.alpha = None
        conv_system.cell.beta = None
        conv_system.cell.gamma = None
        conv_system.cell.atomic_density = None
        conv_system.cell.mass_density = None
        conv_system.cell.volume = None

        conv_system.material_id = material_id_1d(self.conv_atoms)
        add_system(conv_system, topology, subsystem)
        add_system_info(conv_system, topology)

    def _create_subsystem(self, cluster: Cluster) -> Optional[System]:
        '''
        Creates a new subsystem as detected by MatID.
        '''
        try:
            dimensionality = cluster.dimensionality()
            classification = cluster.classification()
        except Exception as e:
            self.logger.error(
                'matid system classification failed', exc_info=e, error=str(e)
            )
            return None
        structural_type_map = {
            Classification.Class3D: 'bulk',
            Classification.Surface: 'surface',
            Classification.Material2D: '2D',
        }
        structural_type = structural_type_map.get(classification)
        if not structural_type:
            return None
        building_block_map = {
            Classification.Surface: 'surface',
            Classification.Material2D: '2D material',
        }
        subsystem = System(
            method='matid',
            label='subsystem',
            description=subsystem_description,
            system_relation=Relation(type='subsystem'),
            indices=[list(cluster.indices)]
        )
        subsystem.structural_type = structural_type
        subsystem.dimensionality = f'{dimensionality}D'
        subsystem.building_block = building_block_map.get(classification)

        return subsystem

    def _create_conv_cell_system(self, cluster: Cluster, structural_type: str):
        '''
        Creates a new topology item for a conventional cell.
        '''
        symmsystem = System(
            method='matid',
            label='conventional cell',
            system_relation=Relation(type='conventional_cell'),
        )
        if structural_type == '2D':
            self._add_conventional_2d(cluster, symmsystem)
        else:
            self._add_conventional_bulk(cluster, symmsystem)
        symmsystem.description = conventional_description

        return symmsystem

    def _add_conventional_bulk(self, cluster: Cluster, subsystem: System) -> None:
        '''
        Creates the subsystem with the symmetry information of the conventional cell
        '''
        cell = cluster.cell()
        symm = SymmetryAnalyzer(cell)
        conv_system = symm.get_conventional_system()
        subsystem.atoms = nomad_atoms_from_ase_atoms(conv_system)
        spg_number = symm.get_space_group_number()
        subsystem.cell = cell_from_ase_atoms(conv_system)
        symmetry = self._create_symmetry(symm)
        wyckoff_sets = symm.get_wyckoff_sets_conventional()
        material_id = material_id_bulk(spg_number, wyckoff_sets)
        subsystem.structural_type = 'bulk'
        subsystem.dimensionality = '3D'
        subsystem.material_id = material_id
        subsystem.symmetry = symmetry

    def _add_conventional_2d(self, cluster: Cluster, subsystem: System) -> None:
        '''
        Creates the subsystem with the symmetry information of the conventional cell.
        '''
        cell = cluster.cell()
        conv_atoms, _, wyckoff_sets, spg_number = structures_2d(cell)
        subsystem.cell = cell_from_ase_atoms(conv_atoms)
        subsystem.atoms = nomad_atoms_from_ase_atoms(conv_atoms)

        # Here we zero out the irrelevant lattice parameters to correctly handle
        # 2D systems with nonzero thickness (e.g. MoS2).
        subsystem.cell.c = None
        subsystem.cell.alpha = None
        subsystem.cell.beta = None
        subsystem.cell.atomic_density = None
        subsystem.cell.mass_density = None
        subsystem.cell.volume = None

        subsystem.structural_type = '2D'
        subsystem.dimensionality = '2D'
        subsystem.building_block = '2D material'
        subsystem.material_id = material_id_2d(spg_number, wyckoff_sets)

    def _create_symmetry(self, symm: SymmetryAnalyzer) -> Symmetry:
        international_short = symm.get_space_group_international_short()
        conv_system = symm.get_conventional_system()

        sec_symmetry = Symmetry()
        sec_symmetry.symmetry_method = 'MatID'
        sec_symmetry.space_group_number = symm.get_space_group_number()
        sec_symmetry.space_group_symbol = international_short
        sec_symmetry.hall_number = symm.get_hall_number()
        sec_symmetry.hall_symbol = symm.get_hall_symbol()
        sec_symmetry.point_group = symm.get_point_group()
        sec_symmetry.crystal_system = symm.get_crystal_system()
        sec_symmetry.bravais_lattice = symm.get_bravais_lattice()
        sec_symmetry.origin_shift = symm._get_spglib_origin_shift()
        sec_symmetry.transformation_matrix = symm._get_spglib_transformation_matrix()
        sec_symmetry.wyckoff_sets = wyckoff_sets_from_matid(symm.get_wyckoff_sets_conventional())

        spg_number = symm.get_space_group_number()
        atom_species = conv_system.get_atomic_numbers()
        if type(conv_system) == Atoms or conv_system.wyckoff_letters is None:
            wyckoffs = symm.get_wyckoff_letters_conventional()
        else:
            wyckoffs = conv_system.wyckoff_letters
        norm_wyckoff = atomutils.get_normalized_wyckoff(atom_species, wyckoffs)
        protoDict = atomutils.search_aflow_prototype(spg_number, norm_wyckoff)

        if protoDict is not None:
            sec_symmetry.prototype_label_aflow = protoDict.get('aflow_prototype_id')
            sec_symmetry.prototype_name = structure_name_map.get(protoDict.get('Notes'))

        return sec_symmetry
