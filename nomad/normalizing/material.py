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

from typing import Union, Dict, List, Tuple, Optional
from nptyping import NDArray
from collections import defaultdict
import re
import ase.data
from ase import Atoms
from ase.geometry.cell import complete_cell
import numpy as np
import matid.geometry
from matid.classification.structureclusterer import StructureClusterer
from matid import Classifier
from matid.classifications import Class0D, Atom, Class1D, Material2D, Surface, Class3D, Unknown

from nomad.datamodel.results import Symmetry, Material, System, Cell, Relation, StructureOriginal
from nomad import atomutils
from nomad.utils import hash


class MaterialNormalizer():
    def __init__(
            self,
            entry_archive,
            repr_system,
            repr_symmetry,
            spg_number,
            conv_atoms,
            wyckoff_sets,
            properties,
            optimade,
            logger):
        self.entry_archive = entry_archive
        self.repr_system = repr_system
        self.repr_method = None
        self.spg_number = spg_number
        self.conv_atoms = conv_atoms
        self.wyckoff_sets = wyckoff_sets
        self.properties = properties
        self.optimade = optimade
        self.repr_symmetry = repr_symmetry
        self.structural_type = None
        self.run = entry_archive.run[0]
        self.logger = logger

    def material(self) -> Material:
        '''Returns a populated Material subsection.'''
        material = Material()

        if self.repr_system:
            self.structural_type = self.repr_system.type
            material.structural_type = self.repr_system.type
            symbols, reduced_counts = atomutils.get_hill_decomposition(self.repr_system.atoms.labels, reduced=True)
            material.chemical_formula_reduced_fragments = [
                '{}{}'.format(n, int(c) if c != 1 else '') for n, c in zip(symbols, reduced_counts)
            ]
            material.chemical_formula_hill = self.repr_system.chemical_composition_hill
            material.chemical_formula_descriptive = material.chemical_formula_hill
        if self.optimade:
            material.elements = self.optimade.elements
            material.chemical_formula_reduced = self.optimade.chemical_formula_reduced
            material.chemical_formula_anonymous = self.optimade.chemical_formula_anonymous

        material.symmetry = self.symmetry()

        if self.structural_type == 'bulk':
            material.material_id = self.material_id_bulk(self.spg_number, self.wyckoff_sets)
            material.material_name = self.material_name(symbols, reduced_counts)
            classes = self.material_classification()
            if classes:
                material.functional_type = classes.get('material_class_springer')
                material.compound_type = classes.get('compound_class_springer')
        if self.structural_type == '2D':
            material.material_id = self.material_id_2d(self.spg_number, self.wyckoff_sets)
        elif self.structural_type == '1D':
            material.material_id = self.material_id_1d(self.conv_atoms)

        material.topology = self.topology(material)

        return material

    def material_id_bulk(self, spg_number: int, wyckoff_sets) -> str:
        if spg_number is None or wyckoff_sets is None:
            return None
        norm_hash_string = atomutils.get_symmetry_string(spg_number, wyckoff_sets)
        return hash(norm_hash_string)

    def material_id_2d(self, spg_number: int, wyckoff_sets) -> str:
        if spg_number is None or wyckoff_sets is None:
            return None
        norm_hash_string = atomutils.get_symmetry_string(spg_number, wyckoff_sets, is_2d=True)
        return hash(norm_hash_string)

    def material_id_1d(self, conv_atoms: Atoms) -> str:
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

    def material_classification(self) -> Dict[str, List[str]]:
        try:
            sec_springer = self.repr_system['springer_material'][0]
        except Exception:
            return None

        classes: Dict[str, List[str]] = {}
        try:
            classifications = sec_springer['classification']
        except KeyError:
            pass
        else:
            classes['material_class_springer'] = classifications
        try:
            compound_classes = sec_springer['compound_class']
        except KeyError:
            pass
        else:
            classes['compound_class_springer'] = compound_classes
        return classes

    def material_name(
            self, symbols: Union[List, NDArray], counts: Union[List, NDArray]) -> str:
        if symbols is None or counts is None:
            return None
        name = None
        # Systems with one element are named after it
        if len(symbols) == 1:
            number = ase.data.atomic_numbers[symbols[0]]
            name = ase.data.atomic_names[number]
        # Binary systems have specific names
        elif len(symbols) == 2:
            atomicnumbers = [ase.data.atomic_numbers[i] for i in symbols]
            names = [ase.data.atomic_names[i] for i in atomicnumbers]

            # Non-metal elements are anions in the binary compounds and receive the -ide suffix
            if names[1] == 'Antimony':
                names[1] = names[1][:-1] + 'ide'
            if names[1] == 'Arsenic':
                names[1] = names[1][:-1] + 'de'
            if names[1] == 'Boron' or names[1] == 'Carbon':
                names[1] = names[1][:-2] + 'ide'
            if names[1] == 'Chlorine' or names[1] == 'Germanium' or names[1] == 'Selenium' or names[1] == 'Bromine' \
               or names[1] == 'Tellurium' or names[1] == 'Iodine' or names[1] == 'Polonium' or names[1] == 'Astatine' or \
               names[1] == 'Fluorine':
                names[1] = names[1][:-2] + 'de'
            if names[1] == 'Silicon' or names[1] == 'Sulfur':
                names[1] = names[1][:-2] + 'ide'
            if names[1] == 'Nitrogen' or names[1] == 'Oxygen' or names[1] == 'Hydrogen' or names[1] == 'Phosphorus':
                names[1] = names[1][:-4] + 'ide'

            name = names[0] + ' ' + names[1]

            if names[1] == 'Fluoride' or names[1] == 'Chloride' or names[1] == 'Bromide' or \
               names[1] == 'Iodide' or names[1] == 'Hydride':

                # Non-metals with elements of variable valence, therefore we remove alkaline and
                # alkaline-earth elements, which have fixed valence
                # Only the most electronegative non-metals are supposed to make ionic compounds
                if names[0] != 'Lithium' and names[0] != 'Sodium' and names[0] != 'Potassium' and \
                   names[0] != 'Rubidium' and names[0] != 'Cesium' and names[0] != 'Francium' and \
                   names[0] != 'Beryllium' and names[0] != 'Magnesium' and names[0] != 'Calcium' and \
                   names[0] != 'Strontium' and names[0] != 'Barium' and names[0] != 'Radium' and \
                   names[0] != 'Aluminum':

                    if counts[1] == 2:
                        name = names[0] + '(II)' + ' ' + names[1]
                    elif counts[1] == 3:
                        name = names[0] + '(III)' + ' ' + names[1]
                    elif counts[1] == 4:
                        name = names[0] + '(IV)' + ' ' + names[1]
                    elif counts[1] == 5:
                        name = names[0] + '(V)' + ' ' + names[1]
                    elif counts[1] == 6:
                        name = names[0] + '(VI)' + ' ' + names[1]
                    elif counts[1] == 7:
                        name = names[0] + '(VII)' + ' ' + names[1]

            if names[1] == 'Oxide' or names[1] == 'Sulfide' or names[1] == 'Selenide':
                if names[0] != 'Lithium' and names[0] != 'Sodium' and names[0] != 'Potassium' and \
                   names[0] != 'Rubidium' and names[0] != 'Cesium' and names[0] != 'Francium' and \
                   names[0] != 'Beryllium' and names[0] != 'Magnesium' and names[0] != 'Calcium' and \
                   names[0] != 'Strontium' and names[0] != 'Barium' and names[0] != 'Radium' and \
                   names[0] != 'Aluminum':

                    if counts[0] == 1 and counts[1] == 1:
                        name = names[0] + '(II)' + ' ' + names[1]
                    elif counts[0] == 2 and counts[1] == 1:
                        name = names[0] + '(I)' + ' ' + names[1]
                    elif counts[0] == 1 and counts[1] == 2:
                        name = names[0] + '(IV)' + ' ' + names[1]
                    elif counts[0] == 2 and counts[1] == 3:
                        name = names[0] + '(III)' + ' ' + names[1]
                    elif counts[0] == 2 and counts[1] == 5:
                        name = names[0] + '(V)' + ' ' + names[1]
                    elif counts[0] == 1 and counts[1] == 3:
                        name = names[0] + '(VI)' + ' ' + names[1]
                    elif counts[0] == 2 and counts[1] == 7:
                        name = names[0] + '(VII)' + ' ' + names[1]

            if names[1] == 'Nitride' or names[1] == 'Phosphide':
                if names[0] != 'Lithium' and names[0] != 'Sodium' and names[0] != 'Potassium' and \
                   names[0] != 'Rubidium' and names[0] != 'Cesium' and names[0] != 'Francium' and \
                   names[0] != 'Beryllium' and names[0] != 'Magnesium' and names[0] != 'Calcium' and \
                   names[0] != 'Strontium' and names[0] != 'Barium' and names[0] != 'Radium' and \
                   names[0] != 'Aluminum':

                    if counts[0] == 1 and counts[1] == 1:
                        name = names[0] + '(III)' + ' ' + names[1]
                    if counts[0] == 1 and counts[1] == 2:
                        name = names[0] + '(VI)' + ' ' + names[1]
                    elif counts[0] == 3 and counts[1] == 2:
                        name = names[0] + '(II)' + ' ' + names[1]
                    elif counts[0] == 3 and counts[1] == 4:
                        name = names[0] + '(IV)' + ' ' + names[1]
                    elif counts[0] == 3 and counts[1] == 5:
                        name = names[0] + '(V)' + ' ' + names[1]
                    elif counts[0] == 3 and counts[1] == 7:
                        name = names[0] + '(VII)' + ' ' + names[1]

            if names[1] == 'Carbide':
                if names[0] != 'Lithium' and names[0] != 'Sodium' and names[0] != 'Potassium' and \
                   names[0] != 'Rubidium' and names[0] != 'Cesium' and names[0] != 'Francium' and \
                   names[0] != 'Beryllium' and names[0] != 'Magnesium' and names[0] != 'Calcium' and \
                   names[0] != 'Strontium' and names[0] != 'Barium' and names[0] != 'Radium' and \
                   names[0] != 'Aluminum':

                    if counts[0] == 1 and counts[1] == 1:
                        name = names[0] + '(IV)' + ' ' + names[1]
                    if counts[0] == 2 and counts[1] == 1:
                        name = names[0] + '(II)' + ' ' + names[1]
                    if counts[0] == 4 and counts[1] == 1:
                        name = names[0] + '(I)' + ' ' + names[1]
                    if counts[0] == 4 and counts[1] == 3:
                        name = names[0] + '(III)' + ' ' + names[1]
                    if counts[0] == 4 and counts[1] == 5:
                        name = names[0] + '(V)' + ' ' + names[1]
                    if counts[0] == 2 and counts[1] == 3:
                        name = names[0] + '(VI)' + ' ' + names[1]
                    if counts[0] == 4 and counts[1] == 7:
                        name = names[0] + '(VII)' + ' ' + names[1]

        return name

    def symmetry(self) -> Symmetry:
        '''Returns a populated Symmetry subsection.'''
        result = Symmetry()
        filled = False

        if self.repr_symmetry:
            result.hall_number = self.repr_symmetry.hall_number
            result.hall_symbol = self.repr_symmetry.hall_symbol
            result.bravais_lattice = self.repr_symmetry.bravais_lattice
            result.crystal_system = self.repr_symmetry.crystal_system
            result.space_group_number = self.repr_symmetry.space_group_number
            result.space_group_symbol = self.repr_symmetry.international_short_symbol
            result.point_group = self.repr_symmetry.point_group
            filled = True

        # Fill in prototype information. SystemNormalizer has cached many of
        # the values during it's own analysis. These cached values are used
        # here.
        proto = self.repr_system.prototype if self.repr_system else None
        proto = proto[0] if proto else None
        if proto:
            # Prototype id and formula
            result.prototype_aflow_id = proto.aflow_id
            result.prototype_formula = proto.m_cache.get('prototype_name')

            # Strukturbericht: replace LaTeX with plain text
            strukturbericht = proto.m_cache.get('strukturbericht_designation')
            if strukturbericht:
                strukturbericht = re.sub('[$_{}]', '', strukturbericht)
                result.strukturbericht_designation = strukturbericht

            # Structure name. Only relevant information hidden in 'notes' is
            # handed over TODO: review and eventually add more -ites which
            # are commonly used (see wurzite)
            note = proto.m_cache.get('prototype_notes')
            if note:
                structure_name_map = {
                    'CaTiO<sub>3</sub> Pnma Perovskite Structure': 'perovskite',
                    'Hypothetical Tetrahedrally Bonded Carbon with 4&ndash;Member Rings': '4-member ring',
                    'In (A6) Structure': 'fct',
                    '$\\alpha$&ndash;Pa (A<sub>a</sub>) Structure': 'bct',
                    'Hypothetical BCT5 Si Structure': 'bct5',
                    'Wurtzite (ZnS, B4) Structure': 'wurtzite',
                    'Hexagonal Close Packed (Mg, A3) Structure': 'hcp',
                    'Half&ndash;Heusler (C1<sub>b</sub>) Structure': 'half-Heusler',
                    'Zincblende (ZnS, B3) Structure': 'zincblende',
                    'Cubic Perovskite (CaTiO<sub>3</sub>, E2<sub>1</sub>) Structure': 'cubic perovskite',
                    '$\\alpha$&ndash;Po (A<sub>h</sub>) Structure': 'simple cubic',
                    'Si<sub>46</sub> Clathrate Structure': 'clathrate',
                    'Cuprite (Cu<sub>2</sub>O, C3) Structure': 'cuprite',
                    'Heusler (L2<sub>1</sub>) Structure': 'Heusler',
                    'Rock Salt (NaCl, B1) Structure': 'rock salt',
                    'Face&ndash;Centered Cubic (Cu, A1) Structure': 'fcc',
                    'Diamond (A4) Structure': 'diamond',
                    'Body&ndash;Centered Cubic (W, A2) Structure': 'bcc',
                }
                result.structure_name = structure_name_map.get(note, None)

            filled = True

        if filled:
            return result
        return None

    def topology(self, material) -> List[System]:
        '''Extracts the system topology if one is available
        '''
        # Use the calculation topology primarily
        topology_calc = self.topology_calculation(material)
        if topology_calc:
            return topology_calc
        topology_matid = self.topology_matid(material)
        if topology_matid:
            return topology_matid

        return []

    def topology_calculation(self, material: Material) -> List[System]:
        '''Extracts the system topology as defined in the original calculation.
        This topology typically comes from e.g. classical force fields that
        require a specific topology for the system.
        '''
        try:
            groups = self.entry_archive.run[0].system[0].atoms_group
            if len(groups) == 0:
                return None
        except Exception:
            return None

        topology, original, structure_original = self._init_orig_topology(material)
        if structure_original is None:
            return []

        top_id = 0
        top_id_original = f'/results/material/topology/{top_id}'
        top_id += 1
        if groups:
            label_to_instances: Dict[str, List] = defaultdict(list)
            label_to_id: Dict[str, str] = {}

            def add_group(groups, parent=None, parent_id=None):
                nonlocal top_id
                if groups:
                    for group in groups:
                        label = group.label
                        if label not in label_to_instances:
                            try:
                                formula = atomutils.Formula(group.composition_formula)
                            except Exception:
                                formula_hill = None
                                formula_anonymous = None
                                formula_reduced = None
                                elements = None
                            else:
                                formula_hill = formula.format('hill')
                                formula_anonymous = formula.format('anonymous')
                                formula_reduced = formula.format('reduce')
                                elements = formula.elements()
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
                            top_id_str = f'/results/material/topology/{top_id}'
                            description = description_map.get(group.type, None)
                            structural_type = structural_type_map.get(group.type, None)

                            system = System(
                                system_id=top_id_str,
                                method='parser',
                                description=description,
                                label=group.label,
                                formula_hill=formula_hill,
                                formula_anonymous=formula_anonymous,
                                formula_reduced=formula_reduced,
                                elements=elements,
                                structural_type=structural_type,
                                n_atoms=group.n_atoms,
                                system_relation=Relation(type='subsystem'),
                                parent_system=None if not parent_id else parent_id
                            )

                            topology[top_id_str] = system
                            label_to_id[label] = top_id_str
                            if parent:
                                parent_children = parent.child_systems if parent.child_systems else []
                                parent_children.append(top_id_str)
                                parent.child_systems = parent_children
                            top_id += 1
                            add_group(group.atoms_group, system, top_id_str)

                        # Add the instance indices to the topology. Only added
                        # if the length matches, otherwise log an error
                        # (instances with the same label should have the same
                        # atoms)
                        old_instances = label_to_instances[label]
                        save = False
                        if len(old_instances) == 0:
                            save = True
                        else:
                            if len(old_instances[0]) == len(group.atom_indices):
                                save = True
                            else:
                                self.logger.warn((
                                    'the topology contains entries with the same '
                                    'label but with different number of atoms'
                                ))
                        if save:
                            label_to_instances[label].append(group.atom_indices)
            add_group(groups, original, top_id_original)

            # Add the gathered instance information and formula information if
            # not yet present
            for label, value in label_to_instances.items():
                top_id_str = label_to_id[label]
                top = topology[top_id_str]
                indices = np.array(value)
                top.indices = indices
                if structure_original and top.formula_hill is None:
                    symbols = ''.join(np.array(structure_original.species_at_sites)[indices.flatten()])
                    try:
                        formula = atomutils.Formula(symbols)
                    except Exception:
                        pass
                    else:
                        top.formula_hill = formula.format('hill')
                        top.formula_anonymous = formula.format('anonymous')
                        top.formula_reduced = formula.format('reduce')
                        top.elements = formula.elements()

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
        topology, original, structure_original = self._init_orig_topology(material)
        if structure_original is None:
            return None
        if structure_original.nperiodic_dimensions == 0:
            return None
        # TODO: It still needs to be decided what this limit should be.
        n_atoms = len(structure_original.species_at_sites)
        if n_atoms > 500:
            return None
        top_id = 1
        top_id_str = f'/results/material/topology/{top_id}'

        # Perform clustering
        clusters = self._perform_matid_clustering(structure_original)

        # Add all meaningful clusters to the topology
        cluster_indices_list = self._filter_clusters(clusters)
        for indices in cluster_indices_list:
            system_type = self._system_type_analysis(structure_original, indices)
            if system_type not in {"2D", "surface"}:
                continue
            subsystem = self._create_subsystem(structure_original, indices, system_type, top_id)

            topology[top_id_str] = subsystem
            parent_children = original.child_systems if original.child_systems else []
            parent_children.append(top_id_str)
            original.child_systems = parent_children
            top_id += 1
            top_id_str = f'/results/material/topology/{top_id}'

            # TODO: Analyze the connected cells
            # (matid.core.linkedunits.LinkedUnitCollection) in the subsystem to
            # figure out the structure type, conventional cell, orientation,
            # etc. This data should be stored as children of the subsystem
            # within the topology.

        # If the returned clusters contain more than one 2D/surface component,
        # the topology is accepted and returned. TODO: This should be extended
        # to also output more complex topologies.
        if len(topology) < 3:
            return None

        return list(topology.values())

    def _init_orig_topology(self, material: Material) -> Tuple[Dict[str, System], System, StructureOriginal]:
        topology: Dict[str, System] = {}

        structure_original = self._check_original_structure()
        if structure_original is None:
            return None, None, None

        top_id = 0
        original = self._add_orig_topology_as_root(material, top_id)

        if structure_original:
            original.cell = self._create_cell(structure_original)
            original.atoms = structure_original
            original.n_atoms = structure_original.n_sites
        topology[str(top_id)] = original

        return topology, original, structure_original

    def _add_orig_topology_as_root(self, material: Material, top_id: int) -> System:
        '''Adds the original system topology as root
        '''
        top_id_str = f'/results/material/topology/{top_id}'
        original = System(
            system_id=top_id_str,
            method='parser',
            label='original',
            description='A representative system chosen from the original simulation.',
            material_id=material.material_id,
            material_name=material.material_name,
            structural_type=material.structural_type,
            functional_type=material.functional_type,
            compound_type=material.compound_type,
            formula_hill=material.chemical_formula_hill,
            formula_anonymous=material.chemical_formula_anonymous,
            formula_reduced=material.chemical_formula_reduced,
            elements=material.elements,
        )

        return original

    def _check_original_structure(self) -> StructureOriginal:
        '''
        Checks if original system is available and if system size is processable. The
        topology is created only if structural_type == unavailable and a meaningful
        topology can be extracted.
        '''
        try:
            structure_original = self.properties.structures.structure_original
        except Exception:
            return None

        return structure_original

    def _perform_matid_clustering(self, structure_original: StructureOriginal) -> list:
        '''
        Creates an ase.atoms and performs the clustering with MatID
        '''
        system = Atoms(
            symbols=structure_original.species_at_sites,
            positions=structure_original.cartesian_site_positions * 1e10,
            cell=complete_cell(structure_original.lattice_vectors * 1e10),
            pbc=np.array(structure_original.dimension_types, dtype=bool)
        )

        # Perform the clustering
        clusterer = StructureClusterer()
        clusters = clusterer.get_clusters(
            system,
            max_cell_size=6,
            pos_tol=0.70,
            angle_tol=20,
            merge_threshold=0.3,
            merge_radius=5
        )
        return clusters

    def _create_cell(self, structure_original: StructureOriginal) -> Cell:
        '''
        Creates a Cell from the given structure.
        '''
        cell = Cell(
            a=structure_original.lattice_parameters.a,
            b=structure_original.lattice_parameters.b,
            c=structure_original.lattice_parameters.c,
            alpha=structure_original.lattice_parameters.alpha,
            beta=structure_original.lattice_parameters.beta,
            gamma=structure_original.lattice_parameters.gamma,
            volume=structure_original.cell_volume,
            atomic_density=structure_original.atomic_density,
            mass_density=structure_original.mass_density,
        )
        return cell

    def _filter_clusters(self, clusters: StructureClusterer) -> List[List[int]]:
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
        for cluster in filtered_cluster:
            indices = list(cluster.indices)
            cluster_indices_list += [indices]

        return cluster_indices_list

    def _system_type_analysis(self, structure_original: StructureOriginal, indices: List[int]) -> matid.classifications:
        '''
        Classifies ase.atoms and returns the MatID system type as a string.
        '''
        system_types = {Class3D: 'bulk',
                        Atom: 'atom',
                        Class0D: 'molecule / cluster',
                        Class1D: '1D',
                        Surface: 'surface',
                        Material2D: '2D',
                        Unknown: 'unavailable'}
        # Create the system as ASE.Atoms
        cluster_atoms = Atoms(
            symbols=np.array(structure_original.species_at_sites)[indices],
            positions=np.array(structure_original.cartesian_site_positions)[indices] * 1e10,
            cell=complete_cell(structure_original.lattice_vectors * 1e10),
            pbc=np.array(structure_original.dimension_types, dtype=bool)
        )
        try:
            classifier = Classifier(radii='covalent', cluster_threshold=1)
            cls = classifier.classify(cluster_atoms)
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

    def _create_subsystem(self, structure_original: StructureOriginal, indices: List[int], system_type: str, top_id: int) -> System:
        '''
        Creates the subsystem with system type
        '''
        subspecies = np.array(structure_original.species_at_sites)[indices]
        formula = atomutils.Formula(''.join(subspecies))
        formula_hill = formula.format('hill')
        formula_anonymous = formula.format('anonymous')
        formula_reduced = formula.format('reduce')
        elements = formula.elements()
        top_id_str = f'/results/material/topology/{top_id}'

        subsystem = System(
            system_id=top_id_str,
            method='matid',
            label='subsystem',
            description='Automatically detected subsystem.',
            structural_type=system_type,
            indices=[indices],
            system_relation=Relation(type='subsystem'),
            formula_hill=formula_hill,
            formula_anonymous=formula_anonymous,
            formula_reduced=formula_reduced,
            elements=elements,
            parent_system='/results/material/topology/0',
            n_atoms=len(indices)
        )
        return subsystem
