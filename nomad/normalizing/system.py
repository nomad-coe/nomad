# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the 'License');
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an'AS IS' BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from typing import Any, Dict
from collections import Counter
import ase
from ase import Atoms
import numpy as np
import json
import re
from matid import SymmetryAnalyzer, Classifier
from matid.classifications import Class0D, Atom, Class1D, Material2D, Surface, Class3D

from nomad import atomutils, archive
from nomad import utils, config
from nomad.datamodel.metainfo.public import section_symmetry, section_std_system, \
    section_primitive_system, section_original_system, section_springer_material, \
    section_prototype, section_system

from .normalizer import SystemBasedNormalizer

# use a regular expression to check atom labels; expression is build from list of
# all labels sorted desc to find Br and not B when searching for Br.
atom_label_re = re.compile('|'.join(
    sorted(ase.data.chemical_symbols, key=lambda x: len(x), reverse=True)))


def normalized_atom_labels(atom_labels):
    '''
    Normalizes the given atom labels: they either are labels right away, or contain
    additional numbers (to distinguish same species but different labels, see meta-info),
    or we replace them with ase placeholder atom for unknown elements 'X'.
    '''
    return [
        ase.data.chemical_symbols[0] if match is None else match.group(0)
        for match in [re.search(atom_label_re, atom_label) for atom_label in atom_labels]]


def formula_normalizer(atoms):
    '''
    Reads the chemical symbols in ase.atoms and returns a normalized formula.
    Formula normalization is on the basis of atom counting,
    e.g., Tc ->  Tc100, SZn -> S50Zn50, Co2Nb -> Co67Nb33
    '''
    #
    chem_symb = atoms.get_chemical_symbols()
    atoms_counter = Counter(chem_symb)  # dictionary
    atoms_total = sum(atoms_counter.values())

    atoms_normed = []
    for key in atoms_counter.keys():
        norm = str(round(100 * atoms_counter[key] / atoms_total))
        atoms_normed.append(key + norm)
    #
    atoms_normed.sort()
    return ''.join(atoms_normed)


class SystemNormalizer(SystemBasedNormalizer):
    '''
    This normalizer performs all system (atoms, cells, etc.) related normalizations
    of the legacy NOMAD-coe *stats* normalizer.
    '''
    @staticmethod
    def atom_label_to_num(atom_label):
        # Take first three characters and make first letter capitalized.
        atom_label = atom_label[:3].title()

        for symbol_length in reversed(range(1, 4)):
            symbol = atom_label[:symbol_length]
            if symbol in ase.data.chemical_symbols:
                return ase.data.chemical_symbols.index(symbol)

        return 0

    def normalize_system(self, system, is_representative) -> bool:
        '''
        The 'main' method of this :class:`SystemBasedNormalizer`.
        Normalizes the section with the given `index`.
        Normalizes geometry, classifies, system_type, and runs symmetry analysis.

        Returns: True, iff the normalization was successful
        '''
        if self.section_run is None:
            self.logger.error('section_run is not present.')
            return False

        def get_value(quantity_def, default: Any = None, numpy: bool = True) -> Any:
            try:
                value = system.m_get(quantity_def)
                if not numpy and type(value).__module__ == np.__name__:
                    value = value.tolist()

                elif numpy and isinstance(value, list):
                    value = np.array(value)

                return value
            except (KeyError, IndexError):
                return default

        if is_representative:
            system.is_representative = is_representative

        # analyze atoms labels
        atom_labels = get_value(section_system.atom_labels, numpy=False)
        if atom_labels is not None:
            atom_labels = normalized_atom_labels(atom_labels)

        atom_species = get_value(section_system.atom_species, numpy=False)
        if atom_labels is None and atom_species is None:
            self.logger.warn('system has neither atom species nor labels')
            return False

        # If there are no atom labels we create them from atom species data.
        if atom_labels is None:
            try:
                atom_labels = list(ase.data.chemical_symbols[species] for species in atom_species)
            except IndexError:
                self.logger.error('system has atom species that are out of range')
                return False

            system.atom_labels = atom_labels

        # At this point we should have atom labels.
        try:
            atoms = ase.Atoms(symbols=atom_labels)
            chemical_symbols = list(atoms.get_chemical_symbols())
            if atom_labels != chemical_symbols:
                self.logger.error('atom labels are ambiguous', atom_labels=atom_labels[:10])
            atom_labels = chemical_symbols
        except Exception as e:
            self.logger.error(
                'cannot build ase atoms from atom labels',
                atom_labels=atom_labels[:10], exc_info=e, error=str(e))
            raise e

        if atom_species is None:
            atom_species = atoms.get_atomic_numbers().tolist()
            system.atom_species = atom_species
        else:
            if not isinstance(atom_species, list):
                atom_species = [atom_species]
            if atom_species != atoms.get_atomic_numbers().tolist():
                self.logger.warning(
                    'atom species do not match labels',
                    atom_labels=atom_labels[:10], atom_species=atom_species[:10])
                atom_species = atoms.get_atomic_numbers().tolist()
            system.atom_species = atom_species

        # periodic boundary conditions
        pbc = get_value(section_system.configuration_periodic_dimensions, numpy=False)
        if pbc is None:
            pbc = [False, False, False]
            self.logger.warning('missing configuration_periodic_dimensions')
            system.configuration_periodic_dimensions = pbc
        try:
            atoms.set_pbc(pbc)
        except Exception as e:
            self.logger.error(
                'cannot use pbc with ase atoms', exc_info=e, pbc=pbc, error=str(e))
            return False

        # formulas
        system.chemical_composition = atoms.get_chemical_formula(mode='all')
        system.chemical_composition_reduced = atoms.get_chemical_formula(mode='reduce')
        system.chemical_composition_bulk_reduced = atoms.get_chemical_formula(mode='hill')

        # positions
        atom_positions = get_value(section_system.atom_positions, numpy=True)
        if atom_positions is None:
            self.logger.warning('no atom positions, skip further system analysis')
            return False
        if len(atom_positions) != len(atoms):
            self.logger.error(
                'len of atom position does not match number of atoms',
                n_atom_positions=len(atom_positions), n_atoms=len(atoms))
            return False
        try:
            atoms.set_positions(1e10 * atom_positions.magnitude)
        except Exception as e:
            self.logger.error(
                'cannot use positions with ase atoms', exc_info=e, error=str(e))
            return False

        # lattice vectors
        lattice_vectors = get_value(section_system.lattice_vectors, numpy=True)
        if lattice_vectors is None:
            lattice_vectors = get_value(section_system.simulation_cell, numpy=True)
            if lattice_vectors is not None:
                system.lattice_vectors = lattice_vectors
        if lattice_vectors is None:
            if any(pbc):
                self.logger.error('no lattice vectors but periodicity', pbc=pbc)
        else:
            try:
                atoms.set_cell(1e10 * lattice_vectors.magnitude)
            except Exception as e:
                self.logger.error(
                    'cannot use lattice_vectors with ase atoms', exc_info=e, error=str(e))
                return False

        # configuration
        configuration = [
            atom_labels, atoms.positions.tolist(),
            atoms.cell.tolist() if atoms.cell is not None else None,
            atoms.pbc.tolist()]
        configuration_id = utils.hash(json.dumps(configuration).encode('utf-8'))
        system.configuration_raw_gid = configuration_id

        if is_representative:
            # Save the Atoms as a temporary variable
            system.m_cache["representative_atoms"] = atoms

            # System type analysis
            if atom_positions is not None:
                with utils.timer(
                        self.logger, 'system classification executed',
                        system_size=len(atoms)):
                    self.system_type_analysis(atoms)

            system_type = system.system_type
            # Symmetry analysis
            if atom_positions is not None and (lattice_vectors is not None or not any(pbc)) and system_type == "bulk":
                with utils.timer(
                        self.logger, 'symmetry analysis executed',
                        system_size=len(atoms)):
                    self.symmetry_analysis(system, atoms)

        return True

    def system_type_analysis(self, atoms: Atoms) -> None:
        '''
        Determine the system type with MatID. Write the system type to the
        entry_archive.

        Args:
            atoms: The structure to analyse
        '''
        system_type = config.services.unavailable_value
        if len(atoms) <= config.normalize.system_classification_with_clusters_threshold:
            try:
                classifier = Classifier(radii="covalent", cluster_threshold=config.normalize.cluster_threshold)
                cls = classifier.classify(atoms)
            except Exception as e:
                self.logger.error(
                    'matid project system classification failed', exc_info=e, error=str(e))
            else:
                classification = type(cls)
                if classification == Class3D:
                    system_type = 'bulk'
                elif classification == Atom:
                    system_type = 'atom'
                elif classification == Class0D:
                    system_type = 'molecule / cluster'
                elif classification == Class1D:
                    system_type = '1D'
                elif classification == Surface:
                    system_type = 'surface'
                elif classification == Material2D:
                    system_type = '2D'
        else:
            self.logger.info("system type analysis not run due to large system size")
        idx = self.section_run.m_cache["representative_system_idx"]
        self.section_run.section_system[idx].system_type = system_type
        self.section_run.section_system[-1].system_type = system_type

    def symmetry_analysis(self, system, atoms: ase.Atoms) -> None:
        '''Analyze the symmetry of the material being simulated. Only performed
        for bulk materials.

        We feed in the parsed values in section_system to the the symmetry
        analyzer. The analysis results are written to the entry_archive.

        Args:
            atoms: The atomistic structure to analyze.

        Returns:
            None: The method should write symmetry variables
            to the entry_archive which is member of this class.
        '''
        # Try to use MatID's symmetry analyzer to analyze the ASE object.
        try:
            symm = SymmetryAnalyzer(atoms, symmetry_tol=config.normalize.symmetry_tolerance)

            space_group_number = symm.get_space_group_number()

            hall_number = symm.get_hall_number()
            hall_symbol = symm.get_hall_symbol()

            crystal_system = symm.get_crystal_system()
            bravais_lattice = symm.get_bravais_lattice()
            point_group = symm.get_point_group()

            orig_wyckoff = symm.get_wyckoff_letters_original()
            prim_wyckoff = symm.get_wyckoff_letters_primitive()
            conv_wyckoff = symm.get_wyckoff_letters_conventional()

            orig_equivalent_atoms = symm.get_equivalent_atoms_original()
            prim_equivalent_atoms = symm.get_equivalent_atoms_primitive()
            conv_equivalent_atoms = symm.get_equivalent_atoms_conventional()
            international_short = symm.get_space_group_international_short()

            conv_sys = symm.get_conventional_system()
            conv_pos = conv_sys.get_scaled_positions()
            conv_cell = conv_sys.get_cell()
            conv_num = conv_sys.get_atomic_numbers()

            prim_sys = symm.get_primitive_system()
            prim_pos = prim_sys.get_scaled_positions()
            prim_cell = prim_sys.get_cell()
            prim_num = prim_sys.get_atomic_numbers()

            transform = symm._get_spglib_transformation_matrix()
            origin_shift = symm._get_spglib_origin_shift()
        except ValueError as e:
            self.logger.debug('symmetry analysis is not available', details=str(e))
            return
        except Exception as e:
            self.logger.error('matid symmetry analysis fails with exception', exc_info=e)
            return

        # Write data extracted from MatID's symmetry analysis to the
        # representative section_system.

        sec_symmetry = system.m_create(section_symmetry)
        sec_symmetry.m_cache["symmetry_analyzer"] = symm

        sec_symmetry.symmetry_method = 'MatID (spg)'
        sec_symmetry.space_group_number = space_group_number
        sec_symmetry.hall_number = hall_number
        sec_symmetry.hall_symbol = hall_symbol
        sec_symmetry.international_short_symbol = international_short
        sec_symmetry.point_group = point_group
        sec_symmetry.crystal_system = crystal_system
        sec_symmetry.bravais_lattice = bravais_lattice
        sec_symmetry.origin_shift = origin_shift
        sec_symmetry.transformation_matrix = transform

        sec_std = sec_symmetry.m_create(section_std_system)
        sec_std.lattice_vectors_std = conv_cell
        sec_std.atom_positions_std = conv_pos
        sec_std.atomic_numbers_std = conv_num
        sec_std.wyckoff_letters_std = conv_wyckoff
        sec_std.equivalent_atoms_std = conv_equivalent_atoms

        sec_prim = sec_symmetry.m_create(section_primitive_system)
        sec_prim.lattice_vectors_primitive = prim_cell
        sec_prim.atom_positions_primitive = prim_pos
        sec_prim.atomic_numbers_primitive = prim_num
        sec_prim.wyckoff_letters_primitive = prim_wyckoff
        sec_prim.equivalent_atoms_primitive = prim_equivalent_atoms

        sec_orig = sec_symmetry.m_create(section_original_system)
        sec_orig.wyckoff_letters_original = orig_wyckoff
        sec_orig.equivalent_atoms_original = orig_equivalent_atoms

        self.springer_classification(atoms, space_group_number)  # Springer Normalizer
        self.prototypes(system, conv_num, conv_wyckoff, space_group_number)

    def springer_classification(self, atoms, space_group_number):
        normalized_formula = formula_normalizer(atoms)
        springer_data = query_springer_data(normalized_formula, space_group_number)
        idx = self.section_run.m_cache["representative_system_idx"]

        for material in springer_data.values():
            sec_springer_mat = self.section_run.section_system[idx].m_create(section_springer_material)

            sec_springer_mat.springer_id = material['spr_id']
            sec_springer_mat.springer_alphabetical_formula = material['spr_aformula']
            sec_springer_mat.springer_url = material['spr_url']

            compound_classes = material['spr_compound']
            if compound_classes is None:
                compound_classes = []
            sec_springer_mat.springer_compound_class = compound_classes

            classifications = material['spr_classification']
            if classifications is None:
                classifications = []
            sec_springer_mat.springer_classification = classifications

        # Check the 'springer_classification' and 'springer_compound_class' information
        # found is the same for all springer_id's
        springer_data_keys = list(springer_data.keys())
        if len(springer_data_keys) != 0:
            class_0 = springer_data[springer_data_keys[0]]['spr_classification']
            comp_0 = springer_data[springer_data_keys[0]]['spr_compound']

            # compare 'class_0' and 'comp_0' against the rest
            for ii in range(1, len(springer_data_keys)):
                class_test = (class_0 == springer_data[springer_data_keys[ii]]['spr_classification'])
                comp_test = (comp_0 == springer_data[springer_data_keys[ii]]['spr_compound'])

                if (class_test or comp_test) is False:
                    self.logger.info('Mismatch in Springer classification or compounds')

    def prototypes(self, system, atom_species: np.array, wyckoffs: np.array, spg_number: int) -> None:
        '''Tries to match the material to an entry in the AFLOW prototype data.
        If a match is found, a section_prototype is added to section_system.

        Args:
            atomic_numbers: Array of atomic numbers.
            wyckoff_letters: Array of Wyckoff letters as strings.
            spg_number: Space group number.
        '''
        norm_wyckoff = atomutils.get_normalized_wyckoff(atom_species, wyckoffs)
        protoDict = atomutils.search_aflow_prototype(spg_number, norm_wyckoff)
        if protoDict is not None:
            aflow_prototype_id = protoDict["aflow_prototype_id"]
            aflow_prototype_url = protoDict["aflow_prototype_url"]
            aflow_prototype_notes = protoDict["Notes"]
            aflow_prototype_name = protoDict["Prototype"]
            aflow_strukturbericht_designation = protoDict["Strukturbericht Designation"]
            prototype_label = '%d-%s-%s' % (
                spg_number,
                aflow_prototype_name,
                protoDict.get("Pearsons Symbol", "-")
            )
            idx = self.section_run.m_cache["representative_system_idx"]
            sec_prototype = self.section_run.section_system[idx].m_create(section_prototype)
            sec_prototype.prototype_label = prototype_label
            sec_prototype.prototype_aflow_id = aflow_prototype_id
            sec_prototype.prototype_aflow_url = aflow_prototype_url
            sec_prototype.prototype_assignment_method = "normalized-wyckoff"
            sec_prototype.m_cache["prototype_notes"] = aflow_prototype_notes
            sec_prototype.m_cache["prototype_name"] = aflow_prototype_name
            if aflow_strukturbericht_designation != "None":
                sec_prototype.m_cache["strukturbericht_designation"] = aflow_strukturbericht_designation


def query_springer_data(normalized_formula: str, space_group_number: int) -> Dict[str, Any]:
    ''' Queries a msgpack database for springer-related quantities. '''
    if config.normalize.springer_db_path is None:
        return {}

    entries = archive.query_archive(config.normalize.springer_db_path, {str(space_group_number): {normalized_formula: '*'}})
    db_dict = {}
    entries = entries.get(str(space_group_number), {}).get(normalized_formula, {})

    for sp_id, entry in entries.items():
        db_dict[sp_id] = {
            'spr_id': sp_id,
            'spr_aformula': entry['aformula'],
            'spr_url': entry['url'],
            'spr_compound': entry['compound'],
            'spr_classification': entry['classification']}
    return db_dict
