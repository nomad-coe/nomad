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

from collections import Counter
from typing import Any
import ase
import numpy as np
import json
import re
import os
import sqlite3
import functools
import fractions

from matid import SymmetryAnalyzer, Classifier
from matid.classifications import Class0D, Atom, Class1D, Material2D, Surface, Class3D

from nomadcore.structure_types import structure_types_by_spacegroup as str_types_by_spg

from nomad import utils, config
from nomad.normalizing.normalizer import SystemBasedNormalizer

# use a regular expression to check atom labels; expression is build from list of
# all labels sorted desc to find Br and not B when searching for Br.
atom_label_re = re.compile('|'.join(
    sorted(ase.data.chemical_symbols, key=lambda x: len(x), reverse=True)))


springer_db_connection = None


def open_springer_database():
    """
    Create a global connection to the Springer database in a way that
    each worker opens the database just once.
    """
    global springer_db_connection
    if springer_db_connection is None:
        # filepath definition in 'nomad-FAIR/nomad/config.py'
        db_file = config.springer_db_path
        if not os.path.exists(db_file):
            utils.get_logger(__name__).error('Springer database not found')
            return None
        springer_db_connection = sqlite3.connect(db_file, check_same_thread=False, uri=True)
        # we lift the thread check because we share the connection among workers
        # 'uri=True': open a database in read-only mode

    return springer_db_connection


def normalized_atom_labels(atom_labels):
    """
    Normalizes the given atom labels: they either are labels right away, or contain
    additional numbers (to distinguish same species but different labels, see meta-info),
    or we replace them with ase placeholder atom for unknown elements 'X'.
    """
    return [
        ase.data.chemical_symbols[0] if match is None else match.group(0)
        for match in [re.search(atom_label_re, atom_label) for atom_label in atom_labels]]


def formula_normalizer(atoms):
    """
    Reads the chemical symbols in ase.atoms and returns a normalized formula.
    Formula normalization is on the basis of atom counting,
    e.g., Tc ->  Tc100, SZn -> S50Zn50, Co2Nb -> Co67Nb33
    """
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
    """
    This normalizer performs all system (atoms, cells, etc.) related normalizations
    of the legacy NOMAD-coe *stats* normalizer.
    """
    @staticmethod
    def atom_label_to_num(atom_label):
        # Take first three characters and make first letter capitalized.
        atom_label = atom_label[:3].title()

        for symbol_length in reversed(range(1, 4)):
            symbol = atom_label[:symbol_length]
            if symbol in ase.data.chemical_symbols:
                return ase.data.chemical_symbols.index(symbol)

        return 0

    def normalize_system(self, index, is_representative) -> bool:
        """
        The 'main' method of this :class:`SystemBasedNormalizer`.
        Normalizes the section with the given `index`.
        Normalizes geometry, classifies, system_type, and runs symmetry analysis.

        Returns: True, iff the normalization was successful
        """

        def get_value(key: str, default: Any = None, numpy: bool = True) -> Any:
            try:
                value = self._backend.get_value(key, index)
                if not numpy and type(value).__module__ == np.__name__:
                    value = value.tolist()

                elif numpy and isinstance(value, list):
                    value = np.array(value)

                return value
            except (KeyError, IndexError):
                return default

        def set_value(key: str, value: Any):
            self._backend.addValue(key, value)

        if is_representative:
            self._backend.addValue('is_representative', is_representative)

        # analyze atoms labels
        atom_labels = get_value('atom_labels', numpy=False)
        if atom_labels is not None:
            atom_labels = normalized_atom_labels(atom_labels)

        atom_species = get_value('atom_species', numpy=False)
        if atom_labels is None and atom_species is None:
            self.logger.error('system has neither atom species nor labels')
            return False

        # If there are no atom labels we create them from atom species data.
        if atom_labels is None:
            try:
                atom_labels = list(ase.data.chemical_symbols[species] for species in atom_species)
            except IndexError:
                self.logger.error('system has atom species that are out of range')
                return False

            self._backend.addArrayValues('atom_labels', atom_labels)

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
            self._backend.addArrayValues('atom_species', atom_species)
        else:
            if not isinstance(atom_species, list):
                atom_species = [atom_species]
            if atom_species != atoms.get_atomic_numbers().tolist():
                self.logger.warning(
                    'atom species do not match labels',
                    atom_labels=atom_labels[:10], atom_species=atom_species[:10])
                atom_species = atoms.get_atomic_numbers().tolist()
            set_value('atom_species', atom_species)

        # periodic boundary conditions
        pbc = get_value('configuration_periodic_dimensions', numpy=False)
        if pbc is None:
            pbc = [False, False, False]
            self.logger.warning('missing configuration_periodic_dimensions')
            set_value('configuration_periodic_dimensions', pbc)
        try:
            atoms.set_pbc(pbc)
        except Exception as e:
            self.logger.error(
                'cannot use pbc with ase atoms', exc_info=e, pbc=pbc, error=str(e))
            return False

        # formulas
        set_value('chemical_composition', atoms.get_chemical_formula(mode='all'))
        set_value('chemical_composition_reduced', atoms.get_chemical_formula(mode='reduce'))
        set_value('chemical_composition_bulk_reduced', atoms.get_chemical_formula(mode='hill'))

        # positions
        atom_positions = get_value('atom_positions', None, numpy=True)
        if atom_positions is None:
            self.logger.warning('no atom positions, skip further system analysis')
            return False
        if len(atom_positions) != atoms.get_number_of_atoms():
            self.logger.error(
                'len of atom position does not match number of atoms',
                n_atom_positions=len(atom_positions), n_atoms=atoms.get_number_of_atoms())
            return False
        try:
            atoms.set_positions(1e10 * atom_positions)
        except Exception as e:
            self.logger.error(
                'cannot use positions with ase atoms', exc_info=e, error=str(e))
            return False

        # lattice vectors
        lattice_vectors = get_value('lattice_vectors', numpy=True)
        if lattice_vectors is None:
            lattice_vectors = get_value('simulation_cell', numpy=True)
            if lattice_vectors is not None:
                set_value('lattice_vectors', lattice_vectors)
        if lattice_vectors is None:
            if any(pbc):
                self.logger.error('no lattice vectors but periodicity', pbc=pbc)
        else:
            try:
                atoms.set_cell(1e10 * np.array(lattice_vectors))
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
        set_value('configuration_raw_gid', configuration_id)

        if is_representative:
            # Save the Atoms as a temporary variable
            self._backend.add_tmp_value("section_system", "representative_atoms", atoms)

            # system type analysis
            if atom_positions is not None:
                with utils.timer(
                        self.logger, 'system classification executed',
                        system_size=atoms.get_number_of_atoms()):
                    self.system_type_analysis(atoms)

            system_type = self._backend.get_value("system_type")

            # symmetry analysis
            if atom_positions is not None and (lattice_vectors is not None or not any(pbc)) and system_type == "bulk":
                with utils.timer(
                        self.logger, 'symmetry analysis executed',
                        system_size=atoms.get_number_of_atoms()):

                    self.symmetry_analysis(atoms)

        return True

    def system_type_analysis(self, atoms: ase.Atoms) -> None:
        """
        Determine the system type with MatID. Write the system type to the
        backend.

        Args:
            atoms
        """
        system_type = config.services.unavailable_value
        if atoms.get_number_of_atoms() <= config.normalize.system_classification_with_clusters_threshold:
            try:
                classifier = Classifier(cluster_threshold=config.normalize.cluster_threshold)
                cls = classifier.classify(atoms)
            except Exception as e:
                self.logger.error(
                    'matid project system classification failed', exc_info=e, error=str(e))
            else:
                if isinstance(cls, Class3D):
                    system_type = 'bulk'
                elif isinstance(cls, Atom):
                    system_type = 'atom'
                elif isinstance(cls, Class0D):
                    system_type = 'molecule / cluster'
                elif isinstance(cls, Class1D):
                    system_type = '1D'
                elif isinstance(cls, Surface):
                    system_type = 'surface'
                elif isinstance(cls, Material2D):
                    system_type = '2D'
        else:
            self.logger.info("System type analysis not run due to large system size.")

        self._backend.addValue('system_type', system_type)

    def symmetry_analysis(self, atoms: ase.Atoms) -> None:
        """Analyze the symmetry of the material being simulated. Only performed
        for bulk materials.

        We feed in the parsed values in section_system to the the symmetry
        analyzer. The analysis results are written to the backend.

        Args:
            atoms: The atomistic structure to analyze.
        """
        # Try to use MatID's symmetry analyzer to analyze the ASE object.
        # TODO: dts, find out what the symmetry_tol does.
        try:
            symm = SymmetryAnalyzer(atoms, config.normalize.symmetry_tolerance)

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

        # Write data extracted from MatID's symmetry analysis to the backend.
        symmetry_gid = self._backend.openSection('section_symmetry')
        self._backend.add_tmp_value("section_symmetry", "symmetry_analyzer", symm)

        # TODO: @dts, should we change the symmetry_method to MATID?
        self._backend.addValue('symmetry_method', 'MatID (spg)')
        self._backend.addValue('space_group_number', space_group_number)
        self._backend.addValue('hall_number', hall_number)
        self._backend.addValue('hall_symbol', hall_symbol)
        self._backend.addValue('international_short_symbol', international_short)
        self._backend.addValue('point_group', point_group)
        self._backend.addValue('crystal_system', crystal_system)
        self._backend.addValue('bravais_lattice', bravais_lattice)
        self._backend.addArrayValues('origin_shift', origin_shift)
        self._backend.addArrayValues('transformation_matrix', transform)

        std_gid = self._backend.openSection('section_std_system')
        self._backend.addArrayValues('lattice_vectors_std', conv_cell)
        self._backend.addArrayValues('atom_positions_std', conv_pos)
        self._backend.addArrayValues('atomic_numbers_std', conv_num)
        self._backend.addArrayValues('wyckoff_letters_std', conv_wyckoff)
        self._backend.addArrayValues('equivalent_atoms_std', conv_equivalent_atoms)
        self._backend.closeSection('section_std_system', std_gid)

        prim_gid = self._backend.openSection('section_primitive_system')
        self._backend.addArrayValues('lattice_vectors_primitive', prim_cell)
        self._backend.addArrayValues('atom_positions_primitive', prim_pos)
        self._backend.addArrayValues('atomic_numbers_primitive', prim_num)
        self._backend.addArrayValues('wyckoff_letters_primitive', prim_wyckoff)
        self._backend.addArrayValues('equivalent_atoms_primitive', prim_equivalent_atoms)
        self._backend.closeSection('section_primitive_system', prim_gid)

        orig_gid = self._backend.openSection('section_original_system')
        self._backend.addArrayValues('wyckoff_letters_original', orig_wyckoff)
        self._backend.addArrayValues('equivalent_atoms_original', orig_equivalent_atoms)
        self._backend.closeSection('section_original_system', orig_gid)
        self._backend.closeSection('section_symmetry', symmetry_gid)

        self.springer_classification(atoms, space_group_number)  # Springer Normalizer

        self.prototypes(prim_num, prim_wyckoff, space_group_number)

        self._backend.closeSection('section_symmetry', symmetry_gid)

    def springer_classification(self, atoms, space_group_number):
        # SPRINGER NORMALIZER
        normalized_formula = formula_normalizer(atoms)
        #
        springer_db_connection = open_springer_database()
        if springer_db_connection is None:
            return

        cur = springer_db_connection.cursor()

        # SQL QUERY
        # (this replaces the four queries done in the old 'classify4me_SM_normalizer.py')
        cur.execute("""
            SELECT
                entry.entry_id,
                entry.alphabetic_formula,
                GROUP_CONCAT(DISTINCT compound_classes.compound_class_name),
                GROUP_CONCAT(DISTINCT classification.classification_name)
            FROM entry
            LEFT JOIN entry_compound_class as ecc ON ecc.entry_nr = entry.entry_nr
            LEFT JOIN compound_classes ON ecc.compound_class_nr = compound_classes.compound_class_nr
            LEFT JOIN entry_classification as ec ON ec.entry_nr = entry.entry_nr
            LEFT JOIN classification ON ec.classification_nr = classification.classification_nr
            LEFT JOIN entry_reference as er ON er.entry_nr = entry.entry_nr
            LEFT JOIN reference ON reference.reference_nr = er.entry_nr
            WHERE entry.normalized_formula = ( %r ) and entry.space_group_number = '%d'
            GROUP BY entry.entry_id;
            """ % (normalized_formula, space_group_number))

        results = cur.fetchall()
        # 'results' is a list of tuples, i.e. '[(a,b,c,d), ..., (a,b,c,d)]'
        # All SQL queries done

        # Storing 'results' in a dictionary
        dbdict = {}
        for ituple in results:
            # 'spr' means 'springer'
            spr_id = ituple[0]
            spr_aformula = ituple[1]  # alphabetical formula
            spr_url = 'http://materials.springer.com/isp/crystallographic/docs/' + spr_id
            spr_compound = ituple[2].split(',')  # split to convert string to list
            spr_classification = ituple[3].split(',')
            #
            spr_compound.sort()
            spr_classification.sort()
            #
            dbdict[spr_id] = {'spr_id': spr_id,
                              'spr_aformula': spr_aformula,
                              'spr_url': spr_url,
                              'spr_compound': spr_compound,
                              'spr_classification': spr_classification}
        # =============

        # SPRINGER's METAINFO UPDATE
        # LAYOUT: Five sections under 'section_springer_material' for each material ID:
        #     id, alphabetical formula, url, compound_class, clasification.
        # As per Markus/Luca's emails, we don't expose Springer bib references (Springer's paywall)
        for material in dbdict.values():
            self._backend.openNonOverlappingSection('section_springer_material')

            self._backend.addValue('springer_id', material['spr_id'])
            self._backend.addValue('springer_alphabetical_formula', material['spr_aformula'])
            self._backend.addValue('springer_url', material['spr_url'])
            self._backend.addArrayValues('springer_compound_class', material['spr_compound'])
            self._backend.addArrayValues('springer_classification', material['spr_classification'])

            self._backend.closeNonOverlappingSection('section_springer_material')

        # Check the 'springer_classification' and 'springer_compound_class' information
        # found is the same for all springer_id's
        dkeys = list(dbdict.keys())
        if len(dkeys) != 0:
            class_0 = dbdict[dkeys[0]]['spr_classification']
            comp_0 = dbdict[spr_id]['spr_compound']

            # compare 'class_0' and 'comp_0' against the rest
            for ii in range(1, len(dkeys)):
                class_test = (class_0 == dbdict[dkeys[ii]]['spr_classification'])
                comp_test = (comp_0 == dbdict[dkeys[ii]]['spr_compound'])

                if (class_test or comp_test) is False:
                    self.logger.warning('Mismatch in Springer classification or compounds')

    def prototypes(self, atomSpecies, wyckoffs, spg_nr):
        try:
            norm_wyckoff = SystemNormalizer.get_normalized_wyckoff(atomSpecies, wyckoffs)
            protoDict = SystemNormalizer.get_structure_type(spg_nr, norm_wyckoff)
            aflow_prototype_notes = None

            if protoDict is None:
                proto = "%d-_" % spg_nr
                labels = dict(prototype_label=proto)
            else:
                proto = '%d-%s-%s' % (spg_nr, protoDict.get("Prototype", "-"),
                                      protoDict.get("Pearsons Symbol", "-"))
                aflow_prototype_id = protoDict.get("aflow_prototype_id", "-")
                aflow_prototype_url = protoDict.get("aflow_prototype_url", "-")
                aflow_prototype_notes = protoDict.get("Notes", None)
                labels = dict(
                    prototype_label=proto,
                    prototype_aflow_id=aflow_prototype_id,
                    prototype_aflow_url=aflow_prototype_url)
        except Exception as e:
            self.logger.error("cannot create AFLOW prototype", exc_info=e)
            return

        pSect = self._backend.openSection("section_prototype")
        if aflow_prototype_notes is not None and aflow_prototype_id != "-":
            self._backend.add_tmp_value("section_prototype", "prototype_notes", aflow_prototype_notes)
        self._backend.addValue(
            "prototype_assignement_method", "normalized-wyckoff")
        self._backend.addValue("prototype_label", labels['prototype_label'])
        aid = labels.get("prototype_aflow_id")
        if aid:
            self._backend.addValue("prototype_aflow_id", aid)
        aurl = labels.get("prototype_aflow_url")
        if aurl:
            self._backend.addValue("prototype_aflow_url", aurl)
        self._backend.closeSection("section_prototype", pSect)

    @staticmethod
    def get_normalized_wyckoff(atomic_number, wyckoff):
        """Returns a normalized Wyckoff sequence."""
        atomCount = {}
        for nr in atomic_number:
            atomCount[nr] = atomCount.get(nr, 0) + 1
        wycDict = {}

        for i, wk in enumerate(wyckoff):
            oldVal = wycDict.get(wk, {})
            nr = atomic_number[i]
            oldVal[nr] = oldVal.get(nr, 0) + 1
            wycDict[wk] = oldVal
        sortedWyc = list(wycDict.keys())
        sortedWyc.sort()

        def cmpp(a, b):
            return ((a < b) - (a > b))

        def compareAtNr(at1, at2):
            c = cmpp(atomCount[at1], atomCount[at2])
            if (c != 0):
                return c
            for wk in sortedWyc:
                p = wycDict[wk]
                c = cmpp(p.get(at1, 0), p.get(at2, 0))
                if c != 0:
                    return c
            return 0

        sortedAt = list(atomCount.keys())
        sortedAt.sort(key=functools.cmp_to_key(compareAtNr))
        standardAtomNames = {}
        for i, at in enumerate(sortedAt):
            standardAtomNames[at] = ("X_%d" % i)
        standardWyc = {}
        for wk, ats in wycDict.items():
            stdAts = {}
            for at, count in ats.items():
                stdAts[standardAtomNames[at]] = count
            standardWyc[wk] = stdAts
        if standardWyc:
            counts = [c for x in standardWyc.values() for c in x.values()]
            gcd = counts[0]
            for c in counts[1:]:
                gcd = fractions.gcd(gcd, c)
            if gcd != 1:
                for wk, d in standardWyc.items():
                    for at, c in d.items():
                        d[at] = c // gcd
        return standardWyc

    @staticmethod
    def get_structure_type(space_group, norm_wyckoff):
        """Returns the information on the prototype.
        """
        structure_type_info = {}

        type_descriptions = str_types_by_spg.get(space_group, [])
        for type_description in type_descriptions:
            current_norm_wyckoffs = type_description.get("normalized_wysytax")
            if current_norm_wyckoffs and current_norm_wyckoffs == norm_wyckoff:
                structure_type_info = type_description
                break
        if structure_type_info:
            return structure_type_info
        else:
            return None
