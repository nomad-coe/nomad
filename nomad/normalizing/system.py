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

from typing import Any
import ase
import numpy as np
import matid
import json

from matid import SymmetryAnalyzer, Classifier

from nomad import utils
from nomad.normalizing.normalizer import SystemBasedNormalizer


class SystemNormalizer(SystemBasedNormalizer):

    """
    This normalizer performs all system (atoms, cells, etc.) related normalizations
    of the legacy NOMAD-coe *stats* normalizer.
    """
    def __init__(self, backend):
        super().__init__(backend, all_sections=True)

    @staticmethod
    def atom_label_to_num(atom_label):
        # Take first three characters and make first letter capitalized.
        atom_label = atom_label[:3].title()

        for symbol_length in reversed(range(1, 4)):
            symbol = atom_label[:symbol_length]
            if symbol in ase.data.chemical_symbols:
                return ase.data.chemical_symbols.index(symbol)

        return 0

    def normalize_system(self, index) -> None:
        """
        The 'main' method of this :class:`SystemBasedNormalizer`.
        Normalizes the section with the given `index`.
        Normalizes geometry, classifies, system_type, and runs symmetry analysis.
        """

        def get_value(key: str, default: Any = None, nonp: bool = False) -> Any:
            try:
                value = self._backend.get_value(key, index)
                if nonp and type(value).__module__ == np.__name__:
                    value = value.tolist()
                return value
            except KeyError:
                return default

        def set_value(key: str, value: Any):
            self._backend.addValue(key, value)

        # analyze atoms labels
        atom_labels = get_value('atom_labels', nonp=True)
        atom_species = get_value('atom_species', nonp=True)
        if atom_labels is None and atom_species is None:
            self.logger.error('calculation has neither atom species nor labels')
            return

        if atom_labels is None:
            atom_labels = list(ase.data.chemical_symbols[species] for species in atom_species)
        else:
            atom_labels = atom_labels

        symbols = ''.join(atom_labels)
        symbols = symbols.replace('1', '')
        try:
            atoms = ase.Atoms(symbols=symbols)
        except Exception as e:
            self.logger.error(
                'cannot build ase atoms from atom labels',
                atom_labels=atom_labels, exc_info=e)
            return
        chemical_symbols = list(atoms.get_chemical_symbols())
        if atom_labels != chemical_symbols:
            self.logger.error('atom labels are ambiguous', atom_labels=atom_labels)
            return

        if atom_species is None:
            atom_species = atoms.get_atomic_numbers().tolist()
            set_value('atom_species', atom_species)
        else:
            if atom_species != atoms.get_atomic_numbers().tolist():
                self.logger.warning(
                    'atom species do not match labels',
                    atom_labels=atom_labels, atom_species=atom_species)
                atom_species = atoms.get_atomic_numbers().tolist()
            set_value('atom_species', atom_species)

        # periodic boundary conditions
        pbc = get_value('configuration_periodic_dimensions', nonp=True)
        if pbc is None:
            pbc = [False, False, False]
            self.logger.warning('missing configuration_periodic_dimensions')
            set_value('configuration_periodic_dimensions', pbc)
        try:
            atoms.set_pbc(pbc)
        except Exception as e:
            self.logger.error('cannot use pbc with ase atoms', exc_info=e, pbc=pbc)
            return

        # formulas
        formula = atoms.get_chemical_formula(mode='all')
        formula_reduced = atoms.get_chemical_formula(mode='reduce')
        if any(atoms.pbc):
            formula_bulk = formula_reduced
        else:
            formula_bulk = formula
        set_value('chemical_composition', formula)
        set_value('chemical_composition_reduced', formula_reduced)
        set_value('chemical_composition_bulk_reduced', formula_bulk)

        # positions
        atom_positions = get_value('atom_positions', None)
        if atom_positions is None:
            self.logger.warning('no atom positions, skip further system analysis')
            return
        if len(atom_positions) != atoms.get_number_of_atoms():
            self.logger.error(
                'len of atom position does not match number of atoms',
                n_atom_positions=len(atom_positions), n_atoms=atoms.get_number_of_atoms())
            return
        try:
            atoms.set_positions(1e10 * atom_positions)
        except Exception as e:
            self.logger.error('cannot use positions with ase atoms', exc_info=e)
            return

        # lattice vectors
        lattice_vectors = get_value('lattice_vectors')
        if lattice_vectors is None:
            lattice_vectors = get_value('simulation_cell')
            if lattice_vectors is not None:
                set_value('lattice_vectors', lattice_vectors)
        if lattice_vectors is None and any(pbc):
            self.logger.error('no lattice vectors but periodicity', pbc=pbc)
        else:
            try:
                atoms.set_cell(1e10 * lattice_vectors)
            except Exception as e:
                self.logger.error('cannot use lattice_vectors with ase atoms', exc_info=e)
                return

        # configuration
        configuration = [
            atom_labels, atoms.positions.tolist(),
            atoms.cell.tolist() if atoms.cell is not None else None,
            atoms.pbc.tolist()]
        configuration_id = utils.hash(json.dumps(configuration).encode('utf-8'))
        set_value('configuration_raw_gid', configuration_id)

        # system type analysis
        if atom_positions is not None:
            try:
                classifier = Classifier()
                system_type = classifier.classify(atoms)
            except Exception:
                self.logger.error('matid project system classification failed')
            else:
                # Convert Matid classification to a Nomad classification.
                system_type = self.map_matid_to_nomad_system_types(atoms, system_type)
                set_value('system_type', system_type)

        # symmetry analysis
        if atom_positions is not None and (lattice_vectors is not None or not any(pbc)):
            self.symmetry_analysis(atoms)

    def symmetry_analysis(self, atoms) -> None:
        """Analyze the symmetry of the material bein simulated.

        We feed in the parsed values in section_system to the
        the symmetry analyzer. We then use the Matid library
        to classify the system as 0D, 1D, 2D or 3D and more specific
        when possible. When lattice vectors or simulation cells are
        not present we skip this analysis.

        Args:
            None: We feed in the bakend and atoms object from the
            SymmetryAndType normalizer.

        Returns:
            None: The method should write symmetry variables
            to the backend which is member of this class.
        """
        # Try to use Matid's symmetry analyzer to anlyze the ASE object.
        # TODO: dts, find out what the symmetry_tol does.
        try:
            symm = SymmetryAnalyzer(atoms, symmetry_tol=0.1)

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
            point_group = symm.get_point_group()

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

        except Exception as e:
            self.logger.error('matid symmetry analysis fails with exception', exc_info=e)
            return

        # Write data extracted from Matid symmetry analysis to the backend.
        symGid = self._backend.openSection('section_symmetry')
        # TODO: @dts, should we change the symmetry_method to MATID?
        self._backend.addValue('symmetry_method', 'Matid (spg)')
        self._backend.addValue('space_group_number', space_group_number)
        self._backend.addValue('hall_number', hall_number)
        self._backend.addValue('hall_symbol', hall_symbol)
        self._backend.addValue('international_short_symbol', international_short)
        self._backend.addValue('point_group', point_group)
        self._backend.addValue('crystal_system', crystal_system)
        self._backend.addValue('bravais_lattice', bravais_lattice)
        self._backend.addArrayValues('origin_shift', origin_shift)
        self._backend.addArrayValues('transformation_matrix', transform)

        stdGid = self._backend.openSection('section_std_system')
        self._backend.addArrayValues('lattice_vectors_std', conv_cell)
        self._backend.addArrayValues('atom_positions_std', conv_pos)
        self._backend.addArrayValues('atomic_numbers_std', conv_num)
        self._backend.addArrayValues('wyckoff_letters_std', conv_wyckoff)
        self._backend.addArrayValues('equivalent_atoms_std', conv_equivalent_atoms)
        self._backend.closeSection('section_std_system', stdGid)

        primGid = self._backend.openSection('section_primitive_system')
        self._backend.addArrayValues('lattice_vectors_primitive', prim_cell)
        self._backend.addArrayValues('atom_positions_primitive', prim_pos)
        self._backend.addArrayValues('atomic_numbers_primitive', prim_num)
        self._backend.addArrayValues('wyckoff_letters_primitive', prim_wyckoff)
        self._backend.addArrayValues('equivalent_atoms_primitive', prim_equivalent_atoms)
        self._backend.closeSection('section_primitive_system', primGid)

        origGid = self._backend.openSection('section_original_system')
        self._backend.addArrayValues('wyckoff_letters_original', orig_wyckoff)
        self._backend.addArrayValues('equivalent_atoms_original', orig_equivalent_atoms)
        self._backend.closeSection('section_original_system', origGid)

        self._backend.closeSection('section_symmetry', symGid)

    # Create a class static dictionary for mapping Matid classifications
    # to Nomad classifications.
    translation_dict = {
        matid.classifications.Class0D: 'Atom',
        matid.classifications.Class1D: '1D',
        matid.classifications.Material2D: '2D',
        matid.classifications.Surface: 'Surface',
        matid.classifications.Class2DWithCell: '2D',
        matid.classifications.Class2D: '2D',
        matid.classifications.Class3D: 'Bulk',
        matid.classifications.Unknown: 'Unknown'
    }

    def map_matid_to_nomad_system_types(self, atoms, system_type):
        """ We map the system type classification from matid to Nomad values.

        Args:
            system_type: Object of a matid class representing a
                material classification.
        Returns:
            nomad_classification: String representing a material
                classification that fits into Nomad's current way
                of naming material classes.
        """
        nomad_classification = None

        for matid_class in SystemNormalizer.translation_dict:
            if isinstance(system_type, matid_class):
                nomad_classification = SystemNormalizer.translation_dict[matid_class]
                break
        # Check to make sure a match was found in translating classes.
        if (nomad_classification is None) or (nomad_classification == 'Unknown'):
            # Then something unexpected has happened with our system_type.
            self.logger.error(
                'Matid classfication has given us an unexpected type: %s' % system_type)

        if nomad_classification == 'Atom' and (atoms.get_number_of_atoms() > 1):
            nomad_classification = 'Molecule / Cluster'

        return nomad_classification
