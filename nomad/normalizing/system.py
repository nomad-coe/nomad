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

import ase
import numpy as np
import sys
import matid

from matid import SymmetryAnalyzer, Classifier

from nomadcore.json_support import addShasOfJson
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
        atom_label = atom_label[:3].title()

        for symbol_length in reversed(range(1, 4)):
            symbol = atom_label[:symbol_length]
            if symbol in ase.data.chemical_symbols:
                return ase.data.chemical_symbols.index(symbol)

        return 0

    def normalize_system(self, section_system) -> None:
        """ Main normalizer that runs system, syste_type and symmetry analysis."""

        self.atom_labels = section_system['atom_labels']
        self.atom_species = section_system['atom_atom_numbers']
        self.atom_positions = section_system['atom_positions']
        self.periodic_dirs = section_system['configuration_periodic_dimensions']
        # Try to first read the cell information from the renamed metainfo
        # lattice_vectors, if this doesn't work try the depreciated name
        # simulation_cell. Otherwise, if neither are present, assign None.
        self.cell = section_system.get(
            'lattice_vectors', section_system.get('simulation_cell', None)
        )
        # Run a system analysis on the system.
        self.system_analysis()

        if self.cell is None:
            # Then the parser hasn't parsed any information about periodicity.
            # We therefore assume we're simulating a single cell without
            # periodicity and don't try to ascertain symmetry information.
            self.atoms = ase.Atoms(
                positions=1e10 * np.asarray(self.atom_positions),
                symbols=np.asarray(self.atom_labels)
            )
            # Classify the material's system type.
            self.system_type_classification()
            if self.nomad_system_type not in ['Atom', 'Molecule / Cluster']:
                self.logger.error(
                    'Matid has classified the sim as more than 1D despite having'
                    'no simulation_cell/lattice_vectors'
                )
            # Return w/out symmetry analysis since we don't have a sim_cell.
            return None

        self.pbc = section_system.get('configuration_periodic_dimensions', None)
        # If no pbc is found assume there is no periodicity.
        if self.pbc is None:
            self.pbc = np.array([False, False, False])
        # The pbc should be defined as a single-dimensional list.
        if len(np.asarray(self.pbc).shape) == 2:
            self.pbc = self.pbc[0, :]

        # Build an ASE atoms object to feed into Matid.
        try:
            self.atoms = ase.Atoms(
                positions=1e10 * np.asarray(self.atom_positions),
                symbols=np.asarray(self.atom_labels),
                cell=1e10 * np.asarray(self.cell),
                pbc=self.pbc
            )
        except Exception:
            self.logger.error(
                'The ASE library is unable to build an object from the member'
                'variables.'
            )

        # Classify the material's system type.
        self.system_type_classification()
        # Analyze the symmetry of the material.
        self.symmetry_analysis()

    def system_analysis(self) -> None:
        """Analyze system properties of a simulation from parsed values."""
        results = dict()
        if self.atom_labels is not None and self.atom_species is None:
            atom_label_to_num = SystemNormalizer.atom_label_to_num
            self.atom_species = [
                atom_label_to_num(atom_label) for atom_label in self.atom_labels
            ]

        formula = None

        if self.atom_species:
            results['atom_species'] = self.atom_species
            atom_symbols = [
                ase.data.chemical_symbols[atom_number] for atom_number in self.atom_species
            ]
            formula = ase.Atoms(atom_symbols).get_chemical_formula(mode='all')
            formula_reduced = ase.Atoms(atom_symbols).get_chemical_formula(mode='reduce')
            if self.periodic_dirs is not None and any(self.periodic_dirs):
                formula_bulk = formula_reduced
            else:
                formula_bulk = formula

        if self.cell is not None:
            results['lattice_vectors'] = self.cell

        if self.atom_positions is not None:
            results['atom_positions'] = self.atom_positions
            if not formula:
                formula = (
                    'X%d' % len(self.atom_positions) if len(self.atom_positions) != 1 else 'X'
                )

        if self.periodic_dirs is not None:
            results['configuration_periodic_dimensions'] = self.periodic_dirs.tolist()
        # TODO: @dts, might be good to clean this up so it is more readable in the
        # future.
        configuration_id = 's' + addShasOfJson(results).b64digests()[0][0:28]

        self._backend.addValue('configuration_raw_gid', configuration_id)
        self._backend.addValue('atom_species', self.atom_species)
        self._backend.addValue('chemical_composition', formula)
        self._backend.addValue('chemical_composition_reduced', formula_reduced)
        self._backend.addValue('chemical_composition_bulk_reduced', formula_bulk)

    def symmetry_analysis(self) -> None:
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
            symm = SymmetryAnalyzer(self.atoms, symmetry_tol=0.1)

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

        except Exception:
            self.logger.error(
                'The matid project symmetry analyzer fails on the ASE'
                ' object from this section.'
            )
            return None  # Without trying to write any symmetry data.
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
        # nomad-xt: context already closed in nomad-xt.
        # backend.closeContext(context)
        sys.stdout.flush()

    def system_type_classification(self) -> None:
        """Try to classify the ASE materials object using Matid's classification."""
        try:
            # Define the classifier as Matid's Classifier that we've imported.
            classifier = Classifier()
            # Perform classification on the atoms ASE object.
            matid_system_type = classifier.classify(self.atoms)
        except Exception:
            self.logger.error(
                'The matid project clsasification fails on the ASE'
                ' object from this section.'
            )
            return None  # Without saving any system type value.
        # Convert Matid classification to a Nomad classification.
        self.nomad_system_type = self.map_matid_to_nomad_system_types(matid_system_type)
        self._backend.addValue('system_type', self.nomad_system_type)

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

    def map_matid_to_nomad_system_types(self, system_type):
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
            self.logger.error('Matid classfication has given us an unexpected type.'
                              ' No match was found for this system type: %s' % system_type)

        if nomad_classification == 'Atom' and (len(self.atom_labels) > 1):
            nomad_classification = 'Molecule / Cluster'
        return nomad_classification
