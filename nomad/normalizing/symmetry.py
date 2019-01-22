# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
""" This file seeks to find the symmetry and classfication of the material.

We use Matid to analyze the symmetry and classify the material. If we don't have
basic parameters such as lattice vectors or cell structure we simply skip this
normalizer. In the past (pre Jan 2019), we had a seperate normalizer for system type
to classify the material. We have since fused the system type and symmetry normalizers
into one file."""

import numpy as np
# import os.path
import logging
import sys
from ase import Atoms
from matid import SymmetryAnalyzer, Classifier

from nomad.normalizing.normalizer import SystemBasedNormalizer
# from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
# from nomadcore.parser_backend import JsonParseEventsWriterBackend
# from nomadcore.parse_streamed_dicts import ParseStreamedDicts

# TODO, one normalier, called system normalizer.


class SymmetryAndTypeNormalizer(SystemBasedNormalizer):
    """ This is basically a copy of the legace NOMAD-coe symmetry normalizer."""
    def __init__(self, backend):
        super().__init__(backend, all_sections=True)

    def normalize_system(self, section_system) -> None:
        print('symmetry normalizer.')
        self.atom_labels = section_system["atom_labels"]
        self.atom_pos = section_system["atom_positions"]

        # Try to first read the cell information from the renamed metainfo
        # lattice_vectors, if this doesn't work try the depreciated name
        # simulation_cell. Otherwise, if neither are present, assign None.
        self.cell = section_system.get(
            'lattice_vectors', section_system.get('simulation_cell', None)
        )

        if self.cell is None:
            # Then the parser hasn't parsed any information about periodicity.
            # We therefore assume we're simulating a single cell without
            # periodicity and don't try to ascertain symmetry information.
            return None
        # TODO: dts@, should we assign pbc = [False, False, False]
        # when we can't find the pbc? Either skip symmetry without pbc or continue
        # with False, False, False values. Let's ask around.
        self.pbc = section_system.get('configuration_periodic_dimensions', None)
        if self.pbc is None:
            # Without a pbc value we cannot build an ASE atoms object.
            return None

        # The pbc should be defined as a single-dimensional list.
        if len(np.asarray(self.pbc).shape) == 2:
            self.pbc = self.pbc[0, :]

        # Build an ASE atoms object to feed into Matid.
        try:
            self.atoms = Atoms(
                positions=1e10 * np.asarray(self.atom_pos),
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
            print("crystal_system is: %s" % crystal_system)
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
            return None # Without trying to write any symmetry data.

        # Write data extracted from Matid symmetry analysis to the backend.
        symGid = self._backend.openSection("section_symmetry")
        # TODO: @dts, should we change the symmetry_method to MATID?
        self._backend.addValue("symmetry_method", "spg_normalized")
        self._backend.addValue("space_group_number", space_group_number)
        self._backend.addValue("hall_number", hall_number)
        self._backend.addValue("hall_symbol", hall_symbol)
        self._backend.addValue("international_short_symbol", international_short)
        self._backend.addValue("point_group", point_group)
        self._backend.addValue("crystal_system", crystal_system)
        self._backend.addValue("bravais_lattice", bravais_lattice)
        self._backend.addArrayValues("origin_shift", origin_shift)
        self._backend.addArrayValues("transformation_matrix", transform)

        stdGid = self._backend.openSection("section_std_system")
        self._backend.addArrayValues("lattice_vectors_std", conv_cell)
        self._backend.addArrayValues("atom_positions_std", conv_pos)
        self._backend.addArrayValues("atomic_numbers_std", conv_num)
        self._backend.addArrayValues("wyckoff_letters_std", conv_wyckoff)
        self._backend.addArrayValues("equivalent_atoms_std", conv_equivalent_atoms)
        self._backend.closeSection("section_std_system", stdGid)

        primGid = self._backend.openSection("section_primitive_system")
        self._backend.addArrayValues("lattice_vectors_primitive", prim_cell)
        self._backend.addArrayValues("atom_positions_primitive", prim_pos)
        self._backend.addArrayValues("atomic_numbers_primitive", prim_num)
        self._backend.addArrayValues("wyckoff_letters_primitive", prim_wyckoff)
        self._backend.addArrayValues("equivalent_atoms_primitive", prim_equivalent_atoms)
        self._backend.closeSection("section_primitive_system", primGid)

        origGid = self._backend.openSection("section_original_system")
        self._backend.addArrayValues("wyckoff_letters_original", orig_wyckoff)
        self._backend.addArrayValues("equivalent_atoms_original", orig_equivalent_atoms)
        self._backend.closeSection("section_original_system", origGid)

        self._backend.closeSection("section_symmetry", symGid)
        # nomad-xt: context already closed in nomad-xt.
        # backend.closeContext(context)
        sys.stdout.flush()

    def system_type_classification(self):
        # Try to classify the ASE materials object using Matid's classification.
        # TODO: add mapping of matid to Nomad classifications.
        try:
            # Define the classifier as Matid's Classifier that we've imported.
            classifier = Classifier()
            # Perform classification on the atoms ASE object.
            system_type = classifier.classify(self.atoms)
            # For debug:
            print(" The classification from Matid is: %s" % system_type)
        except Exception:
            self.logger.error(
                'The matid project clsasification fails on the ASE'
                ' object from this section.'
            )
            return None  # Without saving any system type value.
        self._backend.addValue('system_type', system_type)
