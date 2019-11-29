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

from .normalizer import Normalizer
import numpy as np


class DosNormalizer(Normalizer):

    def normalize(self, logger=None) -> None:
        if logger is not None:
            self.logger = logger.bind(normalizer=self.__class__.__name__)

        try:
            # Not all calculations will have DOS, so we check.
            # section_system will always be, so we read it first
            section_system_indices = self._backend.get_sections('section_system')
            section_dos_indices = self._backend.get_sections('section_dos')
        except LookupError:
            return  # absent section_dos

        print("\n####  SECTION_DOS ")
        print('section_dos_indices: ', section_dos_indices)
        for index in section_dos_indices:
            dos = self._backend.get_value('dos_values', index)  # a numpy.ndarray

        print("\n####  SECTION_SYSTEM ")
        print('section_system_indices: ', section_system_indices)
        for index in section_system_indices:
            atom_positions = self._backend.get_value('atom_positions', index)
            lattice_vectors = self._backend.get_value('lattice_vectors', index)
            # simulation_cell = self._backend.get_value('simulation_cell', index)

            number_of_atoms = np.shape(atom_positions)[0]
            unit_cell_volume = np.linalg.det(lattice_vectors)

        # Final quantities
        dos_per_atom = dos / number_of_atoms
        dos_per_unit_volume = dos / unit_cell_volume

        # Add quantities to NOMAD's Metainfo
        self._backend.addArrayValues('dos_values_per_atoms', dos_per_atom)
        self._backend.addArrayValues('dos_values_per_unit_volume', dos_per_unit_volume)

        # (everything below will be deleted in the final version)
        # MARKUS:
        # 1. Shall we save 'unit_cell_volume' to metainfo?
        #    N.B. the description of 'dos_values_per_unit_volume' tells the factor

        # 2. Metainfo backend has two quantities that refer to the same concept:
        #      lattice_vectors
        #      simulation_cell
        #   can these differ in some scenarios?

        # 3. there is one DOS array per spin channel, hence, possible shapes are
        # (1,N) or (2,N). All operations done on 'dos_values' are scalar,
        # hence its shape is preserved.

        # METAINFO:

        # dos_values_per_atoms:
        # "Values (number of states for a given energy divided by
        # the numer of atoms, the set of discrete energy values is given in dos_energies)
        # of density (electronic-energy or vibrational-energy) of states."
        # "shape":
        #     ["number_of_spin_channels", "number_of_dos_values"]

        # dos_values_per_unit_volume:
        # "Values (number of states for a given energy divided by volume, the set of
        # discrete energy values is given in dos_energies) of density (electronic-energy
        # or vibrational-energy) of states.",
        # shape:
        #   ["number_of_spin_channels", "number_of_dos_values"
