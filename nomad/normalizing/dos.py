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

        # 'scc': single_configuration_calculation
        section_scc_indices = self._backend.get_sections('section_single_configuration_calculation')

        for scc_index in section_scc_indices:
            section_dos_indices = self._backend.get_sections('section_dos', scc_index)

            for dos_index in section_dos_indices:
                dos = self._backend.get_value('dos_values', dos_index)  # a numpy.ndarray

                system_index = self._backend.get_value(
                    'single_configuration_calculation_to_system_ref', scc_index)

                atom_positions = self._backend.get_value('atom_positions', system_index)
                lattice_vectors = self._backend.get_value('lattice_vectors', system_index)

                number_of_atoms = np.shape(atom_positions)[0]
                unit_cell_volume = np.linalg.det(lattice_vectors)

                # Final quantities
                dos_normed = dos / (number_of_atoms * unit_cell_volume)

                # Add quantities to NOMAD's Metainfo
                self._backend.addArrayValues('dos_values_normalized', dos_normed, dos_index)
