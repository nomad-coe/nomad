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

from nomad.normalizing import Normalizer
from symmetrynormalizer.symmetry_analysis import normalize


class SymmetryNormalizer(Normalizer):
    """
    This is basically a copy of the legace NOMAD-coe symmetry normalizer.
    """

    def normalize(self) -> None:
        quantities = [
            'atom_labels',
            'atom_positions',
            'lattice_vectors',
            'simulation_cell',
            'configuration_periodic_dimensions'
        ]

        for g_index in self._backend.get_sections('section_system'):
            input_data = dict()
            for quantity in quantities:
                try:
                    input_data[quantity] = self._backend.get_value(quantity, g_index)
                except KeyError:
                    pass

            normalize(self._backend, input_data)
