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


class DosNormalizer(Normalizer):

    def normalize(self, logger=None) -> None:
        if logger is not None:
            self.logger = logger.bind(normalizer=self.__class__.__name__)

        # loop over all sccs (few)
        #    loop over section_dos's (nothing or one)

        # get values stuff
        # dos_index = for all section dos
        # sccs_index = backend.get_parent_section(dos_index)
        # system_index = backend.get_value("single_config_system_ref", sccs_index)
        # number of atoms = backend.getValue(@number of atoms, szstem.index)

        # Along this lines, try out
        dos_section_indices = self._backend.get_sections('section_dos')
        for index in dos_section_indices:
            values = self._backend.get_value('dos_values', index)
            print(values)
            print("\n\n\n\n# ", index)
