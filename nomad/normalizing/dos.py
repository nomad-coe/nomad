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

import numpy as np

from nomad_dos_fingerprints import DOSFingerprint

from nomad.datamodel.metainfo.public import section_dos_fingerprint
from nomad.atomutils import get_volume

from .normalizer import Normalizer


class DosNormalizer(Normalizer):

    def normalize(self, logger=None) -> None:

        if logger is not None:
            self.logger = logger.bind(normalizer=self.__class__.__name__)

        # Do nothing if section_run is not present
        if self.section_run is None:
            return

        # 'scc': single_configuration_calculation
        section_sccs = self.section_run.section_single_configuration_calculation
        if section_sccs is None:
            return

        for scc in section_sccs:
            section_dos = scc.section_dos
            if section_dos is None:
                continue

            for dos in section_dos:
                dos_values = dos.dos_values
                if dos_values is None:
                    # section dos without dos_values
                    continue

                system = scc.single_configuration_calculation_to_system_ref
                if system is None:
                    self.logger.error('referenced system for dos calculation could not be found')
                    continue

                atom_positions = system.atom_positions
                lattice_vectors = system.lattice_vectors
                if atom_positions is None:
                    self.logger.error('required quantity atom_positions is not available')
                    return
                if lattice_vectors is None:
                    self.logger.error('required quantity lattice_vectors is not available')
                    return

                number_of_atoms = np.shape(atom_positions)[0]
                unit_cell_volume = get_volume(lattice_vectors.magnitude)

                # Final quantities
                dos_normed = dos_values / (number_of_atoms * unit_cell_volume)

                # Data for DOS fingerprint
                dos_fingerprint = None
                try:
                    dos_energies = dos.dos_energies
                    dos_fingerprint = DOSFingerprint().calculate(np.array(dos_energies), dos_values)
                except Exception as e:
                    self.logger.error('could not generate dos fingerprint', exc_info=e)

                # Add quantities to NOMAD's Metainfo
                scc_url = '/section_run/0/section_single_configuration_calculation/%d/section_dos/0' % scc.m_parent_index
                self._backend.openContext(scc_url)
                dos.dos_values_normalized = dos_normed
                if dos_fingerprint is not None:
                    sec_dos_fingerprint = dos.m_create(section_dos_fingerprint)
                    sec_dos_fingerprint.bins = dos_fingerprint.bins
                    sec_dos_fingerprint.indices = dos_fingerprint.indices
                    sec_dos_fingerprint.stepsize = dos_fingerprint.stepsize
                    sec_dos_fingerprint.grid_id = dos_fingerprint.grid_id
                    sec_dos_fingerprint.filling_factor = dos_fingerprint.filling_factor
                self._backend.closeContext(scc_url)
