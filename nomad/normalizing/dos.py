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

from nomad import config
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
                dos_energies = dos.dos_energies
                energy_reference_fermi = scc.energy_reference_fermi
                if dos_energies is None or dos_values is None or energy_reference_fermi is None:
                    continue

                # Normalize DOS values to be 1/J/atom/m^3
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
                dos_values_normalized = dos_values / (number_of_atoms * unit_cell_volume)

                # Normalize energies so that they are normalized to HOMO.
                i_channel = 0
                fermi_energy = energy_reference_fermi[i_channel]
                fermi_idx = (np.abs(dos_energies - fermi_energy)).argmin()
                energy_threshold = config.normalize.band_structure_energy_tolerance
                value_threshold = 1e-8  # The DOS value that is considered to be zero
                homo_found = False
                zero_found = False
                energy_reference = fermi_energy

                # Walk through the energies in descencing direction to see if a
                # gap is nearby (see energy_threshold). If gap found, continue
                # until HOMO found
                idx = fermi_idx
                while True:
                    try:
                        value = dos_values_normalized[i_channel, idx]
                        energy_distance = fermi_energy - dos_energies[idx]
                    except IndexError:
                        break
                    if energy_distance.magnitude > energy_threshold and not zero_found:
                        break
                    if value <= value_threshold:
                        zero_found = True
                    if zero_found and value > value_threshold:
                        energy_reference = dos_energies[idx + 1]
                        homo_found = True
                        break
                    idx -= 1
                # If gap was not found in descending direction, check the
                # ascending direction for a nearby (see energy_threshold) HOMO
                # value
                if not homo_found:
                    idx = fermi_idx + 1
                    while True:
                        try:
                            value = dos_values_normalized[i_channel, idx]
                            energy_distance = dos_energies[idx] - fermi_energy
                        except IndexError:
                            break
                        if energy_distance.magnitude > energy_threshold:
                            break
                        if value <= value_threshold:
                            energy_reference = dos_energies[idx]
                            break
                        idx += 1

                dos_energies_normalized = dos_energies - energy_reference

                # Data for DOS fingerprint
                dos_fingerprint = None
                try:
                    dos_fingerprint = DOSFingerprint().calculate(
                        dos_energies_normalized.magnitude,
                        dos_values_normalized,
                        n_atoms=number_of_atoms
                    )
                except Exception as e:
                    self.logger.error('could not generate dos fingerprint', exc_info=e)

                # Add quantities to NOMAD's Metainfo
                scc_url = '/section_run/0/section_single_configuration_calculation/%d/section_dos/0' % scc.m_parent_index
                self._backend.openContext(scc_url)
                dos.dos_values_normalized = dos_values_normalized
                dos.dos_energies_normalized = dos_energies_normalized
                if dos_fingerprint is not None:
                    sec_dos_fingerprint = dos.m_create(section_dos_fingerprint)
                    sec_dos_fingerprint.bins = dos_fingerprint.bins
                    sec_dos_fingerprint.indices = dos_fingerprint.indices
                    sec_dos_fingerprint.stepsize = dos_fingerprint.stepsize
                    sec_dos_fingerprint.grid_id = dos_fingerprint.grid_id
                    sec_dos_fingerprint.filling_factor = dos_fingerprint.filling_factor
                self._backend.closeContext(scc_url)
