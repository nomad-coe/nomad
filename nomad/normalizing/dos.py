#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import numpy as np
from nptyping import NDArray

from nomad import config
from nomad_dos_fingerprints import DOSFingerprint
from nomad.datamodel.metainfo.public import (
    section_dos_fingerprint,
    Dos,
    ChannelInfo,
)
from nomad.atomutils import get_volume

from .normalizer import Normalizer


class DosNormalizer(Normalizer):
    """Normalizer with the following responsibilities:

      - Determines highest occupied and lowest unocupied energies for a DOS.
      - Adds normalized (intensive) DOS values which are not tied to the
        current simulation cell size.
    """
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

            energy_fermi = scc.energy_reference_fermi
            energy_highest = scc.energy_reference_highest_occupied
            energy_lowest = scc.energy_reference_lowest_unoccupied
            for dos in section_dos:
                dos_values = dos.dos_values
                dos_energies = dos.dos_energies

                # Normalize DOS values to be 1/J/atom/m^3
                if dos_values is not None:
                    system = scc.single_configuration_calculation_to_system_ref
                    if system is None:
                        self.logger.error('referenced system for dos calculation could not be found')
                        continue
                    atom_positions = system.atom_positions
                    lattice_vectors = system.lattice_vectors
                    if atom_positions is None:
                        self.logger.error('required quantity atom_positions is not available')
                        continue
                    if lattice_vectors is None:
                        self.logger.error('required quantity lattice_vectors is not available')
                        continue
                    number_of_atoms = np.shape(atom_positions)[0]
                    unit_cell_volume = get_volume(lattice_vectors.magnitude)
                    dos_values_normalized = dos_values / (number_of_atoms * unit_cell_volume)
                    dos.dos_values_normalized = dos_values_normalized

                # Add energy references
                self.add_energy_references(dos, energy_fermi, energy_highest, energy_lowest)

                # Save normalized values. The normalization is performed
                # against the spin channel with higher highest occupied energy.
                # TODO: The normalized versions should not be stored after all
                # data has been migrated to store the highest occupied energy
                # explicitly.
                normalization_reference = None
                for info in dos.channel_info:
                    energy_highest = info.energy_highest_occupied
                    if energy_highest is not None:
                        if normalization_reference is None:
                            normalization_reference = energy_highest
                        else:
                            normalization_reference = max(normalization_reference, energy_highest)
                if normalization_reference is not None:
                    dos_energies_normalized = dos_energies - normalization_reference
                    dos.dos_energies_normalized = dos_energies_normalized

                    # Data for DOS fingerprint
                    try:
                        dos_fingerprint = DOSFingerprint().calculate(
                            dos_energies_normalized.magnitude,
                            dos_values_normalized,
                            n_atoms=number_of_atoms
                        )
                    except Exception as e:
                        self.logger.error('could not generate dos fingerprint', exc_info=e)
                    else:
                        sec_dos_fingerprint = dos.m_create(section_dos_fingerprint)
                        sec_dos_fingerprint.bins = dos_fingerprint.bins
                        sec_dos_fingerprint.indices = dos_fingerprint.indices
                        sec_dos_fingerprint.stepsize = dos_fingerprint.stepsize
                        sec_dos_fingerprint.grid_id = dos_fingerprint.grid_id
                        sec_dos_fingerprint.filling_factor = dos_fingerprint.filling_factor

    def add_energy_references(
            self,
            dos: Dos,
            energy_fermi: NDArray,
            energy_highest: NDArray,
            energy_lowest: NDArray) -> None:
        """Given the band structure and information about energy references,
        determines the energy references separately for all spin channels.
        """
        # No reference available
        eref = energy_highest if energy_fermi is None else energy_fermi
        if eref is None:
            self.logger.info("could not resolve energy references for dos")
            return

        # Create energy reference sections for each spin channel, add fermi
        # energy if present
        dos_values = dos.dos_values
        dos_values_normalized = dos.dos_values_normalized
        dos_energies = dos.dos_energies
        n_channels = dos_values.shape[0]
        for i_channel in range(n_channels):
            info = ChannelInfo()
            info.index = i_channel
            if energy_highest is not None:
                info.energy_highest_occupied = energy_highest[i_channel]
            if energy_lowest is not None:
                info.energy_lowest_unoccupied = energy_lowest[i_channel]
            if energy_fermi is not None:
                info.energy_fermi = energy_fermi[i_channel]
            dos.m_add_sub_section(Dos.channel_info, info)

        # Use a reference energy (fermi or highest occupied) to determine the
        # energy references from the DOS (discretization will affect the exact
        # location).
        energy_threshold = config.normalize.band_structure_energy_tolerance
        value_threshold = 1e-8  # The DOS value that is considered to be zero
        for i_channel in range(n_channels):
            info = dos.channel_info[i_channel]
            i_eref = eref[i_channel]
            fermi_idx = (np.abs(dos_energies - i_eref)).argmin()

            # First check that the closest dos energy to energy reference
            # is not too far away. If it is very far away, the
            # normalization may be very inaccurate and we do not report it.
            fermi_energy_closest = dos_energies[fermi_idx]
            distance = np.abs(fermi_energy_closest - i_eref)
            if distance.magnitude <= energy_threshold:

                # See if there are zero values close below the energy reference.
                idx = fermi_idx
                idx_descend = fermi_idx
                while True:
                    try:
                        value = dos_values_normalized[i_channel, idx]
                        energy_distance = np.abs(i_eref - dos_energies[idx])
                    except IndexError:
                        break
                    if energy_distance.magnitude > energy_threshold:
                        break
                    if value <= value_threshold:
                        idx_descend = idx
                        break
                    idx -= 1

                # See if there are zero values close above the fermi energy.
                idx = fermi_idx
                idx_ascend = fermi_idx
                while True:
                    try:
                        value = dos_values_normalized[i_channel, idx]
                        energy_distance = np.abs(i_eref - dos_energies[idx])
                    except IndexError:
                        break
                    if energy_distance.magnitude > energy_threshold:
                        break
                    if value <= value_threshold:
                        idx_ascend = idx
                        break
                    idx += 1

                # If there is a single peak at fermi energy, no
                # search needs to be performed.
                if idx_ascend != fermi_idx and idx_descend != fermi_idx:
                    info.energy_highest_occupied = fermi_energy_closest
                    info.energy_lowest_unoccupied = fermi_energy_closest
                    continue

                # Look for highest occupied energy below the descend index
                idx = idx_descend
                while True:
                    try:
                        value = dos_values_normalized[i_channel, idx]
                    except IndexError:
                        break
                    if value > value_threshold:
                        idx = idx if idx == idx_descend else idx + 1
                        info.energy_highest_occupied = dos_energies[idx]
                        break
                    idx -= 1

                # Look for lowest unoccupied energy above idx_ascend
                idx = idx_ascend
                while True:
                    try:
                        value = dos_values_normalized[i_channel, idx]
                    except IndexError:
                        break
                    if value > value_threshold:
                        idx = idx if idx == idx_ascend else idx - 1
                        info.energy_lowest_unoccupied = dos_energies[idx]
                        break
                    idx += 1
