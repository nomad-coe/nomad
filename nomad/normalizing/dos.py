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
from nomad_dos_fingerprints import DOSFingerprint  # pylint: disable=import-error
from nomad.datamodel.metainfo.simulation.calculation import (
    Dos, DosFingerprint, BandGap)

from .normalizer import Normalizer


class DosNormalizer(Normalizer):
    """Normalizer with the following responsibilities:

      - Determines highest occupied and lowest unoccupied energies for both
        spin channels.
      - Adds normalization factor for intensive (=not tied to the imulation
        cell size) DOS values which are
    """
    def normalize(self, logger=None) -> None:

        if logger is not None:
            self.logger = logger.bind(normalizer=self.__class__.__name__)

        # Do nothing if section_run is not present
        if self.section_run is None:
            return

        calculations = self.section_run.calculation
        if calculations is None:
            return

        for calc in calculations:
            # Normalize electronic DOS
            dos_electronic = calc.dos_electronic
            if dos_electronic is not None:
                energy_fermi = calc.energy.fermi if calc.energy is not None else None
                energy_highest = calc.energy.highest_occupied if calc.energy is not None else None
                energy_lowest = calc.energy.lowest_unoccupied if calc.energy is not None else None
                for dos in dos_electronic:

                    # Add normalization factor
                    set_normalization_factor = self.add_electronic_normalization_factor(calc, dos)
                    if not set_normalization_factor:
                        continue

                    # Add energy references
                    dos_values = [dos_total.value.magnitude for dos_total in dos.total]
                    self.add_energy_references(dos, energy_fermi, energy_highest, energy_lowest, dos_values)

                    # Calculate the DOS fingerprint for successfully normalized DOS
                    normalization_reference = None
                    for info in dos.band_gap:
                        energy_highest = info.energy_highest_occupied
                        if energy_highest is not None:
                            if normalization_reference is None:
                                normalization_reference = energy_highest
                            else:
                                normalization_reference = max(normalization_reference, energy_highest)
                    if normalization_reference is not None:
                        dos_energies_normalized = dos.energies - normalization_reference

                        try:
                            dos_fingerprint = DOSFingerprint().calculate(
                                dos_energies_normalized.magnitude,
                                dos_values
                            )
                        except Exception as e:
                            self.logger.error('could not generate dos fingerprint', exc_info=e)
                        else:
                            sec_dos_fingerprint = dos.m_create(DosFingerprint)
                            sec_dos_fingerprint.bins = dos_fingerprint.bins
                            sec_dos_fingerprint.indices = dos_fingerprint.indices
                            sec_dos_fingerprint.stepsize = dos_fingerprint.stepsize
                            sec_dos_fingerprint.grid_id = dos_fingerprint.grid_id
                            sec_dos_fingerprint.filling_factor = dos_fingerprint.filling_factor

            # Normalize phonon DOS
            dos_phonons = calc.dos_phonon
            if dos_phonons is not None:
                for dos_phonon in dos_phonons:
                    self.add_phononic_normalization_factor(calc, dos_phonon)

    def add_electronic_normalization_factor(self, calc, dos):
        """Returns a factor that returns a size intensive electronic DOS.
        The values are divided by integral(DOS, lowest state, Fermi energy), or likewise sum(<atomic numbers>)."""
        atoms = calc.system_ref.atoms
        if not len(dos.total):
            self.logger.error('referenced system for dos calculation could not be found')
            return
        elif not len(atoms.species):
            self.logger.error('referenced system for atomic species could not be found')
            return
        else:
            normalization_factor = 1 / sum(atoms.species)
            for dos_total in dos.total:
                dos_total.normalization_factor = normalization_factor
            return normalization_factor

    def add_phononic_normalization_factor(self, calc, dos):
        """Returns a factor that returns a size intensive phononic DOS.
        The values are divided by integral(DOS, 0, infinity), or likewise <no. degrees of freedom>"""
        atoms = calc.system_ref.atoms
        if not len(dos.total):
            self.logger.error('referenced system for dos calculation could not be found')
            return
        elif not len(atoms.species):
            self.logger.error('referenced system for the number of atoms could not be found')
            return
        else:
            normalization_factor = 1 / (3 * len(atoms.species))
            for dos_total in dos.total:
                dos_total.normalization_factor = normalization_factor
            return normalization_factor

    def add_energy_references(
            self,
            dos: Dos,
            energy_fermi: float,
            energy_highest: float,
            energy_lowest: float,
            dos_values: NDArray) -> None:
        """Given the band structure and information about energy references,
        determines the energy references separately for all spin channels.
        """
        dos.energy_fermi = energy_fermi
        eref = energy_highest if energy_fermi is None else energy_fermi
        # No reference available
        if eref is None:
            self.logger.info("could not resolve energy references for dos")
            return

        # Create channel information for each spin channel and populate with
        # initial values.
        dos_total = dos.total
        n_channels = len(dos_total)
        for i_channel in range(n_channels):
            info = dos.band_gap[i_channel] if len(dos.band_gap) > i_channel else dos.m_create(BandGap)
            info.index = i_channel
            if energy_highest is not None:
                info.energy_highest_occupied = energy_highest
            if energy_lowest is not None:
                info.energy_lowest_unoccupied = energy_lowest

        # Use a reference energy (fermi or highest occupied) to determine the
        # energy references from the DOS (discretization will affect the exact
        # location).
        energy_threshold = config.normalize.band_structure_energy_tolerance
        value_threshold = 1e-8  # The DOS value that is considered to be zero
        dos_energies = dos.energies

        for i_channel in range(n_channels):
            dos_channel = dos_values[i_channel]
            info = dos.band_gap[i_channel]
            fermi_idx = (np.abs(dos.energies - eref)).argmin()

            # First check that the closest dos energy to energy reference
            # is not too far away. If it is very far away, the
            # normalization may be very inaccurate and we do not report it.
            fermi_energy_closest = dos_energies[fermi_idx]
            distance = np.abs(fermi_energy_closest - eref)
            if distance.magnitude <= energy_threshold:

                # See if there are zero values close below the energy reference.
                idx = fermi_idx
                idx_descend = fermi_idx
                while True:
                    try:
                        value = dos_channel[idx]
                        energy_distance = np.abs(eref - dos_energies[idx])
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
                        value = dos_channel[idx]
                        energy_distance = np.abs(eref - dos_energies[idx])
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
                        value = dos_channel[idx]
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
                        value = dos_channel[idx]
                    except IndexError:
                        break
                    if value > value_threshold:
                        idx = idx if idx == idx_ascend else idx - 1
                        info.energy_lowest_unoccupied = dos_energies[idx]
                        break
                    idx += 1
