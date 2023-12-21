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
import re
from collections import defaultdict
from typing import Optional

from nomad import config
from nomad_dos_fingerprints import DOSFingerprint  # pylint: disable=import-error
from nomad.datamodel.metainfo.simulation.calculation import (
    BandGap,
    BandGapDeprecated,
    Calculation,
    Dos,
    DosValues,
    DosFingerprint,
    ElectronicStructureProvenance,
)

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
            if dos_electronic:
                energy_fermi = calc.energy.fermi if calc.energy is not None else None
                energy_highest = (
                    calc.energy.highest_occupied if calc.energy is not None else None
                )
                energy_lowest = (
                    calc.energy.lowest_unoccupied if calc.energy is not None else None
                )
                n_spin_channels = len(dos_electronic)
            for index, dos in enumerate(dos_electronic):
                orbital_projected = dos.orbital_projected
                atom_projected = dos.atom_projected
                species_projected = dos.species_projected
                total_dos = dos.total
                # If the atom_projected does not exist, we resolve the orbital_projected and
                # store the summation of each orbital for each atom in atom_projected
                if not atom_projected:
                    # `atom_data` is a dictionary of unknown keys (resolved inside the loop
                    # with `atom`) and whose values are list where the orbital-projected DOS
                    # has to be appended.
                    atom_data = defaultdict(list)
                    for orbital_pdos in orbital_projected:
                        atom = f'{orbital_pdos.atom_label}{orbital_pdos.atom_index if orbital_pdos.atom_index else ""}'
                        atom_data[atom].append(orbital_pdos.value.magnitude)
                    for atom, data in atom_data.items():
                        atom_value = np.sum(data, axis=0)
                        sec_dos_atom = dos.m_create(DosValues, Dos.atom_projected)
                        atom_label = re.sub(r'\d', '', atom)
                        atom_index = re.sub(r'[a-zA-Z]', '', atom)
                        sec_dos_atom.atom_label = atom_label
                        sec_dos_atom.atom_index = atom_index if atom_index else None
                        sec_dos_atom.value = atom_value
                # If the species_projected does not exist, we resolve the atom_projected and
                # store the summation of each atom for each species in species_projected
                if not species_projected:
                    # `species_data` is a dictionary of unknown keys (resolved inside the loop
                    # with `species`) and whose values are list where the atom-projected DOS
                    # has to be appended.
                    species_data = defaultdict(list)
                    for atom_pdos in atom_projected:
                        species = atom_pdos.atom_label
                        species_data[species].append(atom_pdos.value.magnitude)
                    for species, data in species_data.items():
                        species_value = np.sum(data, axis=0)
                        sec_dos_species = dos.m_create(DosValues, Dos.species_projected)
                        sec_dos_species.atom_label = species
                        sec_dos_species.value = species_value
                # If the total DOS does not exist but the species_projected dos, we store the
                # species_projected summation into the total DOS and show a warning about
                # our procedure
                if not total_dos and species_projected:
                    self.logger.info(
                        'Total DOS not found, but Projected DOS was found. We '
                        'will sum up contributions leading to the final DOS.'
                    )
                    total_data = []
                    for species_pdos in species_projected:
                        total_data.append(species_pdos.value.magnitude)
                    total_value = np.sum(total_data, axis=0)
                    sec_dos_total = dos.m_create(DosValues, Dos.total)
                    sec_dos_total.value = total_value
                # TODO add provenance for the Dos. This will require that the spin_polarized
                # and n_spin_channels are moved to `run.method`, and then add the provenance
                # here for each Dos as a reference to the methodology.

                # Add energy references
                dos_values = [dos_total.value.magnitude for dos_total in dos.total]
                self.add_energy_references(
                    calc,
                    index,
                    dos,
                    energy_fermi,
                    energy_highest,
                    energy_lowest,
                    dos_values,
                )

                # Calculate the DOS fingerprint for successfully normalized DOS for finding
                # the energy_highest_occupied
                normalization_reference = None
                for info in calc.band_gap:
                    energy_highest_occupied = info.energy_highest_occupied
                    if energy_highest_occupied is not None:
                        if normalization_reference is None:
                            normalization_reference = energy_highest_occupied
                        else:
                            normalization_reference = max(
                                normalization_reference, energy_highest_occupied
                            )
                if normalization_reference is not None:
                    dos_energies_normalized = dos.energies - normalization_reference

                    try:
                        dos_fingerprint = DOSFingerprint().calculate(
                            dos_energies_normalized.magnitude, dos_values
                        )
                    except Exception as e:
                        self.logger.error(
                            'could not generate dos fingerprint', exc_info=e
                        )
                    else:
                        sec_dos_fingerprint = dos.m_create(DosFingerprint)
                        sec_dos_fingerprint.bins = dos_fingerprint.bins
                        sec_dos_fingerprint.indices = dos_fingerprint.indices
                        sec_dos_fingerprint.stepsize = dos_fingerprint.stepsize
                        sec_dos_fingerprint.grid_id = dos_fingerprint.grid_id
                        sec_dos_fingerprint.filling_factor = (
                            dos_fingerprint.filling_factor
                        )

                # Add normalization factor
                set_normalization_factor = self.add_electronic_normalization_factor(
                    calc, dos, n_spin_channels
                )
                if not set_normalization_factor:
                    continue

            # Normalize phonon DOS
            dos_phonons = calc.dos_phonon
            for dos_phonon in dos_phonons:
                self.add_phononic_normalization_factor(calc, dos_phonon)

    def add_electronic_normalization_factor(
        self, calc: Calculation, dos: Dos, n_spin_channels: int
    ) -> Optional[np.float64]:
        """Returns a factor that returns a size intensive electronic DOS.
        The values are divided by integral(DOS, lowest state, Fermi energy), or likewise sum(<atomic numbers>)."""
        if not calc.system_ref:
            self.logger.warning(
                'Could not resolve the system reference from calculation.system_ref, '
                'thus electronic normalization factor not reported.'
            )
            return None
        atoms = calc.system_ref.atoms
        if not len(dos.total):
            self.logger.warning(
                'Could not resolve total DOS from calculation.dos.total, '
                'thus DOS electronic normalization factor not reported.'
            )
            return None
        elif not len(atoms.species):
            self.logger.warning(
                'Could not resolve atoms information from calculation.system_ref.atoms, '
                'thus DOS electronic normalization factor not reported.'
            )
            return None
        else:
            normalization_factor = 1 / (n_spin_channels * sum(atoms.species))
            for dos_total in dos.total:
                dos_total.normalization_factor = normalization_factor
            return normalization_factor

    def add_phononic_normalization_factor(
        self, calc: Calculation, dos: Dos
    ) -> Optional[float]:
        """Returns a factor that returns a size intensive phononic DOS.
        The values are divided by integral(DOS, 0, infinity), or likewise <no. degrees of freedom>"""
        atoms = calc.system_ref.atoms
        if not len(dos.total):
            self.logger.warning(
                'Could not resolve total DOS from calculation.dos.total, '
                'thus DOS phononic normalization factor not reported.'
            )
            return None
        elif not len(atoms.species):
            self.logger.warning(
                'Could not resolve atoms information from calculation.system_ref.atoms, '
                'thus DOS phononic normalization factor not reported.'
            )
            return None
        else:
            normalization_factor = 1 / (3 * len(atoms.species))
            for dos_total in dos.total:
                dos_total.normalization_factor = normalization_factor
            return normalization_factor

    def add_energy_references(
        self,
        calc: Calculation,
        i_channel: int,
        dos: Dos,
        energy_fermi: float,
        energy_highest: float,
        energy_lowest: float,
        dos_values: list,
    ) -> None:
        """Given the DOS and information about energy references,
        determines the energy references separately for all spin channels.
        """
        dos.energy_fermi = energy_fermi
        eref = energy_highest if energy_fermi is None else energy_fermi
        # No reference available
        if eref is None:
            self.logger.info('could not resolve energy references for dos')
            return

        # Use a reference energy (fermi or highest occupied) to determine the
        # energy references from the DOS (discretization will affect the exact
        # location).
        energy_threshold = config.normalize.band_structure_energy_tolerance
        value_threshold = 1e-8  # The DOS value that is considered to be zero
        dos_energies = dos.energies

        # Create channel information for each spin channel and populate with
        # initial values.
        info = BandGapDeprecated()
        info.index = i_channel
        if energy_highest is not None:
            info.energy_highest_occupied = energy_highest
        if energy_lowest is not None:
            info.energy_lowest_unoccupied = energy_lowest

        # First check that the closest dos energy to energy reference
        # is not too far away. If it is very far away, the
        # normalization may be very inaccurate and we do not report it.
        dos_channel = dos_values[0]
        fermi_idx = (np.abs(dos_energies - eref)).argmin()
        fermi_energy_closest = dos_energies[fermi_idx]
        distance = np.abs(fermi_energy_closest - eref)
        single_peak_fermi = False
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
                single_peak_fermi = True

            if not single_peak_fermi:
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

        # save energy_ref for dos
        dos.energy_ref = (
            info.energy_highest_occupied
            if info.energy_highest_occupied
            else energy_fermi
        )

        # save band gap value
        if info.energy_lowest_unoccupied and info.energy_highest_occupied:
            gap_value = info.energy_lowest_unoccupied - info.energy_highest_occupied
            gap_value = gap_value if gap_value > 0.0 else 0.0
            info.value = gap_value

            if info.value is not None:
                proper_info = BandGap().m_from_dict(info.m_to_dict())
                proper_info.provenance = ElectronicStructureProvenance(
                    dos=dos.total[0], label='dos'
                )
                calc.m_add_sub_section(Calculation.band_gap, proper_info)
                dos.m_add_sub_section(Dos.band_gap, info)
