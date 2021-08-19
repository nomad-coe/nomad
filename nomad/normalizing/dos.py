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

from nomad import config
from nomad_dos_fingerprints import DOSFingerprint
from nomad.datamodel.metainfo.run.calculation import (
    Dos, DosFingerprint, ElectronicStructureInfo)
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
        section_sccs = self.section_run.calculation
        if section_sccs is None:
            return

        for scc in section_sccs:
            for section_dos in [scc.dos_electronic, scc.dos_phonon]:
                if section_dos is None:
                    continue

                energy_fermi = scc.energy.fermi
                energy_highest = scc.energy.highest_occupied
                energy_lowest = scc.energy.lowest_unoccupied
                for dos in section_dos:
                    # perform normalization only for total dos
                    if dos.total is None:
                        continue

                    # Normalize DOS values to be 1/J/atom/m^3
                    system = scc.system_ref
                    if not system or system[0].value.atoms is None:
                        self.logger.error('referenced system for dos calculation could not be found')
                        continue
                    atom_positions = system[0].value.atoms.positions
                    lattice_vectors = system[0].value.atoms.lattice_vectors
                    if atom_positions is None:
                        self.logger.error('required quantity atom_positions is not available')
                        continue
                    if lattice_vectors is None:
                        self.logger.error('required quantity lattice_vectors is not available')
                        continue
                    number_of_atoms = np.shape(atom_positions)[0]
                    unit_cell_volume = get_volume(lattice_vectors.magnitude)
                    for dos_total in dos.total:
                        dos_total.normalization_factor = number_of_atoms * unit_cell_volume

                    # Add energy references
                    self.add_energy_references(dos, energy_fermi, energy_highest, energy_lowest)

                    if dos.info[-1].energy_highest_occupied is None:
                        continue

                    dos_energies_normalized = dos.energies - dos.info[-1].energy_highest_occupied
                    dos.dos_energies_normalized = dos_energies_normalized
                    dos_values_normalized = [
                        dos_total.value.magnitude / dos_total.normalization_factor for dos_total in dos.total]

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
                        sec_dos_fingerprint = dos.m_create(DosFingerprint)
                        sec_dos_fingerprint.bins = dos_fingerprint.bins
                        sec_dos_fingerprint.indices = dos_fingerprint.indices
                        sec_dos_fingerprint.stepsize = dos_fingerprint.stepsize
                        sec_dos_fingerprint.grid_id = dos_fingerprint.grid_id
                        sec_dos_fingerprint.filling_factor = dos_fingerprint.filling_factor

    def add_energy_references(
            self,
            dos: Dos,
            energy_fermi: float,
            energy_highest: float,
            energy_lowest: float) -> None:
        """Given the band structure and information about energy references,
        determines the energy references separately for all spin channels.
        """
        eref = energy_highest if energy_fermi is None else energy_fermi
        # No reference available
        if eref is None:
            self.logger.info("could not resolve energy references for dos")
            return

        # Use a reference energy (fermi or highest occupied) to determine the
        # energy references from the DOS (discretization will affect the exact
        # location).
        energy_threshold = config.normalize.band_structure_energy_tolerance
        value_threshold = 1e-8  # The DOS value that is considered to be zero
        # Create energy reference sections for each spin channel, add fermi
        # energy if present
        info = dos.info[0] if dos.info else dos.m_create(ElectronicStructureInfo)
        if energy_highest is not None:
            info.energy_highest_occupied = energy_highest
        if energy_lowest is not None:
            info.energy_lowest_unoccupied = energy_lowest
        if energy_fermi is not None:
            info.energy_fermi = energy_fermi

        # energy references is only relevant for the total
        dos_values_normalized = np.sum([
            dos_total.value.magnitude / dos_total.normalization_factor for dos_total in dos.total], axis=0)

        fermi_idx = (np.abs(dos.energies - eref)).argmin()

        # First check that the closest dos energy to energy reference
        # is not too far away. If it is very far away, the
        # normalization may be very inaccurate and we do not report it.
        fermi_energy_closest = dos.energies[fermi_idx]
        distance = np.abs(fermi_energy_closest - eref)
        if distance.magnitude <= energy_threshold:

            # See if there are zero values close below the energy reference.
            idx = fermi_idx
            idx_descend = fermi_idx
            while True:
                try:
                    value = dos_values_normalized[idx]
                    energy_distance = np.abs(eref - dos.energies[idx])
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
                    value = dos_values_normalized[idx]
                    energy_distance = np.abs(eref - dos.energies[idx])
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
                return

            # Look for highest occupied energy below the descend index
            idx = idx_descend
            while True:
                try:
                    value = dos_values_normalized[idx]
                except IndexError:
                    break
                if value > value_threshold:
                    idx = idx if idx == idx_descend else idx + 1
                    info.energy_highest_occupied = dos.energies[idx]
                    break
                idx -= 1

            # Look for lowest unoccupied energy above idx_ascend
            idx = idx_ascend
            while True:
                try:
                    value = dos_values_normalized[idx]
                except IndexError:
                    break
                if value > value_threshold:
                    idx = idx if idx == idx_ascend else idx - 1
                    info.energy_lowest_unoccupied = dos.energies[idx]
                    break
                idx += 1
