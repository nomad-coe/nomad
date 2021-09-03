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

from nomad.datamodel.encyclopedia import (
    Calculation,
    Properties,
    Energies,
)
from nomad.metainfo import Section
from nomad.normalizing.encyclopedia.context import Context


class PropertiesNormalizer():
    """A base class that is used for processing calculated quantities that
    should be extracted to Encyclopedia.
    """
    def __init__(self, entry_archive, logger):
        self.entry_archive = entry_archive
        self.logger = logger

    def electronic_band_structure(self, properties: Properties, calc_type: str, material_type: str, context: Context, sec_system: Section) -> None:
        """Tries to resolve a reference to a representative electronic band
        structure. Will loop through the available band in reversed order until
        a valid band is found.
        """
        try:
            representative_scc = context.representative_scc
            bands = representative_scc.band_structure_electronic
            if bands is None:
                return

            representative_band = None
            valid = True
            for band in reversed(bands):
                if not band.segment:
                    valid = False
                    break
                for segment in band.segment:
                    energies = segment.energies
                    k_points = segment.kpoints
                    labels = segment.endpoints_labels
                    if energies is None or k_points is None or labels is None:
                        valid = False
                        break
                if valid:
                    representative_band = band
                    break

        except Exception:
            return
        if representative_band is not None:
            properties.electronic_band_structure = representative_band.m_path()

            # Add band gap information to metadata if present. The channel with
            # smallest band gap index is chosen as a representative one.
            channel_info = properties.electronic_band_structure.channel_info
            if channel_info is not None and len(channel_info) > 0:
                min_gap_index = 0
                min_gap = float("Inf")
                for i, info in enumerate(channel_info):
                    band_gap = info.band_gap.magnitude
                    if band_gap is not None and band_gap < min_gap:
                        min_gap_index = i
                        min_gap = band_gap
                representative_channel = channel_info[min_gap_index]
                bg_value = representative_channel.band_gap
                if bg_value is not None:
                    properties.band_gap = bg_value
                    properties.band_gap_direct = representative_channel.band_gap_type == "direct"

    def electronic_dos(self, properties: Properties, context: Context) -> None:
        """Tries to resolve a reference to a representative electonic density
        of states. Will not store any data if it cannot be resolved
        unambiguously.
        """
        try:
            representative_scc = context.representative_scc
            doses = representative_scc.dos_electronic
            if doses is None:
                return

            representative_dos = None
            for dos in reversed(doses):
                energies = dos.energies
                values = dos.total
                if energies is not None and len(values) > 0:
                    representative_dos = dos
                    break

        except Exception:
            return
        if representative_dos is not None:
            properties.electronic_dos = representative_dos.m_path()

    def elastic_constants_matrix(self) -> None:
        pass

    def elastic_deformation_energies(self) -> None:
        pass

    def elastic_fitting_parameters(self) -> None:
        pass

    def elastic_moduli(self) -> None:
        pass

    def elastic_properties(self) -> None:
        pass

    def fermi_surface(self) -> None:
        pass

    def has_fermi_surface(self) -> None:
        pass

    def thermodynamical_properties(self, properties: Properties) -> None:
        """Tries to resolve a reference to a representative set of
        thermodynamical properties. Will not store any data if it cannot be
        resolved unambiguously.
        """
        try:
            resolved_section = None
            for workflow in reversed(self.entry_archive.workflow):
                if workflow.thermodynamics is not None:
                    resolved_section = workflow.thermodynamics
            if resolved_section is None:
                self.logger("Could not unambiguously select data to display for specific heat.")
                return
        except Exception:
            return
        if resolved_section is not None:
            properties.thermodynamical_properties = resolved_section.m_path()

    def phonon_band_structure(self, properties: Properties, context: Context) -> None:
        """Tries to resolve a reference to a representative phonon band
        structure. Will not store any data if it cannot be resolved
        unambiguously.
        """
        try:
            representative_scc = context.representative_scc
            bands = representative_scc.band_structure_phonon
            if bands is None:
                return

            representative_phonon_band = None
            for band in reversed(bands):
                valid = True
                for segment in band.segment:
                    energies = segment.energies
                    k_points = segment.kpoints
                    labels = segment.endpoints_labels
                    if energies is None or k_points is None or labels is None or "?" in labels:
                        valid = False
                if valid:
                    representative_phonon_band = band
                    break

        except Exception:
            return
        if representative_phonon_band is not None:
            properties.phonon_band_structure = representative_phonon_band.m_path()

    def phonon_dos(self, properties: Properties, context: Context) -> None:
        """Tries to resolve a reference to a representative phonon density of
        states. Will not store any data if it cannot be resolved unambiguously.
        """
        try:
            representative_scc = context.representative_scc
            doses = representative_scc.dos_phonon
            if doses is None:
                return

            representative_phonon_dos = None
            for dos in reversed(doses):
                energies = dos.energies
                values = dos.total
                if energies is not None and len(values) > 0:
                    if representative_phonon_dos is None:
                        representative_phonon_dos = dos
                        break

        except Exception:
            self.logger("Could not unambiguously select data to display for phonon density of states.")
            return
        if representative_phonon_dos is not None:
            properties.phonon_dos = representative_phonon_dos.m_path()

    def energies(self, properties: Properties, n_atoms: int, representative_scc: Section) -> None:
        if representative_scc is not None:
            energies = Energies()
            energy_found = False
            energy = representative_scc.energy
            if energy is not None:
                if energy.total is not None:
                    energies.energy_total = energy.total.value
                if energy.total_t0 is not None:
                    energies.energy_total_T0 = energy.total_t0.value
                if energy.free is not None:
                    energies.energy_free = energy.free.value
                energy_found = True
            if energy_found:
                properties.m_add_sub_section(Properties.energies, energies)

    def normalize(self, context: Context) -> None:
        # There needs to be a valid SCC in order to extract any properties
        representative_scc = context.representative_scc
        if representative_scc is None:
            return

        # Fetch resources
        sec_enc = self.entry_archive.metadata.encyclopedia
        properties = sec_enc.properties
        calc_type = context.calc_type
        material_type = context.material_type
        sec_system = context.representative_system

        # Save metainfo
        self.electronic_band_structure(properties, calc_type, material_type, context, sec_system)
        self.electronic_dos(properties, context)
        if sec_system.atoms is not None:
            n_atoms = len(sec_system.atoms.labels)
            self.energies(properties, n_atoms, representative_scc)

        # Phonon calculations have a specific set of properties to extract
        if context.calc_type == Calculation.calculation_type.type.phonon_calculation:
            self.thermodynamical_properties(properties)
            self.phonon_band_structure(properties, context)
            self.phonon_dos(properties, context)
