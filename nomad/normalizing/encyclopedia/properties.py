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
            bands = representative_scc.section_k_band
            if bands is None:
                return

            representative_band = None
            for band in reversed(bands):
                kind = band.band_structure_kind
                if kind != "vibrational":
                    valid = True
                    if not band.section_k_band_segment:
                        valid = False
                        break
                    for segment in band.section_k_band_segment:
                        energies = segment.band_energies
                        k_points = segment.band_k_points
                        labels = segment.band_segm_labels
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
            band_gaps = properties.electronic_band_structure.section_band_gap
            if band_gaps is not None and len(band_gaps) > 0:
                min_gap_index = 0
                min_gap = float("Inf")
                for i, gap in enumerate(band_gaps):
                    value = gap.value
                    if value < min_gap:
                        min_gap_index = i
                        min_gap = value
                representative_gap = band_gaps[min_gap_index]
                bg_value = representative_gap.value
                if bg_value is not None:
                    properties.band_gap = representative_gap.value
                    properties.band_gap_direct = representative_gap.type == "direct"

    def electronic_dos(self, properties: Properties, context: Context) -> None:
        """Tries to resolve a reference to a representative electonic density
        of states. Will not store any data if it cannot be resolved
        unambiguously.
        """
        try:
            representative_scc = context.representative_scc
            doses = representative_scc.section_dos
            if doses is None:
                return

            representative_dos = None
            for dos in reversed(doses):
                kind = dos.dos_kind
                energies = dos.dos_energies
                values = dos.dos_values
                if kind != "vibrational" and energies is not None and values is not None:
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
            frame_sequences = self.entry_archive.section_run[0].section_frame_sequence
            for frame_sequence in reversed(frame_sequences):
                thermodynamical_props = frame_sequence.section_thermodynamical_properties
                for thermodynamical_prop in thermodynamical_props:
                    if resolved_section is None:
                        resolved_section = thermodynamical_prop
                    else:
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
            bands = representative_scc.section_k_band
            if bands is None:
                return

            representative_phonon_band = None
            for band in reversed(bands):
                kind = band.band_structure_kind
                if kind == "vibrational":
                    valid = True
                    for segment in band.section_k_band_segment:
                        energies = segment.band_energies
                        k_points = segment.band_k_points
                        labels = segment.band_segm_labels
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
            doses = representative_scc.section_dos
            if doses is None:
                return

            representative_phonon_dos = None
            for dos in reversed(doses):
                kind = dos.dos_kind
                energies = dos.dos_energies
                values = dos.dos_values
                if kind == "vibrational" and energies is not None and values is not None:
                    if representative_phonon_dos is None:
                        representative_phonon_dos = dos
                    else:
                        self.logger("Could not unambiguously select data to display for phonon density of states.")
                        return

        except Exception:
            return
        if representative_phonon_dos is not None:
            properties.phonon_dos = representative_phonon_dos.m_path()

    def energies(self, properties: Properties, n_atoms: int, representative_scc: Section) -> None:
        if representative_scc is not None:
            energies = Energies()
            energy_found = False
            for energy_name in ["energy_total", "energy_total_T0", "energy_free"]:
                energy_value = getattr(representative_scc, energy_name)
                if energy_value is not None:
                    energy_found = True

                    # The energies are normalized to be per atom
                    setattr(energies, energy_name, energy_value.magnitude / n_atoms)
            if energy_found:
                properties.m_add_sub_section(Properties.energies, energies)

    def normalize(self, context: Context) -> None:
        # There needs to be a valid SCC in order to extract any properties
        representative_scc = context.representative_scc
        if representative_scc is None:
            return

        # Fetch resources
        sec_enc = self.entry_archive.section_metadata.encyclopedia
        properties = sec_enc.properties
        calc_type = context.calc_type
        material_type = context.material_type
        sec_system = context.representative_system
        n_atoms = len(sec_system.atom_labels)

        # Save metainfo
        self.electronic_band_structure(properties, calc_type, material_type, context, sec_system)
        self.electronic_dos(properties, context)
        self.energies(properties, n_atoms, representative_scc)

        # Phonon calculations have a specific set of properties to extract
        if context.calc_type == Calculation.calculation_type.type.phonon_calculation:
            self.thermodynamical_properties(properties)
            self.phonon_band_structure(properties, context)
            self.phonon_dos(properties, context)
