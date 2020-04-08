# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the 'License');
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an'AS IS' BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from typing import List
import json
import numpy as np
from ase import Atoms

from nomad.metainfo.encyclopedia import (
    Calculation,
    Properties,
    Material,
    ElectronicBandStructure,
    BandGap,
)
from nomad.parsing.legacy import Backend
from nomad.metainfo import Section
from nomad.normalizing.encyclopedia.context import Context
from nomad import atomutils
from nomad import config

J_to_Ry = 4.587425e+17


class PropertiesNormalizer():
    """A base class that is used for processing calculated quantities that
    should be extracted to Encyclopedia.
    """
    def __init__(self, backend: Backend, logger):
        self.backend = backend
        self.logger = logger

    def get_reciprocal_cell(self, band_structure: ElectronicBandStructure, prim_atoms: Atoms, orig_atoms: Atoms):
        """A reciprocal cell for this calculation. If the original unit cell is
        not a primitive one, then we will use the one given by spglib.

        If the used code is calculating it's own primitive cell, then the
        reciprocal cell used in e.g. band structure calculation might differ
        from the one given by spglib. To overcome this we would need to get the
        exact reciprocal cell used by the calculation or then test for
        different primitive cells
        (https://atztogo.github.io/spglib/definition.html#transformation-to-the-primitive-cell)
        whether the k-point path stays inside the first Brillouin zone.

        Args:
            prim_atoms: The primitive system as detected by MatID.
            orig_atoms: The original simulation cell.

        Returns:
            Reciprocal cell as 3x3 numpy array. Units are in SI (1/m).
        """
        primitive_cell = prim_atoms.get_cell()
        source_cell = orig_atoms.get_cell()

        volume_primitive = primitive_cell.volume
        volume_source = source_cell.volume
        volume_diff = abs(volume_primitive - volume_source)

        if volume_diff > (0.001)**3:
            recip_cell = primitive_cell.reciprocal() * 1e10
        else:
            recip_cell = source_cell.reciprocal() * 1e10

        band_structure.reciprocal_cell = recip_cell

    def get_band_gaps(self, band_structure: ElectronicBandStructure, energies: np.array, path: np.array) -> None:
        """Given the band structure and fermi level, calculates the band gap
        for spin channels and also reports the total band gap as the minum gap
        found.
        """
        # Handle spin channels separately to find gaps for spin up and down
        reciprocal_cell = band_structure.reciprocal_cell.magnitude
        fermi_level = band_structure.fermi_level
        n_channels = energies.shape[0]

        gaps: List[BandGap] = [None, None]
        for channel in range(n_channels):
            channel_energies = energies[channel, :, :]
            lower_defined = False
            upper_defined = False
            num_bands = channel_energies.shape[0]
            band_indices = np.arange(num_bands)
            band_minima_idx = channel_energies.argmin(axis=1)
            band_maxima_idx = channel_energies.argmax(axis=1)
            band_minima = channel_energies[band_indices, band_minima_idx]
            band_maxima = channel_energies[band_indices, band_maxima_idx]

            # Add a tolerance to minima and maxima
            band_minima_tol = band_minima + config.normalize.fermi_level_precision
            band_maxima_tol = band_maxima - config.normalize.fermi_level_precision

            for band_idx in range(num_bands):
                band_min = band_minima[band_idx]
                band_max = band_maxima[band_idx]
                band_min_tol = band_minima_tol[band_idx]
                band_max_tol = band_maxima_tol[band_idx]

                # If any of the bands band crosses the Fermi level, there is no
                # band gap
                if band_min_tol <= fermi_level and band_max_tol >= fermi_level:
                    break
                # Whole band below Fermi level, save the current highest
                # occupied band point
                elif band_min_tol <= fermi_level and band_max_tol <= fermi_level:
                    gap_lower_energy = band_max
                    gap_lower_idx = band_maxima_idx[band_idx]
                    lower_defined = True
                # Whole band above Fermi level, save the current lowest
                # unoccupied band point
                elif band_min_tol >= fermi_level:
                    gap_upper_energy = band_min
                    gap_upper_idx = band_minima_idx[band_idx]
                    upper_defined = True
                    break

            # If a highest point of the valence band and a lowest point of the
            # conduction band are found, and the difference between them is
            # positive, save the information location and value.
            if lower_defined and upper_defined and gap_upper_energy - gap_lower_energy >= 0:

                # See if the gap is direct or indirect by comparing the k-point
                # locations with some tolerance
                k_point_lower = path[gap_lower_idx]
                k_point_upper = path[gap_upper_idx]
                k_point_distance = self.get_k_space_distance(reciprocal_cell, k_point_lower, k_point_upper)
                is_direct_gap = k_point_distance <= config.normalize.k_space_precision

                gap = BandGap()
                gap.type = "direct" if is_direct_gap else "indirect"
                gap.value = float(gap_upper_energy - gap_lower_energy)
                gap.conduction_band_min_k_point = k_point_upper
                gap.conduction_band_min_energy = float(gap_upper_energy)
                gap.valence_band_max_k_point = k_point_lower
                gap.valence_band_max_energy = float(gap_lower_energy)
                gaps[channel] = gap

        # For unpolarized calculations we simply report the gap if it is found.
        if n_channels == 1:
            if gaps[0] is not None:
                band_structure.m_add_sub_section(ElectronicBandStructure.band_gap, gaps[0])
        # For polarized calculations we report the gap separately for both
        # channels. Also we report the smaller gap as the the total gap for the
        # calculation.
        elif n_channels == 2:
            if gaps[0] is not None:
                band_structure.m_add_sub_section(ElectronicBandStructure.band_gap_spin_up, gaps[0])
            if gaps[1] is not None:
                band_structure.m_add_sub_section(ElectronicBandStructure.band_gap_spin_down, gaps[1])
            if gaps[0] is not None and gaps[1] is not None:
                if gaps[0].value <= gaps[1].value:
                    band_structure.m_add_sub_section(ElectronicBandStructure.band_gap, gaps[0])
                else:
                    band_structure.m_add_sub_section(ElectronicBandStructure.band_gap, gaps[1])
            else:
                if gaps[0] is not None:
                    band_structure.m_add_sub_section(ElectronicBandStructure.band_gap, gaps[0])
                elif gaps[1] is not None:
                    band_structure.m_add_sub_section(ElectronicBandStructure.band_gap, gaps[1])

    def get_k_space_distance(self, reciprocal_cell: np.array, point1: np.array, point2: np.array) -> float:
        """Used to calculate the Euclidean distance of two points in k-space,
        given relative positions in the reciprocal cell.

        Args:
            reciprocal_cell: Reciprocal cell.
            point1: The first position in k-space.
            point2: The second position in k-space.

        Returns:
            float: Euclidian distance of the two points in k-space in SI units.
        """
        k_point_displacement = np.dot(reciprocal_cell, point1 - point2)
        k_point_distance = np.linalg.norm(k_point_displacement)

        return k_point_distance

    def get_brillouin_zone(self, band_structure: ElectronicBandStructure) -> None:
        """Returns a dictionary containing the information needed to display
        the Brillouin zone for this material. This functionality could be put
        into the GUI directly, with the Brillouin zone construction performed
        from the reciprocal cell.

        The Brillouin Zone is a Wigner-Seitz cell, and is thus uniquely
        defined. It's shape does not depend on the used primitive cell.
        """
        recip_cell = band_structure.reciprocal_cell.magnitude
        brillouin_zone = atomutils.get_brillouin_zone(recip_cell)
        bz_json = json.dumps(brillouin_zone)
        band_structure.brillouin_zone = bz_json

    def band_structure(self, properties: Properties, calc_type: str, material_type: str, context: Context, sec_system: Section) -> None:
        """Band structure data following arbitrary path.

        Currently this function is only taking into account the normalized band
        structures, and if they have a k-point that is labeled as '?', the
        whole band strcture will be ignored.
        """
        # Band structure data is extracted only from bulk calculations. This is
        # because for now we need the reciprocal cell of the primitive cell
        # that is only available from the symmetry analysis. Once the
        # reciprocal cell is directly reported with the band structure this
        # restriction can go away.
        if calc_type != Calculation.calculation_type.type.single_point or material_type != Material.material_type.type.bulk:
            return

        representative_scc = context.representative_scc
        orig_atoms = sec_system.m_cache["representative_atoms"]
        symmetry_analyzer = sec_system.section_symmetry[0].m_cache["symmetry_analyzer"]
        prim_atoms = symmetry_analyzer.get_primitive_system()

        # Try to find an SCC with band structure data, give priority to
        # normalized data
        for src_name in ["section_k_band_normalized", "section_k_band"]:
            bands = getattr(representative_scc, src_name)
            if bands is None:
                continue
            norm = "_normalized" if src_name == "section_k_band_normalized" else ""

            # Loop over bands
            for band_data in bands:
                band_structure = ElectronicBandStructure()
                band_structure.scc_index = int(context.representative_scc_idx)
                kpoints = []
                energies = []
                segments = band_data['section_k_band_segment' + norm]
                if not segments:
                    return

                # Loop over segments
                for segment_src in segments:
                    try:
                        seg_k_points = segment_src["band_k_points" + norm]
                        seg_energies = segment_src["band_energies" + norm]
                        seg_labels = segment_src['band_segm_labels' + norm]
                    except Exception:
                        return
                    if seg_k_points is None or seg_energies is None or seg_labels is None:
                        return
                    else:
                        seg_energies = seg_energies.magnitude
                    if "?" in seg_labels:
                        return

                    seg_energies = np.swapaxes(seg_energies, 1, 2)
                    kpoints.append(seg_k_points)
                    energies.append(seg_energies)

                # Continue to calculate band gaps and other properties.
                kpoints = np.concatenate(kpoints, axis=0)
                energies = np.concatenate(energies, axis=2)
                self.get_reciprocal_cell(band_structure, prim_atoms, orig_atoms)
                self.get_brillouin_zone(band_structure)

                # If we are using a normalized band structure (or the Fermi level
                # is defined somehow else), we can add band gap information
                fermi_level = None
                if src_name == "section_k_band_normalized":
                    fermi_level = 0.0
                if fermi_level is not None:
                    band_structure.fermi_level = fermi_level
                    self.get_band_gaps(band_structure, energies, kpoints)

                properties.m_add_sub_section(Properties.electronic_band_structure, band_structure)

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

    def has_thermal_properties(self) -> None:
        pass

    def phonon_dispersion(self) -> None:
        pass

    def phonon_dos(self) -> None:
        pass

    def specific_heat_cv(self) -> None:
        pass

    def helmholtz_free_energy(self) -> None:
        pass

    def energies(self, properties: Properties, gcd: int, representative_scc: Section) -> None:
        energy_dict = {}
        if representative_scc is not None:
            energies_entries = {
                "energy_total": "Total E",
                "energy_total_T0": "Total E projected to T=0",
                "energy_free": "Free E",
            }
            for energy_name, label in energies_entries.items():
                result = getattr(representative_scc, energy_name)
                if result is not None:
                    energy_dict[label] = result.magnitude / gcd

            if len(energy_dict) == 0:
                energy_dict = None
        energies = json.dumps(energy_dict)
        properties.energies = energies

    def normalize(self, context: Context) -> None:
        # There needs to be a valid SCC in order to extract any properties
        representative_scc = context.representative_scc
        if representative_scc is None:
            return

        # Fetch resources
        sec_enc = self.backend.entry_archive.section_encyclopedia
        properties = sec_enc.properties
        calc_type = context.calc_type
        material_type = context.material_type
        sec_system = context.representative_system
        gcd = context.greatest_common_divisor

        # Save metainfo
        self.band_structure(properties, calc_type, material_type, context, sec_system)
        self.energies(properties, gcd, representative_scc)
