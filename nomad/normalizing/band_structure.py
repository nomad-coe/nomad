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

import json
import numpy as np
from typing import List

from nomad.datamodel.metainfo.public import section_k_band, section_band_gap, section_system
from nomad.normalizing.normalizer import Normalizer
from nomad import config, atomutils


class BandStructureNormalizer(Normalizer):
    """Normalizer with the following responsibilities:

      - Calculates band gap(s) if present (section_band_gap, section_band_gap_spin_up, section_band_gap_spin_down).
      - TODO: Creates labels for special points within the band path (band_path_labels)
      - TODO: Determines if the path is a standard one or not (is_standard)
    """
    def normalize(self, logger=None) -> None:
        # Setup logger
        if logger is not None:
            self.logger = logger.bind(normalizer=self.__class__.__name__)

        # Do nothing if section_run is not present
        if self.section_run is None:
            return

        # Loop through the bands
        for scc in self.section_run.section_single_configuration_calculation:

            # In order to resolve band gaps, we need a reference to the highest
            # occupied energy.
            valence_band_maximum = scc.energy_reference_highest_occupied

            # In order to resolve the special points and the reciprocal cell,
            # we need informatoin about the system.
            system = scc.single_configuration_calculation_to_system_ref

            for band in scc.section_k_band:
                if band.band_structure_kind != "vibrational":
                    self.add_reciprocal_cell(band, system)
                    self.add_brillouin_zone(band)
                    self.add_band_gaps(band, valence_band_maximum)

    def add_reciprocal_cell(self, band: section_k_band, system: section_system):
        """A reciprocal cell for this calculation. If the original unit cell is
        not a primitive one, then we will use the one given by spglib.

        If the used code is calculating it's own primitive cell, then the
        reciprocal cell used in e.g. band structure calculation might differ
        from the one given by spglib. To overcome this we would need to get the
        exact reciprocal cell used by the calculation or then test for
        different primitive cells
        (https://atztogo.github.io/spglib/definition.html#transformation-to-the-primitive-cell)
        whether the k-point path stays inside the first Brillouin zone.
        """
        try:
            orig_atoms = system.m_cache["representative_atoms"]
            symmetry_analyzer = system.section_symmetry[0].m_cache["symmetry_analyzer"]
            prim_atoms = symmetry_analyzer.get_primitive_system()
        except Exception:
            self.logger.info("Could not resolve reciprocal cell.")
            return

        primitive_cell = prim_atoms.get_cell()
        source_cell = orig_atoms.get_cell()

        volume_primitive = primitive_cell.volume
        volume_source = source_cell.volume
        volume_diff = abs(volume_primitive - volume_source)

        if volume_diff > (0.001)**3:
            recip_cell = primitive_cell.reciprocal() * 1e10
        else:
            recip_cell = source_cell.reciprocal() * 1e10

        band.reciprocal_cell = recip_cell

    def add_brillouin_zone(self, band: section_k_band) -> None:
        """Adds a dictionary containing the information needed to display
        the Brillouin zone for this material. This functionality could be put
        into the GUI directly, with the Brillouin zone construction performed
        from the reciprocal cell.

        The Brillouin Zone is a Wigner-Seitz cell, and is thus uniquely
        defined. It's shape does not depend on the used primitive cell.
        """
        recip_cell = band.reciprocal_cell
        if recip_cell is None:
            self.logger.info("Could not resolve Brillouin zone as reciprocal cell is missing.")
            return

        brillouin_zone = atomutils.get_brillouin_zone(recip_cell.magnitude)
        bz_json = json.dumps(brillouin_zone)
        band.brillouin_zone = bz_json

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

    def add_band_gaps(self, band: section_k_band, valence_band_maximum: np.array) -> None:
        """Given the band structure and fermi level, calculates the band gap
        for spin channels and also reports the total band gap as the minum gap
        found.
        """
        if valence_band_maximum is None:
            self.logger.info("Could not resolve band gaps as the energy reference is missing.")
            return

        reciprocal_cell = band.reciprocal_cell
        if reciprocal_cell is None:
            self.logger.info("Could not resolve band gaps as reciprocal cell is missing.")
            return

        # Gather the energies and k points from each segment into one big
        # array
        reciprocal_cell = reciprocal_cell.magnitude
        valence_band_maximum = valence_band_maximum.magnitude
        path: np.array = []
        energies: np.array = []
        for segment in band.section_k_band_segment:
            try:
                seg_k_points = segment.band_k_points
                seg_energies = segment.band_energies
            except Exception:
                return
            if seg_k_points is None or seg_energies is None:
                self.logger.info("Could not resolve band gaps as energies or k points are missing.")
                return
            else:
                seg_energies = seg_energies.magnitude

            seg_energies = np.swapaxes(seg_energies, 1, 2)
            path.append(seg_k_points)
            energies.append(seg_energies)

        path = np.concatenate(path, axis=0)
        energies = np.concatenate(energies, axis=2)

        # Handle spin channels separately to find gaps for spin up and down
        n_channels = energies.shape[0]  # pylint: disable=E1136  # pylint/issues/3139
        gaps: List[section_band_gap] = [None, None]

        for channel in range(n_channels):
            if n_channels == 1:
                vbm = valence_band_maximum
            else:
                vbm = valence_band_maximum[channel]
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
            band_minima_tol = band_minima + config.normalize.band_structure_energy_tolerance
            band_maxima_tol = band_maxima - config.normalize.band_structure_energy_tolerance

            for band_idx in range(num_bands):
                band_min = band_minima[band_idx]
                band_max = band_maxima[band_idx]
                band_min_tol = band_minima_tol[band_idx]
                band_max_tol = band_maxima_tol[band_idx]

                # If any of the bands band crosses the Fermi level, there is no
                # band gap
                if band_min_tol <= vbm and band_max_tol >= vbm:
                    break
                # Whole band below Fermi level, save the current highest
                # occupied band point
                elif band_min_tol <= vbm and band_max_tol <= vbm:
                    gap_lower_energy = band_max
                    gap_lower_idx = band_maxima_idx[band_idx]
                    lower_defined = True
                # Whole band above Fermi level, save the current lowest
                # unoccupied band point
                elif band_min_tol >= vbm:
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

                gap = section_band_gap()
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
                band.m_add_sub_section(section_k_band.section_band_gap, gaps[0])
        # For polarized calculations we report the gap separately for both
        # channels. Also we report the smaller gap as the the total gap for the
        # calculation.
        elif n_channels == 2:
            if gaps[0] is not None:
                band.m_add_sub_section(section_k_band.section_band_gap_spin_up, gaps[0])
            if gaps[1] is not None:
                band.m_add_sub_section(section_k_band.section_band_gap_spin_down, gaps[1])
            if gaps[0] is not None and gaps[1] is not None:
                if gaps[0].value <= gaps[1].value:
                    band.m_add_sub_section(section_k_band.section_band_gap, gaps[0])
                else:
                    band.m_add_sub_section(section_k_band.section_band_gap, gaps[1])
            else:
                if gaps[0] is not None:
                    band.m_add_sub_section(section_k_band.section_band_gap, gaps[0])
                elif gaps[1] is not None:
                    band.m_add_sub_section(section_k_band.section_band_gap, gaps[1])
