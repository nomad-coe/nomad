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

from nptyping import NDArray
import numpy as np
import ase

from nomad.datamodel.metainfo.simulation.calculation import (
    BandStructure, BandGap
)
from nomad.datamodel.metainfo.simulation.system import System
from nomad.normalizing.normalizer import Normalizer
from nomad import config, atomutils
from nomad.constants import pi


class BandStructureNormalizer(Normalizer):
    """Normalizer with the following responsibilities:

      - Calculates band gaps and energy references.
      - TODO: Creates labels for special points within the band path (band_segm_labels)
      - TODO: Determines if the path is a standard one or not (is_standard)
    """
    def normalize(self, logger=None) -> None:
        # Setup logger
        if logger is not None:
            self.logger = logger.bind(normalizer=self.__class__.__name__)

        # Do nothing if section run is not present
        if self.section_run is None:
            return

        # Loop through the bands
        for scc in self.section_run.calculation:

            # In order to resolve band gaps, we need a reference to the highest
            # occupied energy or the Fermi energy
            energy_fermi = scc.energy.fermi if scc.energy is not None else None
            energy_highest = scc.energy.highest_occupied if scc.energy is not None else None
            energy_lowest = scc.energy.lowest_unoccupied if scc.energy is not None else None

            # In order to resolve the special points and the reciprocal cell,
            # we need information about the system.
            system = scc.system_ref
            for band in scc.band_structure_electronic:
                valid_band = self.validate_band(band)
                if valid_band:
                    self.add_reciprocal_cell(band, system)
                    self.add_band_gap(
                        band,
                        energy_fermi,
                        energy_highest,
                        energy_lowest
                    )
                    self.add_path_labels(band, system)

    def validate_band(self, band: BandStructure) -> bool:
        """Used to check that a band has all required information for normalization.
        """
        if len(band.segment) == 0:
            self.logger.info("Could not normalize band structure as band segments are missing.")
            return False
        for segment in band.segment:
            seg_k_points = segment.kpoints
            seg_energies = segment.energies
            if seg_k_points is None or seg_energies is None:
                self.logger.info("Could not normalize band structure as energies or k points are missing.")
                return False
        return True

    def add_reciprocal_cell(self, band: BandStructure, system: System):
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
        # If reciprocal cell is reported by parser, use it.
        if band.reciprocal_cell is not None:
            return
        try:
            orig_atoms = system.m_cache["representative_atoms"]
            symmetry_analyzer = system.symmetry[0].m_cache["symmetry_analyzer"]
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

    def get_k_space_distance(self, reciprocal_cell: NDArray, point1: NDArray, point2: NDArray) -> float:
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

    def add_band_gap(
            self,
            band: BandStructure,
            energy_fermi: NDArray,
            energy_highest: NDArray,
            energy_lowest: NDArray) -> None:
        """Given the band structure and information about energy references,
        determines the band gap and energy references separately for all spin
        channels.
        """
        band.energy_fermi = energy_fermi

        # No reference data available
        eref = energy_highest if energy_fermi is None else energy_fermi
        if eref is None:
            self.logger.info("could not resolve energy references or band gaps for band structure")
            return
        eref = eref.magnitude

        # Create energy reference sections for each spin channel, add fermi
        # energy if present
        n_channels = band.segment[0].energies.shape[0]
        for i_channel in range(n_channels):
            info = band.band_gap[i_channel] if len(band.band_gap) > i_channel else band.m_create(BandGap)
            info.index = i_channel
            if energy_highest is not None:
                info.energy_highest_occupied = energy_highest
            if energy_lowest is not None:
                info.energy_lowest_unoccupied = energy_lowest

        # Gather the energies and k points from each segment into one big array
        path: NDArray = []
        energies: NDArray = []
        for segment in band.segment:
            seg_k_points = segment.kpoints
            seg_energies = segment.energies
            seg_energies = seg_energies.magnitude
            seg_energies = np.swapaxes(seg_energies, 1, 2)
            path.append(seg_k_points)
            energies.append(seg_energies)

        path = np.concatenate(path, axis=0)
        energies = np.concatenate(energies, axis=2)

        # Use a reference energy (fermi or highest occupied) to determine the
        # energy references from the band structure (discretization will affect
        # the exact location).
        for i_channel in range(n_channels):
            i_energy_highest = None
            i_energy_lowest = None
            channel_energies = energies[i_channel, :, :]
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
                if band_min_tol <= eref and band_max_tol >= eref:
                    break
                # Whole band below Fermi level, save the current highest
                # occupied band point
                elif band_min_tol <= eref and band_max_tol <= eref:
                    i_energy_highest = band_max
                    gap_lower_idx = band_maxima_idx[band_idx]
                # Whole band above Fermi level, save the current lowest
                # unoccupied band point
                elif band_min_tol >= eref:
                    i_energy_lowest = band_min
                    gap_upper_idx = band_minima_idx[band_idx]
                    break

            # Save the found energy references
            if i_energy_highest is not None:
                band.band_gap[i_channel].energy_highest_occupied = i_energy_highest
            if i_energy_lowest is not None:
                band.band_gap[i_channel].energy_lowest_unoccupied = i_energy_lowest

            # If highest occupied energy and a lowest unoccupied energy are
            # found, and the difference between them is positive, save
            # information about the band gap.
            gap_value = 0.0
            info = band.band_gap[i_channel]
            if i_energy_lowest is not None and i_energy_highest is not None:
                gap_value = float(i_energy_lowest - i_energy_highest)
                if gap_value > 0:
                    # See if the gap is direct or indirect by comparing the k-point
                    # locations with some tolerance
                    k_point_lower = path[gap_lower_idx]
                    k_point_upper = path[gap_upper_idx]
                    reciprocal_cell = band.reciprocal_cell
                    if reciprocal_cell is not None:
                        reciprocal_cell = reciprocal_cell.magnitude
                        k_point_distance = self.get_k_space_distance(reciprocal_cell, k_point_lower, k_point_upper)
                        is_direct_gap = k_point_distance <= config.normalize.k_space_precision
                        band_gap_type = "direct" if is_direct_gap else "indirect"
                        info.type = band_gap_type
            info.value = gap_value

    def add_path_labels(self, band: BandStructure, system: System) -> None:
        """Adds special high symmmetry point labels to the band path. Only k
        points that land on the special points defined by Setyawan/Curtarolo
        are automatically labeled.
        """
        # If labels are already set by the parser dot nothing.
        for segment in band.segment:
            labels = segment.endpoints_labels
            if labels is not None:
                self.logger.info("Existing band segment labels detected, skipping label detection.")
                return

        # Try to get the required data. Fail if not found.
        try:
            lattice_vectors = system.atoms.lattice_vectors.to("angstrom").magnitude
            reciprocal_cell_trans = band.reciprocal_cell.magnitude.T
            bravais_lattice = system.symmetry[0].bravais_lattice
        except Exception:
            self.logger.info("Could not resolve path labels as required information is missing.")
            return

        # Find special points for this lattice. If an error occurs, the labels
        # are simply not written.
        special_points = self.get_special_points(bravais_lattice, lattice_vectors)
        if not special_points:
            self.logger.warning("Could not resolve high-symmetry points for the given simulation cell.")
            return

        # Form a contiguous array of k points for faster operations
        special_point_labels = list(special_points.keys())
        special_k_points = np.empty((len(special_points), 3))
        for i, kpt in enumerate(special_points.values()):
            special_k_points[i, :] = kpt
        special_k_points_cartesian = np.dot(special_k_points, reciprocal_cell_trans)

        # Match tolerance in 1/m. Taken from the VASP parser.
        eps = config.normalize.k_space_precision

        # Try to find matches for the special points. We only attempt to match
        # points at the start and end of a segment. Any labels set by the
        # parser are overridden, because one cannot ascertain that those labels
        # are consistent across codes.
        for segment in band.segment:
            start_point_cartesian = np.dot(segment.kpoints[0], reciprocal_cell_trans)
            end_point_cartesian = np.dot(segment.kpoints[-1], reciprocal_cell_trans)

            # Calculate distance in cartesian space
            start_index = atomutils.find_match(start_point_cartesian, special_k_points_cartesian, eps)
            end_index = atomutils.find_match(end_point_cartesian, special_k_points_cartesian, eps)

            if start_index is None:
                start_label = ""
            else:
                start_label = special_point_labels[start_index]
            if end_index is None:
                end_label = ""
            else:
                end_label = special_point_labels[end_index]
            segment.endpoints_labels = [start_label, end_label]

    def get_special_points(self, bravais_lattice, lattice_vectors, eps=3e-3):
        """Return dict of special points.

        The definitions are from a paper by Wahyu Setyawana and Stefano
        Curtarolo:

            http://dx.doi.org/10.1016/j.commatsci.2010.05.010

        bravais_lattice: str
            Bravais lattice in Pearson notation.
        lattice_vectors: 3x3 ndarray
            Lattice vectors of the primitive cell.
        eps: float
            tolerance factor to define the Bravais lattice parameters
        """
        cell = ase.cell.Cell(lattice_vectors)
        lattice = cell.get_bravais_lattice(eps)
        special_points = {}
        try:
            # Non-conventional ordering for certain lattices:
            if bravais_lattice in ['oP', 'oF', 'oI', 'oS']:
                a, b, c = lattice.a, lattice.b, lattice.c
                assert a < b
                if bravais_lattice != 'oS':
                    assert b < c
            elif bravais_lattice in ['mP', 'mS']:
                a, b, c = lattice.a, lattice.b, lattice.c
                alpha = lattice.alpha * pi / 180
                assert a <= c and b <= c  # ordering of the conventional lattice
                assert alpha < pi / 2

            for (keys, values) in lattice.get_special_points().items():
                if keys == 'G':
                    keys = 'Γ'
                if bravais_lattice == 'tI':
                    if keys == 'S':
                        keys = 'Σ'
                    elif keys == 'S1':
                        keys = 'Σ1'
                special_points[keys] = list(values)
        except Exception:
            pass

        return special_points
