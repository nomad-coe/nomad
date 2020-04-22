import ase
from nomad.normalizing.normalizer import Normalizer
from nomad.constants import pi


class BandStructureNormalizer(Normalizer):
    def __init__(self):
        pass

    def normalize():
        pass


def crystal_structure_from_cell(cell, eps=1e-4):
    """Return the crystal structure as a string calculated from the cell.
    """
    cellpar = ase.geometry.cell_to_cellpar(cell=cell)
    abc = cellpar[:3]
    angles = cellpar[3:] / 180 * pi
    a, b, c = abc
    alpha, beta, gamma = angles

    # According to:
    # https://www.physics-in-a-nutshell.com/article/6/symmetry-crystal-systems-and-bravais-lattices#the-seven-crystal-systems
    # If a=b=c and alpha=beta=gamma=90degrees we have cubic.
    if abc.ptp() < eps and abs(angles - pi / 2).max() < eps:
        return 'cubic'
    elif abc.ptp() < eps and abs(angles - pi / 3).max() < eps:
        return 'fcc'
    elif abc.ptp() < eps and abs(angles - np.arccos(-1 / 3)).max() < eps:
        return 'bcc'
    # If a=b!=c, alpha=beta=gamma=90deg, tetragonal.
    elif abs(a - b) < eps and abs(angles - pi / 2).max() < eps:
        return 'tetragonal'
    elif abs(angles - pi / 2).max() < eps:
        return 'orthorhombic'
    # if a = b != c , alpha = beta = 90deg, gamma = 120deg, hexagonal
    elif (abs(a - b) < eps
          and abs(gamma - pi / 3 * 2) < eps
          and abs(angles[:2] - pi / 2).max() < eps):
        return 'hexagonal'
    elif (c >= a and c >= b and alpha < pi / 2 and
          abs(angles[1:] - pi / 2).max() < eps):
        return 'monoclinic'
    else:
        raise ValueError('Cannot find crystal structure')


def get_special_points(cell, eps=1e-4):
    """Return dict of special points.

    The definitions are from a paper by Wahyu Setyawana and Stefano
    Curtarolo::

        http://dx.doi.org/10.1016/j.commatsci.2010.05.010

    lattice: str
        One of the following: cubic, fcc, bcc, orthorhombic, tetragonal,
        hexagonal or monoclinic.
    cell: 3x3 ndarray
        Unit cell.
    eps: float
        Tolerance for cell-check.
    """

    lattice = crystal_structure_from_cell(cell)

    cellpar = ase.geometry.cell_to_cellpar(cell=cell)
    abc = cellpar[:3]
    angles = cellpar[3:] / 180 * pi
    a, b, c = abc
    alpha, beta, gamma = angles

    # Check that the unit-cells are as in the Setyawana-Curtarolo paper:
    if lattice == 'cubic':
        assert abc.ptp() < eps and abs(angles - pi / 2).max() < eps
    elif lattice == 'fcc':
        assert abc.ptp() < eps and abs(angles - pi / 3).max() < eps
    elif lattice == 'bcc':
        angle = np.arccos(-1 / 3)
        assert abc.ptp() < eps and abs(angles - angle).max() < eps
    elif lattice == 'tetragonal':
        assert abs(a - b) < eps and abs(angles - pi / 2).max() < eps
    elif lattice == 'orthorhombic':
        assert abs(angles - pi / 2).max() < eps
    elif lattice == 'hexagonal':
        assert abs(a - b) < eps
        assert abs(gamma - pi / 3 * 2) < eps
        assert abs(angles[:2] - pi / 2).max() < eps
    elif lattice == 'monoclinic':
        assert c >= a and c >= b
        assert alpha < pi / 2 and abs(angles[1:] - pi / 2).max() < eps
    if lattice != 'monoclinic':
        return special_points[lattice]

    # Here, we need the cell:
    eta = (1 - b * np.cos(alpha) / c) / (2 * np.sin(alpha)**2)
    nu = 1 / 2 - eta * c * np.cos(alpha) / b
    return {'Γ': [0, 0, 0],
            'A': [1 / 2, 1 / 2, 0],
            'C': [0, 1 / 2, 1 / 2],
            'D': [1 / 2, 0, 1 / 2],
            'D1': [1 / 2, 0, -1 / 2],
            'E': [1 / 2, 1 / 2, 1 / 2],
            'H': [0, eta, 1 - nu],
            'H1': [0, 1 - eta, nu],
            'H2': [0, eta, -nu],
            'M': [1 / 2, eta, 1 - nu],
            'M1': [1 / 2, 1 - eta, nu],
            'M2': [1 / 2, eta, -nu],
            'X': [0, 1 / 2, 0],
            'Y': [0, 0, 1 / 2],
            'Y1': [0, 0, -1 / 2],
            'Z': [1 / 2, 0, 0]}

specialPoints = {}
try:
    specialPoints = get_special_points(
        convert_unit_function("m", "angstrom")(self.cell))
except Exception as e:
    self.logger.warn("failed to get special points/xtal structure", exc_info=e)


    # Loop over bands
    for band_data in bands:

        # Add reciprocal cell and brillouin zone
        self.get_reciprocal_cell(band_data, prim_atoms, orig_atoms)
        self.get_brillouin_zone(band_data)

        # If the parser has reported a valence band maximum, we can add
        # band gap information.
        if valence_band_maximum is not None:
            kpoints = []
            energies = []
            segments = band_data.section_k_band_segment
            if not segments:
                return

            # Loop over segments
            for segment_src in segments:
                try:
                    seg_k_points = segment_src.band_k_points
                    seg_energies = segment_src.band_energies
                except Exception:
                    return
                if seg_k_points is None or seg_energies is None:
                    return
                else:
                    seg_energies = seg_energies.magnitude

                seg_energies = np.swapaxes(seg_energies, 1, 2)
                kpoints.append(seg_k_points)
                energies.append(seg_energies)

            # Continue to calculate band gaps and other properties.
            kpoints = np.concatenate(kpoints, axis=0)
            energies = np.concatenate(energies, axis=2)

            # If the parser has reported a valence band maximum, we can add
            # band gap information.
            self.get_band_gaps(band_data, energies, kpoints, valence_band_maximum)

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

