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
