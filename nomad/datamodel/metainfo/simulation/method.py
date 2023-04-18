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

import numpy as np            # pylint: disable=unused-import
import typing

from pint.util import SharedRegistryObject                 # pylint: disable=unused-import
from nomad.datamodel.data import ArchiveSection
from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference, MEnum, derived)
from ..common import FastAccess


m_package = Package()


class Mesh(MSection):
    '''
    Contains the settings for a sampling mesh.
    Supports uniformly-spaced meshes and symmetry-reduced representations.
    '''

    m_def = Section(validate=False)

    dimensionality = Quantity(
        type=np.int32,
        shape=[],
        description='''
        Dimensionality of the mesh.
        ''')

    sampling_method = Quantity(
        type=MEnum(
            'Gamma-centered', 'Monkhorst-Pack', 'Gamma-offcenter', 'Line-path',
            'Equidistant', 'Logarithmic', 'Gauss-Legendre', 'Gauss-Laguerre'
            'Clenshaw-Curtis', 'Newton-Cotes', 'Gauss-Hermite'
        ),
        shape=[],
        description='''
        Method used to generate the mesh:

        | Name      | Description                      | Reference             |

        | --------- | -------------------------------- | --------------------- |

        | `'Gamma-centered'` | Regular mesh is centered around Gamma. No offset. |

        | `'Monkhorst-Pack'` | Regular mesh with an offset of half the reciprocal lattice vector. |

        | `'Gamma-offcenter'` | Regular mesh with an offset that is neither `'Gamma-centered'`, nor `'Monkhorst-Pack'`. |

        | `'Line-path'` | Line path along high-symmetry points. Typically employed for simualting band structures. |

        | `'Equidistant'`  | Equidistant 1D grid (also known as 'Newton-Cotes')                      |

        | `'Logarithmic'`  | log distance 1D grid               |

        | `'Gauss-Legendre'` | Quadrature rule for integration using Legendre polynomials |

        | `'Gauss-Laguerre'` | Quadrature rule for integration using Laguerre polynomials |

        | `'Clenshaw-Curtis'`  | Quadrature rule for integration using Chebyshev polynomials using discrete cosine transformations |

        | `'Gauss-Hermite'`  | Quadrature rule for integration using Hermite polynomials |
        ''')

    n_points = Quantity(
        type=np.int32,
        shape=[],
        description='''
        Total number of points in the mesh, accounting for the multiplicities.
        ''')

    grid = Quantity(
        type=np.int32,
        shape=['dimensionality'],
        description='''
        Amount of mesh point sampling along each axis, i.e. [nx, ny, nz].
        ''')

    points = Quantity(
        type=np.complex128,
        shape=['*', 'dimensionality'],
        description='''
        List of all the points in the mesh.
        ''')

    multiplicities = Quantity(
        type=np.float64,
        shape=['*'],
        description='''
        The amount of times the same point reappears. These are accounted for in `n_points`.
        A value larger than 1, typically indicates a symmtery operation that was applied to the mesh.
        ''')

    weights = Quantity(
        type=np.float64,
        shape=['*'],
        description='''
        The frequency of times the same point reappears.
        A value larger than 1, typically indicates a symmtery operation that was applied to the mesh.
        ''')


class LinePathSegment(MSection):
    '''
    Contains the settings for a single line path segment in a mesh.
    '''

    m_def = Section(validate=False)

    start_point = Quantity(
        type=str,
        shape=[],
        description='''
        Name of the hihg-symmetry starting point of the line path segment.
        ''')

    end_point = Quantity(
        type=str,
        shape=[],
        description='''
        Name of the high-symmetry end point of the line path segment.
        ''')

    n_points = Quantity(
        type=np.int32,
        shape=[],
        description='''
        Number of points in the line path segment.
        ''')

    points = Quantity(
        type=np.float64,
        shape=['*', 3],
        description='''
        List of all the points in the line path segment.
        ''')


class KMesh(Mesh):
    '''
    Contains the settings for a sampling mesh in 3D reciprocal space.
    Supports uniformly-spaced meshes, line paths along high-symmetry points,
    as well as symmetry-reduced and full representations.
    '''

    m_def = Section(validate=False)

    offset = Quantity(
        type=np.float64,
        shape=[3],
        description='''
        Offset vector shifting the mesh with respect to a Gamma-centered case.
        '''
    )

    all_points = Quantity(
        type=np.dtype(np.float64),
        shape=['*', 3],
        description='''
        Full list of the mesh points without any symmetry operations.
        ''')

    high_symmetry_points = Quantity(
        type=str,
        shape=['*'],
        description='''
        Named high symmetry points in the mesh.
        ''')

    line_path_segments = SubSection(sub_section=LinePathSegment.m_def, repeats=True)


class FrequencyMesh(Mesh):
    '''
    Contains the settings for a sampling mesh in 1D frequency space, either real or imaginary.
    '''

    m_def = Section(validate=False)

    points = Quantity(
        type=np.complex128,
        shape=['n_points', 'dimensionality'],
        unit='joule',
        description='''
        List of all the points in the mesh in joules.
        ''')

    smearing = Quantity(
        type=np.float64,
        shape=[],
        description='''
        Numerical smearing parameter used for convolutions.
        ''')


class TimeMesh(Mesh):
    '''
    Contains the settings for a sampling mesh in 1D time space, either real or imaginary.
    '''

    m_def = Section(validate=False)

    smearing = Quantity(
        type=np.float64,
        shape=[],
        description='''
        Numerical smearing parameter used for convolutions.
        ''')


class Scf(MSection):
    '''
    Section containing the parameters related to self consistency.
    '''

    m_def = Section(validate=False)

    n_max_iteration = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Specifies the maximum number of allowed self-consistent field (SCF) iterations in
        a calculation.
        ''')

    threshold_energy_change = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Specifies the threshold for the total energy change between two subsequent
        self-consistent field (SCF) iterations. The SCF is considered converged when the
        total-energy change between two SCF cycles is below the threshold (possibly in
        combination with other criteria).
        ''')

    threshold_density_change = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Specifies the threshold for the average  charge density change between two
        subsequent self-consistent field (SCF) iterations. The SCF is considered converged
        when the density change between two SCF cycles is below the threshold (possibly in
        combination with other criteria).
        ''')

    minimization_algorithm = Quantity(
        type=str,
        shape=[],
        description='''
        Specifies the algorithm used for self consistency minimization.
        ''')


class HubbardKanamoriModel(MSection):
    '''
    Setup of the local Hubbard model.
    '''

    m_def = Section(validate=False)

    orbital = Quantity(
        type=str,
        shape=[],
        description='''
        Orbital label corresponding to the Hubbard model. The typical orbitals with strong
        Hubbard interactions have partially filled '3d', '4d' and '4f' orbitals.
        ''')

    n_orbital = Quantity(
        type=np.int32,
        shape=[],
        description='''
        Number of non-degenerated orbitals of the same type (s, p, d, f, ...).
        ''')

    u = Quantity(
        type=np.float64,
        shape=[],
        unit='joule',
        description='''
        Value of the (intraorbital) Hubbard interaction
        ''')

    jh = Quantity(
        type=np.float64,
        shape=[],
        unit='joule',
        description='''
        Value of the (interorbital) Hund's coupling.
        ''')

    up = Quantity(
        type=np.float64,
        shape=[],
        unit='joule',
        description='''
        Value of the (interorbital) Coulomb interaction. In rotational invariant
        systems, up = u - 2 * jh.
        ''')

    j = Quantity(
        type=np.float64,
        shape=[],
        unit='joule',
        description='''
        Value of the exchange interaction. In rotational invariant systems, j = jh.
        ''')

    u_effective = Quantity(
        type=np.float64,
        shape=[],
        unit='joule',
        description='''
        Value of the effective U parameter (u - j).
        ''')

    slater_integrals = Quantity(
        type=np.float64,
        shape=[3],
        unit='joule',
        description='''
        Value of the Slater integrals (F0, F2, F4) in spherical harmonics used to derive
        the local Hubbard interactions:

            u = ((2.0 / 7.0) ** 2) * (F0 + 5.0 * F2 + 9.0 * F4) / (4.0*np.pi)

            up = ((2.0 / 7.0) ** 2) * (F0 - 5.0 * F2 + 3.0 * 0.5 * F4) / (4.0*np.pi)

            jh = ((2.0 / 7.0) ** 2) * (5.0 * F2 + 15.0 * 0.25 * F4) / (4.0*np.pi)

        Ref.: Elbio Dagotto, Nanoscale Phase Separation and Colossal Magnetoresistance,
        Chapter 4, Springer Berlin (2003).
        ''')

    umn = Quantity(
        type=np.float64,
        shape=['n_orbital', 'n_orbital'],
        unit='joule',
        description='''
        Value of the local Coulomb interaction matrix.
        ''')

    double_counting_correction = Quantity(
        type=str,
        shape=[],
        description='''
        Name of the double counting correction algorithm applied.
        ''')


class AtomParameters(MSection):
    '''
    Contains method-related information about a kind of atom identified by label. This
    allows the assignment of an atom-centered basis set or pseudopotential for different
    atoms belonging to the same kind.

    Through this section we use the wording "active" mainly for defining orbital-related
    quantities. Active refers to the relevant orbital parameters in the atom.
    '''

    m_def = Section(validate=False)

    atom_number = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Atomic number (number of protons) of this atom kind, use 0 if not an atom.
        ''')

    n_valence_electrons = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Number of valence electrons.
        ''')

    label = Quantity(
        type=str,
        shape=[],
        description='''
        String used to identify the atoms of this kind. This should correspond to the
        atom labels of the configuration. It is possible for one atom kind to have
        multiple labels (in order to allow two atoms of the same kind to have two
        differently defined sets of atom-centered basis functions or two different pseudo-
        potentials). Atom kind is typically the symbol of the atomic species but it can be
        also a ghost or pseudo-atom.
        ''')

    mass = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='kg',
        description='''
        Mass of the atom.
        ''')

    pseudopotential_name = Quantity(
        type=str,
        shape=[],
        description='''
        Name identifying the pseudopotential used.
        ''')

    n_orbitals = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of active orbitals of the atom.
        ''')

    orbitals = Quantity(
        type=str,
        shape=['n_orbitals'],
        description='''
        Label of the active orbitals of the atoms.
        ''')

    onsite_energies = Quantity(
        type=np.dtype(np.float64),
        shape=['n_orbitals'],
        unit='joule',
        description='''
        Values of the atomic onsite energy corresponding to each orbital.
        ''')

    charge = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='coulomb',
        description='''
        Total charge of the atom.
        ''')

    charges = Quantity(
        type=np.dtype(np.float64),
        shape=['n_orbitals'],
        unit='coulomb',
        description='''
        Values of the charge corresponding to each orbital.
        ''')

    hubbard_kanamori_model = SubSection(sub_section=HubbardKanamoriModel.m_def, repeats=False)


class MoleculeParameters(MSection):
    '''
    Contains method-related information about a kind of atom identified by label. This
    allows the assignment of an atom-centered basis set or pseudopotential for different
    atoms belonging to the same kind.
    '''

    m_def = Section(validate=False)

    label = Quantity(
        type=str,
        shape=[],
        description='''
        String to identify the molecule.
        ''')

    n_atoms = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of atoms in the molecule.
        ''')

    atom_parameters = SubSection(sub_section=AtomParameters.m_def, repeats=True)


class Photon(MSection):
    '''
    Section containing the details of the photon field used for spectrum calculations.
    '''

    m_def = Section(validate=False)

    multipole_type = Quantity(
        type=str,
        description='''
        Type used for the multipolar expansion: dipole, quadrupole, NRIXS, Raman, etc.
        ''')

    polarization = Quantity(
        type=np.float64,
        shape=[3],
        description='''
        Direction of the photon polarization in cartesian coordinates.
        ''')

    energy = Quantity(
        type=np.float64,
        unit='joule',
        description='''
        Photon energy.
        ''')

    momentum_transfer = Quantity(
        type=np.float64,
        shape=[3],
        description='''
        Momentum transfer which would be important for quadrupole or NRIXS or Raman.
        ''')


class GaussianBasisGroup(MSection):
    '''
    Section that describes a group of Gaussian contractions. Groups allow one to calculate
    the primitive Gaussian integrals once for several different linear combinations of
    them. This defines basis functions with radial part $f_i(r) = r^{l_i} \\sum_{j} c_{i j}
    A(l_i, \\alpha_j) exp(-\\alpha_j r^2)$ where $A(l_i, \\alpha_j)$ is a the normalization
    coefficient for primitive Gaussian basis functions. Here, $\\alpha_j$ is defined in
    gaussian_basis_group_exponents, $l_i$ is given in gaussian_basis_group_ls, and $c_{i
    j}$ is given in gaussian_basis_group_contractions, whereas the radial part is given by
    the spherical harmonics $Y_{l m}$.

    This section is defined only if the original basis function uses Gaussian basis
    functions, and the sequence of radial functions $f_i$ across all
    section_gaussian_basis_group in section_basis_set_atom_centered should match the one
    of basis_set_atom_centered_radial_functions.
    '''

    m_def = Section(validate=False)

    n_contractions = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of different contractions, i.e. resulting basis functions in a
        gaussian_basis_group section.
        ''')

    n_exponents = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of different Gaussian exponents in a section_gaussian_basis_group
        section.
        ''')

    contractions = Quantity(
        type=np.dtype(np.float64),
        shape=['n_contractions', 'n_exponents'],
        description='''
        contraction coefficients $c_{i j}$ defining the contracted basis functions with
        respect to *normalized* primitive Gaussian functions. They define the Gaussian
        basis functions as described in section_gaussian_basis_group.
        ''')

    exponents = Quantity(
        type=np.dtype(np.float64),
        shape=['n_exponents'],
        unit='1 / meter ** 2',
        description='''
        Exponents $\\alpha_j$ of the Gaussian functions defining this basis set
        $exp(-\\alpha_j r^2)$. One should be careful about the units of the coefficients.
        ''')

    ls = Quantity(
        type=np.dtype(np.float64),
        shape=['n_contractions'],
        description='''
        Azimuthal quantum number ($l$) values (of the angular part given by the spherical
        harmonic $Y_{l m}$ of the various contracted basis functions).
        ''')


class BasisSetAtomCentered(MSection):
    '''
    This section describes the atom-centered basis set. The main contained information is
    a short, non unique but human-interpretable, name for identifying the basis set
    (short_name), a longer unique name, the atomic number of the atomic species the
    basis set is meant for.
    '''

    m_def = Section(validate=False)

    name = Quantity(
        type=str,
        shape=[],
        description='''
        Code-specific, but explicative, base name for the basis set.
        ''')

    atom_number = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Atomic number (i.e., number of protons) of the atom for which this basis set is
        constructed (0 means unspecified or a pseudo atom).
        ''')

    n_basis_functions = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Gives the number of different basis functions in a basis_set_atom_centered
        section. This equals the number of actual coefficients that are specified when
        using this basis set.
        ''')

    gaussian_basis_group = SubSection(sub_section=GaussianBasisGroup.m_def, repeats=True)


class BasisSetCellDependent(MSection):
    '''
    Section describing a cell-dependent (atom-independent) basis set, e.g. plane waves.
    The contained information is the type of basis set (in basis_set_cell_dependent_kind),
    its parameters (e.g., for plane waves in basis_set_planewave_cutoff), and a name that
    identifies the actually used basis set (a string combining the type and the
    parameter(s), stored in name).
    '''

    m_def = Section(validate=False)

    kind = Quantity(
        type=str,
        shape=[],
        description='''
        A string defining the type of the cell-dependent basis set (i.e., non atom
        centered such as plane-waves). Can be one of plane_waves, realspace_grids or
        wavelets.
        ''')

    name = Quantity(
        type=str,
        shape=[],
        description='''
        A label identifying the cell-dependent basis set (i.e., non atom centered such as
        plane-waves). The following convetion should be followed:
        plane_waves ("PW_" + cutoff in Ry) realspace_grids ("GR_" + grid spacing in fm)
        wavelets (WV_" + smallest wavelet spacing in fm).
        ''')

    planewave_cutoff = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Spherical cutoff  in reciprocal space for a plane-wave basis set. It is the energy
        of the highest plan-ewave ($\\frac{\\hbar^2|k+G|^2}{2m_e}$) included in the basis
        set. Note that normally this basis set is used for the wavefunctions, and the
        density would have 4 times the cutoff, but this actually depends on the use of the
        basis set by the method.
        ''')

    grid_spacing = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter',
        description='''
        Grid spacing used for the realspace representation of the wave functions.
        ''')


class BasisSet(MSection):
    '''
    This section contains all basis sets used to represent the wavefunction or electron
    density.
    '''

    m_def = Section(validate=False)

    kind = Quantity(
        type=str,
        shape=[],
        description='''
        String describing the use of the basis set, i.e, if it used for expanding a
        wavefunction or an electron density.
        ''')

    type = Quantity(
        type=str,
        shape=[],
        description='''
        The type of basis set used by the program. Valid values are: [`Numeric AOs`,
        `Gaussians`, `(L)APW+lo`, `plane waves`, `psinc functions`, `real-space grid`].
        ''')

    name = Quantity(
        type=str,
        shape=[],
        description='''
        Identifies the basis set.
        ''')

    cell_dependent = SubSection(sub_section=BasisSetCellDependent.m_def, repeats=True)

    atom_centered = SubSection(sub_section=BasisSetAtomCentered.m_def, repeats=True)


class Interaction(MSection):
    '''
    Section containing the parameters of a contribution to a force field model.
    '''

    m_def = Section(validate=False)

    type = Quantity(
        type=str,
        shape=[],
        description='''
        Denotes the classification of the potential.
        ''')

    name = Quantity(
        type=str,
        shape=[],
        description='''
        Specifies the name of the potential. Can contain information on the species,
        cut-offs, potential versions.
        ''')

    n_atoms = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of atoms included in the interaction
        '''
    )

    atom_labels = Quantity(
        type=str,
        shape=['n_atoms'],
        description='''
        Labels of the atoms described by the interaction.
        ''')

    atom_indices = Quantity(
        type=np.dtype(np.int32),
        shape=['n_atoms'],
        description='''
        Indices of the atoms in the system described by the interaction.
        ''')

    functional_form = Quantity(
        type=str,
        shape=[],
        description='''
        Specifies the functional form of the interaction potential.
        ''')

    n_parameters = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Specifies the number of parameters in the interaction potential.
        ''')

    parameters = Quantity(
        type=typing.Any,
        shape=[],
        description='''
        Dictionary of label and parameters of the interaction potential.
        ''')


class Model(MSection):
    '''
    Section containing the parameters of a force field model. If specified, the parameters
    corresponding to the individual contributions to the model are given in contributions.
    Otherwise, the parameters can also be found in a reference to the published model.
    '''

    m_def = Section(validate=False)

    name = Quantity(
        type=str,
        shape=[],
        description='''
        Identifies the name of the model.
        ''')

    reference = Quantity(
        type=str,
        shape=[],
        description='''
        Reference to the model e.g. DOI, URL.
        ''')

    contributions = SubSection(sub_section=Interaction.m_def, repeats=True)


class Functional(MSection):
    '''
    Section containing the parameters of an exchange or correlation functional.
    '''

    m_def = Section(validate=False)

    name = Quantity(
        type=str,
        shape=[],
        description='''
        Provides the name of one of the exchange and/or correlation (XC) functional
        following the libbx convention.
        ''')

    parameters = Quantity(
        type=typing.Any,
        shape=[],
        description='''
        Contains an associative list of non-default values of the parameters for the
        functional.

        For example, if a calculations using a hybrid XC functional (e.g., HSE06)
        specifies a user-given value of the mixing parameter between exact and GGA
        exchange, then this non-default value is stored in this metadata.

        The labels and units of these values may be defined in name.

        If this metadata is not given, the default parameter values for the functional
        are assumed.
        ''')

    weight = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Provides the value of the weight for the functional.

        This weight is used in the linear combination of the different functionals. If not
        specified then the default is set to 1.
        ''')


class XCFunctional(Model):
    '''
    Section describing the exchange-correlation functional used in the DFT calculation.
    The name of the exchange-correlation functional is given by name and the reference to
    the published functional is provided by reference. Other contributions to the
    functional not covered by exchange, correlation or hybrid types may be specified in
    contributions.
    '''

    m_def = Section(validate=False)

    exchange = SubSection(sub_section=Functional.m_def, repeats=True)

    correlation = SubSection(sub_section=Functional.m_def, repeats=True)

    hybrid = SubSection(sub_section=Functional.m_def, repeats=True)

    contributions = SubSection(sub_section=Functional.m_def, repeats=True)

    def normalize_hybrid(self):
        for hyb in self.hybrid:
            if 'exact_exchange_mixing_factor' in hyb.parameters.keys() and not hyb.name:
                hyb.name += '+alpha'


class DFT(MSection):
    '''
    Section containing the various parameters that define a DFT calculation. These include
    settings for the exchange correlation functionals, LDA+U, etc.
    '''

    m_def = Section(validate=False)

    self_interaction_correction_method = Quantity(
        type=str,
        shape=[],
        description='''
        Contains the name for the self-interaction correction (SIC) treatment used to
        calculate the final energy and related quantities. If skipped or empty, no special
        correction is applied.

        The following SIC methods are available:

        | SIC method                | Description                       |

        | ------------------------- | --------------------------------  |

        | `""`                      | No correction                     |

        | `"SIC_AD"`                | The average density correction    |

        | `"SIC_SOSEX"`             | Second order screened exchange    |

        | `"SIC_EXPLICIT_ORBITALS"` | (scaled) Perdew-Zunger correction explicitly on a
        set of orbitals |

        | `"SIC_MAURI_SPZ"`         | (scaled) Perdew-Zunger expression on the spin
        density / doublet unpaired orbital |

        | `"SIC_MAURI_US"`          | A (scaled) correction proposed by Mauri and co-
        workers on the spin density / doublet unpaired orbital |
        ''')

    xc_functional = SubSection(sub_section=XCFunctional.m_def)


class Wannier(MSection):
    '''
    Section containing the various parameters that define a Wannier tight-binding method.
    '''

    m_def = Section(validate=False)

    n_projected_orbitals = Quantity(
        type=np.int32,
        description='''
        Number of Wannier orbitals used to fit the DFT band structure
        ''')

    n_bands = Quantity(
        type=np.int32,
        description='''
        Number of input Bloch bands to calculate the projection matrix.
        ''')

    is_maximally_localized = Quantity(
        type=bool,
        description='''
        Are the projected orbitals maximally localized or just a single-shot projection?
        ''')

    convergence_tolerance_max_localization = Quantity(
        type=np.float64,
        description='''
        Convergence tolerance for maximal localization of the projected orbitals.
        ''')

    energy_window_outer = Quantity(
        type=np.float64,
        unit='electron_volt',
        shape=[2],
        description='''
        Bottom and top of the outer energy window used for the projection.
        ''')

    energy_window_inner = Quantity(
        type=np.float64,
        unit='electron_volt',
        shape=[2],
        description='''
        Bottom and top of the inner energy window used for the projection.
        ''')


class Projection(MSection):
    '''
    Section containing the various parameters that define a Wannier90-like projection
    '''

    m_def = Section(validate=False)

    wannier = SubSection(sub_section=Wannier.m_def, repeats=False)


class HoppingMatrix(MSection):
    '''
    Section containing the hopping/overlap matrix elements between N projected
    orbitals.
    '''

    m_def = Section(validate=False)

    n_orbitals = Quantity(
        type=np.int32,
        description='''
        Number of projected orbitals.
        ''')

    n_wigner_seitz_points = Quantity(
        type=np.int32,
        description='''
        Number of Wigner-Seitz real points.
        ''')

    degeneracy_factors = Quantity(
        type=np.int32,
        shape=['n_wigner_seitz_points'],
        description='''
        Degeneracy of each Wigner-Seitz grid point.
        ''')

    value = Quantity(
        type=np.float64,
        shape=['n_wigner_seitz_points', 'n_orbitals * n_orbitals', 7],
        description='''
        Real space hopping matrix for each Wigner-Seitz grid point. The elements are
        defined as follows:

            n_x   n_y   n_z   orb_1   orb_2   real_part + j * imag_part

        where (n_x, n_y, n_z) define the Wigner-Seitz cell vector in fractional coordinates,
        (orb_1, orb_2) indicates the hopping amplitude between orb_1 and orb_2, and the
        real and imaginary parts of the hopping in electron_volt.
        ''')


class LatticeModelHamiltonian(MSection):
    '''
    Section containing the parameters of the non-interacting parts of a lattice model Hamiltonian.
    '''

    m_def = Section(validate=False)

    # I defined t_parameters apart from HoppingMatrix to simplify the parsing (writing
    # everything in the HoppingMatrix basis might be tedious).
    # TODO generalize parsers to write HoppingMatrix as default even for simple models.
    lattice_name = Quantity(
        type=str,
        shape=[],
        description='''
        Name of the lattice to identify the model. E.g., 'Square', 'Honeycomb'.
        ''')

    n_neighbors = Quantity(
        type=np.int32,
        description='''
        Number of direct neighbors considered for the hopping integrals.
        ''')

    t_parameters = Quantity(
        type=np.complex128,
        shape=['n_neighbors'],
        description='''
        Hopping parameters for simple models, with [t, t`, t``, etc].
        ''')

    # TODO determine standard dimensions for projection_matrix:
    #   KMesh.n_points x n_atoms_per_unit_cell x n_orbitals? x n_bands?
    projection_matrix = Quantity(
        type=np.complex128,
        shape=['*', '*', '*', '*'],
        description='''
        Projection matrices from Bloch bands to virtual projected orbitals.
        ''')

    hopping_matrix = SubSection(sub_section=HoppingMatrix.m_def, repeats=False)

    hubbard_kanamori_model = SubSection(sub_section=HubbardKanamoriModel.m_def, repeats=True)


class CoreHole(MSection):
    '''
    Section containing the various parameters that define a core-hole calculation. It can
    be within BSE as a "core" subsection.
    '''

    m_def = Section(validate=False)

    solver = Quantity(
        type=str,
        description='''
        Solver algorithm used for the core-hole spectra.
        ''')

    edge = Quantity(
        type=MEnum(
            'K', 'L1', 'L2', 'L3', 'L23', 'M1', 'M2', 'M3', 'M23', 'M4', 'M5', 'M45',
            'N1', 'N2', 'N3', 'N23', 'N4', 'N5', 'N45'),
        description='''
        Edge to be calculated for the core-hole spectra.
        ''')

    mode = Quantity(
        type=MEnum('absorption', 'emission'),
        description='''
        Type of spectra to be calculated: absorption or emission.
        ''')

    broadening = Quantity(
        type=np.float64,
        unit='joule',
        description='''
        Core-hole lifetime broadening applied to the edge spectra in full-width at half maximum.
        ''')


class ExcitedStateMethodology(MSection):
    '''
    Base class containing the common numerical parameters typical of excited-state
    calculations.
    '''

    m_def = Section(validate=False)

    type = Quantity(
        type=str,
        shape=[],
        description='''
        Type which allows to identify the excited-state calculation with a
        common string.
        ''')

    n_states = Quantity(
        type=np.int32,
        shape=[],
        description='''
        Number of states used to calculate the excitations.
        ''')

    n_empty_states = Quantity(
        type=np.int32,
        shape=[],
        description='''
        Number of empty states used to calculate the excitations. This quantity is
        complementary to `n_states`.
        ''')

    broadening = Quantity(
        type=np.float64,
        unit='joule',
        description='''
        Lifetime broadening applied to the spectra in full-width at half maximum.
        ''')

    # Used to define the meshes of sub_sections that are not base classes yet, e.g., screening
    k_mesh = SubSection(sub_section=KMesh.m_def)

    q_mesh = SubSection(sub_section=KMesh.m_def)

    frequency_mesh = SubSection(sub_section=FrequencyMesh.m_def)


class Screening(ExcitedStateMethodology):
    '''
    Section containing the various parameters that define a screening calculation, as for
    example, in RPA.
    '''

    m_def = Section(validate=False)

    dielectric_infinity = Quantity(
        type=np.int32,
        description='''
        Value of the static dielectric constant at infinite q. For metals, this is infinite
        (or a very large value), while for insulators is finite.
        ''')


class GW(ExcitedStateMethodology):
    '''
    Section containing the various parameters that define a GW calculation.
    '''

    m_def = Section(validate=False)

    type = Quantity(
        type=MEnum('G0W0', 'scGW', 'scGW0', 'scG0W', 'ev-scGW0', 'ev-scGW', 'qp-scGW0', 'qp-scGW'),
        shape=[],
        description='''
        GW Hedin's self-consistency cycle:

        | Name      | Description                      | Reference             |

        | --------- | -------------------------------- | --------------------- |

        | `'G0W0'`  | single-shot                      | PRB 74, 035101 (2006) |

        | `'scGW'`  | self-consistent G and W               | PRB 75, 235102 (2007) |

        | `'scGW0'` | self-consistent G with fixed W0  | PRB 54, 8411 (1996)   |

        | `'scG0W'` | self-consistent W with fixed G0  | -                     |

        | `'ev-scGW0'`  | eigenvalues self-consistent G with fixed W0   | PRB 34, 5390 (1986)   |

        | `'ev-scGW'`  | eigenvalues self-consistent G and W   | PRB 74, 045102 (2006)   |

        | `'qp-scGW0'`  | quasiparticle self-consistent G with fixed W0 | PRL 99, 115109 (2007) |

        | `'qp-scGW'`  | quasiparticle self-consistent G and W | PRL 96, 226402 (2006) |
        ''')

    analytical_continuation = Quantity(
        type=MEnum(
            'pade', 'contour_deformation', 'ppm_GodbyNeeds', 'ppm_HybertsenLouie',
            'ppm_vonderLindenHorsh', 'ppm_FaridEngel', 'multi_pole'),
        shape=[],
        description='''
        Analytical continuation approximations of the GW self-energy:

        | Name           | Description         | Reference                        |

        | -------------- | ------------------- | -------------------------------- |

        | `'pade'` | Pade's approximant  | J. Low Temp. Phys 29, 179 (1977) |

        | `'contour_deformation'` | Contour deformation | PRB 67, 155208 (2003) |

        | `'ppm_GodbyNeeds'` | Godby-Needs plasmon-pole model | PRL 62, 1169 (1989) |

        | `'ppm_HybertsenLouie'` | Hybertsen and Louie plasmon-pole model | PRB 34, 5390 (1986) |

        | `'ppm_vonderLindenHorsh'` | von der Linden and P. Horsh plasmon-pole model | PRB 37, 8351 (1988) |

        | `'ppm_FaridEngel'` | Farid and Engel plasmon-pole model  | PRB 47, 15931 (1993) |

        | `'multi_pole'` | Multi-pole fitting  | PRL 74, 1827 (1995) |
        ''')

    interval_qp_corrections = Quantity(
        type=np.int32,
        shape=[2],
        description='''
        Band indices (in an interval) for which the GW quasiparticle corrections are
        calculated.
        ''')

    screening = SubSection(sub_section=Screening.m_def)


class BSE(ExcitedStateMethodology):
    '''
    Section containing the various parameters that define a BSE calculation.
    '''

    m_def = Section(validate=False)

    screening = SubSection(sub_section=Screening.m_def)

    core_hole = SubSection(sub_section=CoreHole.m_def)


class DMFT(MSection):
    '''
    Section containing the various parameters that define a DMFT calculation
    '''

    m_def = Section(validate=False)

    n_atoms_per_unit_cell = Quantity(
        type=np.int32,
        shape=[],
        description='''
        Number of atoms per unit cell.
        ''')

    n_correlated_orbitals = Quantity(
        type=np.int32,
        shape=['n_atoms_per_unit_cell'],
        description='''
        Number of correlated orbitals per atom in the unit cell.
        ''')

    n_correlated_electrons = Quantity(
        type=np.float64,
        shape=['n_atoms_per_unit_cell'],
        description='''
        Number of correlated electrons per atom in the unit cell.
        ''')

    inverse_temperature = Quantity(
        type=np.float64,
        unit='1/joule',
        shape=[],
        description='''
        Inverse temperature = 1/(kB*T).
        ''')

    magnetic_state = Quantity(
        type=MEnum('paramagnetic', 'ferromagnetic', 'antiferromagnetic'),
        shape=[],
        description='''
        Magnetic state in which the DMFT calculation is done:

        | Name                  | State                   |

        | --------------------- | ----------------------- |

        | `'paramagnetic'`      | paramagnetic state      |

        | `'ferromagnetic'`     | ferromagnetic state     |

        | `'antiferromagnetic'` | antiferromagnetic state |
        ''')

    impurity_solver = Quantity(
        type=MEnum(
            'CT-INT', 'CT-HYB', 'CT-AUX', 'ED', 'NRG', 'MPS', 'IPT', 'NCA', 'OCA',
            'slave_bosons', 'hubbard_I'),
        shape=[],
        description='''
        Impurity solver method used in the DMFT loop:

        | Name              | Reference                            |

        | ----------------- | ------------------------------------ |

        | `'CT-INT'`        | Rubtsov et al., JEPT Lett 80 (2004)  |

        | `'CT-HYB'`        | Werner et al., PRL 97 (2006)         |

        | `'CT-AUX'`        | Gull et al., EPL 82 (2008)           |

        | `'ED'`            | Caffarrel et al, PRL 72 (1994)       |

        | `'NRG'`           | Bulla et al., RMP 80 (2008)          |

        | `'MPS'`           | Ganahl et al., PRB 90 (2014)         |

        | `'IPT'`           | Georges et al., PRB 45 (1992)        |

        | `'NCA'`           | Pruschke et al., PRB 47 (1993)       |

        | `'OCA'`           | Pruschke et al., PRB 47 (1993)       |

        | `'slave_bosons'`  | Kotliar et al., PRL 57 (1986)        |

        | `'hubbard_I'`     | -                                    |
        ''')


class TBModel(Model):
    '''
    Section containing the parameters pertaining to a tight-binding calculation.
    '''

    m_def = Section(validate=False)

    hamiltonian = SubSection(sub_section=Interaction.m_def, repeats=True)

    overlap = SubSection(sub_section=Interaction.m_def, repeats=True)

    repulsion = SubSection(sub_section=Interaction.m_def, repeats=True)

    magnetic = SubSection(sub_section=Interaction.m_def, repeats=True)

    coulomb = SubSection(sub_section=Interaction.m_def, repeats=True)


class TB(MSection):
    '''
    Section containing the parameters pertaining to a tight-binding calculation.
    '''

    m_def = Section(validate=False)

    kind = Quantity(
        type=str,
        shape=[],
        description='''
        The kind of tight-binding model used. Can be orthogonal or non-orthogonal.
        ''')

    model = SubSection(sub_section=TBModel.m_def, repeats=True)


class NeighborSearching(MSection):
    '''
    Section containing the parameters for neighbor searching/lists during a molecular dynamics run.
    '''

    m_def = Section(validate=False)

    neighbor_update_frequency = Quantity(
        type=int,
        shape=[],
        description='''
        Number of timesteps between updating the neighbor list.
        ''')

    neighbor_update_cutoff = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='m',
        description='''
        The distance cutoff for determining the neighbor list.
        ''')


class ForceCalculations(MSection):
    '''
    Section containing the parameters for force calculations according to the referenced force field
    during a molecular dynamics run.
    '''

    m_def = Section(validate=False)

    vdw_cutoff = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='m',
        description='''
        Cutoff for calculating VDW forces.
        ''')

    coulomb_type = Quantity(
        type=MEnum('cutoff', 'ewald', 'multilevel_summation', 'particle_mesh_ewald',
                   'particle_particle_particle_mesh', 'reaction_field'),
        shape=[],
        description='''
        Method used for calculating long-ranged Coulomb forces.

        Allowed values are:

        | Barostat Name          | Description                               |

        | ---------------------- | ----------------------------------------- |

        | `""`                   | No thermostat               |

        | `"Cutoff"`          | Simple cutoff scheme. |

        | `"Ewald"` | Standard Ewald summation as described in any solid-state physics text. |

        | `"Multi-Level Summation"` |  D. Hardy, J.E. Stone, and K. Schulten,
        [Parallel. Comput. **35**, 164](https://doi.org/10.1016/j.parco.2008.12.005)|

        | `"Particle-Mesh-Ewald"`        | T. Darden, D. York, and L. Pedersen,
        [J. Chem. Phys. **98**, 10089 (1993)](https://doi.org/10.1063/1.464397) |

        | `"Particle-Particle Particle-Mesh"` | See e.g. Hockney and Eastwood, Computer Simulation Using Particles,
        Adam Hilger, NY (1989). |

        | `"Reaction-Field"` | J.A. Barker and R.O. Watts,
        [Mol. Phys. **26**, 789 (1973)](https://doi.org/10.1080/00268977300102101)|
        ''')

    coulomb_cutoff = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='m',
        description='''
        Cutoff for calculating short-ranged Coulomb forces.
        ''')

    neighbor_searching = SubSection(sub_section=NeighborSearching.m_def, repeats=False)


class ForceField(MSection):
    '''
    Section containing the parameters pertaining to a force field calculation.
    '''

    m_def = Section(validate=False)

    model = SubSection(sub_section=Model.m_def, repeats=True)

    force_calculations = SubSection(sub_section=ForceCalculations.m_def, repeats=False)


class Smearing(MSection):
    '''
    Section containing the parameters related to the smearing of the electronic density of
    states at the Fermi level.
    '''

    m_def = Section(validate=False)

    kind = Quantity(
        type=str,
        shape=[],
        description='''
        Specifies the kind of smearing on the electron occupation used to calculate the
        free energy (see energy_free)

        Valid values are:

        | Smearing kind             | Description                       |

        | ------------------------- | --------------------------------- |

        | `"empty"`                 | No smearing is applied            |

        | `"gaussian"`              | Gaussian smearing                 |

        | `"fermi"`                 | Fermi smearing                    |

        | `"marzari-vanderbilt"`    | Marzari-Vanderbilt smearing       |

        | `"methfessel-paxton"`     | Methfessel-Paxton smearing        |

        | `"tetrahedra"`            | Interpolation of state energies and occupations
        (ignores smearing_width) |
        ''')

    width = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Specifies the width of the smearing in energy for the electron occupation used to
        calculate the free energy (see energy_free).

        *NOTE:* Not all methods specified in smearing_kind uses this value.
        ''')


class Electronic(MSection):
    '''
    Section containing the parameters related to the electronic structure.
    '''

    m_def = Section(validate=False)

    spin_target = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Stores the target (user-imposed) value of the spin multiplicity $M=2S+1$, where
        $S$ is the total spin. It is an integer number. This value is not necessarily the
        value obtained at the end of the calculation. See spin_S2 for the converged value
        of the spin moment.
        ''')

    charge = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        unit='coulomb',
        description='''
        Stores the total charge of the system.
        ''')

    n_bands = Quantity(
        type=int,
        shape=[],
        description='''
        Specifies the number of bands used in the calculation.
        ''')

    n_spin_channels = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of spin channels.
        ''')

    n_electrons = Quantity(
        type=np.dtype(np.float64),
        shape=['n_spin_channels'],
        description='''
        Number of electrons in system
        ''')

    method = Quantity(
        type=str,
        shape=[],
        description='''
        Non-unique string identifying the used electronic structure method. It is not
        unique in the sense that two calculations with the same
        electronic structure method string may have not been performed with exactly the
        same method.
        ''')

    relativity_method = Quantity(
        type=MEnum(['scalar_relativistic', 'pseudo_scalar_relativistic', 'scalar_relativistic_atomic_ZORA']),
        shape=[],
        description='''
        Describes the relativistic treatment used for the calculation of the final energy
        and related quantities. If skipped or empty, no relativistic treatment is applied.
        ''')

    van_der_waals_method = Quantity(
        type=str,
        shape=[],
        description='''
        Describes the Van der Waals method. If skipped or an empty string is used, it
        means no Van der Waals correction is applied.

        Allowed values are:

        | Van der Waals method  | Description                               |

        | --------------------- | ----------------------------------------- |

        | `""`                  | No Van der Waals correction               |

        | `"TS"`                | A. Tkatchenko and M. Scheffler, [Phys. Rev. Lett.
        **102**, 073005 (2009)](http://dx.doi.org/10.1103/PhysRevLett.102.073005) |

        | `"OBS"`               | F. Ortmann, F. Bechstedt, and W. G. Schmidt, [Phys. Rev.
        B **73**, 205101 (2006)](http://dx.doi.org/10.1103/PhysRevB.73.205101) |

        | `"G06"`               | S. Grimme, [J. Comput. Chem. **27**, 1787
        (2006)](http://dx.doi.org/10.1002/jcc.20495) |

        | `"JCHS"`              | P. Jurečka, J. Černý, P. Hobza, and D. R. Salahub,
        [Journal of Computational Chemistry **28**, 555
        (2007)](http://dx.doi.org/10.1002/jcc.20570) |

        | `"MDB"`               | Many-body dispersion. A. Tkatchenko, R. A. Di Stasio Jr,
        R. Car, and M. Scheffler, [Physical Review Letters **108**, 236402
        (2012)](http://dx.doi.org/10.1103/PhysRevLett.108.236402) and A. Ambrosetti, A. M.
        Reilly, R. A. Di Stasio Jr, and A. Tkatchenko, [The Journal of Chemical Physics
        **140**, 18A508 (2014)](http://dx.doi.org/10.1063/1.4865104) |

        | `"XC"`                | The method to calculate the Van der Waals energy uses a
        non-local functional which is described in section_XC_functionals. |
        ''')

    smearing = SubSection(sub_section=Smearing.m_def)


class Method(ArchiveSection):
    '''
    Section containing the various parameters that define the theory and the
    approximations (convergence, thresholds, etc.) behind the calculation.
    '''

    m_def = Section(validate=False)

    stress_tensor_method = Quantity(
        type=str,
        shape=[],
        description='''
        Specifies the method used to calculate stress_tensor for, e.g., molecular dynamics
        and geometry optimization.

        The allowed values are:

        * numeric

        * analytic
        ''')

    starting_method_ref = Quantity(
        type=Reference(SectionProxy('Method')),
        shape=[],
        description='''
        Links the current section method to a section method containing the starting
        parameters.
        ''',
        categories=[FastAccess])

    core_method_ref = Quantity(
        type=Reference(SectionProxy('Method')),
        shape=[],
        description='''
        Links the current section method to a section method containing the core settings.
        ''',
        categories=[FastAccess])

    n_references = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of references to the current method.
        ''')

    methods_ref = Quantity(
        type=Reference(SectionProxy('Method')),
        shape=['n_references'],
        description='''
        Links the section method to other method sections. For instance, one calculation
        is a perturbation performed using a self-consistent field (SCF) calculation as
        starting point, or a simulated system is partitioned in regions with different but
        connected Hamiltonians (e.g., QM/MM, or a region treated via Kohn-Sham DFT
        embedded into a region treated via orbital-free DFT).
        ''',
        categories=[FastAccess])

    dft = SubSection(sub_section=DFT.m_def)

    projection = SubSection(sub_section=Projection.m_def)

    lattice_model_hamiltonian = SubSection(sub_section=LatticeModelHamiltonian.m_def, repeats=True)

    gw = SubSection(sub_section=GW.m_def)

    bse = SubSection(sub_section=BSE.m_def)

    dmft = SubSection(sub_section=DMFT.m_def)

    tb = SubSection(sub_section=TB.m_def)

    force_field = SubSection(sub_section=ForceField.m_def)

    k_mesh = SubSection(sub_section=KMesh.m_def)

    frequency_mesh = SubSection(sub_section=FrequencyMesh.m_def)

    time_mesh = SubSection(sub_section=TimeMesh.m_def)

    electronic = SubSection(sub_section=Electronic.m_def)

    scf = SubSection(sub_section=Scf.m_def)

    atom_parameters = SubSection(sub_section=AtomParameters.m_def, repeats=True, label_quantity='label')

    molecule_parameters = SubSection(sub_section=MoleculeParameters.m_def, repeats=True, label_quantity='label')

    basis_set = SubSection(sub_section=BasisSet.m_def, repeats=True, label_quantity='type')

    photon = SubSection(sub_section=Photon.m_def, repeats=True)

    core_hole = SubSection(sub_section=CoreHole.m_def, repeats=True)


m_package.__init_metainfo__()
