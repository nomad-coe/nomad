import numpy as np            # pylint: disable=unused-import
import typing                 # pylint: disable=unused-import
from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference
)

from nomad.datamodel.metainfo import public

m_package = Package(name='common', description='None')


class settings_atom_in_molecule(MCategory):
    '''
    Parameters of an atom within a molecule.
    '''


class settings_constraint(MCategory):
    '''
    Some parameters that describe a constraint
    '''


class settings_interaction(MCategory):
    '''
    Some parameters that describe a bonded interaction.
    '''


class soap_parameter(MCategory):
    '''
    A soap parameter
    '''


class response_context(MSection):
    '''
    The top level context containing the reponse to an api query, when using jsonAPI they
    are tipically in the meta part
    '''

    m_def = Section(validate=False)

    shortened_meta_info = Quantity(
        type=str,
        shape=[],
        description='''
        A meta info whose corresponding data has been shortened
        ''')

    section_response_message = SubSection(
        sub_section=SectionProxy('section_response_message'),
        repeats=True)


class section_atom_type(MSection):
    '''
    Section describing a type of atom in the system.
    '''

    m_def = Section(validate=False)

    atom_type_charge = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='coulomb',
        description='''
        Charge of the atom type.
        ''')

    atom_type_mass = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='kilogram',
        description='''
        Mass of the atom type.
        ''')

    atom_type_name = Quantity(
        type=str,
        shape=[],
        description='''
        Name (label) of the atom type.
        ''')


class section_constraint(MSection):
    '''
    Section describing a constraint between arbitrary atoms.
    '''

    m_def = Section(validate=False)

    constraint_atoms = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_constraints', 'number_of_atoms_per_constraint'],
        description='''
        List of the indexes involved in this constraint. The fist atom has index 1, the
        last number_of_topology_atoms.
        ''')

    constraint_kind = Quantity(
        type=str,
        shape=[],
        description='''
        Short and unique name for this constraint type. Valid names are described in the
        [constraint\\_kind wiki page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-
        info/wikis/metainfo/constraint-kind).
        ''')

    constraint_parameters = Quantity(
        type=typing.Any,
        shape=[],
        description='''
        Explicit constraint parameters for this kind of constraint (depending on the
        constraint type, some might be given implicitly through other means).
        ''')

    number_of_atoms_per_constraint = Quantity(
        type=int,
        shape=[],
        description='''
        Number of atoms involved in this constraint.
        ''')

    number_of_constraints = Quantity(
        type=int,
        shape=[],
        description='''
        Number of constraints of this type.
        ''')


class section_dft_plus_u_orbital(MSection):
    '''
    Section for DFT+U-settings of a single orbital
    '''

    m_def = Section(validate=False)

    dft_plus_u_orbital_atom = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        DFT+U-orbital setting: atom index (references index of atom_labels/atom_positions)
        ''')

    dft_plus_u_orbital_J = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        DFT+U-orbital setting: value J (exchange interaction)
        ''',
        categories=[public.energy_value])

    dft_plus_u_orbital_label = Quantity(
        type=str,
        shape=[],
        description='''
        DFT+U-orbital setting: orbital label (normally (n,l)), notation: '3d', '4f', ...
        ''')

    dft_plus_u_orbital_U_effective = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        DFT+U-orbital setting: value U_{effective} (U-J), if implementation uses it
        ''',
        categories=[public.energy_value])

    dft_plus_u_orbital_U = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        DFT+U-orbital setting: value U (on-site Coulomb interaction)
        ''',
        categories=[public.energy_value])


class section_excited_states(MSection):
    '''
    Excited states properties.
    '''

    m_def = Section(validate=False)

    excitation_energies = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_excited_states'],
        description='''
        Excitation energies.
        ''',
        categories=[public.energy_value])

    number_of_excited_states = Quantity(
        type=int,
        shape=[],
        description='''
        Number of excited states.
        ''')

    oscillator_strengths = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_excited_states'],
        description='''
        Excited states oscillator strengths.
        ''')

    transition_dipole_moments = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_excited_states', 3],
        description='''
        Transition dipole moments.
        ''')


class section_interaction(MSection):
    '''
    Section containing the description of a bonded interaction between arbitrary atoms.
    '''

    m_def = Section(validate=False)

    interaction_atoms = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_interactions', 'number_of_atoms_per_interaction'],
        description='''
        List of the indexes involved in this interaction. The fist atom has index 1, the
        last atom index number_of_topology_atoms.
        ''')

    interaction_kind = Quantity(
        type=str,
        shape=[],
        description='''
        Short and unique name for this interaction type. Valid names are described in the
        [interaction\\_kind wiki page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-
        info/wikis/metainfo/interaction-kind).
        ''')

    interaction_parameters = Quantity(
        type=typing.Any,
        shape=[],
        description='''
        Explicit interaction parameters for this kind of interaction (depending on the
        interaction_kind some might be given implicitly through other means).
        ''')

    number_of_atoms_per_interaction = Quantity(
        type=int,
        shape=[],
        description='''
        Number of atoms involved in this interaction.
        ''')

    number_of_interactions = Quantity(
        type=int,
        shape=[],
        description='''
        Number of interactions of this type.
        ''')


class section_method_basis_set(MSection):
    '''
    This section contains the definition of the basis sets that are defined independently
    of the atomic configuration.
    '''

    m_def = Section(validate=False)

    mapping_section_method_basis_set_atom_centered = Quantity(
        type=np.dtype(np.int64),
        shape=['number_of_basis_sets_atom_centered', 2],
        description='''
        Reference to an atom-centered basis set defined in section_basis_set_atom_centered
        and to the atom kind as defined in section_method_atom_kind.
        ''')

    mapping_section_method_basis_set_cell_associated = Quantity(
        type=public.section_basis_set_cell_dependent,
        shape=[],
        description='''
        Reference to a cell-associated basis set.
        ''')

    method_basis_set_kind = Quantity(
        type=str,
        shape=[],
        description='''
        String describing the use of the basis set, i.e, if it used for expanding a
        wavefunction or an electron density. Allowed values are listed in the
        [basis\\_set\\_kind wiki page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-
        info/wikis/metainfo/basis-set-kind).
        ''')

    number_of_basis_sets_atom_centered = Quantity(
        type=int,
        shape=[],
        description='''
        String describing the use of the basis set, i.e, if it used for expanding a
        wavefunction or an electron density. Allowed values are listed in the
        [basis\\_set\\_kind wiki page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-
        info/wikis/metainfo/basis-set-kind).
        ''')


class section_molecule_constraint(MSection):
    '''
    Section describing a constraint between atoms within a molecule.
    '''

    m_def = Section(validate=False)

    molecule_constraint_atoms = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_molecule_constraints', 'number_of_atoms_per_molecule_constraint'],
        description='''
        List of the indexes involved in this constraint. The fist atom has index 1, the
        last index is number_of_atoms_in_molecule.
        ''')

    molecule_constraint_kind = Quantity(
        type=str,
        shape=[],
        description='''
        Short and unique name for this constraint type. Valid names are described in the
        [constraint\\_kind wiki page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-
        info/wikis/metainfo/constraint-kind).
        ''')

    molecule_constraint_parameters = Quantity(
        type=typing.Any,
        shape=[],
        description='''
        Explicit constraint parameters for this kind of constraint (depending on the
        constraint type some might be given implicitly through other means).
        ''')

    number_of_atoms_per_molecule_constraint = Quantity(
        type=int,
        shape=[],
        description='''
        Number of atoms, in this molecule, involved in this constraint.
        ''')

    number_of_molecule_constraints = Quantity(
        type=int,
        shape=[],
        description='''
        Number of constraints of this type.
        ''')


class section_molecule_interaction(MSection):
    '''
    Section describing a bonded interaction between atoms within a molecule.
    '''

    m_def = Section(validate=False)

    molecule_interaction_atoms = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_molecule_interactions', 'number_of_atoms_per_molecule_interaction'],
        description='''
        List of the indexes involved in this bonded interaction within a molecule. The
        first atom has index 1, the last index is number_of_atoms_in_.
        ''')

    molecule_interaction_kind = Quantity(
        type=str,
        shape=[],
        description='''
        Short and unique name for this interaction type, used for bonded interactions for
        atoms in a molecule. Valid names are described in the [interaction\\_kind wiki
        page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-
        info/wikis/metainfo/interaction-kind).
        ''')

    molecule_interaction_parameters = Quantity(
        type=typing.Any,
        shape=[],
        description='''
        Explicit interaction parameters for this kind of interaction (depending on the
        interaction type some might be given implicitly through other means), used for
        bonded interactions for atoms in a molecule.
        ''')

    number_of_atoms_per_molecule_interaction = Quantity(
        type=int,
        shape=[],
        description='''
        Number of atoms, in this molecule, involved in this interaction.
        ''')

    number_of_molecule_interactions = Quantity(
        type=int,
        shape=[],
        description='''
        Number of bonded interactions of this type.
        ''')


class section_molecule_type(MSection):
    '''
    Section describing a type of molecule in the system.
    '''

    m_def = Section(validate=False)

    atom_in_molecule_charge = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms_in_molecule'],
        unit='coulomb',
        description='''
        Charge of each atom in the molecule.
        ''',
        categories=[settings_atom_in_molecule])

    atom_in_molecule_name = Quantity(
        type=str,
        shape=['number_of_atoms_in_molecule'],
        description='''
        Name (label) of each atom in the molecule.
        ''',
        categories=[settings_atom_in_molecule])

    atom_in_molecule_to_atom_type_ref = Quantity(
        type=Reference(SectionProxy('section_atom_type')),
        shape=['number_of_atoms_in_molecule'],
        description='''
        Reference to the atom type of each atom in the molecule.
        ''',
        categories=[settings_atom_in_molecule])

    molecule_type_name = Quantity(
        type=str,
        shape=[],
        description='''
        Name of the molecule.
        ''')

    number_of_atoms_in_molecule = Quantity(
        type=int,
        shape=[],
        description='''
        Number of atoms in this molecule.
        ''')

    section_molecule_constraint = SubSection(
        sub_section=SectionProxy('section_molecule_constraint'),
        repeats=True)

    section_molecule_interaction = SubSection(
        sub_section=SectionProxy('section_molecule_interaction'),
        repeats=True)


class section_response_message(MSection):
    '''
    Messages outputted by the program formatting the data in the current response
    '''

    m_def = Section(validate=False)

    response_message_count = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        How many times this message was repeated
        ''')

    response_message_level = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        level of the message: 0 fatal, 1 error, 2 warning, 3 debug
        ''')

    response_message = Quantity(
        type=str,
        shape=[],
        description='''
        Message outputted by the program formatting the data in the current format
        ''')


class section_soap_coefficients(MSection):
    '''
    Stores the soap coefficients for the pair of atoms given in
    soap_coefficients_atom_pair.
    '''

    m_def = Section(validate=False)

    number_of_soap_coefficients = Quantity(
        type=int,
        shape=[],
        description='''
        number of soap coefficients
        ''')

    soap_coefficients_atom_pair = Quantity(
        type=str,
        shape=[],
        description='''
        Pair of atoms described in the current section
        ''')

    soap_coefficients = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_soap_coefficients'],
        description='''
        Compressed coefficient of the soap descriptor for the atom pair
        soap_coefficients_atom_pair
        ''')


class section_soap(MSection):
    '''
    Stores a soap descriptor for this configuration.
    '''

    m_def = Section(validate=False)

    soap_angular_basis_L = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        angular basis L
        ''',
        categories=[soap_parameter])

    soap_angular_basis_type = Quantity(
        type=str,
        shape=[],
        description='''
        angular basis type
        ''',
        categories=[soap_parameter])

    soap_kernel_adaptor = Quantity(
        type=str,
        shape=[],
        description='''
        kernel adaptor
        ''',
        categories=[soap_parameter])

    soap_parameters_gid = Quantity(
        type=str,
        shape=[],
        description='''
        Unique checksum of all the soap parameters (all those with abstract type
        soap_parameter) with prefix psoap
        ''')

    soap_radial_basis_integration_steps = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        radial basis integration steps
        ''',
        categories=[soap_parameter])

    soap_radial_basis_mode = Quantity(
        type=str,
        shape=[],
        description='''
        radial basis mode
        ''',
        categories=[soap_parameter])

    soap_radial_basis_n = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        radial basis N
        ''',
        categories=[soap_parameter])

    soap_radial_basis_sigma = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        radial basis sigma
        ''',
        categories=[soap_parameter])

    soap_radial_basis_type = Quantity(
        type=str,
        shape=[],
        description='''
        radial basis type
        ''',
        categories=[soap_parameter])

    soap_radial_cutoff_center_weight = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        radial cutoff center weight
        ''',
        categories=[soap_parameter])

    soap_radial_cutoff_rc_width = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        radial cutoff width
        ''',
        categories=[soap_parameter])

    soap_radial_cutoff_rc = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        radial cutoff
        ''',
        categories=[soap_parameter])

    soap_radial_cutoff_type = Quantity(
        type=str,
        shape=[],
        description='''
        radial cutoff type
        ''',
        categories=[soap_parameter])

    soap_spectrum_2l1_norm = Quantity(
        type=bool,
        shape=[],
        description='''
        2l1 norm spectrum
        ''',
        categories=[soap_parameter])

    soap_spectrum_global = Quantity(
        type=bool,
        shape=[],
        description='''
        global spectrum
        ''',
        categories=[soap_parameter])

    soap_spectrum_gradients = Quantity(
        type=bool,
        shape=[],
        description='''
        gradients in specturm
        ''',
        categories=[soap_parameter])

    soap_type_list = Quantity(
        type=str,
        shape=[],
        description='''
        Type list
        ''',
        categories=[soap_parameter])

    section_soap_coefficients = SubSection(
        sub_section=SectionProxy('section_soap_coefficients'),
        repeats=True)


class section_topology(MSection):
    '''
    Section containing the definition of topology (connectivity among atoms in force
    fileds), force field, and constraints of a system.
    '''

    m_def = Section(validate=False)

    atom_to_molecule = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_topology_atoms', 2],
        description='''
        Table mapping atom to molecules: the first column is the index of the molecule and
        the second column the index of the atom, signifying that the atom in the second
        column belongs to the molecule in the first column in the same row.
        ''')

    molecule_to_molecule_type_map = Quantity(
        type=Reference(SectionProxy('section_molecule_type')),
        shape=['number_of_topology_molecules'],
        description='''
        Mapping from molecules to molecule types.
        ''')

    number_of_topology_atoms = Quantity(
        type=int,
        shape=[],
        description='''
        Number of atoms in the system described by this topology.
        ''')

    number_of_topology_molecules = Quantity(
        type=int,
        shape=[],
        description='''
        Number of molecules in the system, as described by this topology.
        ''')

    topology_force_field_name = Quantity(
        type=str,
        shape=[],
        description='''
        A unique string idenfiying the force field defined in this section. Strategies to
        define it are discussed in the
        [topology\\_force\\_field\\_name](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-
        info/wikis/metainfo/topology-force-field-name).
        ''')

    section_atom_type = SubSection(
        sub_section=SectionProxy('section_atom_type'),
        repeats=True)

    section_constraint = SubSection(
        sub_section=SectionProxy('section_constraint'),
        repeats=True)

    section_interaction = SubSection(
        sub_section=SectionProxy('section_interaction'),
        repeats=True)

    section_molecule_type = SubSection(
        sub_section=SectionProxy('section_molecule_type'),
        repeats=True)


class section_method(public.section_method):

    m_def = Section(validate=False, extends_base_section=True)

    dft_plus_u_functional = Quantity(
        type=str,
        shape=[],
        description='''
        Type of DFT+U functional (such as DFT/DFT+U double-counting compensation). Valid
        names are described in the [dft\\_plus\\_u\\_functional wiki
        page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/dft-
        plus-u-functional).
        ''')

    dft_plus_u_projection_type = Quantity(
        type=str,
        shape=[],
        description='''
        DFT+U: Type of orbitals used for projection in order to calculate occupation
        numbers. Valid names are described in the [dft\\_plus\\_u\\_projection\\_type wiki
        page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/dft-
        plus-u-projection-type).
        ''')

    gw_bare_coulomb_cutofftype = Quantity(
        type=str,
        shape=[],
        description='''
        Cutoff type for the calculation of the bare Coulomb potential: none, 0d, 1d, 2d.
        See Rozzi et al., PRB 73, 205119 (2006)
        ''')

    gw_bare_coulomb_gmax = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='1 / meter',
        description='''
        Maximum G for the pw basis for the Coulomb potential.
        ''')

    gw_basis_set = Quantity(
        type=str,
        shape=[],
        description='''
        Auxillary basis set used for non-local operators: mixed - mixed basis set, Kotani
        and Schilfgaarde, Solid State Comm. 121, 461 (2002).
        ''')

    gw_core_treatment = Quantity(
        type=str,
        shape=[],
        description='''
        It specifies whether the core states are treated in the GW calculation: all - All
        electron calculation; val - Valence electron only calculation; vab - Core
        electrons are excluded from the mixed product basis; xal - All electron treatment
        of the exchange self-energy only
        ''')

    gw_frequency_grid_type = Quantity(
        type=str,
        shape=[],
        description='''
        Frequency integration grid type for the correlational self energy: 'eqdis' -
        equidistant frequencies from 0 to freqmax; 'gaulag' - Gauss-Laguerre quadrature
        from 0 to infinity; 'gauleg' - Gauss-Legendre quadrature from 0 to freqmax;
        'gaule2' (default) - double Gauss-Legendre quadrature from 0 to freqmax and from
        freqmax to infinity.
        ''')

    gw_max_frequency = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Maximum frequency for the calculation of the self energy.
        ''')

    gw_mixed_basis_gmax = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='1 / meter',
        description='''
        Cut-off parameter for the truncation of the expansion of the plane waves in the
        interstitial region.
        ''')

    gw_mixed_basis_lmax = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Maximum l value used for the radial functions within the muffin-tin.
        ''')

    gw_mixed_basis_tolerance = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Eigenvalue threshold below which the egenvectors are discarded in the construction
        of the radial basis set.
        ''')

    gw_ngridq = Quantity(
        type=np.dtype(np.int32),
        shape=[3],
        description='''
        k/q-point grid size used in the GW calculation.
        ''')

    gw_frequency_number = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number referring to the frequency used in the calculation of the self energy.
        ''')

    gw_frequency_values = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Values of the frequency used in the calculation of the self energy.
        ''')

    gw_frequency_weights = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Weights of the frequency used in the calculation of the self energy.
        ''')

    gw_number_of_frequencies = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of frequency points used in the calculation of the self energy.
        ''')

    gw_polarizability_number_of_empty_states = Quantity(
        type=int,
        shape=[],
        description='''
        Number of empty states used to compute the polarizability P
        ''')

    gw_qp_equation_treatment = Quantity(
        type=str,
        shape=[],
        description='''
        Methods to solve the quasi-particle equation: 'linearization', 'self-consistent'
        ''')

    gw_screened_coulomb_volume_average = Quantity(
        type=str,
        shape=[],
        description='''
        Type of volume averaging for the dynamically screened Coulomb potential: isotropic
        - Simple averaging along a specified direction using only diagonal components of
        the dielectric tensor; anisotropic - Anisotropic screening by C. Freysoldt et al.,
        CPC 176, 1-13 (2007)
        ''')

    gw_screened_Coulomb = Quantity(
        type=str,
        shape=[],
        description='''
        Model used to calculate the dinamically-screened Coulomb potential: 'rpa' - Full-
        frequency random-phase approximation; 'ppm' - Godby-Needs plasmon-pole model Godby
        and Needs, Phys. Rev. Lett. 62, 1169 (1989); 'ppm_hl' - Hybertsen and Louie, Phys.
        Rev. B 34, 5390 (1986); 'ppm_lh' - von der Linden and P. Horsh, Phys. Rev. B 37,
        8351 (1988); 'ppm_fe' - Farid and Engel, Phys. Rev. B 47,15931 (1993); 'cdm' -
        Contour deformation method, Phys. Rev. B 67, 155208 (2003).)
        ''')

    gw_self_energy_c_analytical_continuation = Quantity(
        type=str,
        shape=[],
        description='''
        Models for the correlation self-energy analytical continuation: 'pade' -  Pade's
        approximant (by H. J. Vidberg and J. W. Serence, J. Low Temp. Phys. 29, 179
        (1977)); 'mpf' -  Multi-Pole Fitting (by H. N Rojas, R. W. Godby and R. J. Needs,
        Phys. Rev. Lett. 74, 1827 (1995)); 'cd' - contour deformation; 'ra' - real axis
        ''')

    gw_self_energy_c_number_of_empty_states = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of empty states to be used to calculate the correlation self energy.
        ''')

    gw_self_energy_c_number_of_poles = Quantity(
        type=int,
        shape=[],
        description='''
        Number of poles used in the analytical continuation.
        ''')

    gw_self_energy_singularity_treatment = Quantity(
        type=str,
        shape=[],
        description='''
        Treatment of the integrable singular terms in the calculation of the self energy.
        Values: 'mpb' - Auxiliary function method by S. Massidda, M. Posternak, and A.
        Baldereschi, PRB 48, 5058 (1993); 'crg' - Auxiliary function method by P. Carrier,
        S. Rohra, and A. Goerling, PRB 75, 205126 (2007).
        ''')

    gw_starting_point = Quantity(
        type=str,
        shape=[],
        description='''
        Exchange-correlation functional of the ground-state calculation. See XC_functional
        list at https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/XC-
        functional
        ''')

    gw_type_test = Quantity(
        type=str,
        shape=[],
        description='''
        GW methodology: exciting test variable
        ''')

    gw_type = Quantity(
        type=str,
        shape=[],
        description='''
        GW methodology: G0W0; ev-scGW: (eigenvalues self-consistent GW) – Phys.Rev.B 34,
        5390 (1986); qp-scGW: (quasi-particle self-consistent GW) – Phys. Rev. Lett. 96,
        226402 (2006)  scGW0: (self-consistent G with fixed W0) – Phys.Rev.B 54, 8411
        (1996); scG0W: (self-consistent W with fixed G0); scGW: (self-consistent GW) –
        Phys. Rev. B 88, 075105 (2013)
        ''')

    method_to_topology_ref = Quantity(
        type=Reference(SectionProxy('section_topology')),
        shape=[],
        description='''
        Reference to the topology and force fields to be used.
        ''')

    section_dft_plus_u_orbital = SubSection(
        sub_section=SectionProxy('section_dft_plus_u_orbital'),
        repeats=True)

    section_method_basis_set = SubSection(
        sub_section=SectionProxy('section_method_basis_set'),
        repeats=True)


class section_single_configuration_calculation(public.section_single_configuration_calculation):

    m_def = Section(validate=False, extends_base_section=True)

    energy_C_mGGA = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Component of the correlation (C) energy at the GGA (or MetaGGA) level using the
        self-consistent density of the target XC functional (full unscaled value, i.e.,
        not scaled due to exact-exchange mixing).
        ''',
        categories=[public.energy_value, public.energy_type_C, public.energy_component])

    energy_reference_fermi = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_spin_channels'],
        unit='joule',
        description='''
        Fermi energy (separates occupied from unoccupied single-particle states in metals)
        ''',
        categories=[public.energy_type_reference, public.energy_value])

    energy_reference_highest_occupied = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_spin_channels'],
        unit='joule',
        description='''
        Highest occupied single-particle state energy (in insulators or HOMO energy in
        finite systems)
        ''',
        categories=[public.energy_type_reference, public.energy_value])

    energy_reference_lowest_unoccupied = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_spin_channels'],
        unit='joule',
        description='''
        Lowest unoccupied single-particle state energy (in insulators or LUMO energy in
        finite systems)
        ''',
        categories=[public.energy_type_reference, public.energy_value])

    energy_X_mGGA_scaled = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Component of the exchange (X) energy at the GGA (or MetaGGA) level, using the self
        consistent density of the target functional, scaled accordingly to the mixing
        parameter.
        ''',
        categories=[public.energy_value, public.energy_component])

    energy_X_mGGA = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Component of the exchange (X) energy at the GGA (or MetaGGA) level using the self
        consistent density of the target functional (full unscaled value, i.e., not scaled
        due to exact-exchange mixing).
        ''',
        categories=[public.energy_type_X, public.energy_value, public.energy_component])

    gw_fermi_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        GW Fermi energy
        ''')

    gw_fundamental_gap = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        GW fundamental band gap
        ''')

    gw_optical_gap = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        GW optical band gap
        ''')

    gw_self_energy_c = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_spin_channels', 'number_of_eigenvalues_kpoints', 'number_of_eigenvalues'],
        unit='joule',
        description='''
        Diagonal matrix elements of the correlation self-energy
        ''')

    gw_self_energy_x = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_spin_channels', 'number_of_eigenvalues_kpoints', 'number_of_eigenvalues'],
        unit='joule',
        description='''
        Diagonal matrix elements of the exchange self-energy
        ''')

    gw_xc_potential = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_spin_channels', 'number_of_eigenvalues_kpoints', 'number_of_eigenvalues'],
        unit='joule',
        description='''
        Diagonal matrix elements of the exchange-correlation potential
        ''')

    section_excited_states = SubSection(
        sub_section=SectionProxy('section_excited_states'),
        repeats=True)


class section_scf_iteration(public.section_scf_iteration):

    m_def = Section(validate=False, extends_base_section=True)

    energy_reference_fermi_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_spin_channels'],
        unit='joule',
        description='''
        Fermi energy (separates occupied from unoccupied single-particle states in metals)
        during the self-consistent field (SCF) iterations.
        ''',
        categories=[public.energy_type_reference, public.energy_value])

    energy_reference_highest_occupied_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_spin_channels'],
        unit='joule',
        description='''
        Highest occupied single-particle state energy (in insulators or HOMO energy in
        finite systems) during the self-consistent field (SCF) iterations.
        ''',
        categories=[public.energy_type_reference, public.energy_value])

    energy_reference_lowest_unoccupied_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_spin_channels'],
        unit='joule',
        description='''
        Lowest unoccupied single-particle state energy (in insulators or LUMO energy in
        finite systems) during the self-consistent field (SCF) iterations.
        ''',
        categories=[public.energy_type_reference, public.energy_value])


class section_eigenvalues(public.section_eigenvalues):

    m_def = Section(validate=False, extends_base_section=True)

    gw_qp_linearization_prefactor = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_spin_channels', 'number_of_eigenvalues_kpoints', 'number_of_eigenvalues'],
        description='''
        Linearization prefactor
        ''')


class section_system(public.section_system):

    m_def = Section(validate=False, extends_base_section=True)

    number_of_electrons = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_spin_channels'],
        description='''
        Number of electrons in system
        ''',
        categories=[public.configuration_core])

    topology_ref = Quantity(
        type=Reference(SectionProxy('section_topology')),
        shape=[],
        description='''
        Reference to the topology used for this system; if not given, the trivial topology
        should be assumed.
        ''')

    is_representative = Quantity(
        type=bool,
        shape=[],
        description='''
        Most systems in a run are only minor variations of each other. Systems marked
        representative where chosen to be representative for all systems in the run.
        ''')

    section_soap = SubSection(
        sub_section=SectionProxy('section_soap'),
        repeats=True)


class section_run(public.section_run):

    m_def = Section(validate=False, extends_base_section=True)

    section_topology = SubSection(
        sub_section=SectionProxy('section_topology'),
        repeats=True)


m_package.__init_metainfo__()
