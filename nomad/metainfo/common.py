import numpy as np
import typing
from nomad.metainfo import MSection, Package, Quantity, SubSection, MProxy

m_package = Package(name='common', description='None')


class section_atom_type(MSection):
    '''
    Reference to the atom type of each atom in the molecule.
    '''

    atom_type_charge = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='coulomb',
        description='''
        Charge of the atom type.
        ''')

    atom_type_mass = Quantity(
        type=np.dtype(np.float32),
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


class section_topology(MSection):
    '''
    Table mapping atom to molecules: the first column is the index of the molecule and the
    second column the index of the atom, signifying that the atom in the second column
    belongs to the molecule in the first column in the same row.
    '''

    atom_to_molecule = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_topology_atoms', 2],
        description='''
        Table mapping atom to molecules: the first column is the index of the molecule and
        the second column the index of the atom, signifying that the atom in the second
        column belongs to the molecule in the first column in the same row.
        ''')

    molecule_to_molecule_type_map = Quantity(
        type=MProxy('section_molecule_type'),
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
        sub_section=MProxy('section_atom_type'),
        repeats=False)

    section_constraint = SubSection(
        sub_section=MProxy('section_constraint'),
        repeats=False)

    section_interaction = SubSection(
        sub_section=MProxy('section_interaction'),
        repeats=False)

    section_molecule_type = SubSection(
        sub_section=MProxy('section_molecule_type'),
        repeats=False)


class section_constraint(MSection):
    '''
    List of the indexes involved in this constraint. The fist atom has index 1, the last
    number_of_topology_atoms.
    '''

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


class section_method(MSection):
    '''
    Type of DFT+U functional (such as DFT/DFT+U double-counting compensation). Valid names
    are described in the [dft\\_plus\\_u\\_functional wiki
    page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/dft-plus-u-
    functional).
    '''

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
        type=np.dtype(np.float32),
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
        type=np.dtype(np.float32),
        shape=[],
        description='''
        Maximum frequency for the calculation of the self energy.
        ''')

    gw_mixed_basis_gmax = Quantity(
        type=np.dtype(np.float32),
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
        type=np.dtype(np.float32),
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
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Values of the frequency used in the calculation of the self energy.
        ''')

    gw_frequency_weights = Quantity(
        type=np.dtype(np.float32),
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
        type=MProxy('section_topology'),
        shape=[],
        description='''
        Reference to the topology and force fields to be used.
        ''')

    calculation_method_current = Quantity(
        type=str,
        shape=[],
        description='''
        String that represents the method used to calculate the energy_current. If the
        method is perturbative, this string does not describe the starting point method,
        the latter being referenced to by section_method_to_method_refs. For self-
        consistent field (SCF) ab initio calculations, for example, this is composed by
        concatenating XC_method_current and basis_set. See [calculation_method_current
        wiki page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-
        info/wikis/metainfo/calculation-method-current) for the details.
        ''')

    calculation_method_kind = Quantity(
        type=str,
        shape=[],
        description='''
        Kind of method in calculation_method_current.

        Accepted values are:

        - absolute

        - perturbative.
        ''')

    calculation_method = Quantity(
        type=str,
        shape=[],
        description='''
        String that uniquely represents the method used to calculate energy_total, If the
        present calculation_method_current is a perturbative method Y that uses method X
        as starting point, this string is automatically created as X@Y, where X is taken
        from calculation_method_current and Y from method_to_method_ref. In order to
        activate this, method_to_method_kind must have the value starting_point (see the
        [method_to_method_kind wiki page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-
        meta-info/wikis/metainfo/method-to-method-kind)).
        ''')

    number_of_spin_channels = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of spin channels, see section_method.
        ''')

    spin_target_multiplicity = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Stores the target (user-imposed) value of the spin multiplicity $M=2S+1$, where
        $S$ is the total spin. It is an integer number. This value is not necessarily the
        value obtained at the end of the calculation. See spin_S2 for the converged value
        of the spin moment.
        ''')

    total_charge = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        unit='coulomb',
        description='''
        Provides the total amount of charge of the system in a run.
        ''')

    section_dft_plus_u_orbital = SubSection(
        sub_section=MProxy('section_dft_plus_u_orbital'),
        repeats=False)

    section_method_basis_set = SubSection(
        sub_section=MProxy('section_method_basis_set'),
        repeats=False)

    section_method_atom_kind = SubSection(
        sub_section=MProxy('section_method_atom_kind'),
        repeats=False)

    section_method_to_method_refs = SubSection(
        sub_section=MProxy('section_method_to_method_refs'),
        repeats=False)


class section_dft_plus_u_orbital(MSection):
    '''
    DFT+U-orbital setting: atom index (references index of atom_labels/atom_positions)
    '''

    dft_plus_u_orbital_atom = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        DFT+U-orbital setting: atom index (references index of atom_labels/atom_positions)
        ''')

    dft_plus_u_orbital_J = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        description='''
        DFT+U-orbital setting: value J (exchange interaction)
        ''')

    dft_plus_u_orbital_label = Quantity(
        type=str,
        shape=[],
        description='''
        DFT+U-orbital setting: orbital label (normally (n,l)), notation: '3d', '4f', ...
        ''')

    dft_plus_u_orbital_U_effective = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        description='''
        DFT+U-orbital setting: value U_{effective} (U-J), if implementation uses it
        ''')

    dft_plus_u_orbital_U = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        description='''
        DFT+U-orbital setting: value U (on-site Coulomb interaction)
        ''')


class section_scf_iteration(MSection):
    '''
    Fermi energy (separates occupied from unoccupied single-particle states in metals)
    during the self-consistent field (SCF) iterations.
    '''

    energy_reference_fermi_iteration = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_spin_channels'],
        unit='joule',
        description='''
        Fermi energy (separates occupied from unoccupied single-particle states in metals)
        during the self-consistent field (SCF) iterations.
        ''')

    energy_reference_highest_occupied_iteration = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_spin_channels'],
        unit='joule',
        description='''
        Highest occupied single-particle state energy (in insulators or HOMO energy in
        finite systems) during the self-consistent field (SCF) iterations.
        ''')

    energy_reference_lowest_unoccupied_iteration = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_spin_channels'],
        unit='joule',
        description='''
        Lowest unoccupied single-particle state energy (in insulators or LUMO energy in
        finite systems) during the self-consistent field (SCF) iterations.
        ''')

    electronic_kinetic_energy_scf_iteration = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Electronic kinetic energy as defined in XC_method during the self-consistent field
        (SCF) iterations.
        ''')

    energy_change_scf_iteration = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Stores the change of total energy with respect to the previous self-consistent
        field (SCF) iteration.
        ''')

    energy_correction_entropy_scf_iteration = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Entropy correction to the potential energy to compensate for the change in
        occupation so that forces at finite T do not need to keep the change of occupation
        in account. The array lists the values of the entropy correction for each self-
        consistent field (SCF) iteration. Defined consistently with XC_method.
        ''')

    energy_correction_hartree_scf_iteration = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Correction to the density-density electrostatic energy in the sum of eigenvalues
        (that uses the mixed density on one side), and the fully consistent density-
        density electrostatic energy during the self-consistent field (SCF) iterations.
        Defined consistently with XC_method.
        ''')

    energy_electrostatic_scf_iteration = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Total electrostatic energy (nuclei + electrons) during each self-consistent field
        (SCF) iteration.
        ''')

    energy_free_per_atom_scf_iteration = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Free energy per atom (whose minimum gives the smeared occupation density
        calculated with smearing_kind) calculated with XC_method during the self-
        consistent field (SCF) iterations.
        ''')

    energy_free_scf_iteration = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Free energy (whose minimum gives the smeared occupation density calculated with
        smearing_kind) calculated with the method described in XC_method during the self-
        consistent field (SCF) iterations.
        ''')

    energy_hartree_error_scf_iteration = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Error in the Hartree (electrostatic) potential energy during each self-consistent
        field (SCF) iteration. Defined consistently with XC_method.
        ''')

    energy_sum_eigenvalues_per_atom_scf_iteration = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Value of the energy per atom, where the energy is defined as the sum of the
        eigenvalues of the Hamiltonian matrix given by XC_method, during each self-
        consistent field (SCF) iteration.
        ''')

    energy_sum_eigenvalues_scf_iteration = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Sum of the eigenvalues of the Hamiltonian matrix defined by XC_method, during each
        self-consistent field (SCF) iteration.
        ''')

    energy_total_scf_iteration = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Value of the total electronic energy calculated with the method described in
        XC_method during each self-consistent field (SCF) iteration.
        ''')

    energy_total_T0_per_atom_scf_iteration = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Value of the total energy, calculated with the method described in XC_method per
        atom extrapolated to $T=0$, based on a free-electron gas argument, during each
        self-consistent field (SCF) iteration.
        ''')

    energy_total_T0_scf_iteration = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Value of the total energy (or equivalently free energy), calculated with the
        method described in XC_method and extrapolated to $T=0$, based on a free-electron
        gas argument, during each self-consistent field (SCF) iteration.
        ''')

    energy_XC_potential_scf_iteration = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Value for exchange-correlation (XC) potential energy: the integral of the first
        order derivative of the functional stored in XC_functional (integral of
        v_xc*electron_density), i.e., the component of XC that is in the sum of the
        eigenvalues. Values are given for each self-consistent field (SCF) iteration
        (i.e., not the converged value, the latter being stored in energy_XC_potential).
        ''')

    energy_XC_scf_iteration = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Value for exchange-correlation (XC) energy obtained during each self-consistent
        field (SCF) iteration, using the method described in XC_method.
        ''')

    spin_S2_scf_iteration = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        description='''
        Stores the value of the total spin moment operator $S^2$ during the self-
        consistent field (SCF) iterations of the XC_method. It can be used to calculate
        the spin contamination in spin-unrestricted calculations.
        ''')

    time_scf_iteration_cpu1_end = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='second',
        description='''
        Stores the end time of a self-consistent field (SCF) iteration on CPU 1.
        ''')

    time_scf_iteration_cpu1_start = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='second',
        description='''
        Stores the start time of a self-consistent field (SCF) iteration on CPU 1.
        ''')

    time_scf_iteration_date_end = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='second',
        description='''
        Stores the end date of a self-consistent field (SCF) iteration as time since the
        *Unix epoch* (00:00:00 UTC on 1 January 1970) in seconds. For date and times
        without a timezone, the default timezone GMT is used.
        ''')

    time_scf_iteration_date_start = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='second',
        description='''
        Stores the start date of a self-consistent field (SCF) iteration as time since the
        *Unix epoch* (00:00:00 UTC on 1 January 1970) in seconds. For date and times
        without a timezone, the default timezone GMT is used.
        ''')

    time_scf_iteration_wall_end = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='second',
        description='''
        Stores the internal wall-clock time at the end of a self-consistent field (SCF)
        iteration.
        ''')

    time_scf_iteration_wall_start = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='second',
        description='''
        Stores the internal wall-clock time from the start of a self-consistent field
        (SCF) iteration.
        ''')


class section_single_configuration_calculation(MSection):
    '''
    Fermi energy (separates occupied from unoccupied single-particle states in metals)
    '''

    energy_reference_fermi = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_spin_channels'],
        unit='joule',
        description='''
        Fermi energy (separates occupied from unoccupied single-particle states in metals)
        ''')

    energy_reference_highest_occupied = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_spin_channels'],
        unit='joule',
        description='''
        Highest occupied single-particle state energy (in insulators or HOMO energy in
        finite systems)
        ''')

    energy_reference_lowest_unoccupied = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_spin_channels'],
        unit='joule',
        description='''
        Lowest unoccupied single-particle state energy (in insulators or LUMO energy in
        finite systems)
        ''')

    energy_X_mGGA_scaled = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Component of the exchange (X) energy at the GGA (or MetaGGA) level, using the self
        consistent density of the target functional, scaled accordingly to the mixing
        parameter.
        ''')

    gw_fermi_energy = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        GW Fermi energy
        ''')

    gw_fundamental_gap = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        GW fundamental band gap
        ''')

    gw_optical_gap = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        GW optical band gap
        ''')

    gw_self_energy_c = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_spin_channels', 'number_of_eigenvalues_kpoints', 'number_of_eigenvalues'],
        unit='joule',
        description='''
        Diagonal matrix elements of the correlation self-energy
        ''')

    gw_self_energy_x = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_spin_channels', 'number_of_eigenvalues_kpoints', 'number_of_eigenvalues'],
        unit='joule',
        description='''
        Diagonal matrix elements of the exchange self-energy
        ''')

    gw_xc_potential = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_spin_channels', 'number_of_eigenvalues_kpoints', 'number_of_eigenvalues'],
        unit='joule',
        description='''
        Diagonal matrix elements of the exchange-correlation potential
        ''')

    electronic_kinetic_energy = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Self-consistent electronic kinetic energy as defined in XC_method.
        ''')

    energy_correction_entropy = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Entropy correction to the potential energy to compensate for the change in
        occupation so that forces at finite T do not need to keep the change of occupation
        in account. Defined consistently with XC_method.
        ''')

    energy_correction_hartree = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Correction to the density-density electrostatic energy in the sum of eigenvalues
        (that uses the mixed density on one side), and the fully consistent density-
        density electrostatic energy. Defined consistently with XC_method.
        ''')

    energy_current = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Value of the energy calculated with calculation_method_current. energy_current is
        equal to energy_total for non-perturbative methods. For perturbative methods,
        energy_current is equal to the correction: energy_total minus energy_total of the
        calculation_to_calculation_ref with calculation_to_calculation_kind =
        starting_point (see the [method_to_method_kind wiki
        page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/method-
        to-method-kind)). See also [energy_current wiki
        page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/energy-
        current).
        ''')

    energy_electrostatic = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Total electrostatic energy (nuclei + electrons), defined consistently with
        calculation_method.
        ''')

    energy_free_per_atom = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Free energy per atom (whose minimum gives the smeared occupation density
        calculated with smearing_kind) calculated with XC_method.
        ''')

    energy_free = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Free energy (nuclei + electrons) (whose minimum gives the smeared occupation
        density calculated with smearing_kind) calculated with the method described in
        XC_method.
        ''')

    energy_hartree_error = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Error in the Hartree (electrostatic) potential energy. Defined consistently with
        XC_method.
        ''')

    energy_hartree_fock_X_scaled = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Scaled exact-exchange energy that depends on the mixing parameter of the
        functional. For example in hybrid functionals, the exchange energy is given as a
        linear combination of exact-energy and exchange energy of an approximate DFT
        functional; the exact exchange energy multiplied by the mixing coefficient of the
        hybrid functional would be stored in this metadata. Defined consistently with
        XC_method.
        ''')

    energy_method_current = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Value of the energy calculated with the method calculation_method_current.
        Depending on calculation_method_kind it might be a total energy or only a
        correction.
        ''')

    energy_sum_eigenvalues_per_atom = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Value of the energy per atom, where the energy is defined as the sum of the
        eigenvalues of the Hamiltonian matrix given by XC_method.
        ''')

    energy_sum_eigenvalues = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Sum of the eigenvalues of the Hamiltonian matrix defined by XC_method.
        ''')

    energy_T0_per_atom = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Value of the total energy per atom, calculated with the method described in
        XC_method and extrapolated to $T=0$, based on a free-electron gas argument.
        ''')

    energy_total_T0_per_atom = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Value of the total energy, calculated with the method described in XC_method per
        atom extrapolated to $T=0$, based on a free-electron gas argument.
        ''')

    energy_total_T0 = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Value of the total energy (or equivalently free energy), calculated with the
        method described in XC_method and extrapolated to $T=0$, based on a free-electron
        gas argument.
        ''')

    energy_total = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Value of the total energy, calculated with the method described in XC_method and
        extrapolated to $T=0$, based on a free-electron gas argument.
        ''')

    energy_XC_potential = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Value of the exchange-correlation (XC) potential energy: the integral of the first
        order derivative of the functional stored in XC_functional (integral of
        v_xc*electron_density), i.e., the component of XC that is in the sum of the
        eigenvalues. Value associated with the configuration, should be the most converged
        value.
        ''')

    energy_zero_point = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        description='''
        Value for the converged zero-point vibrations energy calculated using the method
        described in zero_point_method , and used in energy_current .
        ''')

    hessian_matrix = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_atoms', 'number_of_atoms', 3, 3],
        description='''
        The matrix with the second derivative with respect to atom displacements.
        ''')

    message_debug_evaluation = Quantity(
        type=str,
        shape=[],
        description='''
        A debugging message of the computational program, associated with a *single
        configuration calculation* (see section_single_configuration_calculation).
        ''')

    message_error_evaluation = Quantity(
        type=str,
        shape=[],
        description='''
        An error message of the computational program, associated with a *single
        configuration calculation* (see section_single_configuration_calculation).
        ''')

    message_info_evaluation = Quantity(
        type=str,
        shape=[],
        description='''
        An information message of the computational program, associated with a *single
        configuration calculation* (see section_single_configuration_calculation).
        ''')

    message_warning_evaluation = Quantity(
        type=str,
        shape=[],
        description='''
        A warning message of the computational program.
        ''')

    parsing_message_debug_evaluation = Quantity(
        type=str,
        shape=[],
        description='''
        This field is used for debugging messages of the parsing program associated with a
        run, see section_run.
        ''')

    parsing_message_error_single_configuration = Quantity(
        type=str,
        shape=[],
        description='''
        This field is used for error messages of the parsing program associated with a
        single configuration calculation, see section_single_configuration_calculation.
        ''')

    parsing_message_info_single_configuration = Quantity(
        type=str,
        shape=[],
        description='''
        This field is used for info messages of the parsing program associated with a
        single configuration calculation, see section_single_configuration_calculation.
        ''')

    parsing_message_warning_evaluation = Quantity(
        type=str,
        shape=[],
        description='''
        This field is used for warning messages of the parsing program associated with a
        run, see section_run.
        ''')

    single_configuration_calculation_converged = Quantity(
        type=np.dtype(np.int8),
        shape=[],
        description='''
        Determines whether a *single configuration calculation* in
        section_single_configuration_calculation is converged.
        ''')

    single_configuration_calculation_to_system_ref = Quantity(
        type=MProxy('section_system'),
        shape=[],
        description='''
        Reference to the system (atomic configuration, cell, ...) that is calculated in
        section_single_configuration_calculation.
        ''')

    single_configuration_to_calculation_method_ref = Quantity(
        type=MProxy('section_method'),
        shape=[],
        description='''
        Reference to the method used for the calculation in
        section_single_configuration_calculation.
        ''')

    spin_S2 = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        description='''
        Stores the value of the total spin moment operator $S^2$ for the converged
        wavefunctions calculated with the XC_method. It can be used to calculate the spin
        contamination in spin-unrestricted calculations.
        ''')

    stress_tensor = Quantity(
        type=np.dtype(np.float32),
        shape=[3, 3],
        unit='pascal',
        description='''
        Stores the final value of the default stress tensor consistent with energy_total
        and calculated with the method specified in stress_tensor_method.

        This value is used (if needed) for, e.g., molecular dynamics and geometry
        optimization. Alternative definitions of the stress tensor can be assigned with
        stress_tensor_kind
        ''')

    time_calculation = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='second',
        description='''
        Stores the wall-clock time needed for a calculation using
        calculation_method_current. Basically, it tracks the real time that has been
        elapsed from start to end.
        ''')

    time_single_configuration_calculation_cpu1_end = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='second',
        description='''
        Stores the end time of the *single configuration calculation* (see
        section_single_configuration_calculation) on CPU 1.
        ''')

    time_single_configuration_calculation_cpu1_start = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='second',
        description='''
        Stores the start time of the *single configuration calculation* (see
        section_single_configuration_calculation) on CPU 1.
        ''')

    time_single_configuration_calculation_date_end = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='second',
        description='''
        Stores the end date of the *single configuration calculation* (see
        section_single_configuration_calculation) as time since the *Unix epoch* (00:00:00
        UTC on 1 January 1970) in seconds. For date and times without a timezone, the
        default timezone GMT is used.
        ''')

    time_single_configuration_calculation_date_start = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='second',
        description='''
        Stores the start date of the *single configuration calculation* (see
        section_single_configuration_calculation) as time since the *Unix epoch* (00:00:00
        UTC on 1 January 1970) in seconds. For date and times without a timezone, the
        default timezone GMT is used.
        ''')

    time_single_configuration_calculation_wall_end = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='second',
        description='''
        Stores the internal wall-clock time at the end of the *single configuration
        calculation* (see section_single_configuration_calculation).
        ''')

    time_single_configuration_calculation_wall_start = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='second',
        description='''
        Stores the internal wall-clock time from the start of the *single configuration
        calculation* (see section_single_configuration_calculation).
        ''')

    zero_point_method = Quantity(
        type=str,
        shape=[],
        description='''
        Describes the zero-point vibrations method. If skipped or an empty string is used,
        it means no zero-point vibrations correction is applied.
        ''')

    section_excited_states = SubSection(
        sub_section=MProxy('section_excited_states'),
        repeats=False)

    section_atom_projected_dos = SubSection(
        sub_section=MProxy('section_atom_projected_dos'),
        repeats=True)

    section_atomic_multipoles = SubSection(
        sub_section=MProxy('section_atomic_multipoles'),
        repeats=False)

    section_basis_set = SubSection(
        sub_section=MProxy('section_basis_set'),
        repeats=False)

    section_calculation_to_calculation_refs = SubSection(
        sub_section=MProxy('section_calculation_to_calculation_refs'),
        repeats=False)

    section_calculation_to_folder_refs = SubSection(
        sub_section=MProxy('section_calculation_to_folder_refs'),
        repeats=False)

    section_dos = SubSection(
        sub_section=MProxy('section_dos'),
        repeats=True)

    section_eigenvalues = SubSection(
        sub_section=MProxy('section_eigenvalues'),
        repeats=False)

    section_energy_code_independent = SubSection(
        sub_section=MProxy('section_energy_code_independent'),
        repeats=False)

    section_energy_van_der_Waals = SubSection(
        sub_section=MProxy('section_energy_van_der_Waals'),
        repeats=True)

    section_k_band_normalized = SubSection(
        sub_section=MProxy('section_k_band_normalized'),
        repeats=True)

    section_k_band = SubSection(
        sub_section=MProxy('section_k_band'),
        repeats=True)

    section_scf_iteration = SubSection(
        sub_section=MProxy('section_scf_iteration'),
        repeats=True)

    section_species_projected_dos = SubSection(
        sub_section=MProxy('section_species_projected_dos'),
        repeats=True)

    section_stress_tensor = SubSection(
        sub_section=MProxy('section_stress_tensor'),
        repeats=False)

    section_volumetric_data = SubSection(
        sub_section=MProxy('section_volumetric_data'),
        repeats=False)


class section_excited_states(MSection):
    '''
    Excitation energies.
    '''

    excitation_energies = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_excited_states'],
        description='''
        Excitation energies.
        ''')

    number_of_excited_states = Quantity(
        type=int,
        shape=[],
        description='''
        Number of excited states.
        ''')

    oscillator_strengths = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_excited_states'],
        description='''
        Excited states oscillator strengths.
        ''')

    transition_dipole_moments = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_excited_states', 3],
        description='''
        Transition dipole moments.
        ''')


class section_eigenvalues(MSection):
    '''
    Linearization prefactor
    '''

    gw_qp_linearization_prefactor = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_spin_channels', 'number_of_eigenvalues_kpoints', 'number_of_eigenvalues'],
        description='''
        Linearization prefactor
        ''')

    eigenvalues_kind = Quantity(
        type=str,
        shape=[],
        description='''
        A short string describing the kind of eigenvalues, as defined in the
        [eigenvalues_kind wiki page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-
        info/wikis/metainfo/eigenvalues-kind).
        ''')

    eigenvalues_kpoints_multiplicity = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_eigenvalues_kpoints'],
        description='''
        Multiplicity of the $k$ point (i.e., how many distinct points per cell this
        expands to after applying all symmetries). This defaults to 1. If expansion is
        preformed then each point will have weight
        eigenvalues_kpoints_weights/eigenvalues_kpoints_multiplicity.
        ''')

    eigenvalues_kpoints_weights = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_eigenvalues_kpoints'],
        description='''
        Weights of the $k$ points (in the basis of the reciprocal lattice vectors) used
        for the evaluation of the eigenvalues tabulated in eigenvalues_values, should
        account for symmetry too.
        ''')

    eigenvalues_kpoints = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_eigenvalues_kpoints', 3],
        description='''
        Coordinates of the $k$ points (in the basis of the reciprocal lattice vectors)
        used for the evaluation of the eigenvalues tabulated in eigenvalues_values.
        ''')

    eigenvalues_occupation = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_spin_channels', 'number_of_eigenvalues_kpoints', 'number_of_eigenvalues'],
        description='''
        Occupation of the eigenstates. The corresponding eigenvalues (energy) are given in
        eigenvalues_values. The coordinates in the reciprocal space are defined in
        eigenvalues_kpoints.
        ''')

    eigenvalues_values = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_spin_channels', 'number_of_eigenvalues_kpoints', 'number_of_eigenvalues'],
        unit='joule',
        description='''
        Values of the (electronic-energy) eigenvalues. The coordinates of the
        corresponding eigenstates in the reciprocal space are defined in
        eigenvalues_kpoints and their occupations are given in eigenvalues_occupation.
        ''')

    number_of_band_segment_eigenvalues = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of eigenvalues in a band segment, see band_energies.
        ''')

    number_of_eigenvalues_kpoints = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of $k$ points, see eigenvalues_kpoints. $k$ points are calculated
        within a run and are irreducible if a symmetry is used.
        ''')

    number_of_eigenvalues = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of eigenvalues, see eigenvalues_values.
        ''')

    number_of_normalized_band_segment_eigenvalues = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of normalized eigenvalues in a band segment, see

        band_energies_normalized.
        ''')


class section_interaction(MSection):
    '''
    List of the indexes involved in this interaction. The fist atom has index 1, the last
    atom index number_of_topology_atoms.
    '''

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
    Reference to an atom-centered basis set defined in section_basis_set_atom_centered and
    to the atom kind as defined in section_method_atom_kind.
    '''

    mapping_section_method_basis_set_atom_centered = Quantity(
        type=np.dtype(np.int64),
        shape=['number_of_basis_sets_atom_centered', 2],
        description='''
        Reference to an atom-centered basis set defined in section_basis_set_atom_centered
        and to the atom kind as defined in section_method_atom_kind.
        ''')

    mapping_section_method_basis_set_cell_associated = Quantity(
        type=MProxy('section_basis_set_cell_dependent'),
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


class section_basis_set_cell_dependent(MSection):
    '''
    Reference to a cell-associated basis set.
    '''

    basis_set_cell_dependent_kind = Quantity(
        type=str,
        shape=[],
        description='''
        A string defining the type of the cell-dependent basis set (i.e., non atom
        centered such as plane-waves). Allowed values are listed in the
        [basis_set_cell_dependent_kind wiki page](https://gitlab.mpcdf.mpg.de/nomad-
        lab/nomad-meta-info/wikis/metainfo/basis-set-cell-dependent-kind).
        ''')

    basis_set_cell_dependent_name = Quantity(
        type=str,
        shape=[],
        description='''
        A label identifying the cell-dependent basis set (i.e., non atom centered such as
        plane-waves). Allowed values are listed in the [basis_set_cell_dependent_name wiki
        page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/basis-
        set-cell-dependent-name).
        ''')

    basis_set_planewave_cutoff = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Spherical cutoff  in reciprocal space for a plane-wave basis set. It is the energy
        of the highest plan-ewave ($\\frac{\\hbar^2|k+G|^2}{2m_e}$) included in the basis
        set. Note that normally this basis set is used for the wavefunctions, and the
        density would have 4 times the cutoff, but this actually depends on the use of the
        basis set by the method.
        ''')


class section_molecule_constraint(MSection):
    '''
    List of the indexes involved in this constraint. The fist atom has index 1, the last
    index is number_of_atoms_in_molecule.
    '''

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
    List of the indexes involved in this bonded interaction within a molecule. The first
    atom has index 1, the last index is number_of_atoms_in_.
    '''

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
    Mapping from molecules to molecule types.
    '''

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
        sub_section=MProxy('section_molecule_constraint'),
        repeats=False)

    section_molecule_interaction = SubSection(
        sub_section=MProxy('section_molecule_interaction'),
        repeats=False)


class section_soap_coefficients(MSection):
    '''
    number of soap coefficients
    '''

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
        type=np.dtype(np.float32),
        shape=['number_of_soap_coefficients'],
        description='''
        Compressed coefficient of the soap descriptor for the atom pair
        soap_coefficients_atom_pair
        ''')


class response_context(MSection):
    '''
    The top level context containing the reponse to an api query, when using jsonAPI they
    are tipically in the meta part
    '''

    shortened_meta_info = Quantity(
        type=str,
        shape=[],
        description='''
        A meta info whose corresponding data has been shortened
        ''')

    section_response_message = SubSection(
        sub_section=MProxy('section_response_message'),
        repeats=False)


class section_response_message(MSection):
    '''
    How many times this message was repeated
    '''

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


class section_soap(MSection):
    '''
    Stores the soap coefficients for the pair of atoms given in
    soap_coefficients_atom_pair.
    '''

    soap_angular_basis_L = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        angular basis L
        ''')

    soap_angular_basis_type = Quantity(
        type=str,
        shape=[],
        description='''
        angular basis type
        ''')

    soap_kernel_adaptor = Quantity(
        type=str,
        shape=[],
        description='''
        kernel adaptor
        ''')

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
        ''')

    soap_radial_basis_mode = Quantity(
        type=str,
        shape=[],
        description='''
        radial basis mode
        ''')

    soap_radial_basis_n = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        radial basis N
        ''')

    soap_radial_basis_sigma = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        description='''
        radial basis sigma
        ''')

    soap_radial_basis_type = Quantity(
        type=str,
        shape=[],
        description='''
        radial basis type
        ''')

    soap_radial_cutoff_center_weight = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        description='''
        radial cutoff center weight
        ''')

    soap_radial_cutoff_rc_width = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        radial cutoff width
        ''')

    soap_radial_cutoff_rc = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        description='''
        radial cutoff
        ''')

    soap_radial_cutoff_type = Quantity(
        type=str,
        shape=[],
        description='''
        radial cutoff type
        ''')

    soap_spectrum_2l1_norm = Quantity(
        type=np.dtype(np.int8),
        shape=[],
        description='''
        2l1 norm spectrum
        ''')

    soap_spectrum_global = Quantity(
        type=np.dtype(np.int8),
        shape=[],
        description='''
        global spectrum
        ''')

    soap_spectrum_gradients = Quantity(
        type=np.dtype(np.int8),
        shape=[],
        description='''
        gradients in specturm
        ''')

    soap_type_list = Quantity(
        type=str,
        shape=[],
        description='''
        Type list
        ''')

    section_soap_coefficients = SubSection(
        sub_section=MProxy('section_soap_coefficients'),
        repeats=False)


class section_system(MSection):
    '''
    Stores a soap descriptor for this configuration.
    '''

    topology_ref = Quantity(
        type=MProxy('section_topology'),
        shape=[],
        description='''
        Reference to the topology used for this system; if not given, the trivial topology
        should be assumed.
        ''')

    is_representative = Quantity(
        type=np.dtype(np.int8),
        shape=[],
        description='''
        Most systems in a run are only minor variations of each other. Systems marked
        representative where chosen to be representative for all systems in the run.
        ''')

    atom_atom_number = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_sites'],
        description='''
        Atomic number Z of the atom.
        ''')

    atom_concentrations = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_atoms'],
        description='''
        concentration of the atom species in a variable composition, by default it should
        be considered an array of ones. Summing these should give the number_of_sites
        ''')

    atom_species = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Species of the atom (normally the atomic number Z, 0 or negative for unidentifed
        species or particles that are not atoms.
        ''')

    atom_velocities = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_atoms', 3],
        unit='meter / second',
        description='''
        Velocities of the nuclei, defined as the change in Cartesian coordinates of the
        nuclei with respect to time.
        ''')

    configuration_raw_gid = Quantity(
        type=str,
        shape=[],
        description='''
        checksum of the configuration_core, i.e. the geometry of the system. The values
        are not normalized in any way so equivalent configurations might have different
        values
        ''')

    local_rotations = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_atoms', 3, 3],
        description='''
        A rotation matrix defining the orientation of each atom. If the rotation matrix
        only needs to be specified for some atoms, the remaining atoms should set it to
        the zero matrix (not the identity!)
        ''')

    number_of_atoms = Quantity(
        type=int,
        shape=[],
        description='''
        Stores the total number of atoms used in the calculation. For alloys where
        concentrations of species are given for each site in the unit cell, it stores the
        number of sites.
        ''')

    number_of_sites = Quantity(
        type=int,
        shape=[],
        description='''
        number of sites in a variable composition representation. By default (no variable
        composition) it is the same as number_of_atoms.
        ''')

    SC_matrix = Quantity(
        type=np.dtype(np.int32),
        shape=[3, 3],
        description='''
        Specifies the matrix that transforms the unit-cell into the super-cell in which
        the actual calculation is performed.
        ''')

    symmorphic = Quantity(
        type=np.dtype(np.int8),
        shape=[],
        description='''
        Is the space group symmorphic? Set to True if all translations are zero.
        ''')

    system_composition = Quantity(
        type=str,
        shape=[],
        description='''
        Composition, i.e. cumulative chemical formula with atoms ordered by decreasing
        atomic number Z.
        ''')

    system_configuration_consistent = Quantity(
        type=np.dtype(np.int8),
        shape=[],
        description='''
        Flag set is the configuration is consistent
        ''')

    system_name = Quantity(
        type=str,
        shape=[],
        description='''
        Specifies the name of the system. This information is provided by the user in some
        codes and is stored here for debugging or visualization purposes.
        ''')

    system_reweighted_composition = Quantity(
        type=str,
        shape=[],
        description='''
        Composition, i.e. cumulative chemical with atoms ordered by decreasing atomic
        number Z reweighted so that the sum is close to 100, and values are rounded up,
        and are stable (i.e. it is a fixed point).
        ''')

    system_type = Quantity(
        type=str,
        shape=[],
        description='''
        Type of the system
        ''')

    time_reversal_symmetry = Quantity(
        type=np.dtype(np.int8),
        shape=[],
        description='''
        Is time-reversal symmetry present?
        ''')

    chemical_composition = Quantity(
        type=str,
        shape=[],
        description='''
        The chemical composition as full formula of the system, based on atom species.
        ''')

    chemical_composition_reduced = Quantity(
        type=str,
        shape=[],
        description='''
        The chemical composition as reduced formula of the system, based on atom species.
        ''')

    chemical_composition_bulk_reduced = Quantity(
        type=str,
        shape=[],
        description='''
        The chemical composition as reduced bulk formula of the system, based on atom
        species.
        ''')

    section_soap = SubSection(
        sub_section=MProxy('section_soap'),
        repeats=False)

    section_prototype = SubSection(
        sub_section=MProxy('section_prototype'),
        repeats=False)

    section_symmetry = SubSection(
        sub_section=MProxy('section_symmetry'),
        repeats=False)

    section_system_to_system_refs = SubSection(
        sub_section=MProxy('section_system_to_system_refs'),
        repeats=False)


class section_run(MSection):
    '''
    Section containing the definition of topology (connectivity among atoms in force
    fileds), force field, and constraints of a system.
    '''

    calculation_file_uri = Quantity(
        type=str,
        shape=[],
        description='''
        Contains the nomad uri of a raw the data file connected to the current run. There
        should be an value for the main_file_uri and all ancillary files.
        ''')

    message_debug_run = Quantity(
        type=str,
        shape=[],
        description='''
        A debugging message of the computational program, associated with a run.
        ''')

    message_error_run = Quantity(
        type=str,
        shape=[],
        description='''
        An error message of the computational program, associated with a run.
        ''')

    message_info_run = Quantity(
        type=str,
        shape=[],
        description='''
        An information message of the computational program, associated with a run.
        ''')

    message_warning_run = Quantity(
        type=str,
        shape=[],
        description='''
        A warning message of the computational program, associated with a run.
        ''')

    parsing_message_debug_run = Quantity(
        type=str,
        shape=[],
        description='''
        This field is used for debugging messages of the parsing program associated with a
        single configuration calculation, see section_single_configuration_calculation.
        ''')

    parsing_message_error_run = Quantity(
        type=str,
        shape=[],
        description='''
        This field is used for error messages of the parsing program associated with a
        run, see section_run.
        ''')

    parsing_message_info_run = Quantity(
        type=str,
        shape=[],
        description='''
        This field is used for info messages of the parsing program associated with a run,
        see section_run.
        ''')

    parsing_message_warning_run = Quantity(
        type=str,
        shape=[],
        description='''
        This field is used for warning messages of the parsing program associated with a
        run, see section_run.
        ''')

    program_basis_set_type = Quantity(
        type=str,
        shape=[],
        description='''
        The type of basis set used by the program to represent wave functions.

        Valid values are:

        * Numeric AOs

        * Gaussians

        * (L)APW+lo

        * FLAPW (full-potential linearized augmented planewave)

        * Plane waves

        * Real-space grid

        * Local-orbital minimum-basis
        ''')

    run_clean_end = Quantity(
        type=np.dtype(np.int8),
        shape=[],
        description='''
        Indicates whether this run terminated properly (true), or if it was killed or
        exited with an error code unequal to zero (false).
        ''')

    run_hosts = Quantity(
        type=typing.Any,
        shape=[],
        description='''
        An associative list of host(s) that performed this simulation. This is an
        associative list that contains program-dependent information (*key*) on how the
        host was used (*value*). Useful for debugging purposes.
        ''')

    time_run_cpu1_end = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='second',
        description='''
        Stores the end time of the run on CPU 1.
        ''')

    time_run_cpu1_start = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='second',
        description='''
        Stores the start time of the run on CPU 1.
        ''')

    time_run_date_end = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='second',
        description='''
        Stores the end date of the run as time since the *Unix epoch* (00:00:00 UTC on 1
        January 1970) in seconds. For date and times without a timezone, the default
        timezone GMT is used.
        ''')

    time_run_date_start = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='second',
        description='''
        Stores the start date of the run as time since the *Unix epoch* (00:00:00 UTC on 1
        January 1970) in seconds. For date and times without a timezone, the default
        timezone GMT is used.
        ''')

    time_run_wall_end = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='second',
        description='''
        Stores the internal wall-clock time at the end of the run.
        ''')

    time_run_wall_start = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='second',
        description='''
        Stores the internal wall-clock time from the start of the run.
        ''')

    section_topology = SubSection(
        sub_section=MProxy('section_topology'),
        repeats=False)

    section_frame_sequence = SubSection(
        sub_section=MProxy('section_frame_sequence'),
        repeats=False)

    section_method = SubSection(
        sub_section=MProxy('section_method'),
        repeats=False)

    section_sampling_method = SubSection(
        sub_section=MProxy('section_sampling_method'),
        repeats=False)

    section_single_configuration_calculation = SubSection(
        sub_section=MProxy('section_single_configuration_calculation'),
        repeats=False)

    section_system = SubSection(
        sub_section=MProxy('section_system'),
        repeats=False)


class archive_context(MSection):
    '''
    Contains information relating to an archive.
    '''

    archive_gid = Quantity(
        type=str,
        shape=[],
        description='''
        unique identifier of an archive.
        ''')


class section_primitive_system(MSection):
    '''
    Atom positions in the primitive cell in reduced units.
    '''

    atom_positions_primitive = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_atoms_primitive', 3],
        description='''
        Atom positions in the primitive cell in reduced units.
        ''')

    atomic_numbers_primitive = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_atoms_primitive'],
        description='''
        Atomic numbers in the primitive cell.
        ''')

    equivalent_atoms_primitive = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_atoms_primitive'],
        description='''
        Gives a mapping table of atoms to symmetrically independent atoms in the primitive
        cell. This is used to find symmetrically equivalent atoms.
        ''')

    lattice_vectors_primitive = Quantity(
        type=np.dtype(np.float32),
        shape=[3, 3],
        unit='meter',
        description='''
        Primitive lattice vectors. The vectors are the rows of this matrix.
        ''')

    number_of_atoms_primitive = Quantity(
        type=int,
        shape=[],
        description='''
        Number of atoms in primitive system.
        ''')

    wyckoff_letters_primitive = Quantity(
        type=str,
        shape=['number_of_atoms_primitive'],
        description='''
        Wyckoff letters for atoms in the primitive cell.
        ''')


class section_std_system(MSection):
    '''
    Standardized atom positions in reduced units.
    '''

    atom_positions_std = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_atoms_std', 3],
        description='''
        Standardized atom positions in reduced units.
        ''')

    atomic_numbers_std = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_atoms_std'],
        description='''
        Atomic numbers of the atoms in the standardized cell.
        ''')

    equivalent_atoms_std = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_atoms_std'],
        description='''
        Gives a mapping table of atoms to symmetrically independent atoms in the
        standardized cell. This is used to find symmetrically equivalent atoms.
        ''')

    lattice_vectors_std = Quantity(
        type=np.dtype(np.float32),
        shape=[3, 3],
        unit='meter',
        description='''
        Standardized lattice vectors of the conventional cell. The vectors are the rows of
        this matrix.
        ''')

    number_of_atoms_std = Quantity(
        type=int,
        shape=[],
        description='''
        Number of atoms in standardized system.
        ''')

    wyckoff_letters_std = Quantity(
        type=str,
        shape=['number_of_atoms_std'],
        description='''
        Wyckoff letters for atoms in the standardized cell.
        ''')


class section_atom_projected_dos(MSection):
    '''
    Array containing the set of discrete energy values for the atom-projected density
    (electronic-energy) of states (DOS).
    '''

    atom_projected_dos_energies = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_atom_projected_dos_values'],
        unit='joule',
        description='''
        Array containing the set of discrete energy values for the atom-projected density
        (electronic-energy) of states (DOS).
        ''')

    atom_projected_dos_lm = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_lm_atom_projected_dos', 2],
        description='''
        Tuples of $l$ and $m$ values for which atom_projected_dos_values_lm are given. For
        the quantum number $l$ the conventional meaning of azimuthal quantum number is
        always adopted. For the integer number $m$, besides the conventional use as
        magnetic quantum number ($l+1$ integer values from $-l$ to $l$), a set of
        different conventions is accepted (see the [m_kind wiki
        page](https://gitlab.rzg.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/m-kind).
        The adopted convention is specified by atom_projected_dos_m_kind.
        ''')

    atom_projected_dos_m_kind = Quantity(
        type=str,
        shape=[],
        description='''
        String describing what the integer numbers of $m$ in atom_projected_dos_lm mean.
        The allowed values are listed in the [m_kind wiki
        page](https://gitlab.rzg.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/m-kind).
        ''')

    atom_projected_dos_values_lm = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_lm_atom_projected_dos', 'number_of_spin_channels', 'number_of_atoms', 'number_of_atom_projected_dos_values'],
        description='''
        Values correspond to the number of states for a given energy (the set of discrete
        energy values is given in atom_projected_dos_energies) divided into contributions
        from each $l,m$ channel for the atom-projected density (electronic-energy) of
        states. Here, there are as many atom-projected DOS as the number_of_atoms, the
        list of labels of the atoms and their meanings are in atom_labels.
        ''')

    atom_projected_dos_values_total = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_spin_channels', 'number_of_atoms', 'number_of_atom_projected_dos_values'],
        description='''
        Values correspond to the number of states for a given energy (the set of discrete
        energy values is given in atom_projected_dos_energies) divided into contributions
        summed up over all $l$ channels for the atom-projected density (electronic-energy)
        of states (DOS). Here, there are as many atom-projected DOS as the
        number_of_atoms, the list of labels of the atoms and their meanings are in
        atom_labels.
        ''')

    number_of_atom_projected_dos_values = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of energy values for the atom-projected density of states (DOS)
        based on atom_projected_dos_values_lm and atom_projected_dos_values_total.
        ''')

    number_of_lm_atom_projected_dos = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of $l$, $m$ combinations for the atom projected density of states
        (DOS) defined in section_atom_projected_dos.
        ''')


class section_atomic_multipoles(MSection):
    '''
    String describing the method used to obtain the electrostatic multipoles (including
    the electric charge, dipole, etc.) for each atom. Such multipoles require a charge-
    density partitioning scheme, specified by the value of this metadata. Allowed values
    are listed in the [atomic_multipole_kind wiki page](https://gitlab.mpcdf.mpg.de/nomad-
    lab/nomad-meta-info/wikis/metainfo/atomic-multipole-kind).
    '''

    atomic_multipole_kind = Quantity(
        type=str,
        shape=[],
        description='''
        String describing the method used to obtain the electrostatic multipoles
        (including the electric charge, dipole, etc.) for each atom. Such multipoles
        require a charge-density partitioning scheme, specified by the value of this
        metadata. Allowed values are listed in the [atomic_multipole_kind wiki
        page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/atomic-
        multipole-kind).
        ''')

    atomic_multipole_lm = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_lm_atomic_multipoles', 2],
        description='''
        Tuples of $l$ and $m$ values for which the atomic multipoles (including the
        electric charge, dipole, etc.) are given. The method used to obtain the multipoles
        is specified by atomic_multipole_kind. The meaning of the integer number $l$ is
        monopole/charge for $l=0$, dipole for $l=1$, quadrupole for $l=2$, etc. The
        meaning of the integer numbers $m$ is specified by atomic_multipole_m_kind.
        ''')

    atomic_multipole_m_kind = Quantity(
        type=str,
        shape=[],
        description='''
        String describing the definition for each integer number $m$ in
        atomic_multipole_lm. Allowed values are listed in the [m_kind wiki
        page](https://gitlab.rzg.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/m-kind).
        ''')

    atomic_multipole_values = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_lm_atomic_multipoles', 'number_of_atoms'],
        description='''
        Value of the multipoles (including the monopole/charge for $l$ = 0, the dipole for
        $l$ = 1, etc.) for each atom, calculated as described in atomic_multipole_kind.
        ''')

    number_of_lm_atomic_multipoles = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of $l$, $m$ combinations for atomic multipoles
        atomic_multipole_lm.
        ''')


class section_k_band_segment_normalized(MSection):
    '''
    $k$-dependent energies of the electronic band segment (electronic band structure) with
    respect to the top of the valence band. This is a third-order tensor, with one
    dimension used for the spin channels, one for the $k$ points for each segment, and one
    for the eigenvalue sequence.
    '''

    band_energies_normalized = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_spin_channels', 'number_of_normalized_k_points_per_segment', 'number_of_normalized_band_segment_eigenvalues'],
        unit='joule',
        description='''
        $k$-dependent energies of the electronic band segment (electronic band structure)
        with respect to the top of the valence band. This is a third-order tensor, with
        one dimension used for the spin channels, one for the $k$ points for each segment,
        and one for the eigenvalue sequence.
        ''')

    band_k_points_normalized = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_normalized_k_points_per_segment', 3],
        description='''
        Fractional coordinates of the $k$ points (in the basis of the reciprocal-lattice
        vectors) for which the normalized electronic energies are given.
        ''')

    band_occupations_normalized = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_spin_channels', 'number_of_normalized_k_points_per_segment', 'number_of_normalized_band_segment_eigenvalues'],
        description='''
        Occupation of the $k$-points along the normalized electronic band. The size of the
        dimensions of this third-order tensor are the same as for the tensor in
        band_energies.
        ''')

    band_segm_labels_normalized = Quantity(
        type=str,
        shape=[2],
        description='''
        Start and end labels of the points in the segment (one-dimensional pathways)
        sampled in the $k$-space, using the conventional symbols, e.g., Gamma, K, L. The
        coordinates (fractional, in the reciprocal space) of the start and end points for
        each segment are given in band_segm_start_end_normalized
        ''')

    band_segm_start_end_normalized = Quantity(
        type=np.dtype(np.float32),
        shape=[2, 3],
        description='''
        Fractional coordinates of the start and end point (in the basis of the reciprocal
        lattice vectors) of the segment sampled in the $k$ space. The conventional symbols
        (e.g., Gamma, K, L) of the same points are given in band_segm_labels
        ''')

    number_of_normalized_k_points_per_segment = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of $k$ points in the segment of the normalized band structure
        (see section_k_band_segment_normalized).
        ''')


class section_k_band_segment(MSection):
    '''
    $k$-dependent or $q$-dependent  energies of the electronic or vibrational band segment
    (electronic/vibrational band structure). This is a third-order tensor, with one
    dimension used for the spin channels (1 in case of a vibrational band structure), one
    for the $k$ or $q$ points for each segment, and one for the eigenvalue sequence.
    '''

    band_energies = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_spin_channels', 'number_of_k_points_per_segment', 'number_of_band_segment_eigenvalues'],
        unit='joule',
        description='''
        $k$-dependent or $q$-dependent  energies of the electronic or vibrational band
        segment (electronic/vibrational band structure). This is a third-order tensor,
        with one dimension used for the spin channels (1 in case of a vibrational band
        structure), one for the $k$ or $q$ points for each segment, and one for the
        eigenvalue sequence.
        ''')

    band_k_points = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_k_points_per_segment', 3],
        description='''
        Fractional coordinates of the $k$ or $q$ points (in the basis of the reciprocal-
        lattice vectors) for which the electronic energy are given.
        ''')

    band_occupations = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_spin_channels', 'number_of_k_points_per_segment', 'number_of_band_segment_eigenvalues'],
        description='''
        Occupation of the $k$-points along the electronic band. The size of the dimensions
        of this third-order tensor are the same as for the tensor in band_energies.
        ''')

    band_segm_labels = Quantity(
        type=str,
        shape=[2],
        description='''
        Start and end labels of the points in the segment (one-dimensional pathways)
        sampled in the $k$-space or $q$-space, using the conventional symbols, e.g.,
        Gamma, K, L. The coordinates (fractional, in the reciprocal space) of the start
        and end points for each segment are given in band_segm_start_end
        ''')

    band_segm_start_end = Quantity(
        type=np.dtype(np.float32),
        shape=[2, 3],
        description='''
        Fractional coordinates of the start and end point (in the basis of the reciprocal
        lattice vectors) of the segment sampled in the $k$ space. The conventional symbols
        (e.g., Gamma, K, L) of the same points are given in band_segm_labels
        ''')

    number_of_k_points_per_segment = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of $k$ points in the segment of the band structure, see
        section_k_band_segment.
        ''')


class section_k_band(MSection):
    '''
    String to specify the kind of band structure (either electronic or vibrational).
    '''

    band_structure_kind = Quantity(
        type=str,
        shape=[],
        description='''
        String to specify the kind of band structure (either electronic or vibrational).
        ''')

    section_k_band_segment = SubSection(
        sub_section=MProxy('section_k_band_segment'),
        repeats=True)


class section_basis_set_atom_centered(MSection):
    '''
    Azimuthal quantum number ($l$) values (of the angular part given by the spherical
    harmonic $Y_{lm}$) of the atom-centered basis function defined in the current
    section_basis_set_atom_centered.
    '''

    basis_set_atom_centered_ls = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_kinds_in_basis_set_atom_centered'],
        description='''
        Azimuthal quantum number ($l$) values (of the angular part given by the spherical
        harmonic $Y_{lm}$) of the atom-centered basis function defined in the current
        section_basis_set_atom_centered.
        ''')

    basis_set_atom_centered_radial_functions = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_kinds_in_basis_set_atom_centered', 401, 5],
        description='''
        Values of the radial function of the different basis function kinds. The values
        are numerically tabulated on a default 0.01-nm equally spaced grid from 0 to 4 nm.
        The 5 tabulated values are $r$, $f(r)$, $f'(r)$, $f(r) \\cdot r$,
        $\\frac{d}{dr}(f(r) \\cdot r)$.
        ''')

    basis_set_atom_centered_short_name = Quantity(
        type=str,
        shape=[],
        description='''
        Code-specific, but explicative, base name for the basis set (not unique). Details
        are explained in the [basis_set_atom_centered_short_name wiki
        page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/basis-
        set-atom-centered-short-name), this name should not contain the *atom kind* (to
        simplify the use of a single name for multiple elements).
        ''')

    basis_set_atom_centered_unique_name = Quantity(
        type=str,
        shape=[],
        description='''
        Code-specific, but explicative, base name for the basis set (not unique). This
        string starts with basis_set_atom_centered_short_name. If the basis set defined in
        this section_basis_set_atom_centered is not identical to the default definition
        (stored in a database) of the basis set with the same name stored in a database,
        then the string is extended by 10 identifiable characters as explained in the
        [basis_set_atom_centered_name wiki page](https://gitlab.mpcdf.mpg.de/nomad-
        lab/nomad-meta-info/wikis/metainfo/basis-set-atom-centered-unique-name). The
        reason for this procedure is that often atom-centered basis sets are obtained by
        fine tuning basis sets provided by the code developers or other sources. Each
        basis sets, which has normally a standard name, often reported in publications,
        has also several parameters that can be tuned. This metadata tries to keep track
        of the original basis set and its modifications. This string here defined should
        not contain the *atom kind* for which this basis set is intended for, in order to
        simplify the use of a single name for multiple *atom kinds* (see atom_labels for
        the actual meaning of *atom kind*).
        ''')

    basis_set_atom_number = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Atomic number (i.e., number of protons) of the atom for which this basis set is
        constructed (0 means unspecified or a pseudo atom).
        ''')

    number_of_basis_functions_in_basis_set_atom_centered = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of different basis functions in a section_basis_set_atom_centered
        section. This equals the number of actual coefficients that are specified when
        using this basis set.
        ''')

    number_of_kinds_in_basis_set_atom_centered = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of different *kinds* of radial basis functions in the
        section_basis_set_atom_centered section. Specifically, basis functions with the
        same $n$ and $l$ quantum numbers are grouped in sets. Each set counts as one
        *kind*.
        ''')

    section_basis_functions_atom_centered = SubSection(
        sub_section=MProxy('section_basis_functions_atom_centered'),
        repeats=False)

    section_gaussian_basis_group = SubSection(
        sub_section=MProxy('section_gaussian_basis_group'),
        repeats=False)


class section_basis_set(MSection):
    '''
    String describing the use of the basis set, i.e, if it used for expanding a wave-
    function or an electron density. Allowed values are listed in the [basis_set_kind wiki
    page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/basis-set-
    kind).
    '''

    basis_set_kind = Quantity(
        type=str,
        shape=[],
        description='''
        String describing the use of the basis set, i.e, if it used for expanding a wave-
        function or an electron density. Allowed values are listed in the [basis_set_kind
        wiki page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-
        info/wikis/metainfo/basis-set-kind).
        ''')

    basis_set_name = Quantity(
        type=str,
        shape=[],
        description='''
        String identifying the basis set in an unique way. The rules for building this
        string are specified in the [basis_set_name wiki
        page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/basis-
        set-name).
        ''')

    mapping_section_basis_set_atom_centered = Quantity(
        type=MProxy('section_basis_set_atom_centered'),
        shape=['number_of_atoms'],
        description='''
        An array of the dimension of number_of_atoms where each atom (identified by the
        index in the array) is assigned to an atom-centered basis set, for this
        section_single_configuration_calculation. The actual definition of the atom-
        centered basis set is in the section_basis_set_atom_centered that is referred to
        by this metadata.
        ''')

    mapping_section_basis_set_cell_dependent = Quantity(
        type=MProxy('section_basis_set_cell_dependent'),
        shape=[],
        description='''
        Assignment of the cell-dependent (i.e., non atom centered, e.g., plane-waves)
        parts of the basis set, which is defined (type, parameters) in
        section_basis_set_cell_dependent that is referred to by this metadata.
        ''')

    number_of_basis_functions = Quantity(
        type=int,
        shape=[],
        description='''
        Stores the total number of basis functions in a section_basis_set section.
        ''')


class section_symmetry(MSection):
    '''
    Identifier for the Bravais lattice in Pearson notation. The first lowercase letter
    identifies the crystal family and can be one of the following: a (triclinic), b
    (monoclinic), o (orthorhombic), t (tetragonal), h (hexagonal) or c (cubic). The second
    uppercase letter identifies the centring and can be one of the following: P
    (primitive), S (face centred), I (body centred), R (rhombohedral centring) or F (all
    faces centred).
    '''

    bravais_lattice = Quantity(
        type=str,
        shape=[],
        description='''
        Identifier for the Bravais lattice in Pearson notation. The first lowercase letter
        identifies the crystal family and can be one of the following: a (triclinic), b
        (monoclinic), o (orthorhombic), t (tetragonal), h (hexagonal) or c (cubic). The
        second uppercase letter identifies the centring and can be one of the following: P
        (primitive), S (face centred), I (body centred), R (rhombohedral centring) or F
        (all faces centred).
        ''')

    choice = Quantity(
        type=str,
        shape=[],
        description='''
        String that specifies the centering, origin and basis vector settings of the 3D
        space group that defines the symmetry group of the simulated physical system (see
        section_system). Values are as defined by spglib.
        ''')

    crystal_system = Quantity(
        type=str,
        shape=[],
        description='''
        Name of the crystal system. Can be one of the following: triclinic, monoclinic,
        orthorhombic, tetragonal, trigonal, hexagonal or cubic.
        ''')

    hall_number = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        The Hall number for this system.
        ''')

    hall_symbol = Quantity(
        type=str,
        shape=[],
        description='''
        The Hall symbol for this system.
        ''')

    international_short_symbol = Quantity(
        type=str,
        shape=[],
        description='''
        Specifies the International Union of Crystallography (IUC) short symbol of the 3D
        space group of this system
        ''')

    origin_shift = Quantity(
        type=np.dtype(np.float32),
        shape=[3],
        description='''
        Vector $\\mathbf{p}$ from the origin of the standardized system to the origin of
        the original system. Together with the matrix $\\mathbf{P}$, found in
        space_group_3D_transformation_matrix, the transformation between the standardized
        coordinates $\\mathbf{x}_s$ and original coordinates $\\mathbf{x}$ is then given by
        $\\mathbf{x}_s = \\mathbf{P} \\mathbf{x} + \\mathbf{p}$.
        ''')

    point_group = Quantity(
        type=str,
        shape=[],
        description='''
        Symbol of the crystallographic point group in the Hermann-Mauguin notation.
        ''')

    space_group_number = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Specifies the International Union of Crystallography (IUC) number of the 3D space
        group of this system.
        ''')

    symmetry_method = Quantity(
        type=str,
        shape=[],
        description='''
        Identifies the source of the symmetry information contained within this section.
        If equal to 'spg_normalized' the information comes from a normalization step.
        ''')

    transformation_matrix = Quantity(
        type=np.dtype(np.float32),
        shape=[3, 3],
        description='''
        Matrix $\\mathbf{P}$ that is used to transform the standardized coordinates to the
        original coordinates. Together with the vector $\\mathbf{p}$, found in
        space_group_3D_origin_shift, the transformation between the standardized
        coordinates $\\mathbf{x}_s$ and original coordinates $\\mathbf{x}$ is then given by
        $\\mathbf{x}_s = \\mathbf{P} \\mathbf{x} + \\mathbf{p}$.
        ''')

    section_original_system = SubSection(
        sub_section=MProxy('section_original_system'),
        repeats=False)

    section_primitive_system = SubSection(
        sub_section=MProxy('section_primitive_system'),
        repeats=False)

    section_std_system = SubSection(
        sub_section=MProxy('section_std_system'),
        repeats=False)


class calculation_context(MSection):
    '''
    Contains information relating to a calculation.
    '''

    calculation_gid = Quantity(
        type=str,
        shape=[],
        description='''
        unique identifier of a calculation.
        ''')


class section_restricted_uri(MSection):
    '''
    The number of restricted uris in restricted_uri list.
    '''

    number_of_restricted_uri = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        The number of restricted uris in restricted_uri list.
        ''')

    restricted_uri = Quantity(
        type=str,
        shape=['number_of_restricted_uri'],
        description='''
        The list of nomad uri(s) identifying the restricted info/file corresponding to
        this calculation
        ''')

    restricted_uri_reason = Quantity(
        type=str,
        shape=[],
        description='''
        The reason of restriction for the uri or file. The reason can be 'propriety
        license', 'open-source redistribution restricted license', 'other license', or
        'author restricted'.
        ''')

    restricted_uri_issue_authority = Quantity(
        type=str,
        shape=[],
        description='''
        The issue authority is the restriction owner for the uri or file. This can be
        license owner such as 'VASP' or 'AMBER', 'NOMAD', or the author of the uri. For
        example the repository user name of the author.
        ''')

    restricted_uri_end_date = Quantity(
        type=str,
        shape=[],
        description='''
        The deadline date of the restriction for the uri or file. The end date can be in
        date format string for those restrictions set by authors or NOMAD otherwise it is
        set to 'unlimited' for the restriction related to license.
        ''')

    restricted_uri_restriction = Quantity(
        type=str,
        shape=[],
        description='''
        The type of restriction for the uri or file. The type can be 'any access' or
        'license permitted'.
        ''')

    restricted_uri_license = Quantity(
        type=str,
        shape=[],
        description='''
        The info of the license that is the reason of restriction.
        ''')

    number_of_restricted_uri_files = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        The number of restricted files in restricted_uri_files list.
        ''')


class section_calculation_to_calculation_refs(MSection):
    '''
    URL used to reference an externally stored calculation. The kind of relationship
    between the present and the referenced section_single_configuration_calculation is
    specified by calculation_to_calculation_kind.
    '''

    calculation_to_calculation_external_url = Quantity(
        type=str,
        shape=[],
        description='''
        URL used to reference an externally stored calculation. The kind of relationship
        between the present and the referenced section_single_configuration_calculation is
        specified by calculation_to_calculation_kind.
        ''')

    calculation_to_calculation_kind = Quantity(
        type=str,
        shape=[],
        description='''
        String defining the relationship between the referenced
        section_single_configuration_calculation and the present
        section_single_configuration_calculation. Valid values are described in the
        [calculation_to_calculation_kind wiki page](https://gitlab.mpcdf.mpg.de/nomad-
        lab/nomad-meta-info/wikis/metainfo/calculation-to-calculation-kind). Often
        calculations are connected, for instance, one calculation is a perturbation
        performed using a self-consistent field (SCF) calculation as starting point, or a
        simulated system is partitioned in regions with different but connected
        Hamiltonians (e.g., QM/MM, or a region treated via Kohn-Sham DFT embedded into a
        region treated via orbital-free DFT). Hence, the need of keeping track of these
        connected calculations. The referenced calculation is identified via
        calculation_to_calculation_ref (typically used for a calculation in the same
        section_run) or calculation_to_calculation_external_url.
        ''')

    calculation_to_calculation_ref = Quantity(
        type=MProxy('section_single_configuration_calculation'),
        shape=[],
        description='''
        Reference to another calculation. If both this and
        calculation_to_calculation_external_url are given, then
        calculation_to_calculation_ref is a local copy of the URL given in
        calculation_to_calculation_external_url. The kind of relationship between the
        present and the referenced section_single_configuration_calculation is specified
        by calculation_to_calculation_kind.
        ''')


class section_calculation_to_folder_refs(MSection):
    '''
    URL used to reference a folder containing external calculations. The kind of
    relationship between the present and the referenced
    section_single_configuration_calculation is specified by calculation_to_folder_kind.
    '''

    calculation_to_folder_external_url = Quantity(
        type=str,
        shape=[],
        description='''
        URL used to reference a folder containing external calculations. The kind of
        relationship between the present and the referenced
        section_single_configuration_calculation is specified by
        calculation_to_folder_kind.
        ''')

    calculation_to_folder_kind = Quantity(
        type=str,
        shape=[],
        description='''
        String defining the relationship between the referenced
        section_single_configuration_calculation and a folder containing calculations.
        ''')


class section_dos(MSection):
    '''
    Array containing the set of discrete energy values with respect to the top of the
    valence band for the density (electronic-energy) of states (DOS). This is the total
    DOS, see atom_projected_dos_energies and species_projected_dos_energies for partial
    density of states.
    '''

    dos_energies_normalized = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_dos_values'],
        unit='joule',
        description='''
        Array containing the set of discrete energy values with respect to the top of the
        valence band for the density (electronic-energy) of states (DOS). This is the
        total DOS, see atom_projected_dos_energies and species_projected_dos_energies for
        partial density of states.
        ''')

    dos_energies = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_dos_values'],
        unit='joule',
        description='''
        Array containing the set of discrete energy values for the density (electronic-
        energy or vibrational energy) of states (DOS). This is the total DOS, see
        atom_projected_dos_energies and species_projected_dos_energies for partial density
        of states.
        ''')

    dos_fermi_energy = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        description='''
        Stores the Fermi energy of the density of states.
        ''')

    dos_integrated_values = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_spin_channels', 'number_of_dos_values'],
        description='''
        Integrated density of states (starting at $-\\infty$), pseudo potential
        calculations should start with the number of core electrons if they cover only the
        active electrons
        ''')

    dos_kind = Quantity(
        type=str,
        shape=[],
        description='''
        String to specify the kind of density of states (either electronic or
        vibrational).
        ''')

    dos_lm = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_dos_lms', 2],
        description='''
        Tuples of $l$ and $m$ values for which dos_values_lm are given. For the quantum
        number $l$ the conventional meaning of azimuthal quantum number is always adopted.
        For the integer number $m$, besides the conventional use as magnetic quantum
        number ($l+1$ integer values from $-l$ to $l$), a set of different conventions is
        accepted (see the [m_kind wiki page](https://gitlab.rzg.mpg.de/nomad-lab/nomad-
        meta-info/wikis/metainfo/m-kind). The actual adopted convention is specified by
        dos_m_kind.
        ''')

    dos_m_kind = Quantity(
        type=str,
        shape=[],
        description='''
        String describing what the integer numbers of $m$ in dos_lm mean. The allowed
        values are listed in the [m_kind wiki page](https://gitlab.rzg.mpg.de/nomad-
        lab/nomad-meta-info/wikis/metainfo/m-kind).
        ''')

    dos_values_lm = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_dos_lms', 'number_of_spin_channels', 'number_of_atoms', 'number_of_dos_values'],
        unit='joule',
        description='''
        Array containing the density (electronic-energy) of states values projected on the
        various spherical harmonics (integrated on all atoms), see
        atom_projected_dos_values_lm for atom values.
        ''')

    dos_values_per_atoms = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_spin_channels', 'number_of_dos_values'],
        description='''
        Values (number of states for a given energy divided by the numer of atoms, the set
        of discrete energy values is given in dos_energies) of density (electronic-energy
        or vibrational-energy) of states.
        ''')

    dos_values_per_unit_volume = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_spin_channels', 'number_of_dos_values'],
        description='''
        Values (number of states for a given energy divided by volume, the set of discrete
        energy values is given in dos_energies) of density (electronic-energy or
        vibrational-energy) of states.
        ''')

    dos_values = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_spin_channels', 'number_of_dos_values'],
        description='''
        Values (number of states for a given energy, the set of discrete energy values is
        given in dos_energies) of density (electronic-energy or vibrational-energy) of
        states. This refers to the simulation cell, i.e. integrating over all energies
        will give the number of electrons in the simulation cell.
        ''')

    number_of_dos_lms = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of $l$, $m$ combinations for the given projected density of
        states (DOS) in dos_values and dos_values_lm.
        ''')

    number_of_dos_values = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of energy values for the density of states (DOS), see
        dos_energies.
        ''')


class section_energy_code_independent(MSection):
    '''
    Type of the code-independent total energy (obtained by subtracting a reference energy
    calculated with the same code), created to be comparable among different codes and
    numerical settings. Details can be found on the [energy_code_independent wiki
    page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/energy-
    code-independent).
    '''

    energy_code_independent_kind = Quantity(
        type=str,
        shape=[],
        description='''
        Type of the code-independent total energy (obtained by subtracting a reference
        energy calculated with the same code), created to be comparable among different
        codes and numerical settings. Details can be found on the [energy_code_independent
        wiki page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-
        info/wikis/metainfo/energy-code-independent).
        ''')

    energy_code_independent_value = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Value of the code-independent total energy (obtained by subtracting a reference
        energy calculated with the same code). This value is created to be comparable
        among different codes and numerical settings. Details can be found on the
        [energy_code_independent wiki page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-
        meta-info/wikis/metainfo/energy-code-independent).
        ''')


class section_energy_van_der_Waals(MSection):
    '''
    Method used to compute van der Waals energy stored in energy_van_der_Waals_value. This
    metadata is used when more than one van der Waals method is applied in the same
    *single configuration calculation* (see section_single_configuration_calculation). The
    method used for van der Waals  (the one consistent with energy_current and, e.g., for
    evaluating the forces for a relaxation or dynamics) is defined in
    settings_van_der_Waals.
    '''

    energy_van_der_Waals_kind = Quantity(
        type=str,
        shape=[],
        description='''
        Method used to compute van der Waals energy stored in energy_van_der_Waals_value.
        This metadata is used when more than one van der Waals method is applied in the
        same *single configuration calculation* (see
        section_single_configuration_calculation). The method used for van der Waals  (the
        one consistent with energy_current and, e.g., for evaluating the forces for a
        relaxation or dynamics) is defined in settings_van_der_Waals.
        ''')

    energy_van_der_Waals_value = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='joule',
        description='''
        Value of van der Waals energy, calculated with the method defined in
        energy_van_der_Waals_kind. This metadata is used when more than one van der Waals
        method is applied in the same *single configuration calculation* (see
        section_single_configuration_calculation). The value of the van der Waals energy
        consistent with energy_current and used, e.g., for evaluating the forces for a
        relaxation or dynamics, is given in energy_van_der_Waals and defined in
        settings_van_der_Waals.
        ''')


class section_sampling_method(MSection):
    '''
    Kind of sampled ensemble stored in section_frame_sequence; valid values are defined in
    [ensemble_type wiki page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-
    info/wikis/metainfo/ensemble-type).
    '''

    ensemble_type = Quantity(
        type=str,
        shape=[],
        description='''
        Kind of sampled ensemble stored in section_frame_sequence; valid values are
        defined in [ensemble_type wiki page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-
        meta-info/wikis/metainfo/ensemble-type).
        ''')

    sampling_method_expansion_order = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Order up to which the potential energy surface was expanded in a Taylor series
        (see sampling_method).
        ''')

    sampling_method = Quantity(
        type=str,
        shape=[],
        description='''
        Type of method used to do the sampling.

        Allowed values are:

        | Sampling method                | Description                      |

        | ------------------------------ | -------------------------------- |

        | `"geometry_optimization"`      | Geometry optimization            |

        | `"molecular_dynamics"`         | Molecular dynamics               |

        | `"montecarlo"`                 | (Metropolis) Monte Carlo         |

        | `"steered_molecular_dynamics"` | Steered molecular dynamics (with time dependent
        external forces) |

        | `"meta_dynamics"`              | Biased molecular dynamics with history-
        dependent Hamiltonian |

        | `"wang_landau_montecarlo"`     | Monte Carlo according to the Wang-Landau
        formulation. |

        | `"blue_moon"`                  | Blue Moon sampling               |

        | `"langevin_dynamics"`          | Langevin dynamics                |

        | `"taylor_expansion"`           | Taylor expansion of the potential energy
        surface |
        ''')


class section_original_system(MSection):
    '''
    Gives a mapping table of atoms to symmetrically independent atoms in the original
    cell. This is used to find symmetrically equivalent atoms.
    '''

    equivalent_atoms_original = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_atoms'],
        description='''
        Gives a mapping table of atoms to symmetrically independent atoms in the original
        cell. This is used to find symmetrically equivalent atoms.
        ''')

    wyckoff_letters_original = Quantity(
        type=str,
        shape=['number_of_atoms'],
        description='''
        Wyckoff letters for atoms in the original cell.
        ''')


class section_frame_sequence(MSection):
    '''
    Array containing the strictly increasing indices of the frames the
    frame_sequence_conserved_quantity values refers to. If not given it defaults to the
    trivial mapping 0,1,...
    '''

    frame_sequence_conserved_quantity_frames = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_conserved_quantity_evaluations_in_sequence'],
        description='''
        Array containing the strictly increasing indices of the frames the
        frame_sequence_conserved_quantity values refers to. If not given it defaults to
        the trivial mapping 0,1,...
        ''')

    frame_sequence_conserved_quantity_stats = Quantity(
        type=np.dtype(np.float32),
        shape=[2],
        unit='joule',
        description='''
        Average value of energy-like frame_sequence_conserved_quantity, and its standard
        deviation, over this sequence of frames (i.e., a trajectory, a frame is one
        section_single_configuration_calculation).
        ''')

    frame_sequence_conserved_quantity = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_conserved_quantity_evaluations_in_sequence'],
        unit='joule',
        description='''
        Array containing the values of a quantity that should be conserved,  along a
        sequence of frames (i.e., a trajectory). A frame is one
        section_single_configuration_calculation), for example the total energy in the NVE
        ensemble. If not all frames have a value the indices of the frames that have a
        value are stored in frame_sequence_conserved_quantity_frames.
        ''')

    frame_sequence_continuation_kind = Quantity(
        type=MProxy('section_frame_sequence'),
        shape=[],
        description='''
        Type of continuation that has been performed from the previous sequence of frames
        (i.e., a trajectory, a frame is one section_single_configuration_calculation),
        upon restart.
        ''')

    frame_sequence_external_url = Quantity(
        type=str,
        shape=[],
        description='''
        If the energy, forces, and other quantities for the frames (a frame is one
        section_single_configuration_calculation) in  section_frame_sequence are obtained
        by calling a different code than the code that drives the sequence (e.g., a
        wrapper that drives a molecular dynamics, Monte Carlo, geometry optimization and
        calls an electronic-structure code for energy and forces for each configuration),
        this metadata holds the reference to where the
        section_single_configuration_calculation for each frame are located. The format
        for this reference is described in the [frame_sequence_external_url wiki
        page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/frame-
        sequence-external-url).
        ''')

    frame_sequence_kinetic_energy_frames = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_kinetic_energies_in_sequence'],
        description='''
        Array containing the strictly increasing indices referring to the frames of
        frame_sequence_kinetic_energy. If not given it defaults to the trivial mapping
        0,1,...
        ''')

    frame_sequence_kinetic_energy_stats = Quantity(
        type=np.dtype(np.float32),
        shape=[2],
        unit='joule',
        description='''
        Average kinetic energy and its standard deviation over this sequence of frames
        (i.e., a trajectory, a frame is one section_single_configuration_calculation).
        ''')

    frame_sequence_kinetic_energy = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_kinetic_energies_in_sequence'],
        unit='joule',
        description='''
        Array containing the values of the kinetic energy along this sequence of frames
        (i.e., a trajectory, a frame is one section_single_configuration_calculation). If
        not all frames have a value the indices of the frames that have a value are stored
        in frame_sequence_kinetic_energy_frames.
        ''')

    frame_sequence_local_frames_ref = Quantity(
        type=MProxy('section_single_configuration_calculation'),
        shape=['number_of_frames_in_sequence'],
        description='''
        Reference from each frame (a frame is one
        section_single_configuration_calculation) in this section_frame_sequence to the
        corresponding section_single_configuration_calculation. Each
        section_frame_sequence binds a collection of
        section_single_configuration_calculation, because they all belong to, e.g., a
        molecular dynamics trajectory, or geometry optimization. The full information for
        each frame is stored in section_single_configuration_calculation and this metadata
        establishes the link for each frame.
        ''')

    frame_sequence_potential_energy_frames = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_potential_energies_in_sequence'],
        description='''
        Array containing the strictly increasing indices referring to the frames of
        frame_sequence_potential_energy. If not given it defaults to the trivial mapping
        0,1,...
        ''')

    frame_sequence_potential_energy_stats = Quantity(
        type=np.dtype(np.float32),
        shape=[2],
        unit='joule',
        description='''
        Average potential energy and its standard deviation over this sequence of frames
        (i.e., a trajectory, a frame is one section_single_configuration_calculation).
        ''')

    frame_sequence_potential_energy = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_potential_energies_in_sequence'],
        unit='joule',
        description='''
        Array containing the value of the potential energy along this sequence of frames
        (i.e., a trajectory, a frame is one section_single_configuration_calculation).
        This is equal to energy_total of the corresponding
        section_single_configuration_calculation and repeated here in a summary array for
        easier access. If not all frames have a value the indices of the frames that have
        a value are stored in frame_sequence_potential_energy_frames.
        ''')

    frame_sequence_pressure_frames = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_pressure_evaluations_in_sequence'],
        description='''
        Array containing the strictly increasing indices referring to the frames of
        frame_sequence_pressure. If not given it defaults to the trivial mapping 0,1,...
        ''')

    frame_sequence_pressure_stats = Quantity(
        type=np.dtype(np.float32),
        shape=[2],
        unit='pascal',
        description='''
        Average pressure (one third of the trace of the stress tensor) and standard
        deviation over this sequence of frames (i.e., a trajectory, a frame is one
        section_single_configuration_calculation).
        ''')

    frame_sequence_pressure = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_pressure_evaluations_in_sequence'],
        unit='pascal',
        description='''
        Array containing the values of the pressure (one third of the trace of the stress
        tensor) along this sequence of frames (a frame is one
        section_single_configuration_calculation). If not all frames have a value the
        indices of the frames that have a value are stored in
        frame_sequence_pressure_frames.
        ''')

    frame_sequence_temperature_frames = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_temperatures_in_sequence'],
        description='''
        Array containing the strictly increasing indices referring to the frames of
        frame_sequence_temperature. If not given it defaults to the trivial mapping
        0,1,...
        ''')

    frame_sequence_temperature_stats = Quantity(
        type=np.dtype(np.float32),
        shape=[2],
        unit='kelvin',
        description='''
        Average temperature and its standard deviation over this sequence of frames (i.e.,
        a trajectory, a frame is one section_single_configuration_calculation).
        ''')

    frame_sequence_temperature = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_temperatures_in_sequence'],
        unit='kelvin',
        description='''
        Array containing the values of the instantaneous temperature (a quantity,
        proportional to frame_sequence_kinetic_energy, whose ensemble average equals the
        thermodynamic temperature) along this sequence of frames (i.e., a trajectory, a
        frame is one section_single_configuration_calculation). If not all frames have a
        value the indices of the frames that have a value are stored in
        frame_sequence_temperature_frames.
        ''')

    frame_sequence_time = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_frames_in_sequence'],
        unit='second',
        description='''
        Time along this sequence of frames (i.e., a trajectory, a frame is one
        section_single_configuration_calculation). Time start is arbitrary, but when a
        sequence is a continuation of another time should be continued too.
        ''')

    frame_sequence_to_sampling_ref = Quantity(
        type=MProxy('section_sampling_method'),
        shape=[],
        description='''
        Reference from the present section_frame_sequence to the section_sampling_method,
        that defines the parameters used in this sequence of frames (i.e., a trajectory, a
        frame is one section_single_configuration_calculation).
        ''')

    geometry_optimization_converged = Quantity(
        type=np.dtype(np.int8),
        shape=[],
        description='''
        Arrays specify whether a geometry optimization is converged.
        ''')

    number_of_conserved_quantity_evaluations_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of conserved quantity evaluations in this sequence. A sequence is
        a trajectory, which can have number_of_frames_in_sequence each representing one
        section_single_configuration_calculation section.
        ''')

    number_of_frames_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of frames in a sequence. A sequence is a trajectory, which can
        have number_of_frames_in_sequence each representing one
        section_single_configuration_calculation section.
        ''')

    number_of_kinetic_energies_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of kinetic energy evaluations in this sequence of frames, see
        frame_sequence_kinetic_energy.
        ''')

    number_of_potential_energies_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of potential energies evaluation in this sequence. A sequence is
        a trajectory, which can have number_of_frames_in_sequence each representing one
        section_single_configuration_calculation section.
        ''')

    number_of_pressure_evaluations_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of pressure evaluations in this sequence. A sequence is a
        trajectory, which can have number_of_frames_in_sequence each representing one
        section_single_configuration_calculation section.
        ''')

    number_of_temperatures_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of temperature frames (frame_sequence_temperature) used in the
        section_frame_sequence. A sequence is a trajectory, which can have
        number_of_frames_in_sequence each representing one
        section_single_configuration_calculation section.
        ''')

    previous_sequence_ref = Quantity(
        type=MProxy('section_frame_sequence'),
        shape=[],
        description='''
        Contains a reference to the previous sequence. A sequence is a trajectory, which
        can have number_of_frames_in_sequence each representing one
        section_single_configuration_calculation section. If not given, a start from an
        initial configuration is assumed.
        ''')

    section_frame_sequence_user_quantity = SubSection(
        sub_section=MProxy('section_frame_sequence_user_quantity'),
        repeats=False)

    section_thermodynamical_properties = SubSection(
        sub_section=MProxy('section_thermodynamical_properties'),
        repeats=False)


class section_frame_sequence_user_quantity(MSection):
    '''
    Array containing the strictly increasing indices referring to the frames of
    frame_sequence_user_quantity. If not given it defaults to the trivial mapping 0,1,...
    '''

    frame_sequence_user_quantity_frames = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_user_quantity_evaluations_in_sequence'],
        description='''
        Array containing the strictly increasing indices referring to the frames of
        frame_sequence_user_quantity. If not given it defaults to the trivial mapping
        0,1,...
        ''')

    frame_sequence_user_quantity_name = Quantity(
        type=str,
        shape=[],
        description='''
        Descriptive name of a user-defined quantity, sampled along this sequence of frames
        (i.e., a trajectory, a frame is one section_single_configuration_calculation).
        Dedicated metadata are created for the conserved energy-like quantity
        (frame_sequence_conserved_quantity), the kinetic and potential energies
        (frame_sequence_kinetic_energy and frame_sequence_potential_energy), the
        instantaneous temperature (frame_sequence_temperature) and pressure
        (frame_sequence_pressure). This metadata should be used for other quantities that
        are monitored along a sequence of frames.
        ''')

    frame_sequence_user_quantity_stats = Quantity(
        type=np.dtype(np.float32),
        shape=[2, 'number_of_frame_sequence_user_quantity_components'],
        description='''
        Average of frame_sequence_user_quantity and its standard deviation in this
        sequence of frames (i.e., a trajectory, a frame is one
        section_single_configuration_calculation).
        ''')

    frame_sequence_user_quantity = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_user_quantity_evaluations_in_sequence', 'number_of_frame_sequence_user_quantity_components'],
        description='''
        Array containing the values of the user-defined quantity defined in
        frame_sequence_user_quantity_name, evaluated along this sequence of frames (i.e.,
        trajectory, a frame is one section_single_configuration_calculation). If not all
        frames have a value the indices of the frames that have a value are stored in
        frame_sequence_kinetic_energy_frames. If not all frames have a value the indices
        of the frames that have a value are stored in
        frame_sequence_kinetic_energy_frames.
        ''')

    number_of_frame_sequence_user_quantity_components = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of user-defined quantity defined by
        frame_sequence_user_quantity_name and monitored in a sequence of frames. A
        sequence is a trajectory, which can have number_of_frames_in_sequence each
        representing one section_single_configuration_calculation section.

        Dedicated metadata monitored along a sequence of frames are created for the
        conserved energy-like quantity (frame_sequence_conserved_quantity), the kinetic
        and potential energies ([frame_sequence_kinetic_energy and
        frame_sequence_potential_energy](frame_sequence_kinetic_energy and
        frame_sequence_potential_energy)), the instantaneous temperature
        (frame_sequence_temperature) and the pressure (frame_sequence_pressure).
        ''')

    number_of_user_quantity_evaluations_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of user defined quantity evaluations along a sequence of
        frame_sequence_user_quantity frames. A sequence is a trajectory, which can have
        number_of_frames_in_sequence each representing one
        section_single_configuration_calculation section.
        ''')


class section_gaussian_basis_group(MSection):
    '''
    contraction coefficients $c_{i j}$ defining the contracted basis functions with
    respect to *normalized* primitive Gaussian functions. They define the Gaussian basis
    functions as described in section_gaussian_basis_group.
    '''

    gaussian_basis_group_contractions = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_gaussian_basis_group_contractions', 'number_of_gaussian_basis_group_exponents'],
        description='''
        contraction coefficients $c_{i j}$ defining the contracted basis functions with
        respect to *normalized* primitive Gaussian functions. They define the Gaussian
        basis functions as described in section_gaussian_basis_group.
        ''')

    gaussian_basis_group_exponents = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_gaussian_basis_group_exponents'],
        unit='1 / meter ** 2',
        description='''
        Exponents $\\alpha_j$ of the Gaussian functions defining this basis set
        $exp(-\\alpha_j r^2)$. One should be careful about the units of the coefficients.
        ''')

    gaussian_basis_group_ls = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_gaussian_basis_group_contractions'],
        description='''
        Azimuthal quantum number ($l$) values (of the angular part given by the spherical
        harmonic $Y_{l m}$ of the various contracted basis functions).
        ''')

    number_of_gaussian_basis_group_contractions = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of different contractions, i.e. resulting basis functions in a
        section_gaussian_basis_group section.
        ''')

    number_of_gaussian_basis_group_exponents = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of different Gaussian exponents in a section_gaussian_basis_group
        section.
        ''')


class section_thermodynamical_properties(MSection):
    '''
    Stores the Helmholtz free energy per unit cell at constant volume of a thermodynamic
    calculation.
    '''

    helmholz_free_energy = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_thermodynamical_property_values'],
        unit='joule',
        description='''
        Stores the Helmholtz free energy per unit cell at constant volume of a
        thermodynamic calculation.
        ''')

    number_of_thermodynamical_property_values = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of thermal properties values available in
        section_thermodynamical_properties.
        ''')

    thermodynamical_properties_calculation_method = Quantity(
        type=str,
        shape=[],
        description='''
        Method used to calculate the thermodynamic quantities.

        Valid values:

        * harmonic
        ''')

    thermodynamical_property_heat_capacity_C_v = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_thermodynamical_property_values'],
        unit='joule / kelvin',
        description='''
        Stores the heat capacity per cell unit at constant volume.
        ''')

    thermodynamical_property_temperature = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_thermodynamical_property_values'],
        unit='kelvin',
        description='''
        Specifies the temperatures at which properties such as the Helmholtz free energy
        are calculated.
        ''')

    vibrational_free_energy_at_constant_volume = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_thermodynamical_property_values'],
        unit='joule',
        description='''
        Holds the vibrational free energy per atom at constant volume.
        ''')


class section_k_band_normalized(MSection):
    '''
    If the normalized path is along the default path defined in W. Setyawan and S.
    Curtarolo, [Comput. Mater. Sci. **49**, 299-312
    (2010)](http://dx.doi.org/10.1016/j.commatsci.2010.05.010).
    '''

    k_band_path_normalized_is_standard = Quantity(
        type=np.dtype(np.int8),
        shape=[],
        description='''
        If the normalized path is along the default path defined in W. Setyawan and S.
        Curtarolo, [Comput. Mater. Sci. **49**, 299-312
        (2010)](http://dx.doi.org/10.1016/j.commatsci.2010.05.010).
        ''')

    section_k_band_segment_normalized = SubSection(
        sub_section=MProxy('section_k_band_segment_normalized'),
        repeats=True)


class section_method_atom_kind(MSection):
    '''
    Atomic number (number of protons) of this atom kind, use 0 if not an atom.
    '''

    method_atom_kind_atom_number = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Atomic number (number of protons) of this atom kind, use 0 if not an atom.
        ''')

    method_atom_kind_explicit_electrons = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        description='''
        Number of explicit electrons (often called valence).
        ''')

    method_atom_kind_label = Quantity(
        type=str,
        shape=[],
        description='''
        String used to identify the atoms of this kind. This should correspond to the
        atom_labels of the configuration. It is possible for one atom kind to have
        multiple labels (in order to allow two atoms of the same kind to have two
        differently defined sets of atom-centered basis functions or two different pseudo-
        potentials). Atom kind is typically the symbol of the atomic species but it can be
        also a ghost or pseudo-atom.
        ''')

    method_atom_kind_mass = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        unit='atomic_mass_unit',
        description='''
        Mass of the kind of this kind of atoms.
        ''')

    method_atom_kind_pseudopotential_name = Quantity(
        type=str,
        shape=[],
        description='''
        Name identifying the pseudopotential used.
        ''')


class section_method_to_method_refs(MSection):
    '''
    URL used to reference an externally stored section_method. The kind of relationship
    between the present and the referenced section_method is specified by
    method_to_method_kind.
    '''

    method_to_method_external_url = Quantity(
        type=str,
        shape=[],
        description='''
        URL used to reference an externally stored section_method. The kind of
        relationship between the present and the referenced section_method is specified by
        method_to_method_kind.
        ''')

    method_to_method_kind = Quantity(
        type=str,
        shape=[],
        description='''
        String defining the kind of relationship that the referenced section_method has
        with the present section_method. Valid values are described in the
        [method_to_method_kind wiki page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-
        meta-info/wikis/metainfo/method-to-method-kind). Often calculations are connected,
        for instance, one calculation is a perturbation performed using a self-consistent
        field (SCF) calculation as starting point, or a simulated system is partitioned in
        regions with different but connected Hamiltonians (e.g., QM/MM, or a region
        treated via Kohn-Sham DFT embedded into a region treated via orbital-free DFT).
        Hence, the need of keeping track of these connected calculations. The referenced
        section_method is identified via method_to_method_ref (typically used for a
        section_method in the same section_run) or method_to_method_external_url.
        ''')

    method_to_method_ref = Quantity(
        type=int,
        shape=[],
        description='''
        Reference to a local section_method. If both method_to_method_ref and
        method_to_method_external_url are given, then method_to_method_ref is a local copy
        of the value contained in method_to_method_external_url. The kind of relationship
        between the method defined in the present section_method and the referenced one is
        described by method_to_method_kind.
        ''')


class section_species_projected_dos(MSection):
    '''
    Gives the number of $l$, $m$ combinations for the species-projected density of states
    (DOS) defined in section_species_projected_dos.
    '''

    number_of_lm_species_projected_dos = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of $l$, $m$ combinations for the species-projected density of
        states (DOS) defined in section_species_projected_dos.
        ''')

    number_of_species_projected_dos_values = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of energy values for the species-projected density of states
        (DOS) defined in section_species_projected_dos.
        ''')

    number_of_species = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of species for the species-projected density of states (DOS)
        defined in section_species_projected_dos.
        ''')

    species_projected_dos_energies_normalized = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_species_projected_dos_values'],
        unit='joule',
        description='''
        Contains the set of discrete energy values with respect to the top of the valence
        band for the species-projected density of states (DOS). It is derived from the
        species_projected_dos_energies species field.
        ''')

    species_projected_dos_energies = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_species_projected_dos_values'],
        unit='joule',
        description='''
        Contains the set of discrete energy values for the species-projected density of
        states (DOS).
        ''')

    species_projected_dos_lm = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_lm_species_projected_dos', 2],
        description='''
        Consists of tuples of $l$ and $m$ values for all given values in the
        species_projected_dos_values_lm species field.

        The quantum number $l$ represents the azimuthal quantum number, whereas for the
        quantum number $m$, besides the conventional use as magnetic quantum number ($l+1$
        integer values from $-l$ to $l$), a set of different conventions is accepted. The
        adopted convention is specified by atom_projected_dos_m_kind.
        ''')

    species_projected_dos_m_kind = Quantity(
        type=str,
        shape=[],
        description='''
        Specifies the kind of the integer numbers $m$ used in species_projected_dos_lm.

        Allowed values are listed in the [m_kind wiki
        page](https://gitlab.rzg.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/m-kind)
        and can be (quantum) numbers of

        * spherical

        * polynomial

        * real_orbital

        * integrated

        functions or values.
        ''')

    species_projected_dos_species_label = Quantity(
        type=str,
        shape=['number_of_species'],
        description='''
        Contains labels of the atomic species for the species-projected density of states
        (DOS).

        Differently from atom_labels, which allow more than one label for the same atomic
        species (by adding a number or a string to the label), this list is expected to
        refer to actual atomic species, i.e. belonging to the periodic table of elements.
        Thus, the species-projected DOS are expected to be as many as the different atomic
        species in the system.
        ''')

    species_projected_dos_values_lm = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_lm_species_projected_dos', 'number_of_spin_channels', 'number_of_species', 'number_of_species_projected_dos_values'],
        description='''
        Holds species-projected density of states (DOS) values, divided into contributions
        from each $l,m$ channel.

        Here, there are as many species-projected DOS as the number of species,
        number_of_species. The list of labels of the species is given in
        species_projected_dos_species_label.
        ''')

    species_projected_dos_values_total = Quantity(
        type=np.dtype(np.float32),
        shape=['number_of_spin_channels', 'number_of_species', 'number_of_species_projected_dos_values'],
        description='''
        Holds species-projected density of states (DOS) values, summed up over all
        azimuthal quantum numbers $l$.

        Here, there are as many species-projected DOS as the number of species,
        number_of_species. The list of labels of the species is given in
        species_projected_dos_species_label.
        ''')


class section_processor_info(MSection):
    '''
    Id (name+version) of the processor that generated or added information to the current
    calculation.
    '''

    processor_id = Quantity(
        type=str,
        shape=[],
        description='''
        Id (name+version) of the processor that generated or added information to the
        current calculation.
        ''')

    processor_number_of_evaluated_contexts = Quantity(
        type=np.dtype(np.int64),
        shape=[],
        description='''
        number of contexts evaluated with this processor in the current current
        calculation.
        ''')

    processor_number_of_failed_contexts = Quantity(
        type=np.dtype(np.int64),
        shape=[],
        description='''
        number of contexts in the current current calculation that had failure for this
        processor.
        ''')

    processor_number_of_skipped_contexts = Quantity(
        type=np.dtype(np.int64),
        shape=[],
        description='''
        number of contexts skipped by this processor in the current current calculation.
        ''')

    processor_number_of_successful_contexts = Quantity(
        type=np.dtype(np.int64),
        shape=[],
        description='''
        number of contexts in the current calculation that where successfully handled by
        this processor.
        ''')

    processor_version_details = Quantity(
        type=typing.Any,
        shape=[],
        description='''
        detailed version information on the processor that generated or added information
        to the current calculation.
        ''')


class section_processor_log_event(MSection):
    '''
    Level of the logging, a lower number has more priority. The levels are the same as
    log4j: FATAL -> 100, ERROR -> 200, WARN -> 300, INFO -> 400, DEBUG -> 500, TRACE ->
    600
    '''

    processor_log_event_level = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Level of the logging, a lower number has more priority. The levels are the same as
        log4j: FATAL -> 100, ERROR -> 200, WARN -> 300, INFO -> 400, DEBUG -> 500, TRACE
        -> 600
        ''')

    processor_log_event_message = Quantity(
        type=str,
        shape=[],
        description='''
        The log message
        ''')


class section_processor_log(MSection):
    '''
    The processor id of the processor creating this log
    '''

    processor_log_processor_id = Quantity(
        type=str,
        shape=[],
        description='''
        The processor id of the processor creating this log
        ''')

    processor_log_start = Quantity(
        type=str,
        shape=[],
        description='''
        Start of the log (in ansi notation YYYY-MM-TT...)
        ''')

    section_processor_log_event = SubSection(
        sub_section=MProxy('section_processor_log_event'),
        repeats=False)


class section_prototype(MSection):
    '''
    AFLOW id of the prototype (see
    http://aflowlib.org/CrystalDatabase/prototype_index.html) identified on the basis of
    the space_group and normalized_wyckoff.
    '''

    prototype_aflow_id = Quantity(
        type=str,
        shape=[],
        description='''
        AFLOW id of the prototype (see
        http://aflowlib.org/CrystalDatabase/prototype_index.html) identified on the basis
        of the space_group and normalized_wyckoff.
        ''')

    prototype_aflow_url = Quantity(
        type=str,
        shape=[],
        description='''
        Url to the AFLOW definition of the prototype (see
        http://aflowlib.org/CrystalDatabase/prototype_index.html) identified on the basis
        of the space_group and normalized_wyckoff.
        ''')

    prototype_assignement_method = Quantity(
        type=str,
        shape=[],
        description='''
        Method used to identify the prototype
        ''')

    prototype_label = Quantity(
        type=str,
        shape=[],
        description='''
        Label of the prototype identified on the basis of the space_group and
        normalized_wyckoff. The label is in the same format as in the read_prototypes
        function: <space_group_number>-<prototype_name>-<Pearson's symbol>).
        ''')


class section_basis_functions_atom_centered(MSection):
    '''
    This section contains the description of the basis functions (at least one function)
    of the (atom-centered) basis set defined in section_basis_set_atom_centered.
    '''


class section_springer_classification(MSection):
    '''
    Section_springer_classsification contains a classification tag of a material according
    to Springer Materials
    '''

    springer_classification = Quantity(
        type=str,
        shape=[],
        description='''
        Contains the classification name of the current material according to Springer
        Materials
        ''')

    springer_number_of_classification_reference_per_material = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of publications/references using this classification for the current
        material in the Springer Materials database
        ''')


class section_springer_material(MSection):
    '''
    Section_springer_classsification contains a classification tag of a material according
    to Springer Materials
    '''

    springer_formula = Quantity(
        type=str,
        shape=[],
        description='''
        The formula of current material according to Springer Materials
        ''')

    springer_space_group_number = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Information about the space group number of current material according to Springer
        Materials
        ''')

    section_springer_classification = SubSection(
        sub_section=MProxy('section_springer_classification'),
        repeats=False)

    section_springer_compound_class = SubSection(
        sub_section=MProxy('section_springer_compound_class'),
        repeats=False)

    section_springer_id = SubSection(
        sub_section=MProxy('section_springer_id'),
        repeats=False)

    section_springer_references = SubSection(
        sub_section=MProxy('section_springer_references'),
        repeats=False)


class section_springer_compound_class(MSection):
    '''
    Description of a compound class (according to Springer Materials) of the current
    material. This is a property of the chemical formula of the compound
    '''

    springer_compound_class = Quantity(
        type=str,
        shape=[],
        description='''
        Name of a class of the current compound, as defined in by Springer Materials. This
        is a property of the chemical formula of the compound
        ''')

    springer_number_of_compound_class_reference_per_material = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of publications/references using this compound class for the current
        compound in the Springer Materials database
        ''')


class section_springer_id(MSection):
    '''
    Identifiers used by Springer Materials
    '''

    springer_id = Quantity(
        type=str,
        shape=[],
        description='''
        Id of the classified material according to Springer Materials
        ''')

    springer_url = Quantity(
        type=str,
        shape=[],
        description='''
        Url to the source page in Springer Materials describing the current entry
        ''')


class section_springer_references(MSection):
    '''
    Contains the information about references related to current material according to
    Springer Materials
    '''

    springer_reference = Quantity(
        type=str,
        shape=[],
        description='''
        Contains the information about references related to current material according to
        Springer Materials
        ''')


class section_stress_tensor(MSection):
    '''
    Section collecting alternative values to stress_tensor that have been calculated.

    This section allows the storage of multiple definitions and evaluated values of the
    stress tensor, while only one definition is used for, e.g., molecular dynamics or
    geometry optimization (if needed).
    '''

    stress_tensor_kind = Quantity(
        type=str,
        shape=[],
        description='''
        Specifies the method used to compute the stress tensor stored in
        stress_tensor_value. This is an *alternative* to the stress tensor defined in
        stress_tensor_method, which is stored in stress_tensor.

        This field allows for multiple definitions and evaluated values of the stress
        tensor, while only one definition is used for, e.g., molecular dynamics and
        geometry optimization.
        ''')

    stress_tensor_value = Quantity(
        type=np.dtype(np.float32),
        shape=[3, 3],
        unit='pascal',
        description='''
        Contains the value of the stress tensor of the kind defined in stress_tensor_kind.
        This is an *alternative* to the stress tensor defined in stress_tensor_method.

        This field allows for multiple definitions and evaluated values of the stress
        tensor, while only one definition is used for, e.g., molecular dynamics and
        geometry optimization.
        ''')


class section_system_to_system_refs(MSection):
    '''
    Section that describes the relationship between different section_system sections.

    For instance, if a phonon calculation using a finite difference approach is performed
    the force evaluation is typically done in a larger supercell but the properties such
    as the phonon band structure are still calculated for the primitive cell.

    The kind of relationship between the system defined in this section and the referenced
    one is described by system_to_system_kind. The referenced section_system is identified
    via system_to_system_ref.
    '''

    system_to_system_kind = Quantity(
        type=str,
        shape=[],
        description='''
        String defining the relationship between the referenced section_system and the
        present section_system. Often systems are connected for example if a phonon
        calculation using finite differences is performed the force ealuation is done in a
        larger supercell but properties such as the phonon band structure are still
        calculated for the primitive cell. Hence, the need of keeping track of these
        connected systems. The referenced system is identified via system_to_system_ref.
        ''')

    system_to_system_ref = Quantity(
        type=MProxy('section_system'),
        shape=[],
        description='''
        Reference to another system. The kind of relationship between the present and the
        referenced section_system is specified by system_to_system_kind.
        ''')


class section_volumetric_data(MSection):
    '''
    Section defining a set of volumetric data on a uniform real-space

    grid.

    To store an array (e.g. a density or a potential), define:

    * three grid point displacement vectors ("displacements")

    * number of grid points along each axis ("nx", "ny" and "nz")

    * the origin of the coordinate system, i.e. coordinates of the first grid

    point ("origin")

    * how many spatial functions are represented, e.g., two for a

    normal spin-polarized density ("multiplicity")

    * the values for each grid point ("values")

    * the unit that applies to each value ("units")

    * the kind of array represented by the volumetric data ("kind").

    Allowed kinds are (please add new kinds as necessary): "density",

    "potential_hartree" and "potential_effective".  Densities and

    potentials that are spin-polarized should have multiplicity two.

    Rules for more complex spins are to be decided when necessary.
    '''

    volumetric_data_displacements = Quantity(
        type=np.dtype(np.float32),
        shape=[3, 3],
        unit='meter',
        description='''
        displacement vectors between grid points along each axis; same indexing rules as
        lattice_vectors.  In many cases, displacements and number of points are related to
        lattice_vectors through: [displacement] * [number of points + N] =
        [lattice_vector],where N is 1 for periodic directions and 0 for non-periodic ones
        ''')

    volumetric_data_kind = Quantity(
        type=str,
        shape=[],
        description='''
        The kind of function, e.g. density, potential_hartree, potential_effective.  The
        unit of measurement for "volumetric_data_values" depends on the kind: Densities
        are 1/m^3 and potentials are J/m^3.  See [full specification on the
        wiki](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-
        info/wikis/metainfo/volumetric-data).
        ''')

    volumetric_data_multiplicity = Quantity(
        type=int,
        shape=[],
        description='''
        number of functions stored
        ''')

    volumetric_data_nx = Quantity(
        type=int,
        shape=[],
        description='''
        number of points along x axis
        ''')

    volumetric_data_ny = Quantity(
        type=int,
        shape=[],
        description='''
        number of points along y axis
        ''')

    volumetric_data_nz = Quantity(
        type=int,
        shape=[],
        description='''
        number of points along z axis
        ''')

    volumetric_data_origin = Quantity(
        type=np.dtype(np.float32),
        shape=[3],
        description='''
        location of the first grid point; same coordinate system as atom_positions when
        applicable.
        ''')

    volumetric_data_values = Quantity(
        type=np.dtype(np.float32),
        shape=['volumetric_data_multiplicity', 'volumetric_data_nx', 'volumetric_data_ny', 'volumetric_data_nz'],
        description='''
        Array of shape (multiplicity, nx, ny, nz) containing the values.  The units of
        these values depend on which kind of data the values represent; see
        "volumetric_data_kind".
        ''')


class section_XC_functionals(MSection):
    '''
    Section containing one of the exchange-correlation (XC) functionals for the present
    section_method that are combined to form the XC_functional.
    '''

    XC_functional_name = Quantity(
        type=str,
        shape=[],
        description='''
        Provides the name of one of the exchange and/or correlation (XC) functionals
        combined in XC_functional.

        The valid unique names that can be used are listed in the [XC_functional wiki
        page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/XC-
        functional).

        *NOTE*: This value should refer to a correlation, an exchange or an exchange-
        correlation functional only.
        ''')

    XC_functional_parameters = Quantity(
        type=typing.Any,
        shape=[],
        description='''
        Contains an associative list of non-default values of the parameters for the
        functional declared in XC_functional_name of the section_XC_functionals section.

        For example, if a calculations using a hybrid XC functional (e.g., HSE06)
        specifies a user-given value of the mixing parameter between exact and GGA
        exchange, then this non-default value is stored in this metadata.

        The labels and units of these values are defined in the paragraph dedicated to the
        specified functional declared in XC_functional_name of the [XC_functional wiki
        page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/XC-
        functional).

        If this metadata is not given, the default parameter values for the
        XC_functional_name are assumed.
        ''')

    XC_functional_weight = Quantity(
        type=np.dtype(np.float32),
        shape=[],
        description='''
        Provides the value of the weight for the exchange, correlation, or exchange-
        correlation functional declared in XC_functional_name (see
        section_XC_functionals).

        This weight is used in the linear combination of the different XC functional names
        (XC_functional_name) in different section_XC_functionals sections to form the
        XC_functional used for evaluating energy_XC_functional and related quantities.

        If not specified then the default is set to 1.
        ''')
