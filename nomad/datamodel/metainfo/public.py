import numpy as np            # pylint: disable=unused-import
import typing                 # pylint: disable=unused-import
from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference, MEnum, derived
)
from nomad.metainfo.search_extension import Search
from nomad.metainfo.legacy import LegacyDefinition


m_package = Package(
    name='public_nomadmetainfo_json',
    description='None',
    a_legacy=LegacyDefinition(name='public.nomadmetainfo.json'))


class accessory_info(MCategory):
    '''
    Information that *in theory* should not affect the results of the calculations (e.g.,
    timing).
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='accessory_info'))


class atom_forces_type(MCategory):
    '''
    The types of forces acting on the atoms (i.e., minus derivatives of the specific type
    of energy with respect to the atom position).
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='atom_forces_type'))


class basis_set_description(MCategory):
    '''
    One of the parts building the basis set of the system (e.g., some atom-centered basis
    set, plane-waves or both).
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='basis_set_description'))


class configuration_core(MCategory):
    '''
    Properties defining the current configuration.
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='configuration_core'))


class conserved_quantity(MCategory):
    '''
    A quantity that is preserved during the time propagation (for example,
    kinetic+potential energy during NVE).
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='conserved_quantity'))


class energy_value(MCategory):
    '''
    This metadata stores an energy value.
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='energy_value'))


class error_estimate_contribution(MCategory):
    '''
    An estimate of a partial quantity contributing to the error for a given quantity.
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='error_estimate_contribution'))


class message_debug(MCategory):
    '''
    A debugging message of the computational program.
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='message_debug'))


class parsing_message_debug(MCategory):
    '''
    This field is used for debugging messages of the parsing program.
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='parsing_message_debug'))


class scf_info(MCategory):
    '''
    Contains information on the self-consistent field (SCF) procedure, i.e. the number of
    SCF iterations (number_of_scf_iterations) or a section_scf_iteration section with
    detailed information on the SCF procedure of specified quantities.
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='scf_info'))


class settings_numerical_parameter(MCategory):
    '''
    A parameter that can influence the convergence, but not the physics (unlike
    settings_physical_parameter)
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='settings_numerical_parameter'))


class settings_physical_parameter(MCategory):
    '''
    A parameter that defines the physical model used. Use settings_numerical_parameter for
    parameters that that influence only the convergence/accuracy.
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='settings_physical_parameter'))


class settings_potential_energy_surface(MCategory):
    '''
    Contains parameters that control the potential energy surface.
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='settings_potential_energy_surface'))


class settings_run(MCategory):
    '''
    Contains parameters that control the whole run (but not the *single configuration
    calculation*, see section_single_configuration_calculation).
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='settings_run'))


class settings_sampling(MCategory):
    '''
    Contains parameters controlling the sampling.
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='settings_sampling'))


class settings_scf(MCategory):
    '''
    Contains parameters connected with the convergence of the self-consistent field (SCF)
    iterations.
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='settings_scf'))


class settings_smearing(MCategory):
    '''
    Contain parameters that control the smearing of the orbital occupation at finite
    electronic temperatures.
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='settings_smearing'))


class settings_stress_tensor(MCategory):
    '''
    Settings to calculate the stress tensor (stress_tensor) consistent with the total
    energy of the system given in energy_total.
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='settings_stress_tensor'))


class stress_tensor_type(MCategory):
    '''
    Contains the final value of the default stress tensor (stress_tensor) and/or the value
    of the stress tensor (stress_tensor_value) of the kind defined in stress_tensor_kind.
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='stress_tensor_type'))


class energy_component_per_atom(MCategory):
    '''
    A value of an energy component per atom, concurring in defining the total energy per
    atom.
    '''

    m_def = Category(
        categories=[energy_value],
        a_legacy=LegacyDefinition(name='energy_component_per_atom'))


class energy_component(MCategory):
    '''
    A value of an energy component, expected to be an extensive property.
    '''

    m_def = Category(
        categories=[energy_value],
        a_legacy=LegacyDefinition(name='energy_component'))


class energy_type_reference(MCategory):
    '''
    This metadata stores an energy used as reference point.
    '''

    m_def = Category(
        categories=[energy_value],
        a_legacy=LegacyDefinition(name='energy_type_reference'))


class error_estimate(MCategory):
    '''
    An estimate of the error on the converged (final) value.
    '''

    m_def = Category(
        categories=[error_estimate_contribution],
        a_legacy=LegacyDefinition(name='error_estimate'))


class message_info(MCategory):
    '''
    An information message of the computational program.
    '''

    m_def = Category(
        categories=[message_debug],
        a_legacy=LegacyDefinition(name='message_info'))


class parallelization_info(MCategory):
    '''
    Contains information on the parallelization of the program, i.e. which parallel
    programming language was used and its version, how many cores had been working on the
    calculation and the flags and parameters needed to run the parallelization of the
    code.
    '''

    m_def = Category(
        categories=[accessory_info],
        a_legacy=LegacyDefinition(name='parallelization_info'))


class parsing_message_info(MCategory):
    '''
    This field is used for info messages of the parsing program.
    '''

    m_def = Category(
        categories=[parsing_message_debug],
        a_legacy=LegacyDefinition(name='parsing_message_info'))


class program_info(MCategory):
    '''
    Contains information on the program that generated the data, i.e. the program_name,
    program_version, program_compilation_host and program_compilation_datetime as direct
    children of this field.
    '''

    m_def = Category(
        categories=[accessory_info],
        a_legacy=LegacyDefinition(name='program_info'))


class settings_geometry_optimization(MCategory):
    '''
    Contains parameters controlling the geometry optimization.
    '''

    m_def = Category(
        categories=[settings_sampling],
        a_legacy=LegacyDefinition(name='settings_geometry_optimization'))


class settings_k_points(MCategory):
    '''
    Contains parameters that control the $k$-point mesh.
    '''

    m_def = Category(
        categories=[settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='settings_k_points'))


class settings_metadynamics(MCategory):
    '''
    Contains parameters that control the metadynamics sampling.
    '''

    m_def = Category(
        categories=[settings_sampling],
        a_legacy=LegacyDefinition(name='settings_metadynamics'))


class settings_molecular_dynamics(MCategory):
    '''
    Contains parameters that control the molecular dynamics sampling.
    '''

    m_def = Category(
        categories=[settings_sampling],
        a_legacy=LegacyDefinition(name='settings_molecular_dynamics'))


class settings_Monte_Carlo(MCategory):
    '''
    Contains parameters that control the Monte-Carlo sampling.
    '''

    m_def = Category(
        categories=[settings_sampling],
        a_legacy=LegacyDefinition(name='settings_Monte_Carlo'))


class settings_XC(MCategory):
    '''
    Contains parameters connected with the definition of the exchange-correlation (XC)
    *method*. Here, the term *method* is a more general concept than just *functionals*
    and include, e.g., post Hartree-Fock methods, too.
    '''

    m_def = Category(
        categories=[settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='settings_XC'))


class time_info(MCategory):
    '''
    Stores information on the date and timings of the calculation. They are useful for,
    e.g., debugging or visualization purposes.
    '''

    m_def = Category(
        categories=[accessory_info],
        a_legacy=LegacyDefinition(name='time_info'))


class energy_total_potential_per_atom(MCategory):
    '''
    A value of the total potential energy per atom. Note that a direct comparison may not
    be possible because of a difference in the methods for computing total energies and
    numerical implementations of various codes might leads to different energy zeros (see
    section_energy_code_independent for a code-independent definition of the energy).
    '''

    m_def = Category(
        categories=[energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_total_potential_per_atom'))


class energy_total_potential(MCategory):
    '''
    A value of the total potential energy. Note that a direct comparison may not be
    possible because of a difference in the methods for computing total energies and
    numerical implementations of various codes might leads to different energy zeros (see
    section_energy_code_independent for a code-independent definition of the energy).
    '''

    m_def = Category(
        categories=[energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_total_potential'))


class energy_type_C(MCategory):
    '''
    This metadata stores the correlation (C) energy.
    '''

    m_def = Category(
        categories=[energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_type_C'))


class energy_type_van_der_Waals(MCategory):
    '''
    This metadata stores the converged van der Waals energy.
    '''

    m_def = Category(
        categories=[energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_type_van_der_Waals'))


class energy_type_XC(MCategory):
    '''
    This metadata stores the exchange-correlation (XC) energy.
    '''

    m_def = Category(
        categories=[energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_type_XC'))


class energy_type_X(MCategory):
    '''
    This metadata stores the exchange (X) energy.
    '''

    m_def = Category(
        categories=[energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_type_X'))


class message_warning(MCategory):
    '''
    A warning message of the computational program.
    '''

    m_def = Category(
        categories=[message_info, message_debug],
        a_legacy=LegacyDefinition(name='message_warning'))


class parsing_message_warning(MCategory):
    '''
    This field is used for warning messages of the parsing program.
    '''

    m_def = Category(
        categories=[parsing_message_info, parsing_message_debug],
        a_legacy=LegacyDefinition(name='parsing_message_warning'))


class settings_barostat(MCategory):
    '''
    Contains parameters controlling the barostat in a molecular dynamics calculation.
    '''

    m_def = Category(
        categories=[settings_sampling, settings_molecular_dynamics],
        a_legacy=LegacyDefinition(name='settings_barostat'))


class settings_integrator(MCategory):
    '''
    Contains parameters that control the molecular dynamics (MD) integrator.
    '''

    m_def = Category(
        categories=[settings_sampling, settings_molecular_dynamics],
        a_legacy=LegacyDefinition(name='settings_integrator'))


class settings_post_hartree_fock(MCategory):
    '''
    Contains parameters for the post Hartree-Fock method.
    '''

    m_def = Category(
        categories=[settings_XC, settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='settings_post_hartree_fock'))


class settings_relativity(MCategory):
    '''
    Contains parameters and information connected with the relativistic treatment used in
    the calculation.
    '''

    m_def = Category(
        categories=[settings_XC, settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='settings_relativity'))


class settings_self_interaction_correction(MCategory):
    '''
    Contains parameters and information connected with the self-interaction correction
    (SIC) method being used in self_interaction_correction_method.
    '''

    m_def = Category(
        categories=[settings_XC, settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='settings_self_interaction_correction'))


class settings_thermostat(MCategory):
    '''
    Contains parameters that control the thermostat in the molecular dynamics (MD)
    calculations.
    '''

    m_def = Category(
        categories=[settings_sampling, settings_molecular_dynamics],
        a_legacy=LegacyDefinition(name='settings_thermostat'))


class settings_van_der_Waals(MCategory):
    '''
    Contain parameters and information connected with the Van der Waals treatment used in
    the calculation to compute the Van der Waals energy (energy_van_der_Waals).
    '''

    m_def = Category(
        categories=[settings_XC, settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='settings_van_der_Waals'))


class settings_XC_functional(MCategory):
    '''
    Contain parameters connected with the definition of the exchange-correlation (XC)
    functional (see section_XC_functionals and XC_functional).
    '''

    m_def = Category(
        categories=[settings_XC, settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='settings_XC_functional'))


class message_error(MCategory):
    '''
    An error message of the computational program.
    '''

    m_def = Category(
        categories=[message_info, message_debug, message_warning],
        a_legacy=LegacyDefinition(name='message_error'))


class parsing_message_error(MCategory):
    '''
    This field is used for error messages of the parsing program.
    '''

    m_def = Category(
        categories=[parsing_message_info, parsing_message_warning, parsing_message_debug],
        a_legacy=LegacyDefinition(name='parsing_message_error'))


class settings_coupled_cluster(MCategory):
    '''
    Contains parameters for the coupled-cluster method (CC) in the post Hartree-Fock step.
    '''

    m_def = Category(
        categories=[settings_post_hartree_fock, settings_XC, settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='settings_coupled_cluster'))


class settings_GW(MCategory):
    '''
    Contains parameters for the GW-method in the post Hartree-Fock step, that expands the
    self-energy in terms of the single particle Green's function $G$ and the screened
    Coulomb interaction $W$.
    '''

    m_def = Category(
        categories=[settings_post_hartree_fock, settings_XC, settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='settings_GW'))


class settings_MCSCF(MCategory):
    '''
    Contains parameters for the multi-configurational self-consistent-field (MCSCF)
    method.
    '''

    m_def = Category(
        categories=[settings_post_hartree_fock, settings_XC, settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='settings_MCSCF'))


class settings_moller_plesset_perturbation_theory(MCategory):
    '''
    Contains parameters for Møller–Plesset perturbation theory.
    '''

    m_def = Category(
        categories=[settings_post_hartree_fock, settings_XC, settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='settings_moller_plesset_perturbation_theory'))


class settings_multi_reference(MCategory):
    '''
    Contains parameters for the multi-reference single and double configuration
    interaction method.
    '''

    m_def = Category(
        categories=[settings_post_hartree_fock, settings_XC, settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='settings_multi_reference'))


class archive_context(MSection):
    '''
    Contains information relating to an archive.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='archive_context'))

    archive_gid = Quantity(
        type=str,
        shape=[],
        description='''
        unique identifier of an archive.
        ''',
        a_legacy=LegacyDefinition(name='archive_gid'))


class calculation_context(MSection):
    '''
    Contains information relating to a calculation.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='calculation_context'))

    calculation_gid = Quantity(
        type=str,
        shape=[],
        description='''
        unique identifier of a calculation.
        ''',
        a_legacy=LegacyDefinition(name='calculation_gid'))


class section_atom_projected_dos(MSection):
    '''
    Section collecting the information on an atom projected density of states (DOS)
    evaluation.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_atom_projected_dos'))

    atom_projected_dos_energies = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atom_projected_dos_values'],
        unit='joule',
        description='''
        Array containing the set of discrete energy values for the atom-projected density
        (electronic-energy) of states (DOS).
        ''',
        a_legacy=LegacyDefinition(name='atom_projected_dos_energies'))

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
        ''',
        a_legacy=LegacyDefinition(name='atom_projected_dos_lm'))

    atom_projected_dos_m_kind = Quantity(
        type=str,
        shape=[],
        description='''
        String describing what the integer numbers of $m$ in atom_projected_dos_lm mean.
        The allowed values are listed in the [m_kind wiki
        page](https://gitlab.rzg.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/m-kind).
        ''',
        a_legacy=LegacyDefinition(name='atom_projected_dos_m_kind'))

    atom_projected_dos_values_lm = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_lm_atom_projected_dos', 'number_of_spin_channels', 'number_of_atoms', 'number_of_atom_projected_dos_values'],
        description='''
        Values correspond to the number of states for a given energy (the set of discrete
        energy values is given in atom_projected_dos_energies) divided into contributions
        from each $l,m$ channel for the atom-projected density (electronic-energy) of
        states. Here, there are as many atom-projected DOS as the number_of_atoms, the
        list of labels of the atoms and their meanings are in atom_labels.
        ''',
        a_legacy=LegacyDefinition(name='atom_projected_dos_values_lm'))

    atom_projected_dos_values_total = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_spin_channels', 'number_of_atoms', 'number_of_atom_projected_dos_values'],
        description='''
        Values correspond to the number of states for a given energy (the set of discrete
        energy values is given in atom_projected_dos_energies) divided into contributions
        summed up over all $l$ channels for the atom-projected density (electronic-energy)
        of states (DOS). Here, there are as many atom-projected DOS as the
        number_of_atoms, the list of labels of the atoms and their meanings are in
        atom_labels.
        ''',
        a_legacy=LegacyDefinition(name='atom_projected_dos_values_total'))

    number_of_atom_projected_dos_values = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of energy values for the atom-projected density of states (DOS)
        based on atom_projected_dos_values_lm and atom_projected_dos_values_total.
        ''',
        a_legacy=LegacyDefinition(name='number_of_atom_projected_dos_values'))

    number_of_lm_atom_projected_dos = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of $l$, $m$ combinations for the atom projected density of states
        (DOS) defined in section_atom_projected_dos.
        ''',
        a_legacy=LegacyDefinition(name='number_of_lm_atom_projected_dos'))


class section_atomic_multipoles(MSection):
    '''
    Section describing multipoles (charges/monopoles, dipoles, quadrupoles, ...) for each
    atom.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_atomic_multipoles'))

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
        ''',
        a_legacy=LegacyDefinition(name='atomic_multipole_kind'))

    atomic_multipole_lm = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_lm_atomic_multipoles', 2],
        description='''
        Tuples of $l$ and $m$ values for which the atomic multipoles (including the
        electric charge, dipole, etc.) are given. The method used to obtain the multipoles
        is specified by atomic_multipole_kind. The meaning of the integer number $l$ is
        monopole/charge for $l=0$, dipole for $l=1$, quadrupole for $l=2$, etc. The
        meaning of the integer numbers $m$ is specified by atomic_multipole_m_kind.
        ''',
        a_legacy=LegacyDefinition(name='atomic_multipole_lm'))

    atomic_multipole_m_kind = Quantity(
        type=str,
        shape=[],
        description='''
        String describing the definition for each integer number $m$ in
        atomic_multipole_lm. Allowed values are listed in the [m_kind wiki
        page](https://gitlab.rzg.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/m-kind).
        ''',
        a_legacy=LegacyDefinition(name='atomic_multipole_m_kind'))

    atomic_multipole_values = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_lm_atomic_multipoles', 'number_of_atoms'],
        description='''
        Value of the multipoles (including the monopole/charge for $l$ = 0, the dipole for
        $l$ = 1, etc.) for each atom, calculated as described in atomic_multipole_kind.
        ''',
        a_legacy=LegacyDefinition(name='atomic_multipole_values'))

    number_of_lm_atomic_multipoles = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of $l$, $m$ combinations for atomic multipoles
        atomic_multipole_lm.
        ''',
        a_legacy=LegacyDefinition(name='number_of_lm_atomic_multipoles'))


class section_basis_functions_atom_centered(MSection):
    '''
    This section contains the description of the basis functions (at least one function)
    of the (atom-centered) basis set defined in section_basis_set_atom_centered.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_basis_functions_atom_centered'))


class section_basis_set_atom_centered(MSection):
    '''
    This section describes the atom-centered basis set. The main contained information is
    a short, non unique but human-interpretable, name for identifying the basis set
    (basis_set_atom_centered_short_name), a longer, unique name
    (basis_set_atom_centered_unique_name), the atomic number of the atomic species the
    basis set is meant for (basis_set_atom_number), and a list of actual basis functions
    in the section_basis_functions_atom_centered section.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_basis_set_atom_centered'))

    basis_set_atom_centered_ls = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_kinds_in_basis_set_atom_centered'],
        description='''
        Azimuthal quantum number ($l$) values (of the angular part given by the spherical
        harmonic $Y_{lm}$) of the atom-centered basis function defined in the current
        section_basis_set_atom_centered.
        ''',
        a_legacy=LegacyDefinition(name='basis_set_atom_centered_ls'))

    basis_set_atom_centered_radial_functions = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_kinds_in_basis_set_atom_centered', 401, 5],
        description='''
        Values of the radial function of the different basis function kinds. The values
        are numerically tabulated on a default 0.01-nm equally spaced grid from 0 to 4 nm.
        The 5 tabulated values are $r$, $f(r)$, $f'(r)$, $f(r) \\cdot r$,
        $\\frac{d}{dr}(f(r) \\cdot r)$.
        ''',
        a_legacy=LegacyDefinition(name='basis_set_atom_centered_radial_functions'))

    basis_set_atom_centered_short_name = Quantity(
        type=str,
        shape=[],
        description='''
        Code-specific, but explicative, base name for the basis set (not unique). Details
        are explained in the [basis_set_atom_centered_short_name wiki
        page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/basis-
        set-atom-centered-short-name), this name should not contain the *atom kind* (to
        simplify the use of a single name for multiple elements).
        ''',
        a_legacy=LegacyDefinition(name='basis_set_atom_centered_short_name'))

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
        ''',
        a_legacy=LegacyDefinition(name='basis_set_atom_centered_unique_name'))

    basis_set_atom_number = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Atomic number (i.e., number of protons) of the atom for which this basis set is
        constructed (0 means unspecified or a pseudo atom).
        ''',
        a_legacy=LegacyDefinition(name='basis_set_atom_number'))

    number_of_basis_functions_in_basis_set_atom_centered = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of different basis functions in a section_basis_set_atom_centered
        section. This equals the number of actual coefficients that are specified when
        using this basis set.
        ''',
        a_legacy=LegacyDefinition(name='number_of_basis_functions_in_basis_set_atom_centered'))

    number_of_kinds_in_basis_set_atom_centered = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of different *kinds* of radial basis functions in the
        section_basis_set_atom_centered section. Specifically, basis functions with the
        same $n$ and $l$ quantum numbers are grouped in sets. Each set counts as one
        *kind*.
        ''',
        a_legacy=LegacyDefinition(name='number_of_kinds_in_basis_set_atom_centered'))

    section_basis_functions_atom_centered = SubSection(
        sub_section=SectionProxy('section_basis_functions_atom_centered'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_basis_functions_atom_centered'))

    section_gaussian_basis_group = SubSection(
        sub_section=SectionProxy('section_gaussian_basis_group'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_gaussian_basis_group'))


class section_basis_set_cell_dependent(MSection):
    '''
    Section describing a cell-dependent (atom-independent) basis set, e.g. plane waves.
    The contained information is the type of basis set (in basis_set_cell_dependent_kind),
    its parameters (e.g., for plane waves in basis_set_planewave_cutoff), and a name that
    identifies the actually used basis set (a string combining the type and the
    parameter(s), stored in basis_set_cell_dependent_name).
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_basis_set_cell_dependent'))

    basis_set_cell_dependent_kind = Quantity(
        type=str,
        shape=[],
        description='''
        A string defining the type of the cell-dependent basis set (i.e., non atom
        centered such as plane-waves). Allowed values are listed in the
        [basis_set_cell_dependent_kind wiki page](https://gitlab.mpcdf.mpg.de/nomad-
        lab/nomad-meta-info/wikis/metainfo/basis-set-cell-dependent-kind).
        ''',
        a_legacy=LegacyDefinition(name='basis_set_cell_dependent_kind'))

    basis_set_cell_dependent_name = Quantity(
        type=str,
        shape=[],
        description='''
        A label identifying the cell-dependent basis set (i.e., non atom centered such as
        plane-waves). Allowed values are listed in the [basis_set_cell_dependent_name wiki
        page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/basis-
        set-cell-dependent-name).
        ''',
        a_legacy=LegacyDefinition(name='basis_set_cell_dependent_name'))

    basis_set_planewave_cutoff = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Spherical cutoff  in reciprocal space for a plane-wave basis set. It is the energy
        of the highest plan-ewave ($\\frac{\\hbar^2|k+G|^2}{2m_e}$) included in the basis
        set. Note that normally this basis set is used for the wavefunctions, and the
        density would have 4 times the cutoff, but this actually depends on the use of the
        basis set by the method.
        ''',
        a_legacy=LegacyDefinition(name='basis_set_planewave_cutoff'))


class section_basis_set(MSection):
    '''
    This section contains references to *all* basis sets used in this
    section_single_configuration_calculation. More than one basis set instance per *single
    configuration calculation* (see section_single_configuration_calculation) may be
    needed. This is true for example, for codes that implement adaptive basis sets along
    the self-consistent field (SCF) convergence (e.g., exciting). In such cases, there is
    a section_basis_set instance per SCF iteration, if necessary. Another example is
    having a basis set for wavefunctions, a different one for the density, an auxiliary
    basis set for resolution of identity (RI), etc.

    Supported are the two broad classes of basis sets: *atom-centered* (e.g., Gaussian-
    type, numerical atomic orbitals) and *cell-dependent* (like plane waves or real-space
    grids, so named because they are typically used for periodic-system calculations and
    dependent to the simulated cell as a whole).

    Basis sets used in this section_single_configuration_calculation, belonging to either
    class, are defined in the dedicated section: [section_basis_set_cell_dependent
    ](section_basis_set_cell_dependent) or section_basis_set_atom_centered. The
    correspondence between the basis sets listed in this section and the definition given
    in the dedicated sessions is given by the two concrete metadata:
    mapping_section_basis_set_cell_dependent and mapping_section_basis_set_atom_centered.
    The latter metadata is a list that connects each atom in the system with its basis
    set, where the same basis set can be assigned to more than one atom.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_basis_set'))

    basis_set_kind = Quantity(
        type=str,
        shape=[],
        description='''
        String describing the use of the basis set, i.e, if it used for expanding a wave-
        function or an electron density. Allowed values are listed in the [basis_set_kind
        wiki page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-
        info/wikis/metainfo/basis-set-kind).
        ''',
        a_legacy=LegacyDefinition(name='basis_set_kind'))

    basis_set_name = Quantity(
        type=str,
        shape=[],
        description='''
        String identifying the basis set in an unique way. The rules for building this
        string are specified in the [basis_set_name wiki
        page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/basis-
        set-name).
        ''',
        a_legacy=LegacyDefinition(name='basis_set_name'))

    mapping_section_basis_set_atom_centered = Quantity(
        type=Reference(SectionProxy('section_basis_set_atom_centered')),
        shape=['number_of_atoms'],
        description='''
        An array of the dimension of number_of_atoms where each atom (identified by the
        index in the array) is assigned to an atom-centered basis set, for this
        section_single_configuration_calculation. The actual definition of the atom-
        centered basis set is in the section_basis_set_atom_centered that is referred to
        by this metadata.
        ''',
        a_legacy=LegacyDefinition(name='mapping_section_basis_set_atom_centered'))

    mapping_section_basis_set_cell_dependent = Quantity(
        type=Reference(SectionProxy('section_basis_set_cell_dependent')),
        shape=[],
        description='''
        Assignment of the cell-dependent (i.e., non atom centered, e.g., plane-waves)
        parts of the basis set, which is defined (type, parameters) in
        section_basis_set_cell_dependent that is referred to by this metadata.
        ''',
        a_legacy=LegacyDefinition(name='mapping_section_basis_set_cell_dependent'))

    number_of_basis_functions = Quantity(
        type=int,
        shape=[],
        description='''
        Stores the total number of basis functions in a section_basis_set section.
        ''',
        a_legacy=LegacyDefinition(name='number_of_basis_functions'))


class section_restricted_uri(MSection):
    '''
    Restricted URIs on this calculation (Coverage: any info or files that are related with
    this calculation can be subject to restriction)
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_restricted_uri'))

    number_of_restricted_uri = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        The number of restricted uris in restricted_uri list.
        ''',
        a_legacy=LegacyDefinition(name='number_of_restricted_uri'))

    restricted_uri = Quantity(
        type=str,
        shape=['number_of_restricted_uri'],
        description='''
        The list of nomad uri(s) identifying the restricted info/file corresponding to
        this calculation
        ''',
        a_legacy=LegacyDefinition(name='restricted_uri'))

    restricted_uri_reason = Quantity(
        type=str,
        shape=[],
        description='''
        The reason of restriction for the uri or file. The reason can be 'propriety
        license', 'open-source redistribution restricted license', 'other license', or
        'author restricted'.
        ''',
        a_legacy=LegacyDefinition(name='restricted_uri_reason'))

    restricted_uri_issue_authority = Quantity(
        type=str,
        shape=[],
        description='''
        The issue authority is the restriction owner for the uri or file. This can be
        license owner such as 'VASP' or 'AMBER', 'NOMAD', or the author of the uri. For
        example the repository user name of the author.
        ''',
        a_legacy=LegacyDefinition(name='restricted_uri_issue_authority'))

    restricted_uri_end_date = Quantity(
        type=str,
        shape=[],
        description='''
        The deadline date of the restriction for the uri or file. The end date can be in
        date format string for those restrictions set by authors or NOMAD otherwise it is
        set to 'unlimited' for the restriction related to license.
        ''',
        a_legacy=LegacyDefinition(name='restricted_uri_end_date'))

    restricted_uri_restriction = Quantity(
        type=str,
        shape=[],
        description='''
        The type of restriction for the uri or file. The type can be 'any access' or
        'license permitted'.
        ''',
        a_legacy=LegacyDefinition(name='restricted_uri_restriction'))

    restricted_uri_license = Quantity(
        type=str,
        shape=[],
        description='''
        The info of the license that is the reason of restriction.
        ''',
        a_legacy=LegacyDefinition(name='restricted_uri_license'))

    number_of_restricted_uri_files = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        The number of restricted files in restricted_uri_files list.
        ''',
        a_legacy=LegacyDefinition(name='number_of_restricted_uri_files'))


class section_calculation_to_calculation_refs(MSection):
    '''
    Section that describes the relationship between different
    section_single_configuration_calculation sections.

    For instance, one calculation is a perturbation performed using a self-consistent
    field (SCF) calculation as starting point, or a simulated system is partitioned in
    regions with different but connected Hamiltonians (e.g., QM/MM, or a region treated
    via Kohn-Sham DFT embedded into a region treated via orbital-free DFT).

    The kind of relationship between the calculation defined in this section and the
    referenced one is described by calculation_to_calculation_kind. The referenced
    section_single_configuration_calculation is identified via
    calculation_to_calculation_ref (typically used for a
    section_single_configuration_calculation in the same section_run) or
    calculation_to_calculation_external_url.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_calculation_to_calculation_refs'))

    calculation_to_calculation_external_url = Quantity(
        type=str,
        shape=[],
        description='''
        URL used to reference an externally stored calculation. The kind of relationship
        between the present and the referenced section_single_configuration_calculation is
        specified by calculation_to_calculation_kind.
        ''',
        a_legacy=LegacyDefinition(name='calculation_to_calculation_external_url'))

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
        ''',
        a_legacy=LegacyDefinition(name='calculation_to_calculation_kind'))

    calculation_to_calculation_ref = Quantity(
        type=Reference(SectionProxy('section_single_configuration_calculation')),
        shape=[],
        description='''
        Reference to another calculation. If both this and
        calculation_to_calculation_external_url are given, then
        calculation_to_calculation_ref is a local copy of the URL given in
        calculation_to_calculation_external_url. The kind of relationship between the
        present and the referenced section_single_configuration_calculation is specified
        by calculation_to_calculation_kind.
        ''',
        a_legacy=LegacyDefinition(name='calculation_to_calculation_ref'))


class section_calculation_to_folder_refs(MSection):
    '''
    Section that describes the relationship between
    section_single_configuration_calculationa and the folder containing the original
    calulations
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_calculation_to_folder_refs'))

    calculation_to_folder_external_url = Quantity(
        type=str,
        shape=[],
        description='''
        URL used to reference a folder containing external calculations. The kind of
        relationship between the present and the referenced
        section_single_configuration_calculation is specified by
        calculation_to_folder_kind.
        ''',
        a_legacy=LegacyDefinition(name='calculation_to_folder_external_url'))

    calculation_to_folder_kind = Quantity(
        type=str,
        shape=[],
        description='''
        String defining the relationship between the referenced
        section_single_configuration_calculation and a folder containing calculations.
        ''',
        a_legacy=LegacyDefinition(name='calculation_to_folder_kind'))


class section_dos_fingerprint(MSection):
    '''
    Section for the fingerprint of the electronic density-of-states (DOS).
    DOS fingerprints are a modification of the D-Fingerprints reported in Chem. Mater. 2015, 27, 3, 735–743
    (doi:10.1021/cm503507h). The fingerprint consists of a binary representation of the DOS,
    that is used to evaluate the similarity of materials based on their electronic structure.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_dos_fingerprint'))

    bins = Quantity(
        type=str,
        description='''
        Byte representation of the DOS fingerprint.
        ''',
        a_legacy=LegacyDefinition(name='bins'))

    indices = Quantity(
        type=np.dtype(np.int16),
        shape=[2],
        description='''
        Indices used to compare DOS fingerprints of different energy ranges.
        ''',
        a_legacy=LegacyDefinition(name='indices'))

    stepsize = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Stepsize of interpolation in the first step of the generation of DOS fingerprints.
        ''',
        a_legacy=LegacyDefinition(name='stepsize'))

    filling_factor = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Proportion of 1 bins in the DOS fingerprint.
        ''',
        a_legacy=LegacyDefinition(name='filling_factor'))

    grid_id = Quantity(
        type=str,
        description='''
        Identifier of the DOS grid that was used for the creation of the fingerprint.
        Similarity can only be calculated if the same grid was used for both fingerprints.
        ''',
        a_legacy=LegacyDefinition(name='grid_id'))


class section_dos(MSection):
    '''
    Section collecting information of a (electronic-energy or vibrational-energy) density
    of states (DOS) evaluation.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_dos'))

    dos_energies_normalized = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_dos_values'],
        unit='joule',
        description='''
        Array containing the set of discrete energy values with respect to the
        highest occupied energy level. This is the total DOS, see
        atom_projected_dos_energies and species_projected_dos_energies for
        partial density of states.

        If not available through energy_reference_highest_occupied, the highest
        occupied energy level is detected by searching for a non-zero DOS value
        below (or nearby) the reported energy_reference_fermi. In case the
        highest occupied energy level cannot be detected accurately, the
        normalized values are not reported. For calculations with multiple
        spin-channels, the normalization is determined by the first channel.
        ''',
        a_legacy=LegacyDefinition(name='dos_energies_normalized'))

    dos_energies = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_dos_values'],
        unit='joule',
        description='''
        Array containing the set of discrete energy values for the density (electronic-
        energy or vibrational energy) of states (DOS). This is the total DOS, see
        atom_projected_dos_energies and species_projected_dos_energies for partial density
        of states.
        ''',
        a_legacy=LegacyDefinition(name='dos_energies'))

    dos_integrated_values = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_spin_channels', 'number_of_dos_values'],
        description='''
        Integrated density of states (starting at $-\\infty$), pseudo potential
        calculations should start with the number of core electrons if they cover only the
        active electrons
        ''',
        a_legacy=LegacyDefinition(name='dos_integrated_values'))

    dos_kind = Quantity(
        type=str,
        shape=[],
        description='''
        String to specify the kind of density of states (either electronic or
        vibrational).
        ''',
        a_legacy=LegacyDefinition(name='dos_kind'))

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
        ''',
        a_legacy=LegacyDefinition(name='dos_lm'))

    dos_m_kind = Quantity(
        type=str,
        shape=[],
        description='''
        String describing what the integer numbers of $m$ in dos_lm mean. The allowed
        values are listed in the [m_kind wiki page](https://gitlab.rzg.mpg.de/nomad-
        lab/nomad-meta-info/wikis/metainfo/m-kind).
        ''',
        a_legacy=LegacyDefinition(name='dos_m_kind'))

    dos_values_lm = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_dos_lms', 'number_of_spin_channels', 'number_of_atoms', 'number_of_dos_values'],
        unit='joule',
        description='''
        Array containing the density (electronic-energy) of states values projected on the
        various spherical harmonics (integrated on all atoms), see
        atom_projected_dos_values_lm for atom values.
        ''',
        a_legacy=LegacyDefinition(name='dos_values_lm'))

    dos_values_normalized = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_spin_channels', 'number_of_dos_values'],
        description='''
        Density of states (DOS) values divided by the unit cell volume and by the number
        of atoms.
        ''',
        a_legacy=LegacyDefinition(name='dos_values_normalized'))

    dos_values = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_spin_channels', 'number_of_dos_values'],
        description='''
        Values (number of states for a given energy, the set of discrete energy values is
        given in dos_energies) of density (electronic-energy or vibrational-energy) of
        states. This refers to the simulation cell, i.e. integrating over all energies
        will give the number of electrons in the simulation cell.
        ''',
        a_legacy=LegacyDefinition(name='dos_values'))

    number_of_dos_lms = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of $l$, $m$ combinations for the given projected density of
        states (DOS) in dos_values and dos_values_lm.
        ''',
        a_legacy=LegacyDefinition(name='number_of_dos_lms'))

    number_of_dos_values = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of energy values for the density of states (DOS), see
        dos_energies.
        ''',
        a_legacy=LegacyDefinition(name='number_of_dos_values'))

    section_dos_fingerprint = SubSection(
        sub_section=SectionProxy('section_dos_fingerprint'),
        repeats=False,
        a_legacy=LegacyDefinition(name='section_dos_fingerprint'))


class section_eigenvalues(MSection):
    '''
    Section containing (electronic-energy) eigenvalues for one spin channel. If, for
    example, the eigenvalues of the Kohn-Sham operator are to be stored, a string
    identifying this kind of eigenvalues is put in eigenvalues_kind, the coordinates of
    the $k$-points at which the eigenvalues are evaluated is stored in
    eigenvalues_kpoints, and the energy values of the eigenstates and their occupation is
    stored in eigenvalues_values and eigenvalues_occupation, respectively.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_eigenvalues'))

    eigenvalues_kind = Quantity(
        type=str,
        shape=[],
        description='''
        A short string describing the kind of eigenvalues, as defined in the
        [eigenvalues_kind wiki page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-
        info/wikis/metainfo/eigenvalues-kind).
        ''',
        a_legacy=LegacyDefinition(name='eigenvalues_kind'))

    eigenvalues_kpoints_multiplicity = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_eigenvalues_kpoints'],
        description='''
        Multiplicity of the $k$ point (i.e., how many distinct points per cell this
        expands to after applying all symmetries). This defaults to 1. If expansion is
        preformed then each point will have weight
        eigenvalues_kpoints_weights/eigenvalues_kpoints_multiplicity.
        ''',
        a_legacy=LegacyDefinition(name='eigenvalues_kpoints_multiplicity'))

    eigenvalues_kpoints_weights = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_eigenvalues_kpoints'],
        description='''
        Weights of the $k$ points (in the basis of the reciprocal lattice vectors) used
        for the evaluation of the eigenvalues tabulated in eigenvalues_values, should
        account for symmetry too.
        ''',
        a_legacy=LegacyDefinition(name='eigenvalues_kpoints_weights'))

    eigenvalues_kpoints = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_eigenvalues_kpoints', 3],
        description='''
        Coordinates of the $k$ points (in the basis of the reciprocal lattice vectors)
        used for the evaluation of the eigenvalues tabulated in eigenvalues_values.
        ''',
        a_legacy=LegacyDefinition(name='eigenvalues_kpoints'))

    eigenvalues_occupation = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_spin_channels', 'number_of_eigenvalues_kpoints', 'number_of_eigenvalues'],
        description='''
        Occupation of the eigenstates. The corresponding eigenvalues (energy) are given in
        eigenvalues_values. The coordinates in the reciprocal space are defined in
        eigenvalues_kpoints.
        ''',
        a_legacy=LegacyDefinition(name='eigenvalues_occupation'))

    eigenvalues_values = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_spin_channels', 'number_of_eigenvalues_kpoints', 'number_of_eigenvalues'],
        unit='joule',
        description='''
        Values of the (electronic-energy) eigenvalues. The coordinates of the
        corresponding eigenstates in the reciprocal space are defined in
        eigenvalues_kpoints and their occupations are given in eigenvalues_occupation.
        ''',
        a_legacy=LegacyDefinition(name='eigenvalues_values'))

    number_of_band_segment_eigenvalues = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of eigenvalues in a band segment, see band_energies.
        ''',
        a_legacy=LegacyDefinition(name='number_of_band_segment_eigenvalues'))

    number_of_eigenvalues_kpoints = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of $k$ points, see eigenvalues_kpoints. $k$ points are calculated
        within a run and are irreducible if a symmetry is used.
        ''',
        a_legacy=LegacyDefinition(name='number_of_eigenvalues_kpoints'))

    number_of_eigenvalues = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of eigenvalues, see eigenvalues_values.
        ''',
        a_legacy=LegacyDefinition(name='number_of_eigenvalues'))

    number_of_normalized_band_segment_eigenvalues = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of normalized eigenvalues in a band segment, see

        band_energies_normalized.
        ''',
        a_legacy=LegacyDefinition(name='number_of_normalized_band_segment_eigenvalues'))


class section_energy_code_independent(MSection):
    '''
    Section describing a code-independent total energy obtained by subtracting some
    reference energy calculated with the same code. It contains the type in
    energy_code_independent_kind and the computed code-independent total energy in
    energy_code_independent_value. The computed energy allows for comparisons among
    different codes and numerical settings.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_energy_code_independent'))

    energy_code_independent_kind = Quantity(
        type=str,
        shape=[],
        description='''
        Type of the code-independent total energy (obtained by subtracting a reference
        energy calculated with the same code), created to be comparable among different
        codes and numerical settings. Details can be found on the [energy_code_independent
        wiki page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-
        info/wikis/metainfo/energy-code-independent).
        ''',
        a_legacy=LegacyDefinition(name='energy_code_independent_kind'))

    energy_code_independent_value = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Value of the code-independent total energy (obtained by subtracting a reference
        energy calculated with the same code). This value is created to be comparable
        among different codes and numerical settings. Details can be found on the
        [energy_code_independent wiki page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-
        meta-info/wikis/metainfo/energy-code-independent).
        ''',
        categories=[energy_component, energy_value, energy_total_potential],
        a_legacy=LegacyDefinition(name='energy_code_independent_value'))


class section_energy_van_der_Waals(MSection):
    '''
    Section containing the Van der Waals energy value (energy_van_der_Waals_value) of type
    van_der_Waals_kind. This is used when more than one Van der Waals methods are applied
    in the same *single configuration calculation*, see
    section_single_configuration_calculation. The main Van der Waals method (the one
    concurring to energy_current, and used, e.g., for evaluating the forces for a
    relaxation or dynamics) is given in energy_van_der_Waals and is defined in
    settings_van_der_Waals.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_energy_van_der_Waals'))

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
        ''',
        a_legacy=LegacyDefinition(name='energy_van_der_Waals_kind'))

    energy_van_der_Waals_value = Quantity(
        type=np.dtype(np.float64),
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
        ''',
        categories=[energy_type_van_der_Waals, energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_van_der_Waals_value'))

    energy_van_der_Waals = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Value for the converged van der Waals energy calculated using the method described
        in van_der_Waals_method, and used in energy_current. This is the van der Waals
        method consistent with, e.g., forces used for relaxation or dynamics. Alternative
        methods are listed in section_energy_van_der_Waals.
        ''',
        categories=[energy_type_van_der_Waals, energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_van_der_Waals'))


class section_frame_sequence_user_quantity(MSection):
    '''
    Section collecting some user-defined quantities evaluated along a sequence of frame.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_frame_sequence_user_quantity'))

    frame_sequence_user_quantity_frames = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_user_quantity_evaluations_in_sequence'],
        description='''
        Array containing the strictly increasing indices referring to the frames of
        frame_sequence_user_quantity. If not given it defaults to the trivial mapping
        0,1,...
        ''',
        a_legacy=LegacyDefinition(name='frame_sequence_user_quantity_frames'))

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
        ''',
        a_legacy=LegacyDefinition(name='frame_sequence_user_quantity_name'))

    frame_sequence_user_quantity_stats = Quantity(
        type=np.dtype(np.float64),
        shape=[2, 'number_of_frame_sequence_user_quantity_components'],
        description='''
        Average of frame_sequence_user_quantity and its standard deviation in this
        sequence of frames (i.e., a trajectory, a frame is one
        section_single_configuration_calculation).
        ''',
        a_legacy=LegacyDefinition(name='frame_sequence_user_quantity_stats'))

    frame_sequence_user_quantity = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_user_quantity_evaluations_in_sequence', 'number_of_frame_sequence_user_quantity_components'],
        description='''
        Array containing the values of the user-defined quantity defined in
        frame_sequence_user_quantity_name, evaluated along this sequence of frames (i.e.,
        trajectory, a frame is one section_single_configuration_calculation). If not all
        frames have a value the indices of the frames that have a value are stored in
        frame_sequence_kinetic_energy_frames. If not all frames have a value the indices
        of the frames that have a value are stored in
        frame_sequence_kinetic_energy_frames.
        ''',
        a_legacy=LegacyDefinition(name='frame_sequence_user_quantity'))

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
        ''',
        a_legacy=LegacyDefinition(name='number_of_frame_sequence_user_quantity_components'))

    number_of_user_quantity_evaluations_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of user defined quantity evaluations along a sequence of
        frame_sequence_user_quantity frames. A sequence is a trajectory, which can have
        number_of_frames_in_sequence each representing one
        section_single_configuration_calculation section.
        ''',
        a_legacy=LegacyDefinition(name='number_of_user_quantity_evaluations_in_sequence'))


class section_frame_sequence(MSection):
    '''
    Section containing a sequence of frames, i.e. a trajectory which can have
    number_of_frames_in_sequence each representing one
    section_single_configuration_calculation section evaluated with a sampling method
    (e.g, molecular dynamics, Monte Carlo, geometry optimization). The sampling method
    might be a subset of the whole trajectory.

    Information on the method used for the sampling can be found in the
    section_sampling_method section and information of each frame of the sequence are
    found in the section_single_configuration_calculation section.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_frame_sequence'))

    frame_sequence_conserved_quantity_frames = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_conserved_quantity_evaluations_in_sequence'],
        description='''
        Array containing the strictly increasing indices of the frames the
        frame_sequence_conserved_quantity values refers to. If not given it defaults to
        the trivial mapping 0,1,...
        ''',
        a_legacy=LegacyDefinition(name='frame_sequence_conserved_quantity_frames'))

    frame_sequence_conserved_quantity_stats = Quantity(
        type=np.dtype(np.float64),
        shape=[2],
        unit='joule',
        description='''
        Average value of energy-like frame_sequence_conserved_quantity, and its standard
        deviation, over this sequence of frames (i.e., a trajectory, a frame is one
        section_single_configuration_calculation).
        ''',
        a_legacy=LegacyDefinition(name='frame_sequence_conserved_quantity_stats'))

    frame_sequence_conserved_quantity = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_conserved_quantity_evaluations_in_sequence'],
        unit='joule',
        description='''
        Array containing the values of a quantity that should be conserved,  along a
        sequence of frames (i.e., a trajectory). A frame is one
        section_single_configuration_calculation), for example the total energy in the NVE
        ensemble. If not all frames have a value the indices of the frames that have a
        value are stored in frame_sequence_conserved_quantity_frames.
        ''',
        a_legacy=LegacyDefinition(name='frame_sequence_conserved_quantity'))

    frame_sequence_continuation_kind = Quantity(
        type=Reference(SectionProxy('section_frame_sequence')),
        shape=[],
        description='''
        Type of continuation that has been performed from the previous sequence of frames
        (i.e., a trajectory, a frame is one section_single_configuration_calculation),
        upon restart.
        ''',
        a_legacy=LegacyDefinition(name='frame_sequence_continuation_kind'))

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
        ''',
        a_legacy=LegacyDefinition(name='frame_sequence_external_url'))

    frame_sequence_kinetic_energy_frames = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_kinetic_energies_in_sequence'],
        description='''
        Array containing the strictly increasing indices referring to the frames of
        frame_sequence_kinetic_energy. If not given it defaults to the trivial mapping
        0,1,...
        ''',
        a_legacy=LegacyDefinition(name='frame_sequence_kinetic_energy_frames'))

    frame_sequence_kinetic_energy_stats = Quantity(
        type=np.dtype(np.float64),
        shape=[2],
        unit='joule',
        description='''
        Average kinetic energy and its standard deviation over this sequence of frames
        (i.e., a trajectory, a frame is one section_single_configuration_calculation).
        ''',
        a_legacy=LegacyDefinition(name='frame_sequence_kinetic_energy_stats'))

    frame_sequence_kinetic_energy = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_kinetic_energies_in_sequence'],
        unit='joule',
        description='''
        Array containing the values of the kinetic energy along this sequence of frames
        (i.e., a trajectory, a frame is one section_single_configuration_calculation). If
        not all frames have a value the indices of the frames that have a value are stored
        in frame_sequence_kinetic_energy_frames.
        ''',
        a_legacy=LegacyDefinition(name='frame_sequence_kinetic_energy'))

    frame_sequence_local_frames_ref = Quantity(
        type=Reference(SectionProxy('section_single_configuration_calculation')),
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
        ''',
        a_legacy=LegacyDefinition(name='frame_sequence_local_frames_ref'))

    frame_sequence_potential_energy_frames = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_potential_energies_in_sequence'],
        description='''
        Array containing the strictly increasing indices referring to the frames of
        frame_sequence_potential_energy. If not given it defaults to the trivial mapping
        0,1,...
        ''',
        a_legacy=LegacyDefinition(name='frame_sequence_potential_energy_frames'))

    frame_sequence_potential_energy_stats = Quantity(
        type=np.dtype(np.float64),
        shape=[2],
        unit='joule',
        description='''
        Average potential energy and its standard deviation over this sequence of frames
        (i.e., a trajectory, a frame is one section_single_configuration_calculation).
        ''',
        a_legacy=LegacyDefinition(name='frame_sequence_potential_energy_stats'))

    frame_sequence_potential_energy = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_potential_energies_in_sequence'],
        unit='joule',
        description='''
        Array containing the value of the potential energy along this sequence of frames
        (i.e., a trajectory, a frame is one section_single_configuration_calculation).
        This is equal to energy_total of the corresponding
        section_single_configuration_calculation and repeated here in a summary array for
        easier access. If not all frames have a value the indices of the frames that have
        a value are stored in frame_sequence_potential_energy_frames.
        ''',
        a_legacy=LegacyDefinition(name='frame_sequence_potential_energy'))

    frame_sequence_pressure_frames = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_pressure_evaluations_in_sequence'],
        description='''
        Array containing the strictly increasing indices referring to the frames of
        frame_sequence_pressure. If not given it defaults to the trivial mapping 0,1,...
        ''',
        a_legacy=LegacyDefinition(name='frame_sequence_pressure_frames'))

    frame_sequence_pressure_stats = Quantity(
        type=np.dtype(np.float64),
        shape=[2],
        unit='pascal',
        description='''
        Average pressure (one third of the trace of the stress tensor) and standard
        deviation over this sequence of frames (i.e., a trajectory, a frame is one
        section_single_configuration_calculation).
        ''',
        a_legacy=LegacyDefinition(name='frame_sequence_pressure_stats'))

    frame_sequence_pressure = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_pressure_evaluations_in_sequence'],
        unit='pascal',
        description='''
        Array containing the values of the pressure (one third of the trace of the stress
        tensor) along this sequence of frames (a frame is one
        section_single_configuration_calculation). If not all frames have a value the
        indices of the frames that have a value are stored in
        frame_sequence_pressure_frames.
        ''',
        a_legacy=LegacyDefinition(name='frame_sequence_pressure'))

    frame_sequence_temperature_frames = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_temperatures_in_sequence'],
        description='''
        Array containing the strictly increasing indices referring to the frames of
        frame_sequence_temperature. If not given it defaults to the trivial mapping
        0,1,...
        ''',
        a_legacy=LegacyDefinition(name='frame_sequence_temperature_frames'))

    frame_sequence_temperature_stats = Quantity(
        type=np.dtype(np.float64),
        shape=[2],
        unit='kelvin',
        description='''
        Average temperature and its standard deviation over this sequence of frames (i.e.,
        a trajectory, a frame is one section_single_configuration_calculation).
        ''',
        a_legacy=LegacyDefinition(name='frame_sequence_temperature_stats'))

    frame_sequence_temperature = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_temperatures_in_sequence'],
        unit='kelvin',
        description='''
        Array containing the values of the instantaneous temperature (a quantity,
        proportional to frame_sequence_kinetic_energy, whose ensemble average equals the
        thermodynamic temperature) along this sequence of frames (i.e., a trajectory, a
        frame is one section_single_configuration_calculation). If not all frames have a
        value the indices of the frames that have a value are stored in
        frame_sequence_temperature_frames.
        ''',
        a_legacy=LegacyDefinition(name='frame_sequence_temperature'))

    frame_sequence_time = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_frames_in_sequence'],
        unit='second',
        description='''
        Time along this sequence of frames (i.e., a trajectory, a frame is one
        section_single_configuration_calculation). Time start is arbitrary, but when a
        sequence is a continuation of another time should be continued too.
        ''',
        a_legacy=LegacyDefinition(name='frame_sequence_time'))

    frame_sequence_to_sampling_ref = Quantity(
        type=Reference(SectionProxy('section_sampling_method')),
        shape=[],
        description='''
        Reference from the present section_frame_sequence to the section_sampling_method,
        that defines the parameters used in this sequence of frames (i.e., a trajectory, a
        frame is one section_single_configuration_calculation).
        ''',
        a_legacy=LegacyDefinition(name='frame_sequence_to_sampling_ref'))

    geometry_optimization_converged = Quantity(
        type=bool,
        shape=[],
        description='''
        Arrays specify whether a geometry optimization is converged.
        ''',
        a_legacy=LegacyDefinition(name='geometry_optimization_converged'))

    number_of_conserved_quantity_evaluations_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of conserved quantity evaluations in this sequence. A sequence is
        a trajectory, which can have number_of_frames_in_sequence each representing one
        section_single_configuration_calculation section.
        ''',
        a_legacy=LegacyDefinition(name='number_of_conserved_quantity_evaluations_in_sequence'))

    number_of_frames_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of frames in a sequence. A sequence is a trajectory, which can
        have number_of_frames_in_sequence each representing one
        section_single_configuration_calculation section.
        ''',
        a_legacy=LegacyDefinition(name='number_of_frames_in_sequence'))

    number_of_kinetic_energies_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of kinetic energy evaluations in this sequence of frames, see
        frame_sequence_kinetic_energy.
        ''',
        a_legacy=LegacyDefinition(name='number_of_kinetic_energies_in_sequence'))

    number_of_potential_energies_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of potential energies evaluation in this sequence. A sequence is
        a trajectory, which can have number_of_frames_in_sequence each representing one
        section_single_configuration_calculation section.
        ''',
        a_legacy=LegacyDefinition(name='number_of_potential_energies_in_sequence'))

    number_of_pressure_evaluations_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of pressure evaluations in this sequence. A sequence is a
        trajectory, which can have number_of_frames_in_sequence each representing one
        section_single_configuration_calculation section.
        ''',
        a_legacy=LegacyDefinition(name='number_of_pressure_evaluations_in_sequence'))

    number_of_temperatures_in_sequence = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of temperature frames (frame_sequence_temperature) used in the
        section_frame_sequence. A sequence is a trajectory, which can have
        number_of_frames_in_sequence each representing one
        section_single_configuration_calculation section.
        ''',
        a_legacy=LegacyDefinition(name='number_of_temperatures_in_sequence'))

    previous_sequence_ref = Quantity(
        type=Reference(SectionProxy('section_frame_sequence')),
        shape=[],
        description='''
        Contains a reference to the previous sequence. A sequence is a trajectory, which
        can have number_of_frames_in_sequence each representing one
        section_single_configuration_calculation section. If not given, a start from an
        initial configuration is assumed.
        ''',
        a_legacy=LegacyDefinition(name='previous_sequence_ref'))

    section_frame_sequence_user_quantity = SubSection(
        sub_section=SectionProxy('section_frame_sequence_user_quantity'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_frame_sequence_user_quantity'))

    section_thermodynamical_properties = SubSection(
        sub_section=SectionProxy('section_thermodynamical_properties'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_thermodynamical_properties'))


class section_gaussian_basis_group(MSection):
    '''
    Section that describes a group of Gaussian contractions. Groups allow one to calculate
    the primitive Gaussian integrals once for several different linear combinations of
    them. This defines basis functions with radial part $f_i(r) = r^{l_i} \\sum_{j} c_{i
    j} A(l_i, \\alpha_j) exp(-\\alpha_j r^2)$ where $A(l_i, \\alpha_j)$ is a the
    normalization coefficient for primitive Gaussian basis functions. Here, $\\alpha_j$ is
    defined in gaussian_basis_group_exponents, $l_i$ is given in gaussian_basis_group_ls,
    and $c_{i j}$ is given in gaussian_basis_group_contractions, whereas the radial part
    is given by the spherical harmonics $Y_{l m}$.

    This section is defined only if the original basis function uses Gaussian basis
    functions, and the sequence of radial functions $f_i$ across all
    section_gaussian_basis_group in section_basis_set_atom_centered should match the one
    of basis_set_atom_centered_radial_functions.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_gaussian_basis_group'))

    gaussian_basis_group_contractions = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_gaussian_basis_group_contractions', 'number_of_gaussian_basis_group_exponents'],
        description='''
        contraction coefficients $c_{i j}$ defining the contracted basis functions with
        respect to *normalized* primitive Gaussian functions. They define the Gaussian
        basis functions as described in section_gaussian_basis_group.
        ''',
        a_legacy=LegacyDefinition(name='gaussian_basis_group_contractions'))

    gaussian_basis_group_exponents = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_gaussian_basis_group_exponents'],
        unit='1 / meter ** 2',
        description='''
        Exponents $\\alpha_j$ of the Gaussian functions defining this basis set
        $exp(-\\alpha_j r^2)$. One should be careful about the units of the coefficients.
        ''',
        a_legacy=LegacyDefinition(name='gaussian_basis_group_exponents'))

    gaussian_basis_group_ls = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_gaussian_basis_group_contractions'],
        description='''
        Azimuthal quantum number ($l$) values (of the angular part given by the spherical
        harmonic $Y_{l m}$ of the various contracted basis functions).
        ''',
        a_legacy=LegacyDefinition(name='gaussian_basis_group_ls'))

    number_of_gaussian_basis_group_contractions = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of different contractions, i.e. resulting basis functions in a
        section_gaussian_basis_group section.
        ''',
        a_legacy=LegacyDefinition(name='number_of_gaussian_basis_group_contractions'))

    number_of_gaussian_basis_group_exponents = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of different Gaussian exponents in a section_gaussian_basis_group
        section.
        ''',
        a_legacy=LegacyDefinition(name='number_of_gaussian_basis_group_exponents'))


class section_k_band_normalized(MSection):
    '''
    This section stores information on a normalized $k$-band (electronic band structure)
    evaluation along one-dimensional pathways in the $k$ (reciprocal) space given in
    section_k_band_segment. Eigenvalues calculated at the actual $k$-mesh used for
    energy_total evaluations, can be found in the section_eigenvalues section.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_k_band_normalized'))

    k_band_path_normalized_is_standard = Quantity(
        type=bool,
        shape=[],
        description='''
        If the normalized path is along the default path defined in W. Setyawan and S.
        Curtarolo, [Comput. Mater. Sci. **49**, 299-312
        (2010)](http://dx.doi.org/10.1016/j.commatsci.2010.05.010).
        ''',
        a_legacy=LegacyDefinition(name='k_band_path_normalized_is_standard'))

    section_k_band_segment_normalized = SubSection(
        sub_section=SectionProxy('section_k_band_segment_normalized'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_k_band_segment_normalized'))


class section_k_band_segment_normalized(MSection):
    '''
    Section collecting the information on a normalized $k$-band segment. This section
    stores band structures along a one-dimensional pathway in the $k$ (reciprocal) space.

    Eigenvalues calculated at the actual $k$-mesh used for energy_total evaluations are
    defined in section_eigenvalues and the band structures are represented as third-order
    tensors: one dimension for the spin channels, one for the sequence of $k$ points for
    the segment (given in number_of_k_points_per_segment), and one for the sequence of
    eigenvalues at a given $k$ point. The values of the $k$ points in each segment are
    stored in band_k_points. The energies and occupation for each eigenstate, at each $k$
    point, segment, and spin channel are stored in band_energies and band_occupations,
    respectively. The labels for the segment are specified in band_segm_labels.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_k_band_segment_normalized'))

    band_energies_normalized = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_spin_channels', 'number_of_normalized_k_points_per_segment', 'number_of_normalized_band_segment_eigenvalues'],
        unit='joule',
        description='''
        $k$-dependent energies of the electronic band segment (electronic band structure)
        with respect to the top of the valence band. This is a third-order tensor, with
        one dimension used for the spin channels, one for the $k$ points for each segment,
        and one for the eigenvalue sequence.
        ''',
        a_legacy=LegacyDefinition(name='band_energies_normalized'))

    band_k_points_normalized = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_normalized_k_points_per_segment', 3],
        description='''
        Fractional coordinates of the $k$ points (in the basis of the reciprocal-lattice
        vectors) for which the normalized electronic energies are given.
        ''',
        a_legacy=LegacyDefinition(name='band_k_points_normalized'))

    band_occupations_normalized = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_spin_channels', 'number_of_normalized_k_points_per_segment', 'number_of_normalized_band_segment_eigenvalues'],
        description='''
        Occupation of the $k$-points along the normalized electronic band. The size of the
        dimensions of this third-order tensor are the same as for the tensor in
        band_energies.
        ''',
        a_legacy=LegacyDefinition(name='band_occupations_normalized'))

    band_segm_labels_normalized = Quantity(
        type=str,
        shape=[2],
        description='''
        Start and end labels of the points in the segment (one-dimensional pathways)
        sampled in the $k$-space, using the conventional symbols, e.g., Gamma, K, L. The
        coordinates (fractional, in the reciprocal space) of the start and end points for
        each segment are given in band_segm_start_end_normalized
        ''',
        a_legacy=LegacyDefinition(name='band_segm_labels_normalized'))

    band_segm_start_end_normalized = Quantity(
        type=np.dtype(np.float64),
        shape=[2, 3],
        description='''
        Fractional coordinates of the start and end point (in the basis of the reciprocal
        lattice vectors) of the segment sampled in the $k$ space. The conventional symbols
        (e.g., Gamma, K, L) of the same points are given in band_segm_labels
        ''',
        a_legacy=LegacyDefinition(name='band_segm_start_end_normalized'))

    number_of_normalized_k_points_per_segment = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of $k$ points in the segment of the normalized band structure
        (see section_k_band_segment_normalized).
        ''',
        a_legacy=LegacyDefinition(name='number_of_normalized_k_points_per_segment'))


class section_k_band_segment(MSection):
    '''
    Section collecting the information on a $k$-band or $q$-band segment. This section
    stores band structures along a one-dimensional pathway in the $k$ or $q$ (reciprocal)
    space.

    Eigenvalues calculated at the actual $k$-mesh used for energy_total evaluations are
    defined in section_eigenvalues and the band structures are represented as third-order
    tensors: one dimension for the spin channels, one for the sequence of $k$ or $q$
    points for the segment (given in number_of_k_points_per_segment), and one for the
    sequence of eigenvalues at a given $k$ or $q$ point. The values of the $k$ or $q$
    points in each segment are stored in band_k_points. The energies and occupation for
    each eigenstate, at each $k$ or $q$ point, segment, and spin channel are stored in
    band_energies and band_occupations, respectively. The labels for the segment are
    specified in band_segm_labels.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_k_band_segment'))

    band_energies = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_spin_channels', 'number_of_k_points_per_segment', 'number_of_band_segment_eigenvalues'],
        unit='joule',
        description='''
        $k$-dependent or $q$-dependent  energies of the electronic or vibrational band
        segment (electronic/vibrational band structure). This is a third-order tensor,
        with one dimension used for the spin channels (1 in case of a vibrational band
        structure), one for the $k$ or $q$ points for each segment, and one for the
        eigenvalue sequence.
        ''',
        a_legacy=LegacyDefinition(name='band_energies'))

    band_k_points = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_k_points_per_segment', 3],
        description='''
        Fractional coordinates of the $k$ or $q$ points (in the basis of the reciprocal-
        lattice vectors) for which the electronic energy are given.
        ''',
        a_legacy=LegacyDefinition(name='band_k_points'))

    band_occupations = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_spin_channels', 'number_of_k_points_per_segment', 'number_of_band_segment_eigenvalues'],
        description='''
        Occupation of the $k$-points along the electronic band. The size of the dimensions
        of this third-order tensor are the same as for the tensor in band_energies.
        ''',
        a_legacy=LegacyDefinition(name='band_occupations'))

    band_segm_labels = Quantity(
        type=str,
        shape=[2],
        description='''
        Start and end labels of the points in the segment (one-dimensional pathways)
        sampled in the $k$-space or $q$-space, using the conventional symbols, e.g.,
        Gamma, K, L. The coordinates (fractional, in the reciprocal space) of the start
        and end points for each segment are given in band_segm_start_end
        ''',
        a_legacy=LegacyDefinition(name='band_segm_labels'))

    band_segm_start_end = Quantity(
        type=np.dtype(np.float64),
        shape=[2, 3],
        description='''
        Fractional coordinates of the start and end point (in the basis of the reciprocal
        lattice vectors) of the segment sampled in the $k$ space. The conventional symbols
        (e.g., Gamma, K, L) of the same points are given in band_segm_labels
        ''',
        a_legacy=LegacyDefinition(name='band_segm_start_end'))

    number_of_k_points_per_segment = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of $k$ points in the segment of the band structure, see
        section_k_band_segment.
        ''',
        a_legacy=LegacyDefinition(name='number_of_k_points_per_segment'))


class section_band_gap(MSection):
    '''
    This section stores information for a band gap within a band structure.
    '''
    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_band_gap'))

    value = Quantity(
        type=float,
        unit="joule",
        description="""
        Band gap energy. Value of zero corresponds to a band structure without
        a band gap.
        """,
        a_legacy=LegacyDefinition(name='value')
    )
    type = Quantity(
        type=MEnum("direct", "indirect"),
        description="""
        Type of band gap.
        """,
        a_legacy=LegacyDefinition(name='type')
    )
    conduction_band_min_energy = Quantity(
        type=float,
        unit="joule",
        description="""
        Conduction band minimum energy.
        """,
        a_legacy=LegacyDefinition(name='conduction_band_min_energy')
    )
    valence_band_max_energy = Quantity(
        type=float,
        unit="joule",
        description="""
        Valence band maximum energy.
        """,
        a_legacy=LegacyDefinition(name='valence_band_max_energy')
    )
    conduction_band_min_k_point = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        unit="1 / meter",
        description="""
        Coordinate of the conduction band minimum in k-space.
        """,
        a_legacy=LegacyDefinition(name='conduction_band_min_k_point')
    )
    valence_band_max_k_point = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        unit="1 / meter",
        description="""
        Coordinate of the valence band minimum in k-space.
        """,
        a_legacy=LegacyDefinition(name='valence_band_max_k_point')
    )


class section_brillouin_zone(MSection):
    '''Defines a polyhedra for the Brillouin zone in reciprocal space.
    '''
    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_brillouin_zone'))

    vertices = Quantity(
        type=np.dtype(np.float64),
        shape=[3, "1..*"],
        description='''
        The vertices of the Brillouin zone corners as 3D coordinates in reciprocal space.
        ''',
        a_legacy=LegacyDefinition(name='vertices'))
    faces = Quantity(
        type=np.dtype(np.int32),
        shape=["1..*", "3..*"],
        description='''
        The faces of the Brillouin zone polyhedron as vertex indices. The
        surface normal is determined by a right-hand ordering of the points.
        ''',
        a_legacy=LegacyDefinition(name='faces'))


class section_k_band(MSection):
    '''
    This section stores information on a $k$-band (electronic or vibrational band
    structure) evaluation along one-dimensional pathways in the $k$ or $q$ (reciprocal)
    space given in section_k_band_segment. Eigenvalues calculated at the actual $k$-mesh
    used for energy_total evaluations, can be found in the section_eigenvalues section.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_k_band'))

    band_structure_kind = Quantity(
        type=str,
        shape=[],
        description='''
        String to specify the kind of band structure (either electronic or vibrational).
        ''',
        a_legacy=LegacyDefinition(name='band_structure_kind'))

    reciprocal_cell = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        unit="1 / meter",
        description="""
        The reciprocal cell within which the band structure is calculated.
        """,
        a_legacy=LegacyDefinition(name='reciprocal_cell')
    )

    brillouin_zone = SubSection(
        sub_section=SectionProxy('section_brillouin_zone'),
        repeats=False,
        a_legacy=LegacyDefinition(name='brillouin_zone'))

    section_band_gap = SubSection(
        sub_section=section_band_gap.m_def,
        repeats=True,
        description=""",
        Contains information for band gaps detected in the band structure.
        Contains a section for each spin channel in the same order as reported
        for the band energies. For channels without a band gap, a band gap
        value of zero is reported.
        """,
        a_legacy=LegacyDefinition(name='section_band_gap')
    )

    is_standard_path = Quantity(
        type=bool,
        description="""
        Boolean indicating whether the path follows the standard path for this
        bravais lattice. The AFLOW standard by Setyawan and Curtarolo is used
        (https://doi.org/10.1016/j.commatsci.2010.05.010).
        """,
        a_legacy=LegacyDefinition(name='is_standard_path')
    )

    section_k_band_segment = SubSection(
        sub_section=SectionProxy('section_k_band_segment'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_k_band_segment'))


class section_method_atom_kind(MSection):
    '''
    Every section_method_atom_kind section contains method-related information about a
    kind of atom, and is identified by one or more strings stored in
    method_atom_kind_label.

    This categorization into atom kinds is more flexible than just atomic species, because
    to different atoms of the same species different atom-centered basis sets or pseudo-
    potentials may be assigned. For instance, if two different oxygen atoms are assigned
    to different basis sets or pseudo-potentials, they have to distinguished into two
    different *kinds* of O atoms, by creating two distinct section_method_atom_kind
    sections.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_method_atom_kind'))

    method_atom_kind_atom_number = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Atomic number (number of protons) of this atom kind, use 0 if not an atom.
        ''',
        a_legacy=LegacyDefinition(name='method_atom_kind_atom_number'))

    method_atom_kind_explicit_electrons = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Number of explicit electrons (often called valence).
        ''',
        a_legacy=LegacyDefinition(name='method_atom_kind_explicit_electrons'))

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
        ''',
        a_legacy=LegacyDefinition(name='method_atom_kind_label'))

    method_atom_kind_mass = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='unified_atomic_mass_unit',
        description='''
        Mass of the kind of this kind of atoms.
        ''',
        a_legacy=LegacyDefinition(name='method_atom_kind_mass'))

    method_atom_kind_pseudopotential_name = Quantity(
        type=str,
        shape=[],
        description='''
        Name identifying the pseudopotential used.
        ''',
        a_legacy=LegacyDefinition(name='method_atom_kind_pseudopotential_name'))


class section_method_to_method_refs(MSection):
    '''
    Section that describes the relationship between different section_method sections.

    For instance, one calculation is a perturbation performed using a self-consistent
    field (SCF) calculation as starting point, or a simulated system is partitioned in
    regions with different but connected Hamiltonians (e.g., QM/MM, or a region treated
    via Kohn-Sham DFT embedded into a region treated via orbital-free DFT).

    The kind of relationship between the method defined in this section and the referenced
    one is described by method_to_method_kind. The referenced section section_method is
    identified via method_to_method_ref (typically used for a section_method section in
    the same section_run) or method_to_method_external_url.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_method_to_method_refs'))

    method_to_method_external_url = Quantity(
        type=str,
        shape=[],
        description='''
        URL used to reference an externally stored section_method. The kind of
        relationship between the present and the referenced section_method is specified by
        method_to_method_kind.
        ''',
        a_legacy=LegacyDefinition(name='method_to_method_external_url'))

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
        ''',
        a_legacy=LegacyDefinition(name='method_to_method_kind'))

    method_to_method_ref = Quantity(
        type=Reference(SectionProxy('section_method')),
        shape=[],
        description='''
        Reference to a local section_method. If both method_to_method_ref and
        method_to_method_external_url are given, then method_to_method_ref is a local copy
        of the value contained in method_to_method_external_url. The kind of relationship
        between the method defined in the present section_method and the referenced one is
        described by method_to_method_kind.
        ''',
        a_legacy=LegacyDefinition(name='method_to_method_ref'))


class section_method(MSection):
    '''
    Section containing the various parameters that define the theory and the
    approximations (convergence, thresholds,...) to perform a *single configuration
    calculation*, see section_single_configuration_calculation.

    *NOTE*: This section does not contain settings for molecular dynamics, geometry
    optimization etc. See section frame_sequence for these other settings instead.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_method'))

    basis_set = Quantity(
        type=str,
        shape=[],
        description='''
        Unique string identifying the basis set used for the final wavefunctions
        calculated with XC_method. It might identify a class of basis sets, often matches
        one of the strings given in any of basis_set_name.
        ''',
        categories=[settings_numerical_parameter, settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='basis_set'))

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
        ''',
        a_legacy=LegacyDefinition(name='calculation_method_current'))

    calculation_method_kind = Quantity(
        type=str,
        shape=[],
        description='''
        Kind of method in calculation_method_current.

        Accepted values are:

        - absolute

        - perturbative.
        ''',
        a_legacy=LegacyDefinition(name='calculation_method_kind'))

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
        ''',
        a_legacy=LegacyDefinition(name='calculation_method'))

    electronic_structure_method = Quantity(
        type=str,
        shape=[],
        description='''
        Non-unique string identifying the used electronic structure method. It is not
        unique in the sense that two calculations with the same
        electronic_structure_method string may have not been performed with exactly the
        same method. The allowed strings are given in the [electronic structure method
        wiki page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-
        info/wikis/metainfo/electronic-structure-method).
        ''',
        categories=[settings_XC, settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='electronic_structure_method'))

    k_mesh_points = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_k_mesh_points', 3],
        description='''
        List of all the k points in the $k$-point mesh. These are the k point used to
        evaluate energy_total, and are in fractional coordinates (in the basis of the
        reciprocal-lattice vectors).
        ''',
        categories=[settings_k_points, settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='k_mesh_points'))

    k_mesh_weights = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_k_mesh_points'],
        description='''
        Weights of all the k points in the $k$-point mesh. These are the weights for
        k_mesh_points (i.e. the k point used to evaluate energy_total).
        ''',
        categories=[settings_k_points, settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='k_mesh_weights'))

    number_of_k_mesh_points = Quantity(
        type=int,
        shape=[],
        description='''
        number of k points in the mesh (i.e. the k points used to evaluate energy_total).
        ''',
        categories=[settings_k_points, settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='number_of_k_mesh_points'))

    number_of_spin_channels = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of spin channels, see section_method.
        ''',
        a_legacy=LegacyDefinition(name='number_of_spin_channels'))

    relativity_method = Quantity(
        type=str,
        shape=[],
        description='''
        Describes the relativistic treatment used for the calculation of the final energy
        and related quantities. If skipped or empty, no relativistic treatment is applied.
        ''',
        categories=[settings_relativity, settings_XC, settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='relativity_method'))

    scf_max_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Specifies the maximum number of allowed self-consistent field (SCF) iterations in
        a calculation run, see section_run.
        ''',
        categories=[settings_scf],
        a_legacy=LegacyDefinition(name='scf_max_iteration'))

    scf_threshold_energy_change = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Specifies the threshold for the energy_total_scf_iteration change between two
        subsequent self-consistent field (SCF) iterations. The SCF is considered converged
        when the total-energy change between two SCF cycles is below the threshold
        (possibly in combination with other criteria).
        ''',
        categories=[settings_scf],
        a_legacy=LegacyDefinition(name='scf_threshold_energy_change'))

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
        ''',
        categories=[settings_self_interaction_correction, settings_XC, settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='self_interaction_correction_method'))

    smearing_kind = Quantity(
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
        ''',
        categories=[settings_smearing],
        a_legacy=LegacyDefinition(name='smearing_kind'))

    smearing_width = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Specifies the width of the smearing in energy for the electron occupation used to
        calculate the free energy (see energy_free).

        *NOTE:* Not all methods specified in smearing_kind uses this value.
        ''',
        categories=[settings_smearing],
        a_legacy=LegacyDefinition(name='smearing_width'))

    spin_target_multiplicity = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Stores the target (user-imposed) value of the spin multiplicity $M=2S+1$, where
        $S$ is the total spin. It is an integer number. This value is not necessarily the
        value obtained at the end of the calculation. See spin_S2 for the converged value
        of the spin moment.
        ''',
        a_legacy=LegacyDefinition(name='spin_target_multiplicity'))

    stress_tensor_method = Quantity(
        type=str,
        shape=[],
        description='''
        Specifies the method used to calculate stress_tensor for, e.g., molecular dynamics
        and geometry optimization.

        The allowed values are:

        * numeric

        * analytic
        ''',
        categories=[settings_stress_tensor],
        a_legacy=LegacyDefinition(name='stress_tensor_method'))

    total_charge = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        unit='coulomb',
        description='''
        Provides the total amount of charge of the system in a run.
        ''',
        a_legacy=LegacyDefinition(name='total_charge'))

    van_der_Waals_method = Quantity(
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
        ''',
        categories=[settings_van_der_Waals, settings_XC, settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='van_der_Waals_method'))

    XC_functional = Quantity(
        type=str,
        shape=[],
        description='''
        This value describes a DFT exchange-correlation (XC) functional used for
        evaluating the energy value stored in energy_XC_functional and related quantities
        (e.g., forces).

        It is a unique short name obtained by combining the data stored in
        section_XC_functionals, more specifically by combining different
        XC_functional_name as described in the [XC_functional wiki
        page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/XC-
        functional).
        ''',
        categories=[settings_potential_energy_surface, settings_physical_parameter, settings_XC_functional, settings_XC],
        a_legacy=LegacyDefinition(name='XC_functional'))

    XC_method = Quantity(
        type=str,
        shape=[],
        description='''
        Describes the exchange correlation (XC) method used for evaluating the XC energy
        (energy_XC). Differently from XC_functional, perturbative treatments are also
        accounted for, where the string contains the reference to both the perturbative
        (e.g., MP2) and the starting point (e.g, Hartree-Fock) XC method defined in the
        section section_method.

        The value consists of XC_method_current concatenated with the `@` character and
        the XC method (XC_method) defined in section_method that is referred to by
        method_to_method_ref where method_to_method_kind = "starting_point_method".
        ''',
        categories=[settings_XC, settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='XC_method'))

    XC_method_current = Quantity(
        type=str,
        shape=[],
        description='''
        Identifies the exchange correlation (XC) method used for energy_XC and related
        quantities in a standardized short form as a string.

        It is built by joining the values in the following order using the underscore `_`
        character: electronic_structure_method, XC_functional,
        self_interaction_correction_method, van_der_Waals_method and relativity_method.

        If any of the methods listed in the string contain non-standard settings, then the
        first 10 characters of the Base64 URL encoding of SHA 512 checksum of a normalized
        JSON with all non-redundant non-derived settings_XC are appended to the the string
        preceded by an underscore.

        With empty strings, the underscore `_` character is skipped.

        If the method defined in the section_method section is perturbative, the
        XC_method_current contains only the perturbative method, not the starting point
        (e.g. the DFT XC functional used as a starting point for a RPA perturbative
        calculation). In this case, the string that contains both the perturbative and
        starting point method is stored in XC_method.
        ''',
        categories=[settings_XC, settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='XC_method_current'))

    section_method_atom_kind = SubSection(
        sub_section=SectionProxy('section_method_atom_kind'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_method_atom_kind'))

    section_method_to_method_refs = SubSection(
        sub_section=SectionProxy('section_method_to_method_refs'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_method_to_method_refs'))

    section_XC_functionals = SubSection(
        sub_section=SectionProxy('section_XC_functionals'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_XC_functionals'))


class section_original_system(MSection):
    '''
    Section containing symmetry information that is specific to the original system.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_original_system'))

    equivalent_atoms_original = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_atoms'],
        description='''
        Gives a mapping table of atoms to symmetrically independent atoms in the original
        cell. This is used to find symmetrically equivalent atoms.
        ''',
        a_legacy=LegacyDefinition(name='equivalent_atoms_original'))

    wyckoff_letters_original = Quantity(
        type=str,
        shape=['number_of_atoms'],
        description='''
        Wyckoff letters for atoms in the original cell.
        ''',
        a_legacy=LegacyDefinition(name='wyckoff_letters_original'))


class section_primitive_system(MSection):
    '''
    Section containing symmetry information that is specific to the primitive system. The
    primitive system is derived from the standardized system with a transformation that is
    specific to the centring. The transformation matrices can be found e.g. from here:
    https://atztogo.github.io/spglib/definition.html#transformation-to-the-primitive-cell
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_primitive_system'))

    atom_positions_primitive = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms_primitive', 3],
        description='''
        Atom positions in the primitive cell in reduced units.
        ''',
        a_legacy=LegacyDefinition(name='atom_positions_primitive'))

    atomic_numbers_primitive = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_atoms_primitive'],
        description='''
        Atomic numbers in the primitive cell.
        ''',
        a_legacy=LegacyDefinition(name='atomic_numbers_primitive'))

    equivalent_atoms_primitive = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_atoms_primitive'],
        description='''
        Gives a mapping table of atoms to symmetrically independent atoms in the primitive
        cell. This is used to find symmetrically equivalent atoms.
        ''',
        a_legacy=LegacyDefinition(name='equivalent_atoms_primitive'))

    lattice_vectors_primitive = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        unit='meter',
        description='''
        Primitive lattice vectors. The vectors are the rows of this matrix.
        ''',
        a_legacy=LegacyDefinition(name='lattice_vectors_primitive'))

    number_of_atoms_primitive = Quantity(
        type=int,
        shape=[],
        description='''
        Number of atoms in primitive system.
        ''',
        a_legacy=LegacyDefinition(name='number_of_atoms_primitive'))

    wyckoff_letters_primitive = Quantity(
        type=str,
        shape=['number_of_atoms_primitive'],
        description='''
        Wyckoff letters for atoms in the primitive cell.
        ''',
        a_legacy=LegacyDefinition(name='wyckoff_letters_primitive'))


class section_processor_info(MSection):
    '''
    Section with information about a processor that generated or added information to the
    current calculation.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_processor_info'))

    processor_id = Quantity(
        type=str,
        shape=[],
        description='''
        Id (name+version) of the processor that generated or added information to the
        current calculation.
        ''',
        a_legacy=LegacyDefinition(name='processor_id'))

    processor_number_of_evaluated_contexts = Quantity(
        type=np.dtype(np.int64),
        shape=[],
        description='''
        number of contexts evaluated with this processor in the current current
        calculation.
        ''',
        a_legacy=LegacyDefinition(name='processor_number_of_evaluated_contexts'))

    processor_number_of_failed_contexts = Quantity(
        type=np.dtype(np.int64),
        shape=[],
        description='''
        number of contexts in the current current calculation that had failure for this
        processor.
        ''',
        a_legacy=LegacyDefinition(name='processor_number_of_failed_contexts'))

    processor_number_of_skipped_contexts = Quantity(
        type=np.dtype(np.int64),
        shape=[],
        description='''
        number of contexts skipped by this processor in the current current calculation.
        ''',
        a_legacy=LegacyDefinition(name='processor_number_of_skipped_contexts'))

    processor_number_of_successful_contexts = Quantity(
        type=np.dtype(np.int64),
        shape=[],
        description='''
        number of contexts in the current calculation that where successfully handled by
        this processor.
        ''',
        a_legacy=LegacyDefinition(name='processor_number_of_successful_contexts'))

    processor_version_details = Quantity(
        type=typing.Any,
        shape=[],
        description='''
        detailed version information on the processor that generated or added information
        to the current calculation.
        ''',
        a_legacy=LegacyDefinition(name='processor_version_details'))


class section_processor_log_event(MSection):
    '''
    A log event
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_processor_log_event'))

    processor_log_event_level = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Level of the logging, a lower number has more priority. The levels are the same as
        log4j: FATAL -> 100, ERROR -> 200, WARN -> 300, INFO -> 400, DEBUG -> 500, TRACE
        -> 600
        ''',
        a_legacy=LegacyDefinition(name='processor_log_event_level'))

    processor_log_event_message = Quantity(
        type=str,
        shape=[],
        description='''
        The log message
        ''',
        a_legacy=LegacyDefinition(name='processor_log_event_message'))


class section_processor_log(MSection):
    '''
    log of a processor
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_processor_log'))

    processor_log_processor_id = Quantity(
        type=str,
        shape=[],
        description='''
        The processor id of the processor creating this log
        ''',
        a_legacy=LegacyDefinition(name='processor_log_processor_id'))

    processor_log_start = Quantity(
        type=str,
        shape=[],
        description='''
        Start of the log (in ansi notation YYYY-MM-TT...)
        ''',
        a_legacy=LegacyDefinition(name='processor_log_start'))

    section_processor_log_event = SubSection(
        sub_section=SectionProxy('section_processor_log_event'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_processor_log_event'))


class section_prototype(MSection):
    '''
    Information on the prototype corresponding to the current section.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_prototype'))

    prototype_aflow_id = Quantity(
        type=str,
        shape=[],
        description='''
        AFLOW id of the prototype (see
        http://aflowlib.org/CrystalDatabase/prototype_index.html) identified on the basis
        of the space_group and normalized_wyckoff.
        ''',
        a_legacy=LegacyDefinition(name='prototype_aflow_id'))

    prototype_aflow_url = Quantity(
        type=str,
        shape=[],
        description='''
        Url to the AFLOW definition of the prototype (see
        http://aflowlib.org/CrystalDatabase/prototype_index.html) identified on the basis
        of the space_group and normalized_wyckoff.
        ''',
        a_legacy=LegacyDefinition(name='prototype_aflow_url'))

    prototype_assignment_method = Quantity(
        type=str,
        shape=[],
        description='''
        Method used to identify the prototype.
        ''',
        a_legacy=LegacyDefinition(name='prototype_assignment_method'))

    prototype_label = Quantity(
        type=str,
        shape=[],
        description='''
        Label of the prototype identified on the basis of the space_group and
        normalized_wyckoff. The label is in the same format as in the read_prototypes
        function: <space_group_number>-<prototype_name>-<Pearson's symbol>).
        ''',
        a_legacy=LegacyDefinition(name='prototype_label'))


class section_run(MSection):
    '''
    Every section_run represents a single call of a program. What exactly is contained in
    a run depends on the run type (see for example section_method and
    section_single_configuration_calculation) and the program (see [program_info
    ](program_info)).
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_run'))

    calculation_file_uri = Quantity(
        type=str,
        shape=[],
        description='''
        Contains the nomad uri of a raw the data file connected to the current run. There
        should be an value for the main_file_uri and all ancillary files.
        ''',
        a_legacy=LegacyDefinition(name='calculation_file_uri'))

    message_debug_run = Quantity(
        type=str,
        shape=[],
        description='''
        A debugging message of the computational program, associated with a run.
        ''',
        categories=[message_debug],
        a_legacy=LegacyDefinition(name='message_debug_run'))

    message_error_run = Quantity(
        type=str,
        shape=[],
        description='''
        An error message of the computational program, associated with a run.
        ''',
        categories=[message_info, message_debug, message_error, message_warning],
        a_legacy=LegacyDefinition(name='message_error_run'))

    message_info_run = Quantity(
        type=str,
        shape=[],
        description='''
        An information message of the computational program, associated with a run.
        ''',
        categories=[message_info, message_debug],
        a_legacy=LegacyDefinition(name='message_info_run'))

    message_warning_run = Quantity(
        type=str,
        shape=[],
        description='''
        A warning message of the computational program, associated with a run.
        ''',
        categories=[message_info, message_debug, message_warning],
        a_legacy=LegacyDefinition(name='message_warning_run'))

    parsing_message_debug_run = Quantity(
        type=str,
        shape=[],
        description='''
        This field is used for debugging messages of the parsing program associated with a
        single configuration calculation, see section_single_configuration_calculation.
        ''',
        categories=[parsing_message_debug],
        a_legacy=LegacyDefinition(name='parsing_message_debug_run'))

    parsing_message_error_run = Quantity(
        type=str,
        shape=[],
        description='''
        This field is used for error messages of the parsing program associated with a
        run, see section_run.
        ''',
        categories=[parsing_message_info, parsing_message_error, parsing_message_warning, parsing_message_debug],
        a_legacy=LegacyDefinition(name='parsing_message_error_run'))

    parsing_message_info_run = Quantity(
        type=str,
        shape=[],
        description='''
        This field is used for info messages of the parsing program associated with a run,
        see section_run.
        ''',
        categories=[parsing_message_info, parsing_message_debug],
        a_legacy=LegacyDefinition(name='parsing_message_info_run'))

    parsing_message_warning_run = Quantity(
        type=str,
        shape=[],
        description='''
        This field is used for warning messages of the parsing program associated with a
        run, see section_run.
        ''',
        categories=[parsing_message_info, parsing_message_warning, parsing_message_debug],
        a_legacy=LegacyDefinition(name='parsing_message_warning_run'))

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
        ''',
        a_legacy=LegacyDefinition(name='program_basis_set_type'))

    program_compilation_datetime = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        Contains the program compilation date and time from *Unix epoch* (00:00:00 UTC on
        1 January 1970) in seconds. For date and times without a timezone, the default
        timezone GMT is used.
        ''',
        categories=[accessory_info, program_info],
        a_legacy=LegacyDefinition(name='program_compilation_datetime'))

    program_compilation_host = Quantity(
        type=str,
        shape=[],
        description='''
        Specifies the host on which the program was compiled.
        ''',
        categories=[accessory_info, program_info],
        a_legacy=LegacyDefinition(name='program_compilation_host'))

    program_name = Quantity(
        type=str,
        shape=[],
        description='''
        Specifies the name of the program that generated the data.
        ''',
        categories=[accessory_info, program_info],
        a_legacy=LegacyDefinition(name='program_name'))

    program_version = Quantity(
        type=str,
        shape=[],
        description='''
        Specifies the version of the program that was used. This should be the version
        number of an official release, the version tag or a commit id as well as the
        location of the repository.
        ''',
        categories=[accessory_info, program_info],
        a_legacy=LegacyDefinition(name='program_version'))

    run_clean_end = Quantity(
        type=bool,
        shape=[],
        description='''
        Indicates whether this run terminated properly (true), or if it was killed or
        exited with an error code unequal to zero (false).
        ''',
        a_legacy=LegacyDefinition(name='run_clean_end'))

    run_hosts = Quantity(
        type=typing.Any,
        shape=[],
        description='''
        An associative list of host(s) that performed this simulation. This is an
        associative list that contains program-dependent information (*key*) on how the
        host was used (*value*). Useful for debugging purposes.
        ''',
        categories=[parallelization_info, accessory_info],
        a_legacy=LegacyDefinition(name='run_hosts'))

    raw_id = Quantity(
        type=str,
        shape=[],
        description='''
        An optional calculation id, if one is found in the code input/output files.
        ''',
        a_legacy=LegacyDefinition(name='raw_id'))

    time_run_cpu1_end = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        Stores the end time of the run on CPU 1.
        ''',
        categories=[time_info, accessory_info],
        a_legacy=LegacyDefinition(name='time_run_cpu1_end'))

    time_run_cpu1_start = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        Stores the start time of the run on CPU 1.
        ''',
        categories=[time_info, accessory_info],
        a_legacy=LegacyDefinition(name='time_run_cpu1_start'))

    time_run_date_end = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        Stores the end date of the run as time since the *Unix epoch* (00:00:00 UTC on 1
        January 1970) in seconds. For date and times without a timezone, the default
        timezone GMT is used.
        ''',
        categories=[time_info, accessory_info],
        a_legacy=LegacyDefinition(name='time_run_date_end'))

    time_run_date_start = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        Stores the start date of the run as time since the *Unix epoch* (00:00:00 UTC on 1
        January 1970) in seconds. For date and times without a timezone, the default
        timezone GMT is used.
        ''',
        categories=[time_info, accessory_info],
        a_legacy=LegacyDefinition(name='time_run_date_start'))

    time_run_wall_end = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        Stores the internal wall-clock time at the end of the run.
        ''',
        categories=[time_info, accessory_info],
        a_legacy=LegacyDefinition(name='time_run_wall_end'))

    time_run_wall_start = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        Stores the internal wall-clock time from the start of the run.
        ''',
        categories=[time_info, accessory_info],
        a_legacy=LegacyDefinition(name='time_run_wall_start'))

    section_basis_set_atom_centered = SubSection(
        sub_section=SectionProxy('section_basis_set_atom_centered'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_basis_set_atom_centered'))

    section_basis_set_cell_dependent = SubSection(
        sub_section=SectionProxy('section_basis_set_cell_dependent'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_basis_set_cell_dependent'))

    section_frame_sequence = SubSection(
        sub_section=SectionProxy('section_frame_sequence'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_frame_sequence'))

    section_method = SubSection(
        sub_section=SectionProxy('section_method'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_method'))

    section_sampling_method = SubSection(
        sub_section=SectionProxy('section_sampling_method'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_sampling_method'))

    section_single_configuration_calculation = SubSection(
        sub_section=SectionProxy('section_single_configuration_calculation'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_single_configuration_calculation'))

    section_system = SubSection(
        sub_section=SectionProxy('section_system'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_system'))


class section_sampling_method(MSection):
    '''
    Section containing the settings describing a (potential-energy surface) sampling
    method.

    Results and monitored quantities of such sampling are collected in a sequence of
    frames, section_frame_sequence.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_sampling_method'))

    ensemble_type = Quantity(
        type=str,
        shape=[],
        description='''
        Kind of sampled ensemble stored in section_frame_sequence; valid values are
        defined in [ensemble_type wiki page](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-
        meta-info/wikis/metainfo/ensemble-type).
        ''',
        a_legacy=LegacyDefinition(name='ensemble_type'))

    geometry_optimization_energy_change = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Value of threshold for the energy_total change between two geometry optimization
        steps, as convergence criterion of the geometry_optimization_method. A geometry
        optimization is considered converged when the energy_total change between two
        geometry optimization steps is below the threshold (possibly in combination with
        other criteria)
        ''',
        categories=[settings_geometry_optimization, settings_sampling],
        a_legacy=LegacyDefinition(name='geometry_optimization_energy_change'))

    geometry_optimization_geometry_change = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter',
        description='''
        Value of threshold for the displacement of the nuclei between two geometry
        optimization steps as convergence criterion of the geometry_optimization_method. A
        geometry optimization is considered converged when the maximum among the
        displacements of the nuclei between two geometry optimization steps is below the
        threshold (possibly in combination with other criteria)
        ''',
        categories=[settings_geometry_optimization, settings_sampling],
        a_legacy=LegacyDefinition(name='geometry_optimization_geometry_change'))

    geometry_optimization_method = Quantity(
        type=str,
        shape=[],
        description='''
        Algorithm for the geometry optimization. Allowed values are listed in the
        [geometry_optimization_method wiki page](https://gitlab.mpcdf.mpg.de/nomad-
        lab/nomad-meta-info/wikis/metainfo/geometry-optimization-method).
        ''',
        categories=[settings_geometry_optimization, settings_sampling],
        a_legacy=LegacyDefinition(name='geometry_optimization_method'))

    geometry_optimization_threshold_force = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='newton',
        description='''
        Value of threshold for the force modulus as convergence criterion of the
        geometry_optimization_method. A geometry optimization is considered converged when
        the maximum of the moduli of the force on each of the atoms is below this
        threshold (possibly in combination with other criteria)
        ''',
        categories=[settings_geometry_optimization, settings_sampling],
        a_legacy=LegacyDefinition(name='geometry_optimization_threshold_force'))

    sampling_method_expansion_order = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Order up to which the potential energy surface was expanded in a Taylor series
        (see sampling_method).
        ''',
        a_legacy=LegacyDefinition(name='sampling_method_expansion_order'))

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
        ''',
        a_legacy=LegacyDefinition(name='sampling_method'))


class section_scf_iteration(MSection):
    '''
    Every section_scf_iteration represents a self-consistent field (SCF) iteration, see
    scf_info, and gives detailed information on the SCF procedure of the specified
    quantities.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_scf_iteration'))

    electronic_kinetic_energy_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Electronic kinetic energy as defined in XC_method during the self-consistent field
        (SCF) iterations.
        ''',
        categories=[energy_component, energy_value],
        a_legacy=LegacyDefinition(name='electronic_kinetic_energy_scf_iteration'))

    energy_change_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Stores the change of total energy with respect to the previous self-consistent
        field (SCF) iteration.
        ''',
        categories=[error_estimate_contribution, energy_value],
        a_legacy=LegacyDefinition(name='energy_change_scf_iteration'))

    energy_correction_entropy_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Entropy correction to the potential energy to compensate for the change in
        occupation so that forces at finite T do not need to keep the change of occupation
        in account. The array lists the values of the entropy correction for each self-
        consistent field (SCF) iteration. Defined consistently with XC_method.
        ''',
        categories=[energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_correction_entropy_scf_iteration'))

    energy_correction_hartree_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Correction to the density-density electrostatic energy in the sum of eigenvalues
        (that uses the mixed density on one side), and the fully consistent density-
        density electrostatic energy during the self-consistent field (SCF) iterations.
        Defined consistently with XC_method.
        ''',
        categories=[energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_correction_hartree_scf_iteration'))

    energy_electrostatic_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Total electrostatic energy (nuclei + electrons) during each self-consistent field
        (SCF) iteration.
        ''',
        categories=[energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_electrostatic_scf_iteration'))

    energy_free_per_atom_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Free energy per atom (whose minimum gives the smeared occupation density
        calculated with smearing_kind) calculated with XC_method during the self-
        consistent field (SCF) iterations.
        ''',
        categories=[energy_component_per_atom, energy_value],
        a_legacy=LegacyDefinition(name='energy_free_per_atom_scf_iteration'))

    energy_free_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Free energy (whose minimum gives the smeared occupation density calculated with
        smearing_kind) calculated with the method described in XC_method during the self-
        consistent field (SCF) iterations.
        ''',
        categories=[energy_component, energy_value, energy_total_potential],
        a_legacy=LegacyDefinition(name='energy_free_scf_iteration'))

    energy_hartree_error_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Error in the Hartree (electrostatic) potential energy during each self-consistent
        field (SCF) iteration. Defined consistently with XC_method.
        ''',
        categories=[error_estimate_contribution, energy_value],
        a_legacy=LegacyDefinition(name='energy_hartree_error_scf_iteration'))

    energy_sum_eigenvalues_per_atom_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Value of the energy per atom, where the energy is defined as the sum of the
        eigenvalues of the Hamiltonian matrix given by XC_method, during each self-
        consistent field (SCF) iteration.
        ''',
        categories=[energy_component_per_atom, energy_value],
        a_legacy=LegacyDefinition(name='energy_sum_eigenvalues_per_atom_scf_iteration'))

    energy_sum_eigenvalues_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Sum of the eigenvalues of the Hamiltonian matrix defined by XC_method, during each
        self-consistent field (SCF) iteration.
        ''',
        categories=[energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_sum_eigenvalues_scf_iteration'))

    energy_total_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Value of the total electronic energy calculated with the method described in
        XC_method during each self-consistent field (SCF) iteration.
        ''',
        categories=[energy_component, energy_value, energy_total_potential],
        a_legacy=LegacyDefinition(name='energy_total_scf_iteration'))

    energy_total_T0_per_atom_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Value of the total energy, calculated with the method described in XC_method per
        atom extrapolated to $T=0$, based on a free-electron gas argument, during each
        self-consistent field (SCF) iteration.
        ''',
        categories=[energy_total_potential_per_atom, energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_total_T0_per_atom_scf_iteration'))

    energy_total_T0_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Value of the total energy (or equivalently free energy), calculated with the
        method described in XC_method and extrapolated to $T=0$, based on a free-electron
        gas argument, during each self-consistent field (SCF) iteration.
        ''',
        categories=[energy_component, energy_value, energy_total_potential],
        a_legacy=LegacyDefinition(name='energy_total_T0_scf_iteration'))

    energy_XC_potential_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Value for exchange-correlation (XC) potential energy: the integral of the first
        order derivative of the functional stored in XC_functional (integral of
        v_xc*electron_density), i.e., the component of XC that is in the sum of the
        eigenvalues. Values are given for each self-consistent field (SCF) iteration
        (i.e., not the converged value, the latter being stored in energy_XC_potential).
        ''',
        categories=[energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_XC_potential_scf_iteration'))

    energy_XC_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Value for exchange-correlation (XC) energy obtained during each self-consistent
        field (SCF) iteration, using the method described in XC_method.
        ''',
        categories=[energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_XC_scf_iteration'))

    spin_S2_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Stores the value of the total spin moment operator $S^2$ during the self-
        consistent field (SCF) iterations of the XC_method. It can be used to calculate
        the spin contamination in spin-unrestricted calculations.
        ''',
        a_legacy=LegacyDefinition(name='spin_S2_scf_iteration'))

    time_scf_iteration_cpu1_end = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        Stores the end time of a self-consistent field (SCF) iteration on CPU 1.
        ''',
        categories=[time_info, accessory_info],
        a_legacy=LegacyDefinition(name='time_scf_iteration_cpu1_end'))

    time_scf_iteration_cpu1_start = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        Stores the start time of a self-consistent field (SCF) iteration on CPU 1.
        ''',
        categories=[time_info, accessory_info],
        a_legacy=LegacyDefinition(name='time_scf_iteration_cpu1_start'))

    time_scf_iteration_date_end = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        Stores the end date of a self-consistent field (SCF) iteration as time since the
        *Unix epoch* (00:00:00 UTC on 1 January 1970) in seconds. For date and times
        without a timezone, the default timezone GMT is used.
        ''',
        categories=[time_info, accessory_info],
        a_legacy=LegacyDefinition(name='time_scf_iteration_date_end'))

    time_scf_iteration_date_start = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        Stores the start date of a self-consistent field (SCF) iteration as time since the
        *Unix epoch* (00:00:00 UTC on 1 January 1970) in seconds. For date and times
        without a timezone, the default timezone GMT is used.
        ''',
        categories=[time_info, accessory_info],
        a_legacy=LegacyDefinition(name='time_scf_iteration_date_start'))

    time_scf_iteration_wall_end = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        Stores the internal wall-clock time at the end of a self-consistent field (SCF)
        iteration.
        ''',
        categories=[time_info, accessory_info],
        a_legacy=LegacyDefinition(name='time_scf_iteration_wall_end'))

    time_scf_iteration_wall_start = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        Stores the internal wall-clock time from the start of a self-consistent field
        (SCF) iteration.
        ''',
        categories=[time_info, accessory_info],
        a_legacy=LegacyDefinition(name='time_scf_iteration_wall_start'))


class section_single_configuration_calculation(MSection):
    '''
    Every section_single_configuration_calculation section contains the values computed
    during a *single configuration calculation*, i.e. a calculation performed on a given
    configuration of the system (as defined in section_system) and a given computational
    method (e.g., exchange-correlation method, basis sets, as defined in section_method).

    The link between the current section_single_configuration_calculation and the related
    section_system and section_method sections is established by the values stored in
    single_configuration_calculation_to_system_ref and
    single_configuration_to_calculation_method_ref, respectively.

    The reason why information on the system configuration and computational method is
    stored separately is that several *single configuration calculations* can be performed
    on the same system configuration, viz. several system configurations can be evaluated
    with the same computational method. This storage strategy avoids redundancies.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_single_configuration_calculation'))

    atom_forces_free_raw = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        unit='newton',
        description='''
        Forces acting on the atoms, calculated as minus gradient of energy_free,
        **without** constraints. The derivatives with respect to displacements of nuclei
        are evaluated in Cartesian coordinates. The (electronic) energy_free contains the
        change in (fractional) occupation of the electronic eigenstates, which are
        accounted for in the derivatives, yielding a truly energy-conserved quantity.
        These forces may contain unitary transformations (center-of-mass translations and
        rigid rotations for non-periodic systems) that are normally filtered separately
        (see atom_forces_free for the filtered counterpart). Forces due to constraints
        such as fixed atoms, distances, angles, dihedrals, etc. are also considered
        separately (see atom_forces_free for the filtered counterpart).
        ''',
        categories=[atom_forces_type],
        a_legacy=LegacyDefinition(name='atom_forces_free_raw'))

    atom_forces_free = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        unit='newton',
        description='''
        Forces acting on the atoms, calculated as minus gradient of energy_free,
        **including** constraints, if present. The derivatives with respect to
        displacements of the nuclei are evaluated in Cartesian coordinates. The
        (electronic) energy_free contains the information on the change in (fractional)
        occupation of the electronic eigenstates, which are accounted for in the
        derivatives, yielding a truly energy-conserved quantity. In addition, these forces
        are obtained by filtering out the unitary transformations (center-of-mass
        translations and rigid rotations for non-periodic systems, see
        atom_forces_free_raw for the unfiltered counterpart). Forces due to constraints
        such as fixed atoms, distances, angles, dihedrals, etc. are included (see
        atom_forces_free_raw for the unfiltered counterpart).
        ''',
        categories=[atom_forces_type],
        a_legacy=LegacyDefinition(name='atom_forces_free'))

    atom_forces_raw = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        unit='newton',
        description='''
        Forces acting on the atoms, calculated as minus gradient of energy_total,
        **without** constraints. The derivatives with respect to displacements of the
        nuclei are evaluated in Cartesian coordinates. These forces may contain unitary
        transformations (center-of-mass translations and rigid rotations for non-periodic
        systems) that are normally filtered separately (see atom_forces for the filtered
        counterpart). Forces due to constraints such as fixed atoms, distances, angles,
        dihedrals, etc. are also considered separately (see atom_forces for the filtered
        counterpart).
        ''',
        categories=[atom_forces_type],
        a_legacy=LegacyDefinition(name='atom_forces_raw'))

    atom_forces_T0_raw = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        unit='newton',
        description='''
        Forces acting on the atoms, calculated as minus gradient of energy_total_T0,
        **without** constraints. The derivatives with respect to displacements of the
        nuclei are evaluated in Cartesian coordinates. These forces may contain unitary
        transformations (center-of-mass translations and rigid rotations for non-periodic
        systems) that are normally filtered separately (see atom_forces_T0 for the
        filtered counterpart). Forces due to constraints such as fixed atoms, distances,
        angles, dihedrals, etc. are also considered separately (see atom_forces_T0 for the
        filtered counterpart).
        ''',
        categories=[atom_forces_type],
        a_legacy=LegacyDefinition(name='atom_forces_T0_raw'))

    atom_forces_T0 = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        unit='newton',
        description='''
        Forces acting on the atoms, calculated as minus gradient of energy_total_T0,
        **including** constraints, if present. The derivatives with respect to
        displacements of the nuclei are evaluated in Cartesian coordinates. In addition,
        these forces are obtained by filtering out the unitary transformations (center-of-
        mass translations and rigid rotations for non-periodic systems, see
        atom_forces_free_T0_raw for the unfiltered counterpart). Forces due to constraints
        such as fixed atoms, distances, angles, dihedrals, etc. are also included (see
        atom_forces_free_T0_raw for the unfiltered counterpart).
        ''',
        categories=[atom_forces_type],
        a_legacy=LegacyDefinition(name='atom_forces_T0'))

    atom_forces = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        unit='newton',
        description='''
        Forces acting on the atoms, calculated as minus gradient of energy_total,
        **including** constraints, if present. The derivatives with respect to
        displacements of nuclei are evaluated in Cartesian coordinates. In addition, these
        forces are obtained by filtering out the unitary transformations (center-of-mass
        translations and rigid rotations for non-periodic systems, see
        atom_forces_free_raw for the unfiltered counterpart). Forces due to constraints
        such as fixed atoms, distances, angles, dihedrals, etc. are included (see
        atom_forces_raw for the unfiltered counterpart).
        ''',
        categories=[atom_forces_type],
        a_legacy=LegacyDefinition(name='atom_forces'))

    electronic_kinetic_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Self-consistent electronic kinetic energy as defined in XC_method.
        ''',
        categories=[energy_component, energy_value],
        a_legacy=LegacyDefinition(name='electronic_kinetic_energy'))

    energy_C = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Correlation (C) energy calculated with the method described in XC_functional.
        ''',
        categories=[energy_component, energy_value, energy_type_C],
        a_legacy=LegacyDefinition(name='energy_C'))

    energy_correction_entropy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Entropy correction to the potential energy to compensate for the change in
        occupation so that forces at finite T do not need to keep the change of occupation
        in account. Defined consistently with XC_method.
        ''',
        categories=[energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_correction_entropy'))

    energy_correction_hartree = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Correction to the density-density electrostatic energy in the sum of eigenvalues
        (that uses the mixed density on one side), and the fully consistent density-
        density electrostatic energy. Defined consistently with XC_method.
        ''',
        categories=[energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_correction_hartree'))

    energy_current = Quantity(
        type=np.dtype(np.float64),
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
        ''',
        categories=[energy_component, energy_value, energy_total_potential],
        a_legacy=LegacyDefinition(name='energy_current'))

    energy_electrostatic = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Total electrostatic energy (nuclei + electrons), defined consistently with
        calculation_method.
        ''',
        categories=[energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_electrostatic'))

    energy_free_per_atom = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Free energy per atom (whose minimum gives the smeared occupation density
        calculated with smearing_kind) calculated with XC_method.
        ''',
        categories=[energy_component_per_atom, energy_value],
        a_legacy=LegacyDefinition(name='energy_free_per_atom'))

    energy_free = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Free energy (nuclei + electrons) (whose minimum gives the smeared occupation
        density calculated with smearing_kind) calculated with the method described in
        XC_method.
        ''',
        categories=[energy_component, energy_value, energy_total_potential],
        a_legacy=LegacyDefinition(name='energy_free'))

    energy_hartree_error = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Error in the Hartree (electrostatic) potential energy. Defined consistently with
        XC_method.
        ''',
        categories=[error_estimate_contribution, energy_value],
        a_legacy=LegacyDefinition(name='energy_hartree_error'))

    energy_hartree_fock_X_scaled = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Scaled exact-exchange energy that depends on the mixing parameter of the
        functional. For example in hybrid functionals, the exchange energy is given as a
        linear combination of exact-energy and exchange energy of an approximate DFT
        functional; the exact exchange energy multiplied by the mixing coefficient of the
        hybrid functional would be stored in this metadata. Defined consistently with
        XC_method.
        ''',
        categories=[energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_hartree_fock_X_scaled'))

    energy_hartree_fock_X = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Converged exact-exchange (Hartree-Fock) energy. Defined consistently with
        XC_method.
        ''',
        categories=[energy_type_X, energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_hartree_fock_X'))

    energy_method_current = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Value of the energy calculated with the method calculation_method_current.
        Depending on calculation_method_kind it might be a total energy or only a
        correction.
        ''',
        categories=[energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_method_current'))

    energy_sum_eigenvalues_per_atom = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Value of the energy per atom, where the energy is defined as the sum of the
        eigenvalues of the Hamiltonian matrix given by XC_method.
        ''',
        categories=[energy_component_per_atom, energy_value],
        a_legacy=LegacyDefinition(name='energy_sum_eigenvalues_per_atom'))

    energy_sum_eigenvalues = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Sum of the eigenvalues of the Hamiltonian matrix defined by XC_method.
        ''',
        categories=[energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_sum_eigenvalues'))

    energy_T0_per_atom = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Value of the total energy per atom, calculated with the method described in
        XC_method and extrapolated to $T=0$, based on a free-electron gas argument.
        ''',
        categories=[energy_total_potential_per_atom, energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_T0_per_atom'))

    energy_total_T0_per_atom = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Value of the total energy, calculated with the method described in XC_method per
        atom extrapolated to $T=0$, based on a free-electron gas argument.
        ''',
        categories=[energy_total_potential_per_atom, energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_total_T0_per_atom'))

    energy_total_T0 = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Value of the total energy (or equivalently free energy), calculated with the
        method described in XC_method and extrapolated to $T=0$, based on a free-electron
        gas argument.
        ''',
        categories=[energy_component, energy_value, energy_total_potential],
        a_legacy=LegacyDefinition(name='energy_total_T0'))

    energy_total = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Value of the total energy, calculated with the method described in XC_method and
        extrapolated to $T=0$, based on a free-electron gas argument.
        ''',
        categories=[energy_component, energy_value, energy_total_potential],
        a_legacy=LegacyDefinition(name='energy_total'))

    energy_XC_functional = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Value of the exchange-correlation (XC) energy calculated with the functional
        stored in XC_functional.
        ''',
        categories=[energy_type_XC, energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_XC_functional'))

    energy_XC_potential = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Value of the exchange-correlation (XC) potential energy: the integral of the first
        order derivative of the functional stored in XC_functional (integral of
        v_xc*electron_density), i.e., the component of XC that is in the sum of the
        eigenvalues. Value associated with the configuration, should be the most converged
        value.
        ''',
        categories=[energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_XC_potential'))

    energy_XC = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Value of the exchange-correlation (XC) energy calculated with the method described
        in XC_method.
        ''',
        categories=[energy_type_XC, energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_XC'))

    energy_X = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Value fo the exchange (X) energy calculated with the method described in
        XC_method.
        ''',
        categories=[energy_type_X, energy_component, energy_value],
        a_legacy=LegacyDefinition(name='energy_X'))

    energy_zero_point = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Value for the converged zero-point vibrations energy calculated using the method
        described in zero_point_method , and used in energy_current .
        ''',
        a_legacy=LegacyDefinition(name='energy_zero_point'))

    hessian_matrix = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 'number_of_atoms', 3, 3],
        description='''
        The matrix with the second derivative with respect to atom displacements.
        ''',
        a_legacy=LegacyDefinition(name='hessian_matrix'))

    message_debug_evaluation = Quantity(
        type=str,
        shape=[],
        description='''
        A debugging message of the computational program, associated with a *single
        configuration calculation* (see section_single_configuration_calculation).
        ''',
        categories=[message_debug],
        a_legacy=LegacyDefinition(name='message_debug_evaluation'))

    message_error_evaluation = Quantity(
        type=str,
        shape=[],
        description='''
        An error message of the computational program, associated with a *single
        configuration calculation* (see section_single_configuration_calculation).
        ''',
        categories=[message_info, message_debug, message_error, message_warning],
        a_legacy=LegacyDefinition(name='message_error_evaluation'))

    message_info_evaluation = Quantity(
        type=str,
        shape=[],
        description='''
        An information message of the computational program, associated with a *single
        configuration calculation* (see section_single_configuration_calculation).
        ''',
        categories=[message_info, message_debug],
        a_legacy=LegacyDefinition(name='message_info_evaluation'))

    message_warning_evaluation = Quantity(
        type=str,
        shape=[],
        description='''
        A warning message of the computational program.
        ''',
        categories=[message_info, message_debug, message_warning],
        a_legacy=LegacyDefinition(name='message_warning_evaluation'))

    number_of_scf_iterations = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of performed self-consistent field (SCF) iterations at a specfied
        level of theory.
        ''',
        categories=[scf_info],
        a_legacy=LegacyDefinition(name='number_of_scf_iterations'))

    parsing_message_debug_evaluation = Quantity(
        type=str,
        shape=[],
        description='''
        This field is used for debugging messages of the parsing program associated with a
        run, see section_run.
        ''',
        categories=[parsing_message_debug],
        a_legacy=LegacyDefinition(name='parsing_message_debug_evaluation'))

    parsing_message_error_single_configuration = Quantity(
        type=str,
        shape=[],
        description='''
        This field is used for error messages of the parsing program associated with a
        single configuration calculation, see section_single_configuration_calculation.
        ''',
        categories=[parsing_message_info, parsing_message_error, parsing_message_warning, parsing_message_debug],
        a_legacy=LegacyDefinition(name='parsing_message_error_single_configuration'))

    parsing_message_info_single_configuration = Quantity(
        type=str,
        shape=[],
        description='''
        This field is used for info messages of the parsing program associated with a
        single configuration calculation, see section_single_configuration_calculation.
        ''',
        categories=[parsing_message_info, parsing_message_debug],
        a_legacy=LegacyDefinition(name='parsing_message_info_single_configuration'))

    parsing_message_warning_evaluation = Quantity(
        type=str,
        shape=[],
        description='''
        This field is used for warning messages of the parsing program associated with a
        run, see section_run.
        ''',
        categories=[parsing_message_info, parsing_message_warning, parsing_message_debug],
        a_legacy=LegacyDefinition(name='parsing_message_warning_evaluation'))

    single_configuration_calculation_converged = Quantity(
        type=bool,
        shape=[],
        description='''
        Determines whether a *single configuration calculation* in
        section_single_configuration_calculation is converged.
        ''',
        a_legacy=LegacyDefinition(name='single_configuration_calculation_converged'))

    single_configuration_calculation_to_system_ref = Quantity(
        type=Reference(SectionProxy('section_system')),
        shape=[],
        description='''
        Reference to the system (atomic configuration, cell, ...) that is calculated in
        section_single_configuration_calculation.
        ''',
        a_legacy=LegacyDefinition(name='single_configuration_calculation_to_system_ref'))

    single_configuration_to_calculation_method_ref = Quantity(
        type=Reference(SectionProxy('section_method')),
        shape=[],
        description='''
        Reference to the method used for the calculation in
        section_single_configuration_calculation.
        ''',
        a_legacy=LegacyDefinition(name='single_configuration_to_calculation_method_ref'))

    spin_S2 = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Stores the value of the total spin moment operator $S^2$ for the converged
        wavefunctions calculated with the XC_method. It can be used to calculate the spin
        contamination in spin-unrestricted calculations.
        ''',
        a_legacy=LegacyDefinition(name='spin_S2'))

    stress_tensor = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        unit='pascal',
        description='''
        Stores the final value of the default stress tensor consistent with energy_total
        and calculated with the method specified in stress_tensor_method.

        This value is used (if needed) for, e.g., molecular dynamics and geometry
        optimization. Alternative definitions of the stress tensor can be assigned with
        stress_tensor_kind
        ''',
        categories=[stress_tensor_type],
        a_legacy=LegacyDefinition(name='stress_tensor'))

    time_calculation = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        Stores the wall-clock time needed for a calculation using
        calculation_method_current. Basically, it tracks the real time that has been
        elapsed from start to end.
        ''',
        categories=[time_info, accessory_info],
        a_legacy=LegacyDefinition(name='time_calculation'))

    time_single_configuration_calculation_cpu1_end = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        Stores the end time of the *single configuration calculation* (see
        section_single_configuration_calculation) on CPU 1.
        ''',
        categories=[time_info, accessory_info],
        a_legacy=LegacyDefinition(name='time_single_configuration_calculation_cpu1_end'))

    time_single_configuration_calculation_cpu1_start = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        Stores the start time of the *single configuration calculation* (see
        section_single_configuration_calculation) on CPU 1.
        ''',
        categories=[time_info, accessory_info],
        a_legacy=LegacyDefinition(name='time_single_configuration_calculation_cpu1_start'))

    time_single_configuration_calculation_date_end = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        Stores the end date of the *single configuration calculation* (see
        section_single_configuration_calculation) as time since the *Unix epoch* (00:00:00
        UTC on 1 January 1970) in seconds. For date and times without a timezone, the
        default timezone GMT is used.
        ''',
        categories=[time_info, accessory_info],
        a_legacy=LegacyDefinition(name='time_single_configuration_calculation_date_end'))

    time_single_configuration_calculation_date_start = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        Stores the start date of the *single configuration calculation* (see
        section_single_configuration_calculation) as time since the *Unix epoch* (00:00:00
        UTC on 1 January 1970) in seconds. For date and times without a timezone, the
        default timezone GMT is used.
        ''',
        categories=[time_info, accessory_info],
        a_legacy=LegacyDefinition(name='time_single_configuration_calculation_date_start'))

    time_single_configuration_calculation_wall_end = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        Stores the internal wall-clock time at the end of the *single configuration
        calculation* (see section_single_configuration_calculation).
        ''',
        categories=[time_info, accessory_info],
        a_legacy=LegacyDefinition(name='time_single_configuration_calculation_wall_end'))

    time_single_configuration_calculation_wall_start = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        Stores the internal wall-clock time from the start of the *single configuration
        calculation* (see section_single_configuration_calculation).
        ''',
        categories=[time_info, accessory_info],
        a_legacy=LegacyDefinition(name='time_single_configuration_calculation_wall_start'))

    zero_point_method = Quantity(
        type=str,
        shape=[],
        description='''
        Describes the zero-point vibrations method. If skipped or an empty string is used,
        it means no zero-point vibrations correction is applied.
        ''',
        a_legacy=LegacyDefinition(name='zero_point_method'))

    section_atom_projected_dos = SubSection(
        sub_section=SectionProxy('section_atom_projected_dos'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_atom_projected_dos'))

    section_atomic_multipoles = SubSection(
        sub_section=SectionProxy('section_atomic_multipoles'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_atomic_multipoles'))

    section_basis_set = SubSection(
        sub_section=SectionProxy('section_basis_set'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_basis_set'))

    section_calculation_to_calculation_refs = SubSection(
        sub_section=SectionProxy('section_calculation_to_calculation_refs'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_calculation_to_calculation_refs'))

    section_calculation_to_folder_refs = SubSection(
        sub_section=SectionProxy('section_calculation_to_folder_refs'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_calculation_to_folder_refs'))

    section_dos = SubSection(
        sub_section=SectionProxy('section_dos'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_dos'))

    section_eigenvalues = SubSection(
        sub_section=SectionProxy('section_eigenvalues'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_eigenvalues'))

    section_energy_code_independent = SubSection(
        sub_section=SectionProxy('section_energy_code_independent'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_energy_code_independent'))

    section_energy_van_der_Waals = SubSection(
        sub_section=SectionProxy('section_energy_van_der_Waals'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_energy_van_der_Waals'))

    section_k_band_normalized = SubSection(
        sub_section=SectionProxy('section_k_band_normalized'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_k_band_normalized'))

    section_k_band = SubSection(
        sub_section=SectionProxy('section_k_band'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_k_band'))

    section_scf_iteration = SubSection(
        sub_section=SectionProxy('section_scf_iteration'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_scf_iteration'))

    section_species_projected_dos = SubSection(
        sub_section=SectionProxy('section_species_projected_dos'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_species_projected_dos'))

    section_stress_tensor = SubSection(
        sub_section=SectionProxy('section_stress_tensor'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_stress_tensor'))

    section_volumetric_data = SubSection(
        sub_section=SectionProxy('section_volumetric_data'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_volumetric_data'))


class section_species_projected_dos(MSection):
    '''
    Section collecting the information on a species-projected density of states (DOS)
    evaluation.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_species_projected_dos'))

    number_of_lm_species_projected_dos = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of $l$, $m$ combinations for the species-projected density of
        states (DOS) defined in section_species_projected_dos.
        ''',
        a_legacy=LegacyDefinition(name='number_of_lm_species_projected_dos'))

    number_of_species_projected_dos_values = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of energy values for the species-projected density of states
        (DOS) defined in section_species_projected_dos.
        ''',
        a_legacy=LegacyDefinition(name='number_of_species_projected_dos_values'))

    number_of_species = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of species for the species-projected density of states (DOS)
        defined in section_species_projected_dos.
        ''',
        a_legacy=LegacyDefinition(name='number_of_species'))

    species_projected_dos_energies_normalized = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_species_projected_dos_values'],
        unit='joule',
        description='''
        Contains the set of discrete energy values with respect to the top of the valence
        band for the species-projected density of states (DOS). It is derived from the
        species_projected_dos_energies species field.
        ''',
        a_legacy=LegacyDefinition(name='species_projected_dos_energies_normalized'))

    species_projected_dos_energies = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_species_projected_dos_values'],
        unit='joule',
        description='''
        Contains the set of discrete energy values for the species-projected density of
        states (DOS).
        ''',
        a_legacy=LegacyDefinition(name='species_projected_dos_energies'))

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
        ''',
        a_legacy=LegacyDefinition(name='species_projected_dos_lm'))

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
        ''',
        a_legacy=LegacyDefinition(name='species_projected_dos_m_kind'))

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
        ''',
        a_legacy=LegacyDefinition(name='species_projected_dos_species_label'))

    species_projected_dos_values_lm = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_lm_species_projected_dos', 'number_of_spin_channels', 'number_of_species', 'number_of_species_projected_dos_values'],
        description='''
        Holds species-projected density of states (DOS) values, divided into contributions
        from each $l,m$ channel.

        Here, there are as many species-projected DOS as the number of species,
        number_of_species. The list of labels of the species is given in
        species_projected_dos_species_label.
        ''',
        a_legacy=LegacyDefinition(name='species_projected_dos_values_lm'))

    species_projected_dos_values_total = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_spin_channels', 'number_of_species', 'number_of_species_projected_dos_values'],
        description='''
        Holds species-projected density of states (DOS) values, summed up over all
        azimuthal quantum numbers $l$.

        Here, there are as many species-projected DOS as the number of species,
        number_of_species. The list of labels of the species is given in
        species_projected_dos_species_label.
        ''',
        a_legacy=LegacyDefinition(name='species_projected_dos_values_total'))


class section_springer_material(MSection):
    '''
    Every section_springer_material contains results of classification of materials with
    the same formula according to Springer Materials - it contains
    section_springer_classsification, section_springer_compound, section_springer_id,
    section_springer_references
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_springer_material'))

    springer_id = Quantity(
        type=str,
        shape=[],
        description='''
        Id of the classified material according to Springer Materials
        ''',
        a_legacy=LegacyDefinition(name='springer_id'))

    springer_alphabetical_formula = Quantity(
        type=str,
        shape=[],
        description='''
        The alphabetical formula of the material according to Springer Materials Database
        ''',
        a_legacy=LegacyDefinition(name='springer_alphabetical_formula'))

    springer_url = Quantity(
        type=str,
        shape=[],
        description='''
        Url to the source page in Springer Materials describing the current entry
        ''',
        a_legacy=LegacyDefinition(name='springer_url'))

    springer_compound_class = Quantity(
        type=str,
        shape=['N'],
        description='''
        Name of a class of the current compound, as defined in by Springer Materials. This
        is a property of the chemical formula of the compound
        ''',
        a_legacy=LegacyDefinition(name='springer_compound_class'))

    springer_classification = Quantity(
        type=str,
        shape=['N'],
        description='''
        Contains the classification name of the current material according to Springer
        Materials
        ''',
        a_legacy=LegacyDefinition(name='springer_classification'))

    section_springer_id = SubSection(
        sub_section=SectionProxy('section_springer_id'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_springer_id'))


class section_springer_id(MSection):
    '''
    Identifiers used by Springer Materials
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_springer_id'))


class section_std_system(MSection):
    '''
    Section containing symmetry information that is specific to the standardized system.
    The standardized system is defined as given by spglib and the details can be found
    from https://arxiv.org/abs/1506.01455
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_std_system'))

    atom_positions_std = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms_std', 3],
        description='''
        Standardized atom positions in reduced units.
        ''',
        a_legacy=LegacyDefinition(name='atom_positions_std'))

    atomic_numbers_std = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_atoms_std'],
        description='''
        Atomic numbers of the atoms in the standardized cell.
        ''',
        a_legacy=LegacyDefinition(name='atomic_numbers_std'))

    equivalent_atoms_std = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_atoms_std'],
        description='''
        Gives a mapping table of atoms to symmetrically independent atoms in the
        standardized cell. This is used to find symmetrically equivalent atoms.
        ''',
        a_legacy=LegacyDefinition(name='equivalent_atoms_std'))

    lattice_vectors_std = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        unit='meter',
        description='''
        Standardized lattice vectors of the conventional cell. The vectors are the rows of
        this matrix.
        ''',
        a_legacy=LegacyDefinition(name='lattice_vectors_std'))

    number_of_atoms_std = Quantity(
        type=int,
        shape=[],
        description='''
        Number of atoms in standardized system.
        ''',
        a_legacy=LegacyDefinition(name='number_of_atoms_std'))

    wyckoff_letters_std = Quantity(
        type=str,
        shape=['number_of_atoms_std'],
        description='''
        Wyckoff letters for atoms in the standardized cell.
        ''',
        a_legacy=LegacyDefinition(name='wyckoff_letters_std'))


class section_stress_tensor(MSection):
    '''
    Section collecting alternative values to stress_tensor that have been calculated.

    This section allows the storage of multiple definitions and evaluated values of the
    stress tensor, while only one definition is used for, e.g., molecular dynamics or
    geometry optimization (if needed).
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_stress_tensor'))

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
        ''',
        a_legacy=LegacyDefinition(name='stress_tensor_kind'))

    stress_tensor_value = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        unit='pascal',
        description='''
        Contains the value of the stress tensor of the kind defined in stress_tensor_kind.
        This is an *alternative* to the stress tensor defined in stress_tensor_method.

        This field allows for multiple definitions and evaluated values of the stress
        tensor, while only one definition is used for, e.g., molecular dynamics and
        geometry optimization.
        ''',
        categories=[stress_tensor_type],
        a_legacy=LegacyDefinition(name='stress_tensor_value'))


class section_symmetry(MSection):
    '''
    Section containing information about the symmetry properties of the system.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_symmetry'))

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
        ''',
        a_legacy=LegacyDefinition(name='bravais_lattice'))

    choice = Quantity(
        type=str,
        shape=[],
        description='''
        String that specifies the centering, origin and basis vector settings of the 3D
        space group that defines the symmetry group of the simulated physical system (see
        section_system). Values are as defined by spglib.
        ''',
        a_legacy=LegacyDefinition(name='choice'))

    crystal_system = Quantity(
        type=str,
        shape=[],
        description='''
        Name of the crystal system. Can be one of the following: triclinic, monoclinic,
        orthorhombic, tetragonal, trigonal, hexagonal or cubic.
        ''',
        a_legacy=LegacyDefinition(name='crystal_system'))

    hall_number = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        The Hall number for this system.
        ''',
        a_legacy=LegacyDefinition(name='hall_number'))

    hall_symbol = Quantity(
        type=str,
        shape=[],
        description='''
        The Hall symbol for this system.
        ''',
        a_legacy=LegacyDefinition(name='hall_symbol'))

    international_short_symbol = Quantity(
        type=str,
        shape=[],
        description='''
        Specifies the International Union of Crystallography (IUC) short symbol of the 3D
        space group of this system
        ''',
        a_legacy=LegacyDefinition(name='international_short_symbol'))

    origin_shift = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        description='''
        Vector $\\mathbf{p}$ from the origin of the standardized system to the origin of
        the original system. Together with the matrix $\\mathbf{P}$, found in
        space_group_3D_transformation_matrix, the transformation between the standardized
        coordinates $\\mathbf{x}_s$ and original coordinates $\\mathbf{x}$ is then given by
        $\\mathbf{x}_s = \\mathbf{P} \\mathbf{x} + \\mathbf{p}$.
        ''',
        a_legacy=LegacyDefinition(name='origin_shift'))

    point_group = Quantity(
        type=str,
        shape=[],
        description='''
        Symbol of the crystallographic point group in the Hermann-Mauguin notation.
        ''',
        a_legacy=LegacyDefinition(name='point_group'))

    space_group_number = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Specifies the International Union of Crystallography (IUC) number of the 3D space
        group of this system.
        ''',
        a_legacy=LegacyDefinition(name='space_group_number'))

    symmetry_method = Quantity(
        type=str,
        shape=[],
        description='''
        Identifies the source of the symmetry information contained within this section.
        If equal to 'spg_normalized' the information comes from a normalization step.
        ''',
        a_legacy=LegacyDefinition(name='symmetry_method'))

    transformation_matrix = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        description='''
        Matrix $\\mathbf{P}$ that is used to transform the standardized coordinates to the
        original coordinates. Together with the vector $\\mathbf{p}$, found in
        space_group_3D_origin_shift, the transformation between the standardized
        coordinates $\\mathbf{x}_s$ and original coordinates $\\mathbf{x}$ is then given by
        $\\mathbf{x}_s = \\mathbf{P} \\mathbf{x} + \\mathbf{p}$.
        ''',
        a_legacy=LegacyDefinition(name='transformation_matrix'))

    section_original_system = SubSection(
        sub_section=SectionProxy('section_original_system'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_original_system'))

    section_primitive_system = SubSection(
        sub_section=SectionProxy('section_primitive_system'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_primitive_system'))

    section_std_system = SubSection(
        sub_section=SectionProxy('section_std_system'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_std_system'))


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

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_system_to_system_refs'))

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
        ''',
        a_legacy=LegacyDefinition(name='system_to_system_kind'))

    system_to_system_ref = Quantity(
        type=Reference(SectionProxy('section_system')),
        shape=[],
        description='''
        Reference to another system. The kind of relationship between the present and the
        referenced section_system is specified by system_to_system_kind.
        ''',
        a_legacy=LegacyDefinition(name='system_to_system_ref'))


class section_system(MSection):
    '''
    Every section_system contains all needed properties required to describe the simulated
    physical system, e.g. the given atomic configuration, the definition of periodic cell
    (if present), the external potentials and other parameters.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_system'))

    atom_atom_number = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_sites'],
        description='''
        Atomic number Z of the atom.
        ''',
        a_legacy=LegacyDefinition(name='atom_atom_number'))

    atom_concentrations = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms'],
        description='''
        concentration of the atom species in a variable composition, by default it should
        be considered an array of ones. Summing these should give the number_of_sites
        ''',
        a_legacy=LegacyDefinition(name='atom_concentrations'))

    atom_labels = Quantity(
        type=str,
        shape=['number_of_atoms'],
        description='''
        Labels of the atoms. These strings identify the atom kind and conventionally start
        with the symbol of the atomic species, possibly followed by the atomic number. The
        same atomic species can be labeled with more than one atom_labels in order to
        distinguish, e.g., atoms of the same species assigned to different atom-centered
        basis sets or pseudo-potentials, or simply atoms in different locations in the
        structure (e.g., bulk and surface). These labels can also be used for *particles*
        that do not correspond to physical atoms (e.g., ghost atoms in some codes using
        atom-centered basis sets). This metadata defines a configuration and is therefore
        required.
        ''',
        categories=[configuration_core],
        a_legacy=LegacyDefinition(name='atom_labels'))

    atom_positions = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        unit='meter',
        description='''
        Positions of all the atoms, in Cartesian coordinates. This metadata defines a
        configuration and is therefore required. For alloys where concentrations of
        species are given for each site in the unit cell, it stores the position of the
        sites.
        ''',
        categories=[configuration_core],
        a_legacy=LegacyDefinition(name='atom_positions'))

    atom_species = Quantity(
        type=np.dtype(np.int32),
        shape=['number_of_atoms'],
        description='''
        Species of the atom (normally the atomic number Z, 0 or negative for unidentifed
        species or particles that are not atoms.
        ''',
        a_legacy=LegacyDefinition(name='atom_species'))

    atom_velocities = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        unit='meter / second',
        description='''
        Velocities of the nuclei, defined as the change in Cartesian coordinates of the
        nuclei with respect to time.
        ''',
        a_legacy=LegacyDefinition(name='atom_velocities'))

    configuration_periodic_dimensions = Quantity(
        type=bool,
        shape=[3],
        description='''
        Array labeling which of the lattice vectors use periodic boundary conditions. Note
        for the parser developers: This value is not expected to be given for each
        section_single_configuration_calculation. It is assumed to be valid from the
        section_single_configuration_calculation where it is defined for all subsequent
        section_single_configuration_calculation in section_run, until redefined.
        ''',
        categories=[configuration_core],
        a_legacy=LegacyDefinition(name='configuration_periodic_dimensions'))

    configuration_raw_gid = Quantity(
        type=str,
        shape=[],
        description='''
        checksum of the configuration_core, i.e. the geometry of the system. The values
        are not normalized in any way so equivalent configurations might have different
        values
        ''',
        a_legacy=LegacyDefinition(name='configuration_raw_gid'))

    embedded_system = Quantity(
        type=bool,
        shape=[],
        description='''
        Is the system embedded into a host geometry?.
        ''',
        categories=[configuration_core],
        a_legacy=LegacyDefinition(name='embedded_system'))

    lattice_vectors = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        unit='meter',
        description='''
        Holds the lattice vectors (in Cartesian coordinates) of the simulation cell. The
        last (fastest) index runs over the $x,y,z$ Cartesian coordinates, and the first
        index runs over the 3 lattice vectors.
        ''',
        categories=[configuration_core],
        a_legacy=LegacyDefinition(name='lattice_vectors'))

    local_rotations = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3, 3],
        description='''
        A rotation matrix defining the orientation of each atom. If the rotation matrix
        only needs to be specified for some atoms, the remaining atoms should set it to
        the zero matrix (not the identity!)
        ''',
        a_legacy=LegacyDefinition(name='local_rotations'))

    number_of_atoms = Quantity(
        type=int,
        shape=[],
        description='''
        Stores the total number of atoms used in the calculation. For alloys where
        concentrations of species are given for each site in the unit cell, it stores the
        number of sites.
        ''',
        a_legacy=LegacyDefinition(name='number_of_atoms'))

    number_of_sites = Quantity(
        type=int,
        shape=[],
        description='''
        number of sites in a variable composition representation. By default (no variable
        composition) it is the same as number_of_atoms.
        ''',
        a_legacy=LegacyDefinition(name='number_of_sites'))

    SC_matrix = Quantity(
        type=np.dtype(np.int32),
        shape=[3, 3],
        description='''
        Specifies the matrix that transforms the unit-cell into the super-cell in which
        the actual calculation is performed.
        ''',
        a_legacy=LegacyDefinition(name='SC_matrix'))

    simulation_cell = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        unit='meter',
        description='''
        DEPRECATED, use lattice_vectors instead. Holds the lattice vectors (in Cartesian
        coordinates) of the simulation cell. The last (fastest) index runs over the
        $x,y,z$ Cartesian coordinates, and the first index runs over the 3 lattice
        vectors.
        ''',
        categories=[configuration_core],
        a_legacy=LegacyDefinition(name='simulation_cell'))

    symmorphic = Quantity(
        type=bool,
        shape=[],
        description='''
        Is the space group symmorphic? Set to True if all translations are zero.
        ''',
        a_legacy=LegacyDefinition(name='symmorphic'))

    system_composition = Quantity(
        type=str,
        shape=[],
        description='''
        Composition, i.e. cumulative chemical formula with atoms ordered by decreasing
        atomic number Z.
        ''',
        a_legacy=LegacyDefinition(name='system_composition'))

    system_configuration_consistent = Quantity(
        type=bool,
        shape=[],
        description='''
        Flag set is the configuration is consistent
        ''',
        a_legacy=LegacyDefinition(name='system_configuration_consistent'))

    system_name = Quantity(
        type=str,
        shape=[],
        description='''
        Specifies the name of the system. This information is provided by the user in some
        codes and is stored here for debugging or visualization purposes.
        ''',
        a_legacy=LegacyDefinition(name='system_name'))

    system_reweighted_composition = Quantity(
        type=str,
        shape=[],
        description='''
        Composition, i.e. cumulative chemical with atoms ordered by decreasing atomic
        number Z reweighted so that the sum is close to 100, and values are rounded up,
        and are stable (i.e. it is a fixed point).
        ''',
        a_legacy=LegacyDefinition(name='system_reweighted_composition'))

    system_type = Quantity(
        type=str,
        shape=[],
        description='''
        Type of the system
        ''',
        a_legacy=LegacyDefinition(name='system_type'))

    time_reversal_symmetry = Quantity(
        type=bool,
        shape=[],
        description='''
        Is time-reversal symmetry present?
        ''',
        a_legacy=LegacyDefinition(name='time_reversal_symmetry'))

    chemical_composition = Quantity(
        type=str,
        shape=[],
        description='''
        The chemical composition as full formula of the system, based on atom species.
        ''',
        a_legacy=LegacyDefinition(name='chemical_composition'))

    chemical_composition_reduced = Quantity(
        type=str,
        shape=[],
        description='''
        The chemical composition as reduced formula of the system, based on atom species.
        ''',
        a_legacy=LegacyDefinition(name='chemical_composition_reduced'))

    chemical_composition_bulk_reduced = Quantity(
        type=str,
        shape=[],
        description='''
        The chemical composition as reduced bulk formula of the system, based on atom
        species.
        ''',
        a_legacy=LegacyDefinition(name='chemical_composition_bulk_reduced'))

    section_prototype = SubSection(
        sub_section=SectionProxy('section_prototype'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_prototype'))

    section_springer_material = SubSection(
        sub_section=SectionProxy('section_springer_material'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_springer_material'))

    section_symmetry = SubSection(
        sub_section=SectionProxy('section_symmetry'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_symmetry'))

    section_system_to_system_refs = SubSection(
        sub_section=SectionProxy('section_system_to_system_refs'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_system_to_system_refs'))


class section_thermodynamical_properties(MSection):
    '''
    Section that defines thermodynamical properties about the system in a
    section_frame_sequence.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_thermodynamical_properties'))

    helmholz_free_energy = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_thermodynamical_property_values'],
        unit='joule',
        description='''
        Stores the Helmholtz free energy per unit cell at constant volume of a
        thermodynamic calculation.
        ''',
        a_legacy=LegacyDefinition(name='helmholz_free_energy'))

    number_of_thermodynamical_property_values = Quantity(
        type=int,
        shape=[],
        description='''
        Gives the number of thermal properties values available in
        section_thermodynamical_properties.
        ''',
        a_legacy=LegacyDefinition(name='number_of_thermodynamical_property_values'))

    thermodynamical_properties_calculation_method = Quantity(
        type=str,
        shape=[],
        description='''
        Method used to calculate the thermodynamic quantities.

        Valid values:

        * harmonic
        ''',
        a_legacy=LegacyDefinition(name='thermodynamical_properties_calculation_method'))

    thermodynamical_property_heat_capacity_C_v = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_thermodynamical_property_values'],
        unit='joule / kelvin',
        description='''
        Stores the heat capacity per cell unit at constant volume.
        ''',
        a_legacy=LegacyDefinition(name='thermodynamical_property_heat_capacity_C_v'))

    @derived(
        type=np.dtype(np.float64),
        shape=['number_of_thermodynamical_property_values'],
        unit='joule / kelvin * kilogram',
        description='''
        Stores the specific heat capacity at constant volume.
        ''',
        a_legacy=LegacyDefinition(name='specific_heat_capacity'),
        cached=True
    )
    def specific_heat_capacity(self) -> np.array:
        """Returns the specific heat capacity by dividing the heat capacity per
        cell with the mass of the atoms in the cell.
        """
        import nomad.atomutils
        s_frame_sequence = self.m_parent
        first_frame = s_frame_sequence.frame_sequence_local_frames_ref[0]
        system = first_frame.single_configuration_calculation_to_system_ref
        atomic_numbers = system.atom_species
        mass_per_unit_cell = nomad.atomutils.get_summed_atomic_mass(atomic_numbers)
        heat_capacity = self.thermodynamical_property_heat_capacity_C_v
        specific_heat_capacity = heat_capacity / mass_per_unit_cell

        return specific_heat_capacity

    thermodynamical_property_temperature = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_thermodynamical_property_values'],
        unit='kelvin',
        description='''
        Specifies the temperatures at which properties such as the Helmholtz free energy
        are calculated.
        ''',
        a_legacy=LegacyDefinition(name='thermodynamical_property_temperature'))

    vibrational_free_energy_at_constant_volume = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_thermodynamical_property_values'],
        unit='joule',
        description='''
        Holds the vibrational free energy per atom at constant volume.
        ''',
        a_legacy=LegacyDefinition(name='vibrational_free_energy_at_constant_volume'))

    @derived(
        type=np.dtype(np.float64),
        shape=['number_of_thermodynamical_property_values'],
        unit='joule / kilogram',
        description='''
        Stores the specific vibrational free energy at constant volume.
        ''',
        a_legacy=LegacyDefinition(name='specific_vibrational_free_energy_at_constant_volume'),
        cached=True
    )
    def specific_vibrational_free_energy_at_constant_volume(self) -> np.array:
        """Returns the specific vibrational free energy by dividing the vibrational free energy per
        cell with the mass of the atoms in the cell.
        """
        import nomad.atomutils
        s_frame_sequence = self.m_parent
        first_frame = s_frame_sequence.frame_sequence_local_frames_ref[0]
        system = first_frame.single_configuration_calculation_to_system_ref
        atomic_numbers = system.atom_species
        n_atoms = len(atomic_numbers)
        mass_per_atom = nomad.atomutils.get_summed_atomic_mass(atomic_numbers) / n_atoms
        free_energy = self.vibrational_free_energy_at_constant_volume
        specific_vibrational_free_energy_at_constant_volume = free_energy / mass_per_atom

        return specific_vibrational_free_energy_at_constant_volume


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

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_volumetric_data'))

    volumetric_data_displacements = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        unit='meter',
        description='''
        displacement vectors between grid points along each axis; same indexing rules as
        lattice_vectors.  In many cases, displacements and number of points are related to
        lattice_vectors through: [displacement] * [number of points + N] =
        [lattice_vector],where N is 1 for periodic directions and 0 for non-periodic ones
        ''',
        a_legacy=LegacyDefinition(name='volumetric_data_displacements'))

    volumetric_data_kind = Quantity(
        type=str,
        shape=[],
        description='''
        The kind of function, e.g. density, potential_hartree, potential_effective.  The
        unit of measurement for "volumetric_data_values" depends on the kind: Densities
        are 1/m^3 and potentials are J/m^3.  See [full specification on the
        wiki](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-
        info/wikis/metainfo/volumetric-data).
        ''',
        a_legacy=LegacyDefinition(name='volumetric_data_kind'))

    volumetric_data_multiplicity = Quantity(
        type=int,
        shape=[],
        description='''
        number of functions stored
        ''',
        a_legacy=LegacyDefinition(name='volumetric_data_multiplicity'))

    volumetric_data_nx = Quantity(
        type=int,
        shape=[],
        description='''
        number of points along x axis
        ''',
        a_legacy=LegacyDefinition(name='volumetric_data_nx'))

    volumetric_data_ny = Quantity(
        type=int,
        shape=[],
        description='''
        number of points along y axis
        ''',
        a_legacy=LegacyDefinition(name='volumetric_data_ny'))

    volumetric_data_nz = Quantity(
        type=int,
        shape=[],
        description='''
        number of points along z axis
        ''',
        a_legacy=LegacyDefinition(name='volumetric_data_nz'))

    volumetric_data_origin = Quantity(
        type=np.dtype(np.float64),
        shape=[3],
        description='''
        location of the first grid point; same coordinate system as atom_positions when
        applicable.
        ''',
        a_legacy=LegacyDefinition(name='volumetric_data_origin'))

    volumetric_data_values = Quantity(
        type=np.dtype(np.float64),
        shape=['volumetric_data_multiplicity', 'volumetric_data_nx', 'volumetric_data_ny', 'volumetric_data_nz'],
        description='''
        Array of shape (multiplicity, nx, ny, nz) containing the values.  The units of
        these values depend on which kind of data the values represent; see
        "volumetric_data_kind".
        ''',
        a_legacy=LegacyDefinition(name='volumetric_data_values'))


class section_XC_functionals(MSection):
    '''
    Section containing one of the exchange-correlation (XC) functionals for the present
    section_method that are combined to form the XC_functional.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_XC_functionals'))

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
        ''',
        categories=[settings_physical_parameter],
        a_legacy=LegacyDefinition(name='XC_functional_name'))

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
        ''',
        categories=[settings_physical_parameter],
        a_legacy=LegacyDefinition(name='XC_functional_parameters'))

    XC_functional_weight = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Provides the value of the weight for the exchange, correlation, or exchange-
        correlation functional declared in XC_functional_name (see
        section_XC_functionals).

        This weight is used in the linear combination of the different XC functional names
        (XC_functional_name) in different section_XC_functionals sections to form the
        XC_functional used for evaluating energy_XC_functional and related quantities.

        If not specified then the default is set to 1.
        ''',
        categories=[settings_physical_parameter],
        a_legacy=LegacyDefinition(name='XC_functional_weight'))


class Relaxation(MSection):
    '''
    Section containing the results of a relaxation workflow.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_relaxation'))

    relaxation_type = Quantity(
        type=str,
        shape=[],
        description='''
        The type of relaxation ionic, cell_shape, cell_volume.
        ''',
        a_legacy=LegacyDefinition(name='relaxation_type')
    )

    input_energy_difference_tolerance = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        The input energy difference tolerance criterion.
        ''',
        a_legacy=LegacyDefinition(name='input_energy_difference_tolerance'))

    input_force_maximum_tolerance = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='newton',
        description='''
        The input maximum net force tolerance criterion.
        ''',
        a_legacy=LegacyDefinition(name='input_force_maximum_tolerance'))

    final_energy_difference = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        The difference in the energy between the last two steps during relaxation.
        ''',
        a_legacy=LegacyDefinition(name='final_energy_difference'),
        a_search=Search())

    final_force_maximum = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='newton',
        description='''
        The maximum net force in the last relaxation step.
        ''',
        a_legacy=LegacyDefinition(name='final_force_maximum')
    )

    final_calculation_ref = Quantity(
        type=Reference(SectionProxy('section_single_configuration_calculation')),
        shape=[],
        description='''
        Reference to last calculation step.
        ''',
        a_legacy=LegacyDefinition(name='final_calculation_ref'))


class Workflow(MSection):
    '''
    Section containing the  results of a workflow.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_workflow'))

    workflow_type = Quantity(
        type=str,
        shape=[],
        description='''
        The type of calculation workflow. Can be one of relaxation, elastic, phonon,
        molecular dynamics.
        ''',
        a_legacy=LegacyDefinition(name='workflow_type'),
        a_search=Search())

    section_relaxation = SubSection(
        sub_section=SectionProxy('Relaxation'),
        a_legacy=LegacyDefinition(name='section_relaxation'))


m_package.__init_metainfo__()
