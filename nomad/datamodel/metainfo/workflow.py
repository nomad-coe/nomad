import numpy as np            # pylint: disable=unused-import
import typing                 # pylint: disable=unused-import
from nptyping import NDArray
from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference, MEnum, derived)
from nomad.metainfo.search_extension import Search
from nomad.datamodel.metainfo.run.calculation import Calculation
from nomad.datamodel.metainfo.run.run import RunReference


class FastAccess(MCategory):
    '''
    Used to mark archive objects that need to be stored in a fast 2nd-tier storage medium,
    because they are frequently accessed via archive API.

    If applied to a sub_section, the section will be added to the fast storage. Currently
    this only works for *root* sections that are sub_sections of `EntryArchive`.

    If applied to a reference types quantity, the referenced section will also be added to
    the fast storage, regardless if the referenced section has the category or not.
    '''

    m_def = Category(aliases=['fast_access'])


class Raman(MSection):
    '''
    Section containing results of a Raman workflow.
    '''

    m_def = Section(validate=False)

    n_modes = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of evaluated vibrational modes.
        ''')

    n_atoms = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of atoms in the simulation cell.
        ''')

    frequencies = Quantity(
        type=np.dtype(np.float64),
        shape=['n_modes'],
        unit='1 / m',
        description='''
        Calculated value of the Raman frequencies.
        ''')


class MagneticOrdering(MSection):
    '''
    Section containing results of a magnetic ordering workflow.
    '''

    m_def = Section(validate=False)

    n_structures = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of magnetic structures evaluated.
        ''')

    n_atoms = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of atoms in the simulation cell.
        ''')

    labels = Quantity(
        type=str,
        shape=['n_structures'],
        description='''
        Labels corresponding to each magnetic structure.
        ''')

    energies = Quantity(
        type=np.dtype(np.float64),
        shape=['n_structures'],
        unit='joule',
        description='''
        Calculated value of the energies corresponding to each magnetic structure.
        ''')

    magnetic_moments = Quantity(
        type=np.dtype(np.float64),
        shape=['n_structures', 'n_atoms'],
        unit='bohr_magneton',
        description='''
        Resulting atomic magnetic moments corresponding to each magnetic structure.
        ''')

    magnetic_deformations = Quantity(
        type=np.dtype(np.float64),
        shape=['n_structures'],
        unit='m',
        description='''
        Average atomic displacements after relaxation with respect to the non-magnetic
        case for each magnetic structure.
        ''')


class Adsorption(MSection):
    '''
    Section containing results of a surface adsorption workflow.
    '''

    m_def = Section(
        validate=False)

    n_sites = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of sites for which the adsorption energy is evaluated.
        ''')

    slab_miller_index = Quantity(
        type=np.dtype(np.int32),
        shape=[3],
        description='''
        Miller index of the slab.
        ''')

    slab = Quantity(
        type=str,
        shape=[],
        description='''
        Chemical formula of the slab.
        ''')

    adsorbate = Quantity(
        type=str,
        shape=[],
        description='''
        Chemical formula of the adsorbate molecule.
        ''')

    adsorption_sites = Quantity(
        type=np.dtype(np.float64),
        shape=['n_sites'],
        description='''
        Coordinates of the adsorption sites corresponding to a minimum energy.
        ''')

    adsorption_energies = Quantity(
        type=np.dtype(np.float64),
        shape=['n_sites'],
        unit='joule',
        description='''
        Calculated value of the adsorption energy corresponding to each site.
        ''')


class ConvexHull(MSection):
    '''
    Section containing results of a convex hull workflow.
    '''

    m_def = Section(
        validate=False)

    n_elements = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of elements for which the thermal stability is evaluated. This represents
        the dimensionality of the convex hull.
        ''')

    n_points = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of points for which the energies are evaluated.
        ''')

    compositions = Quantity(
        type=np.dtype(np.float64),
        shape=['n_points', 'n_elements'],
        description='''
        Normalized composition of the elements corresponding to each point for which the
        energies are evaluated.
        ''')

    references = Quantity(
        type=str,
        shape=['n_elements'],
        description='''
        Specifies the reference structure for each element.
        ''')

    heat_of_formation = Quantity(
        type=np.dtype(np.float64),
        shape=['n_points'],
        unit='joule',
        description='''
        Values of the heat of formation corresponding to each point.
        ''')

    energy_hulll = Quantity(
        type=np.dtype(np.float64),
        shape=['n_points'],
        unit='joule',
        description='''
        Values of the energy above the convex hull corresponding to each point.
        ''')


class NudgedElasticBand(MSection):
    '''
    Section containing results of a nudged-elastic band workflow.
    '''

    m_def = Section(
        validate=False)

    method = Quantity(
        type=str,
        shape=[],
        description='''
        Specifies the method used in calculating the minumum energy path. Can be one of
        standard, improved_tangeant, full_spring_force, spline_interpolation, string.
        ''')

    climbing_image = Quantity(
        type=bool,
        shape=[],
        description='''
        Indicates if climbing image is used.
        ''')

    solid_state = Quantity(
        type=bool,
        shape=[],
        description='''
        Indicates if solid state nudged-elastic band calculation is performed.
        ''')

    optimizer = Quantity(
        type=str,
        shape=[],
        description='''
        Specifies the method used in energy minimization.
        ''')

    n_images = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of images used in the calculation.
        ''')

    spring_constants = Quantity(
        type=np.dtype(np.float64),
        shape=['n_images'],
        unit='newton',
        description='''
        Spring constants for each spring.
        ''')

    energy_barrier = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Maximum value of the calculated energy barrier.
        ''')

    force_maximum = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='newton',
        description='''
        Maximum force along the minimum energy path.
        ''')


class EquationOfState(MSection):
    '''
    Section containing results of an equation of state workflow.
    '''

    m_def = Section(
        validate=False)

    fit_function = Quantity(
        type=str,
        shape=[],
        description='''
        Specifies the function used to perform the fitting of the volume-energy data. Value
        can be one of birch_euler, birch_lagrange, birch_murnaghan, mie_gruneisen,
        murnaghan, pack_evans_james, poirier_tarantola, tait, vinet.
        ''')

    n_points = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of volume-energy pairs in data.
        ''')

    volumes = Quantity(
        type=np.dtype(np.float64),
        shape=['n_points'],
        unit='m ** 3',
        description='''
        Array of volumes for which the energies are evaluated.
        ''')

    energies = Quantity(
        type=np.dtype(np.float64),
        shape=['n_points'],
        unit='joule',
        description='''
        Array of energies corresponding to each volume.
        ''')

    bulk_modulus = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='pascal',
        description='''
        Calculated value of the bulk modulus by fitting the volume-energy data.
        ''')

    bulk_modulus_derivative = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Calculated value of the pressure derivative of the bulk modulus.
        ''')

    equilibrium_volume = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='m ** 3',
        description='''
        Calculated value of the equilibrium volume by fitting the volume-energy data.
        ''')

    equilibrium_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Calculated value of the equilibrium energy by fitting the volume-energy data.
        ''')

    rms_error = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Root-mean squared value of the error in the fitting.
        ''')


class DebyeModel(MSection):
    '''
    Section containing results of an debye-model workflow.
    '''

    m_def = Section(
        validate=False)

    n_temperatures = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of temperature sampling points
        ''')

    temperatures = Quantity(
        type=np.dtype(np.float64),
        shape=['n_temperatures'],
        unit='kelvin',
        description='''
        Calculated value of the thermal conductity.
        ''')

    thermal_conductivity = Quantity(
        type=np.dtype(np.float64),
        shape=['n_temperatures'],
        unit='watt / meter * kelvin',
        description='''
        Calculated value of the thermal conductity.
        ''')

    debye_temperature = Quantity(
        type=np.dtype(np.float64),
        shape=['n_temperatures'],
        unit='kelvin',
        description='''
        Calculated value of the Debye temperature.
        ''')

    gruneisen_parameter = Quantity(
        type=np.dtype(np.float64),
        shape=['n_temperatures'],
        description='''
        Calculated value of the Gruneisen parameter.
        ''')

    thermal_expansion = Quantity(
        type=np.dtype(np.float64),
        shape=['n_temperatures'],
        unit='1 / kelvin',
        description='''
        Calculated value of the thermal expansion coefficient.
        ''')

    bulk_modulus_static = Quantity(
        type=np.dtype(np.float64),
        shape=['n_temperatures'],
        unit='pascal',
        description='''
        Calculated value of the static bulk modulus.
        ''')

    bulk_modulus_isothermal = Quantity(
        type=np.dtype(np.float64),
        shape=['n_temperatures'],
        unit='pascal',
        description='''
        Calculated value of the static bulk modulus.
        ''')

    free_energy_gibbs = Quantity(
        type=np.dtype(np.float64),
        shape=['n_temperatures'],
        unit='joule',
        description='''
        Calculated value of the Gibbs free energy, G.
        ''')

    free_energy_vibrational = Quantity(
        type=np.dtype(np.float64),
        shape=['n_temperatures'],
        unit='joule',
        description='''
        Calculated value of the vibrational free energy, F_vib.
        ''')

    internal_energy_vibrational = Quantity(
        type=np.dtype(np.float64),
        shape=['n_temperatures'],
        unit='joule',
        description='''
        Calculated value of the vibrational internal energy, U_vib.
        ''')

    entropy_vibrational = Quantity(
        type=np.dtype(np.float64),
        shape=['n_temperatures'],
        unit='joule / kelvin',
        description='''
        Calculated value of the vibrational entropy, S.
        ''')

    heat_capacity_C_v = Quantity(
        type=np.dtype(np.float64),
        shape=['n_temperatures'],
        unit='joule / kelvin',
        description='''
        Stores the heat capacity per cell unit at constant volume.
        ''')

    heat_capacity_C_p = Quantity(
        type=np.dtype(np.float64),
        shape=['n_temperatures'],
        unit='joule / kelvin',
        description='''
        Stores the heat capacity per cell unit at constant pressure.
        ''')


class GeometryOptimization(MSection):
    '''
    Section containing the results of a geometry_optimization workflow.
    '''

    m_def = Section(
        validate=False)

    type = Quantity(
        type=str,
        shape=[],
        description='''
        The type of geometry optimization can either be ionic, cell_shape, cell_volume.
        ''')

    method = Quantity(
        type=str,
        shape=[],
        description='''
        The method used for geometry optimization.
        ''')

    convergence_tolerance_energy_difference = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        The input energy difference tolerance criterion.
        ''')

    convergence_tolerance_force_maximum = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='newton',
        description='''
        The input maximum net force tolerance criterion.
        ''')

    convergence_tolerance_displacement_maximum = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter',
        description='''
        The input maximum displacement tolerance criterion.
        ''')

    final_energy_difference = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        The difference in the energy_total between the last two steps during
        optimization.
        ''',
        a_search=Search())

    final_force_maximum = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='newton',
        description='''
        The maximum net force in the last optimization step.
        ''')

    final_displacement_maximum = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter',
        description='''
        The maximum displacement in the last optimization step with respect to previous.
        ''')

    optimization_steps = Quantity(
        type=int,
        shape=[],
        description='''
        Number of optimization steps.
        ''')

    energies = Quantity(
        type=np.dtype(np.float64),
        unit='joule',
        shape=["optimization_steps"],
        description='''
        List of energy_total values gathered from the single configuration
        calculations that are a part of the optimization trajectory.
        ''')
    is_converged_geometry = Quantity(
        type=bool,
        shape=[],
        description='''
        Indicates if the geometry convergence criteria were fulfilled.
        ''')


class Phonon(MSection):
    '''
    Section containing the results of a phonon workflow.
    '''

    m_def = Section(
        validate=False)

    force_calculator = Quantity(
        type=str,
        shape=[],
        description='''
        Name of the program used to calculate the forces.
        ''')

    mesh_density = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='1 / meter ** 3',
        description='''
        Density of the k-mesh for sampling.
        ''',
        a_search=Search())

    n_imaginary_frequencies = Quantity(
        type=int,
        shape=[],
        description='''
        Number of modes with imaginary frequencies.
        ''',
        a_search=Search())

    random_displacements = Quantity(
        type=bool,
        shape=[],
        description='''
        Identifies if displacements are made randomly.
        ''')

    with_non_analytic_correction = Quantity(
        type=bool,
        shape=[],
        description='''
        Identifies if non-analytical term corrections are applied to dynamical matrix.
        ''',
        a_search=Search())

    with_grueneisen_parameters = Quantity(
        type=bool,
        shape=[],
        description='''
        Identifies if Grueneisen parameters are calculated.
        ''',
        a_search=Search())


class StrainDiagrams(MSection):
    '''
    Section containing the information regarding the elastic strains.
    '''

    m_def = Section(
        validate=False)

    type = Quantity(
        type=str,
        shape=[],
        description='''
        Kind of strain diagram. Possible values are: energy; cross-validation (cross-
        validation error); d2E (second derivative of the energy wrt the strain)
        ''')

    n_eta = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of strain values used in the strain diagram
        ''')

    n_deformations = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of deformations.
        ''')

    value = Quantity(
        type=np.dtype(np.float64),
        shape=['n_deformations', 'n_eta'],
        description='''
        Values of the energy(units:J)/d2E(units:Pa)/cross-validation (depending on the
        value of strain_diagram_type)
        ''')

    eta = Quantity(
        type=np.dtype(np.float64),
        shape=['n_deformations', 'n_eta'],
        description='''
        eta values used the strain diagrams
        ''')

    stress_voigt_component = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Voigt component corresponding to the strain diagram
        ''')

    polynomial_fit_order = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Order of the polynomial fit
        ''')


class Elastic(MSection):
    '''
    Section containing the results of an elastic workflow.
    '''

    m_def = Section(
        validate=False)

    energy_stress_calculator = Quantity(
        type=str,
        shape=[],
        description='''
        Name of program used to calculate energy or stress.
        ''')

    calculation_method = Quantity(
        type=str,
        shape=[],
        description='''
        Method used to calculate elastic constants, can either be energy or stress.
        ''')

    elastic_constants_order = Quantity(
        type=int,
        shape=[],
        description='''
        Order of the calculated elastic constants.
        ''',
        a_search=Search())

    n_deformations = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of deformed structures used to calculate the elastic constants. This is
        determined by the symmetry of the crystal.
        ''')

    deformation_types = Quantity(
        type=np.dtype('U'),
        shape=['n_deformations', 6],
        description='''
        deformation types
        ''')

    n_strains = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        number of equally spaced strains applied to each deformed structure, which are
        generated between the maximum negative strain and the maximum positive one.
        ''')

    is_mechanically_stable = Quantity(
        type=bool,
        shape=[],
        description='''
        Indicates if structure is mechanically stable from the calculated values of the
        elastic constants.
        ''',
        a_search=Search())

    fitting_error_maximum = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Maximum error in polynomial fit.
        ''')

    strain_maximum = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Maximum strain applied to crystal.
        ''')

    elastic_constants_notation_matrix_second_order = Quantity(
        type=np.dtype('U'),
        shape=[6, 6],
        description='''
        Symmetry of the second-order elastic constant matrix in Voigt notation
        ''')

    elastic_constants_matrix_second_order = Quantity(
        type=np.dtype(np.float64),
        shape=[6, 6],
        unit='pascal',
        description='''
        2nd order elastic constant (stiffness) matrix in GPa
        ''')

    elastic_constants_matrix_third_order = Quantity(
        type=np.dtype(np.float64),
        shape=[6, 6, 6],
        unit='pascal',
        description='''
        3rd order elastic constant (stiffness) matrix in GPa
        ''')

    compliance_matrix_second_order = Quantity(
        type=np.dtype(np.float64),
        shape=[6, 6],
        unit='1 / pascal',
        description='''
        Elastic compliance matrix in 1/GPa
        ''')

    bulk_modulus_voigt = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='pascal',
        description='''
        Voigt bulk modulus
        ''')

    shear_modulus_voigt = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='pascal',
        description='''
        Voigt shear modulus
        ''')

    bulk_modulus_reuss = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='pascal',
        description='''
        Reuss bulk modulus
        ''')

    shear_modulus_reuss = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='pascal',
        description='''
        Reuss shear modulus
        ''')

    bulk_modulus_hill = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='pascal',
        description='''
        Hill bulk modulus
        ''')

    shear_modulus_hill = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='pascal',
        description='''
        Hill shear modulus
        ''')

    young_modulus_voigt = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='pascal',
        description='''
        Voigt Young modulus
        ''')

    poisson_ratio_voigt = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Voigt Poisson ratio
        ''')

    young_modulus_reuss = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='pascal',
        description='''
        Reuss Young modulus
        ''')

    poisson_ratio_reuss = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Reuss Poisson ratio
        ''')

    young_modulus_hill = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='pascal',
        description='''
        Hill Young modulus
        ''')

    poisson_ratio_hill = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Hill Poisson ratio
        ''')

    elastic_anisotropy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Elastic anisotropy
        ''')

    pugh_ratio_hill = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Pugh ratio defined as the ratio between the shear modulus and bulk modulus
        ''')

    debye_temperature = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='kelvin',
        description='''
        Debye temperature
        ''')

    speed_sound_transverse = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter / second',
        description='''
        Speed of sound along the transverse direction
        ''')

    speed_sound_longitudinal = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter / second',
        description='''
        Speed of sound along the longitudinal direction
        ''')

    speed_sound_average = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter / second',
        description='''
        Average speed of sound
        ''')

    eigenvalues_elastic = Quantity(
        type=np.dtype(np.float64),
        shape=[6],
        unit='pascal',
        description='''
        Eigenvalues of the stiffness matrix
        ''')

    strain_diagrams = SubSection(
        sub_section=StrainDiagrams.m_def,
        repeats=True)


class Thermodynamics(MSection):
    '''
    Section containing the results of a thermodynamics workflow.
    '''

    m_def = Section(validate=False)

    n_values = Quantity(
        type=int,
        shape=[],
        description='''
        Number of thermodynamics property evaluations.
        ''')

    temperature = Quantity(
        type=np.dtype(np.float64),
        shape=['n_values'],
        unit='kelvin',
        description='''
        Specifies the temperatures at which properties such as the Helmholtz free energy
        are calculated.
        ''')

    pressure = Quantity(
        type=np.dtype(np.float64),
        shape=['n_values'],
        unit='pascal',
        description='''
        Array containing the values of the pressure (one third of the trace of the stress
        tensor) corresponding to each property evaluation.
        ''')

    helmholz_free_energy = Quantity(
        type=np.dtype(np.float64),
        shape=['n_values'],
        unit='joule',
        description='''
        Helmholtz free energy per unit cell at constant volume.
        ''')

    heat_capacity_c_v = Quantity(
        type=np.dtype(np.float64),
        shape=['n_values'],
        unit='joule / kelvin',
        description='''
        Heat capacity per cell unit at constant volume.
        ''')

    @derived(
        type=np.dtype(np.float64),
        shape=['n_values'],
        unit='joule / (kelvin * kilogram)',
        description='''
        Specific heat capacity at constant volume.
        ''',
        cached=True
    )
    def specific_heat_capacity(self) -> NDArray:
        """Returns the specific heat capacity by dividing the heat capacity per
        cell with the mass of the atoms in the cell.
        """
        import nomad.atomutils
        workflow = self.m_parent
        system = workflow.calculations_ref[0].system_ref[0].value
        atomic_numbers = system.atoms.species
        mass_per_unit_cell = nomad.atomutils.get_summed_atomic_mass(atomic_numbers)
        heat_capacity = self.heat_capacity_c_v
        specific_heat_capacity = heat_capacity / mass_per_unit_cell
        return specific_heat_capacity

    vibrational_free_energy_at_constant_volume = Quantity(
        type=np.dtype(np.float64),
        shape=['n_values'],
        unit='joule',
        description='''
        Holds the vibrational free energy per cell unit at constant volume.
        ''')

    @derived(
        type=np.dtype(np.float64),
        shape=['n_values'],
        unit='joule / kilogram',
        description='''
        Stores the specific vibrational free energy at constant volume.
        ''',
        cached=True
    )
    def specific_vibrational_free_energy_at_constant_volume(self) -> NDArray:
        import nomad.atomutils
        workflow = self.m_parent
        system = workflow.calculations_ref[0].system_ref[0].value
        atomic_numbers = system.atoms.species
        mass_per_unit_cell = nomad.atomutils.get_summed_atomic_mass(atomic_numbers)
        free_energy = self.vibrational_free_energy_at_constant_volume
        specific_free_energy = free_energy / mass_per_unit_cell
        return specific_free_energy


class MolecularDynamics(MSection):
    '''
    Section containing results of molecular dynamics workflow.
    '''

    m_def = Section(
        validate=False)

    ensemble_type = Quantity(
        type=str,
        shape=[],
        description='''
        Ensemble type used.
        ''')

    timestep = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='s',
        description='''
        Time step in the simulation.
        ''')

    finished_normally = Quantity(
        type=bool,
        shape=[],
        description='''
        Indicates if calculation terminated normally.
        ''')

    with_trajectory = Quantity(
        type=bool,
        shape=[],
        description='''
        Indicates if calculation includes trajectory data.
        ''',
        a_search=Search())

    with_thermodynamics = Quantity(
        type=bool,
        shape=[],
        description='''
        Indicates if calculation contains thermodynamic data.
        ''',
        a_search=Search())


class SinglePoint(MSection):
    '''
    Section containing results of a single point workflow.
    '''

    m_def = Section(
        validate=False)

    method = Quantity(
        type=str,
        shape=[],
        description='''
        Calculation method used.
        ''')

    n_scf_steps = Quantity(
        type=int,
        shape=[],
        description='''
        Number of self-consistent steps in the calculation
        ''')

    final_scf_energy_difference = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        The difference in the energy between the last two scf steps.
        ''',
        a_search=Search())

    is_converged = Quantity(
        type=bool,
        shape=[],
        description='''
        Indicates if the convergence criteria were fullfilled
        ''')

    with_density_of_states = Quantity(
        type=bool,
        shape=[],
        description='''
        Indicates if the calculation contains density of states data
        ''')

    with_bandstructure = Quantity(
        type=bool,
        shape=[],
        description='''
        Indicates if the calculation contains bandstructure data
        ''')

    with_eigenvalues = Quantity(
        type=bool,
        shape=[],
        description='''
        Indicates if the calculation contains eigenvalue data
        ''')

    with_volumetric_data = Quantity(
        type=bool,
        shape=[],
        description='''
        Indicates if the calculation contains volumetric data
        ''')

    with_excited_states = Quantity(
        type=bool,
        shape=[],
        description='''
        Indicates if the calculation contains excited states data
        ''')


class WorkflowReference(MSection):
    '''
    Section that provides a link between a section to a section workflow. Such a link is
    necessary for example between an Debye model that uses a the poisson ratio calculated
    from an elastic workflow. The relationship should be described by kind and the
    referred section workflow is given by value. An external url can also be provided in
    place of value.
    '''

    m_def = Section(validate=False)

    kind = Quantity(
        type=str,
        shape=[],
        description='''
        Defines the relationship between the referenced section workflow and the present
        section.
        ''')

    external_url = Quantity(
        type=str,
        shape=[],
        description='''
        URL used to reference an externally stored section workflow.
        ''')

    value = Quantity(
        type=Reference(SectionProxy('Workflow')),
        shape=[],
        description='''
        Value of the referenced section wofklow.
        ''')


class Workflow(MSection):
    '''
    Section containing the  results of a workflow.
    '''

    m_def = Section(
        validate=False)

    type = Quantity(
        type=str,
        shape=[],
        description='''
        The type of calculation workflow. Can be one of geometry_optimization, elastic,
        phonon, molecular_dynamics, single_point, debye_model.
        ''',
        a_search=Search())

    calculator = Quantity(
        type=str,
        shape=[],
        description='''
        Energy and force calculator.
        ''')

    calculation_result_ref = Quantity(
        type=Reference(Calculation.m_def),
        shape=[],
        description='''
        Reference to calculation result. In the case of geometry_optimization and
        molecular dynamics, this corresponds to the final step in the simulation. For the
        rest of the workflow types, it refers to the original system.
        ''',
        categories=[FastAccess])

    calculations_ref = Quantity(
        type=Reference(Calculation.m_def),
        shape=['optimization_steps'],
        description='''
        List of references to each section single_configuration_calculation in the
        simulation.
        ''')

    single_point = SubSection(
        sub_section=SinglePoint.m_def,
        # TODO determine if there is a need for this to be a repeating section
        # such as in the context of fhi-vibes single_point
        repeats=False,
        categories=[FastAccess])

    geometry_optimization = SubSection(
        sub_section=GeometryOptimization.m_def,
        repeats=False,
        categories=[FastAccess])

    phonon = SubSection(
        sub_section=Phonon.m_def,
        repeats=False,
        categories=[FastAccess])

    elastic = SubSection(
        sub_section=Elastic.m_def,
        repeats=False)

    molecular_dynamics = SubSection(
        sub_section=MolecularDynamics.m_def,
        repeats=False)

    debye_model = SubSection(
        sub_section=DebyeModel.m_def,
        repeats=False)

    equation_of_state = SubSection(
        sub_section=EquationOfState.m_def,
        repeats=False)

    nudged_elastic_band = SubSection(
        sub_section=NudgedElasticBand.m_def,
        repeats=False)

    convex_hull = SubSection(
        sub_section=ConvexHull.m_def,
        repeats=False)

    adsorption = SubSection(
        sub_section=Adsorption.m_def,
        repeats=False)

    magnetic_ordering = SubSection(
        sub_section=MagneticOrdering.m_def,
        repeats=False)

    raman = SubSection(
        sub_section=Raman.m_def,
        repeats=False)

    thermodynamics = SubSection(
        sub_section=Thermodynamics.m_def,
        repeats=False)

    workflow_ref = SubSection(
        sub_section=WorkflowReference.m_def, repeats=True)

    run_ref = SubSection(
        sub_section=RunReference.m_def, repeats=True)
