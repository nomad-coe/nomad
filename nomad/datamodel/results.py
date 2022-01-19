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

import numpy as np
from elasticsearch_dsl import Text

from ase.data import chemical_symbols
from ase.formula import Formula

from nomad import config
from nomad import atomutils
from nomad.datamodel.metainfo.measurements import Spectrum
from nomad.datamodel.metainfo.workflow import EquationOfState, EOSFit
from nomad.metainfo.elasticsearch_extension import (
    Elasticsearch,
    material_type,
    material_entry_type,
    get_tokenizer
)

from nomad.metainfo import (
    MSection,
    Section,
    SubSection,
    Quantity,
    MEnum,
    Package,
    Datetime,
)

m_package = Package()

from nomad.datamodel.optimade import Species  # noqa
from nomad.datamodel.metainfo.simulation.calculation import (
    Dos,
    BandStructure as BandStructureCalculation,
    BandEnergies,
    DosValues,
    Calculation)  # noqa
from nomad.datamodel.metainfo.simulation.method import (
    BasisSet, Scf, Electronic, Smearing, GW as GWMethod
)  # noqa
from nomad.datamodel.metainfo.workflow import (
    GeometryOptimization, Phonon, Elastic, Thermodynamics
)  # noqa


unavailable = 'unavailable'
not_processed = 'not processed'
structure_classes = [
    'bulk',
    'surface',
    '2D',
    '1D',
    'molecule / cluster',
    'atom',
    unavailable,
    not_processed,
]
bravais_lattices = [
    'aP',
    'mP',
    'mS',
    'oP',
    'oS',
    'oF',
    'oI',
    'tP',
    'tI',
    'hP',
    'hR',
    'cP',
    'cF',
    'cI',
]
crystal_systems = [
    'triclinic',
    'monoclinic',
    'orthorhombic',
    'tetragonal',
    'trigonal',
    'hexagonal',
    'cubic',
]
xc_treatments = {
    'gga': 'GGA',
    'hf_': 'HF',
    'oep': 'OEP',
    'hyb': 'hybrid',
    'mgg': 'meta-GGA',
    'vdw': 'vdW',
    'lda': 'LDA',
}
basis_set_types = [
    '(L)APW+lo',
    'gaussians',
    'numeric AOs',
    'plane waves',
    'psinc functions',
    'real-space grid',
    unavailable,
    not_processed,
]
core_electron_treatments = [
    'full all electron',
    'all electron frozen core',
    'pseudopotential',
    unavailable,
]


def variants_formula(value):
    ''' Creates several common variants for the given formula.'''
    formula = Formula(value)
    formats = ['hill', 'metal', 'abc']
    formulas = [value] + [formula.format(f) for f in formats]
    return list(set(formulas))


# Tokenizes a formula by splitting it by formula fragments.
tokenizer_formula = get_tokenizer(r'[A-Z][a-z]?\d*')


class BandGap(MSection):
    m_def = Section(
        description='''
        Band gap information for each spin channel.
        '''
    )
    index = Quantity(
        type=np.dtype(np.int64),
        description='''
        Spin channel index.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    value = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Band gap value. Value of zero corresponds to not having a band gap.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    type = Quantity(
        type=MEnum('direct', 'indirect'),
        shape=[],
        description='''
        Band gap type.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='simple')
        ],
    )
    energy_highest_occupied = Quantity(
        type=np.dtype(np.float64),
        unit='joule',
        shape=[],
        description='''
        The highest occupied energy.
        ''',
    )
    energy_lowest_unoccupied = Quantity(
        type=np.dtype(np.float64),
        unit='joule',
        shape=[],
        description='''
        The lowest unoccupied energy.
        ''',
    )


class WyckoffSet(MSection):
    m_def = Section(
        description='''
        Section for storing Wyckoff set information. Only available for
        conventional cells that have undergone symmetry analysis.
        '''
    )
    wyckoff_letter = Quantity(
        type=str,
        description='''
        The Wyckoff letter for this set.
        '''
    )
    indices = Quantity(
        type=np.dtype('i4'),
        shape=['1..*'],
        description='''
        Indices of the atoms belonging to this group.
        '''
    )
    element = Quantity(
        type=str,
        description='''
        Chemical element at this Wyckoff position.
        '''
    )
    x = Quantity(
        type=np.dtype(np.float64),
        description='''
        The free parameter x if present.
        '''
    )
    y = Quantity(
        type=np.dtype(np.float64),
        description='''
        The free parameter y if present.
        '''
    )
    z = Quantity(
        type=np.dtype(np.float64),
        description='''
        The free parameter z if present.
        '''
    )


class LatticeParameters(MSection):
    m_def = Section(
        description='''
        Lattice parameters of a cell.
        ''',
    )
    a = Quantity(
        type=np.dtype(np.float64),
        unit='m',
        description='''
        Length of the first basis vector.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    b = Quantity(
        type=np.dtype(np.float64),
        unit='m',
        description='''
        Length of the second basis vector.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    c = Quantity(
        type=np.dtype(np.float64),
        unit='m',
        description='''
        Length of the third basis vector.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    alpha = Quantity(
        type=np.dtype(np.float64),
        unit='radian',
        description='''
        Angle between second and third basis vector.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    beta = Quantity(
        type=np.dtype(np.float64),
        unit='radian',
        description='''
        Angle between first and third basis vector.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    gamma = Quantity(
        type=np.dtype(np.float64),
        unit='radian',
        description='''
        Angle between first and second basis vector.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )


class Structure(MSection):
    m_def = Section(
        description='''
        Describes an atomistic structure.
        '''
    )
    dimension_types = Quantity(
        type=int,
        shape=[3],
        default=[0, 0, 0],
        description='''
        List of three integers. For each of the three directions indicated by
        the three lattice vectors (see property lattice_vectors). This list
        indicates if the direction is periodic (value 1) or non-periodic (value
        0). Note: the elements in this list each refer to the direction of the
        corresponding entry in lattice_vectors and not the Cartesian x, y, z
        directions.
        '''
    )
    nperiodic_dimensions = Quantity(
        type=int,
        derived=lambda a: sum(a.dimension_types),
        description='''
        An integer specifying the number of periodic dimensions in the
        structure, equivalent to the number of non-zero entries in
        dimension_types.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    lattice_vectors = Quantity(
        type=np.dtype('float64'),
        shape=[3, 3],
        unit='m',
        description='''
        The three lattice vectors in Cartesian coordinates.
        '''
    )
    cartesian_site_positions = Quantity(
        type=np.dtype('float64'),
        shape=['n_sites', 3],
        unit='m',
        description='''
        Cartesian positions of each site. A site is an atom, a site potentially
        occupied by an atom, or a placeholder for a virtual mixture of atoms
        (e.g., in a virtual crystal approximation).
        '''
    )
    n_sites = Quantity(
        type=int,
        default=0,
        derived=lambda a: len(a.cartesian_site_positions) if a.cartesian_site_positions is not None else 0,
        description='''
        An integer specifying the length of the cartesian_site_positions property.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    species_at_sites = Quantity(
        type=str,
        shape=['n_sites'],
        description='''
        Name of the species at each site (where values for sites are specified with the same
        order of the cartesian_site_positions property). The properties of the species are
        found in the species property.
        '''
    )
    cell_volume = Quantity(
        type=np.dtype(np.float64),
        unit='m ** 3',
        description='''
        Volume of the cell.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    atomic_density = Quantity(
        type=np.dtype(np.float64),
        unit='1 / m ** 3',
        description='''
        Atomic density of the material (atoms/volume).'
        '''
    )
    mass_density = Quantity(
        type=np.dtype(np.float64),
        unit='kg / m ** 3',
        description='''
        Mass density of the material.
        '''
    )
    species = SubSection(sub_section=Species.m_def, repeats=True)
    lattice_parameters = SubSection(sub_section=LatticeParameters.m_def)


class StructureOriginal(Structure):
    m_def = Section(
        description='''
        Contains a selected representative structure from the the original
        data.
        '''
    )


class StructurePrimitive(Structure):
    m_def = Section(
        description='''
        Contains the primitive structure that is derived from
        structure_original. This primitive stucture has been idealized and the
        conventions employed by spglib are used.
        '''
    )


class StructureConventional(Structure):
    m_def = Section(
        description='''
        Contains the conventional structure that is derived from
        structure_original. This conventional stucture has been idealized and
        the conventions employed by spglib are used.
        '''
    )
    wyckoff_sets = SubSection(sub_section=WyckoffSet.m_def, repeats=True)


class StructureOptimized(Structure):
    m_def = Section(
        description='''
        Contains a structure that is the result of a geometry optimization.
        '''
    )


class Structures(MSection):
    m_def = Section(
        description='''
        Contains full atomistic representations of the material in different
        forms.
        ''',
    )
    structure_original = SubSection(
        sub_section=StructureOriginal.m_def,
        repeats=False,
    )
    structure_conventional = SubSection(
        sub_section=StructureConventional.m_def,
        repeats=False,
    )
    structure_primitive = SubSection(
        sub_section=StructurePrimitive.m_def,
        repeats=False,
    )


class Symmetry(MSection):
    m_def = Section(
        description='''
        Section containing information about the symmetry of the material. All
        of these properties are derived by running a symmetry analysis on a
        representative geometry from the original data. This original geometry
        is stored in results.properties together with the primitive and
        conventional structures.
        '''
    )
    bravais_lattice = Quantity(
        type=MEnum(bravais_lattices),
        shape=[],
        description='''
        Identifier for the Bravais lattice in Pearson notation. The first lowercase letter
        identifies the crystal family and can be one of the following: a (triclinic), b
        (monoclinic), o (orthorhombic), t (tetragonal), h (hexagonal) or c (cubic). The
        second uppercase letter identifies the centring and can be one of the following: P
        (primitive), S (face centred), I (body centred), R (rhombohedral centring) or F
        (all faces centred).
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion='simple')
        ],
    )
    crystal_system = Quantity(
        type=MEnum(crystal_systems),
        shape=[],
        description='''
        Name of the crystal system.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion='simple')
        ],
    )
    hall_number = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        The Hall number for this system.
        ''',
        a_elasticsearch=Elasticsearch(material_type),
    )
    hall_symbol = Quantity(
        type=str,
        shape=[],
        description='''
        The Hall symbol for this system.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion='simple')
        ],
    )
    point_group = Quantity(
        type=str,
        shape=[],
        description='''
        Symbol of the crystallographic point group in the Hermann-Mauguin notation.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion='simple')
        ],
    )
    space_group_number = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Specifies the International Union of Crystallography (IUC) number of the 3D space
        group of this system.
        ''',
        a_elasticsearch=Elasticsearch(material_type),
    )
    space_group_symbol = Quantity(
        type=str,
        shape=[],
        description='''
        The International Union of Crystallography (IUC) short symbol of the 3D
        space group of this system.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion='simple')
        ],
    )
    prototype_formula = Quantity(
        type=str,
        description='''
        The formula of the prototypical material for this structure.
        ''',
        a_elasticsearch=Elasticsearch(material_type),
    )
    prototype_aflow_id = Quantity(
        type=str,
        description='''
        The identifier of this structure in the AFLOW encyclopedia of
        crystallographic prototypes:
        http://www.aflowlib.org/prototype-encyclopedia/index.html
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion='simple')
        ],
    )
    structure_name = Quantity(
        type=str,
        description='''
        A common name for this structure, e.g. fcc, bcc.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion='default')
        ],
    )
    strukturbericht_designation = Quantity(
        type=str,
        description='''
        Classification of the material according to the historically grown
        'strukturbericht'.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion='simple')
        ],
    )


class Material(MSection):
    m_def = Section(
        description='''
        Contains information that is specific to bulk crystalline materials.
        '''
    )
    material_id = Quantity(
        type=str,
        description='''
        A fixed length, unique material identifier in the form of a hash
        digest.
        ''',
        a_elasticsearch=Elasticsearch(material_type, metrics=dict(n_materials='cardinality'))
    )
    material_name = Quantity(
        type=str,
        description='''
        Meaningful names for this a material if any can be assigned.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion='default')
        ],
    )
    structural_type = Quantity(
        type=MEnum(structure_classes), default='not processed',
        description='''
        Classification based on structural features.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion='default')
        ],
    )
    functional_type = Quantity(
        type=str,
        shape=['0..*'],
        description='''
        Classification based on the functional properties.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type, default_aggregation_size=20),
            Elasticsearch(suggestion='default')
        ],
    )
    compound_type = Quantity(
        type=str,
        shape=['0..*'],
        description='''
        Classification based on the chemical formula.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type, default_aggregation_size=20),
            Elasticsearch(suggestion='default')
        ],
    )
    elements = Quantity(
        type=MEnum(chemical_symbols),
        shape=['0..*'],
        default=[],
        description='''
        Names of the different elements present in the structure.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type, many_all=True),
            Elasticsearch(suggestion='simple')
        ]
    )
    n_elements = Quantity(
        type=int,
        default=0,
        derived=lambda s: len(s.elements),
        description='''
        Number of different elements in the structure as an integer.
        ''',
        a_elasticsearch=Elasticsearch(material_type),
    )
    elements_exclusive = Quantity(
        type=str,
        derived=lambda s: ' '.join(sorted(s.elements)),
        description='''
        String containing the chemical elements in alphabetical order and
        separated by a single whitespace. This quantity can be used for
        exclusive element searches where you want to find entries/materials
        with only certain given elements.
        ''',
        a_elasticsearch=Elasticsearch(material_type)
    )
    chemical_formula_descriptive = Quantity(
        type=str,
        description='''
            The chemical formula for a structure as a string in a form chosen by the API
            implementation.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion=tokenizer_formula)
        ],
    )
    chemical_formula_reduced = Quantity(
        type=str,
        description='''
            The reduced chemical formula for a structure as a string with element symbols and
            integer chemical proportion numbers. The proportion number MUST be omitted if it is 1.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion=tokenizer_formula)
        ],
    )
    chemical_formula_hill = Quantity(
        type=str,
        description='''
            The chemical formula for a structure in Hill form with element symbols followed by
            integer chemical proportion numbers. The proportion number MUST be omitted if it is 1.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type, normalizer=atomutils.get_formula_hill),
            Elasticsearch(suggestion=tokenizer_formula, variants=variants_formula)
        ],
    )
    chemical_formula_anonymous = Quantity(
        type=str,
        description='''
            The anonymous formula is the chemical_formula_reduced, but where the elements are
            instead first ordered by their chemical proportion number, and then, in order left to
            right, replaced by anonymous symbols A, B, C, ..., Z, Aa, Ba, ..., Za, Ab, Bb, ... and
            so on.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion=tokenizer_formula)
        ],
    )
    chemical_formula_reduced_fragments = Quantity(
        type=str,
        shape=['*'],
        description='''
        The reduced formula separated into individual terms containing both the atom
        type and count. Used for searching parts of a formula.
        ''',
        a_elasticsearch=Elasticsearch(material_type, mapping=Text(multi=True)),
    )
    symmetry = SubSection(sub_section=Symmetry.m_def, repeats=False)


class DFT(MSection):
    m_def = Section(
        description='''
        Methodology for a DFT calculation.
        '''
    )
    basis_set_type = Quantity(
        type=MEnum(basis_set_types),
        default=unavailable,
        description='The used basis set functions.',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )
    basis_set_name = BasisSet.name.m_copy()
    basis_set_name.m_annotations['elasticsearch'] = [
        Elasticsearch(material_entry_type),
        Elasticsearch(suggestion='default')
    ]
    core_electron_treatment = Quantity(
        type=MEnum(core_electron_treatments),
        default=unavailable,
        description='''
        How the core electrons are described.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )
    spin_polarized = Quantity(
        type=bool,
        description='''
        Whether the calculation is spin-polarized.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    scf_threshold_energy_change = Scf.threshold_energy_change.m_copy()
    scf_threshold_energy_change.m_annotations['elasticsearch'] = Elasticsearch(material_entry_type)
    van_der_Waals_method = Electronic.van_der_waals_method.m_copy()
    van_der_Waals_method.description = 'The used van der Waals method.'
    van_der_Waals_method.m_annotations['elasticsearch'] = [
        Elasticsearch(material_entry_type),
        Elasticsearch(suggestion='default')
    ]

    relativity_method = Electronic.relativity_method.m_copy()
    relativity_method.m_annotations['elasticsearch'] = [
        Elasticsearch(material_entry_type),
        Elasticsearch(suggestion='default')
    ]

    smearing_kind = Smearing.kind.m_copy()
    smearing_kind.m_annotations['elasticsearch'] = [
        Elasticsearch(material_entry_type),
        Elasticsearch(suggestion='default')
    ]

    smearing_width = Smearing.width.m_copy()
    smearing_width.m_annotations['elasticsearch'] = Elasticsearch(material_entry_type)
    xc_functional_type = Quantity(
        type=MEnum(list(xc_treatments.values()) + [unavailable, not_processed]),
        default=not_processed,
        description='The libXC based xc functional classification used in the simulation.',
        a_elasticsearch=Elasticsearch(material_entry_type, default_aggregation_size=100)
    )
    xc_functional_names = Quantity(
        type=str,
        default=[],
        shape=['*'],
        description='The list of libXC functional names that where used in this entry.',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ]
    )


class GW(MSection):
    m_def = Section(
        description='''
        Methodology for a GW calculation.
        '''
    )
    type = GWMethod.type.m_copy()
    type.m_annotations['elasticsearch'] = [
        Elasticsearch(material_entry_type),
        Elasticsearch(suggestion='default')
    ]
    starting_point = Quantity(
        type=str,
        default=[],
        shape=['*'],
        description='The list of libXC functional names that were used for the ground state calculation.',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )


class GeometryOptimizationMethod(MSection):
    m_def = Section(
        description='''
        Geometry optimization methodology. This methodology applies to the
        properties presented in results.properties.geometry_optimization.
        ''',
    )
    type = GeometryOptimization.type.m_copy()
    convergence_tolerance_energy_difference = GeometryOptimization.convergence_tolerance_energy_difference.m_copy()
    convergence_tolerance_energy_difference.m_annotations['elasticsearch'] = Elasticsearch(material_entry_type)
    convergence_tolerance_force_maximum = GeometryOptimization.convergence_tolerance_force_maximum.m_copy()
    convergence_tolerance_force_maximum.m_annotations['elasticsearch'] = Elasticsearch(material_entry_type)


class MolecularDynamicsMethod(MSection):
    m_def = Section(
        description='''
        Molecular dynamics methodology. This methodology applies to the
        properties presented in results.properties.molecular_dynamics.
        ''',
    )


class VibrationalMethod(MSection):
    m_def = Section(
        description='''
        Vibrational properties methodology. This methodology applies to the
        properties presented in results.properties.vibrational.
        ''',
    )
    force_calculator = Phonon.force_calculator.m_copy()
    mesh_density = Phonon.mesh_density.m_copy()
    random_displacements = Phonon.random_displacements.m_copy()
    with_non_analytic_correction = Phonon.with_non_analytic_correction.m_copy()


class ElasticMethod(MSection):
    m_def = Section(
        description='''
        Elastic properties methodology. This methodology applies to the
        properties presented in results.properties.elastic.
        ''',
    )
    energy_stress_calculator = Elastic.energy_stress_calculator.m_copy()
    elastic_calculation_method = Elastic.calculation_method.m_copy()
    elastic_constants_order = Elastic.elastic_constants_order.m_copy()
    fitting_error_maximum = Elastic.fitting_error_maximum.m_copy()
    strain_maximum = Elastic.strain_maximum.m_copy()


class QuantumCircuit(MSection):
    processors = Quantity(type=str, shape=['0..*'])
    number_of_registers = Quantity(type=int)
    simulated = Quantity(type=bool)


class QuantumCMS(MSection):
    transformation = Quantity(type=str)
    quantum_computer_system = Quantity(type=str)
    quantum_computing_libraries = Quantity(type=str, shape=['0..*'])
    computation_datetime = Quantity(type=Datetime)
    number_of_shots = Quantity(type=int)
    quantum_volume = Quantity(type=int)
    quantum_circuit = SubSection(sub_section=QuantumCircuit)


class Simulation(MSection):
    m_def = Section(
        description='''
        Contains method details for a simulation entry.
        '''
    )
    program_name = Quantity(
        type=str,
        default='not processed',
        description='The name of the used program.',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )
    program_version = Quantity(
        type=str,
        default='not processed',
        description='The version of the used program.',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )
    dft = SubSection(sub_section=DFT.m_def, repeats=False)
    gw = SubSection(sub_section=GW.m_def, repeats=False)
    geometry_optimization = SubSection(sub_section=GeometryOptimizationMethod.m_def, repeats=False)
    molecular_dynamics = SubSection(sub_section=MolecularDynamicsMethod.m_def, repeats=False)
    vibrational = SubSection(sub_section=VibrationalMethod.m_def, repeats=False)
    elastic = SubSection(sub_section=ElasticMethod.m_def, repeats=False)
    quantum_cms = SubSection(sub_section=QuantumCMS.m_def, repeats=False)


class Method(MSection):
    m_def = Section(
        description='''
        Contains a summary of the methodology that has been used in this entry.
        This methodology applies to all of the reported properties and
        determines the result of a single energy evalution. The individual
        properties may be further methodological details affect e.g. the
        sampling.
        '''
    )
    method_id = Quantity(
        type=str,
        description='''
        Identifier for the used method. Only available for a subset of entries
        for which the methodology has been identified with precision.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    equation_of_state_id = Quantity(
        type=str,
        description='''
        Identifier that can be used to group entries within an equation of
        state calculation. Only available for a subset of entries for which the
        structure and methodology have been identified with precision.
        ''',
    )
    parameter_variation_id = Quantity(
        type=str,
        description='''
        Identifier that can be used to group entries that target the same
        structure but with varying parameter settings. Only available for a
        subset of entries for which the structure and methodology have been
        identified with precision.
        ''',
    )
    method_name = Quantity(
        type=MEnum(['DFT', 'GW', 'EELS', 'XPS', config.services.unavailable_value]),
        description='''
        Common name for the used method.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )
    simulation = SubSection(sub_section=Simulation.m_def, repeats=False)


class DOS(MSection):
    m_def = Section(
        description='''
        Base class for density of states information.
        ''',
    )
    energies = Quantity(
        type=Dos.energies,
        description='''
        Array containing the set of discrete energy values for the density of
        states (DOS).
        ''',
    )
    total = Quantity(
        type=DosValues,
        shape=['*'],
        description='''
        Density of states (DOS) values normalized with unit cell volume and
        number of atoms.
        ''',
    )


class DOSPhonon(DOS):
    m_def = Section(
        description='''
        Contains the total phonon density of states.
        ''',
    )


class DOSElectronic(DOS):
    m_def = Section(
        description='''
        Contains the total electronic density of states.
        ''',
    )
    spin_polarized = Quantity(
        type=bool,
        description='''
        Whether the DOS is spin-polarized, i.e. is contains channels for both
        spin values.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    band_gap = SubSection(
        sub_section=BandGap.m_def,
        repeats=True,
        a_elasticsearch=Elasticsearch(material_entry_type, nested=True)
    )
    energy_fermi = Quantity(
        type=np.dtype(np.float64),
        unit='joule',
        shape=[],
        description='''
        Fermi energy.
        '''
    )


class BandStructure(MSection):
    m_def = Section(
        description='''
        Base class for band structure information.
        ''',
    )
    reciprocal_cell = Quantity(
        type=BandStructureCalculation.reciprocal_cell,
        description='''
        The reciprocal cell within which the band structure is calculated.
        ''',
    )
    segment = Quantity(
        type=BandEnergies,
        shape=['*'],
        description='''
        Collection of linear path segments in the reciprocal space. The
        segments are represented as third-order tensors: one dimension for the
        spin channels, one for the sequence of reciprocal space points for the
        segment, and one for the sequence of eigenvalues at a given point.
        ''',
    )
    path_standard = Quantity(
        type=str,
        shape=[],
        description='''
        String that identifies the possible standard used in sampling the
        reciprocal space.
        ''',
    )


class BandStructurePhonon(BandStructure):
    m_def = Section(
        description='''
        This section stores information on a vibrational band structure
        evaluation along one-dimensional pathways in the reciprocal space.
        '''
    )


class BandStructureElectronic(BandStructure):
    m_def = Section(
        description='''
        This section stores information on a electonic band structure
        evaluation along one-dimensional pathways in the reciprocal space.
        '''
    )
    spin_polarized = Quantity(
        type=bool,
        description='''
        Whether the band structure is spin-polarized, i.e. is contains channels
        for both spin values.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    band_gap = SubSection(
        sub_section=BandGap.m_def,
        repeats=True,
        a_elasticsearch=Elasticsearch(material_entry_type, nested=True)
    )
    energy_fermi = Quantity(
        type=np.dtype(np.float64),
        unit='joule',
        shape=[],
        description='''
        Fermi energy.
        '''
    )


class HeatCapacityConstantVolume(MSection):
    m_def = Section(
        description='''
        Contains the values of the specific (per mass) and isochoric (constant
        volume) heat capacity at different temperatures.
        '''
    )
    heat_capacities = Quantity(
        type=Thermodynamics.heat_capacity_c_v,
        shape=[],
        description='''
        Specific heat capacity values at constant volume.
        ''',
    )
    temperatures = Quantity(
        type=Thermodynamics.temperature,
        description='''
        The temperatures at which heat capacities are calculated.
        ''',
    )


class EnergyFreeHelmholtz(MSection):
    m_def = Section(
        description='''
        Contains the values of the Helmholtz free energy per atom at constant
        volume and at different temperatures.
        '''
    )
    energies = Quantity(
        type=Thermodynamics.vibrational_free_energy_at_constant_volume,
        shape=[],
        description='''
        The Helmholtz free energies per atom at constant volume.
        ''',
    )
    temperatures = Quantity(
        type=Thermodynamics.temperature,
        description='''
        The temperatures at which Helmholtz free energies are calculated.
        ''',
    )


class GeometryOptimizationProperties(MSection):
    m_def = Section(
        description='''
        Properties from a geometry optimization.
        ''',
    )
    trajectory = Quantity(
        type=Calculation,
        shape=['0..*'],
        description='''
        List of references to each section_single_configuration_calculation in
        the optimization trajectory.
        ''',
    )
    energies = Quantity(
        type=GeometryOptimization.energies,
        description='''
        List of energy_total values gathered from the single configuration
        calculations that are a part of the optimization trajectory.
        ''',
    )
    structure_optimized = SubSection(
        sub_section=StructureOptimized.m_def,
        repeats=False,
    )
    final_force_maximum = GeometryOptimization.final_force_maximum.m_copy()
    final_force_maximum.m_annotations['elasticsearch'] = Elasticsearch(material_entry_type)
    final_energy_difference = GeometryOptimization.final_energy_difference.m_copy()
    final_energy_difference.m_annotations['elasticsearch'] = Elasticsearch(material_entry_type)


class MolecularDynamicsProperties(MSection):
    m_def = Section(
        description='''
        Properties from molecular_dynamics.
        ''',
    )
    trajectory = Quantity(
        type=Calculation,
        shape=['0..*'],
        description='''
        List of references to each section_single_configuration_calculation in
        the molecular dynamics trajectory.
        ''',
    )


class VibrationalProperties(MSection):
    m_def = Section(
        description='''
        Vibrational properties.
        ''',
    )
    band_structure_phonon = SubSection(sub_section=BandStructurePhonon.m_def, repeats=False)
    dos_phonon = SubSection(sub_section=DOSPhonon.m_def, repeats=False)
    heat_capacity_constant_volume = SubSection(sub_section=HeatCapacityConstantVolume.m_def, repeats=False)
    energy_free_helmholtz = SubSection(sub_section=EnergyFreeHelmholtz.m_def, repeats=False)


class EnergyVolumeCurve(MSection):
    m_def = Section(
        description='''
        Energy volume curve.
        ''',
    )
    type = Quantity(
        type=MEnum(
            'raw',
            'mie_gruneisen',
            'pack_evans_james',
            'vinet',
            'tait',
            'birch_euler',
            'pourier_tarantola',
            'birch_lagrange',
            'murnaghan',
        ),
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )
    volumes = Quantity(type=EquationOfState.volumes)
    energies_raw = Quantity(type=EquationOfState.energies)
    energies_fit = Quantity(type=EOSFit.fitted_energies)


class BulkModulus(MSection):
    m_def = Section(
        description='''
        Contains bulk modulus values calculated with different methodologies.
        '''
    )
    type = Quantity(
        type=MEnum(
            'mie_gruneisen',
            'pack_evans_james',
            'vinet',
            'tait',
            'birch_euler',
            'pourier_tarantola',
            'birch_lagrange',
            'murnaghan',
            'voigt_average',
            'reuss_average',
            'voigt_reuss_hill_average',
        ),
        description='Describes the methodology for obtaining the value.',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )
    value = Quantity(
        type=np.dtype(np.float64),
        description='Bulk modulus value.',
        unit='pascal',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )


class ShearModulus(MSection):
    m_def = Section(
        description='''
        Contains shear modulus values calculated with different methodologies.
        ''',
    )
    type = Quantity(
        type=MEnum(
            'voigt_average',
            'reuss_average',
            'voigt_reuss_hill_average',
        ),
        description='Describes the methodology for obtaining the value.',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )
    value = Quantity(
        type=np.dtype(np.float64),
        description='Shear modulus value.',
        unit='pascal',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )


class MechanicalProperties(MSection):
    m_def = Section(
        description='''
        Mechanical properties.
        ''',
    )
    energy_volume_curve = SubSection(sub_section=EnergyVolumeCurve.m_def, repeats=True)
    bulk_modulus = SubSection(
        sub_section=BulkModulus.m_def,
        repeats=True,
        a_elasticsearch=Elasticsearch(material_entry_type, nested=True)
    )
    shear_modulus = SubSection(
        sub_section=ShearModulus.m_def,
        repeats=True,
        a_elasticsearch=Elasticsearch(material_entry_type, nested=True)
    )


class ElectronicProperties(MSection):
    m_def = Section(
        description='''
        Electronic properties.
        ''',
    )
    band_structure_electronic = SubSection(sub_section=BandStructureElectronic.m_def, repeats=False)
    dos_electronic = SubSection(sub_section=DOSElectronic.m_def, repeats=False)


class SpectroscopyProperties(MSection):
    m_def = Section(
        description='''
        Spectroscopic properties.
        ''',
    )
    spectrum = Quantity(type=Spectrum)


class Properties(MSection):
    m_def = Section(
        description='''
        Contains the physical properties that have been calculated or used in
        this entry.
        '''
    )
    structures = SubSection(sub_section=Structures.m_def, repeats=False)
    geometry_optimization = SubSection(sub_section=GeometryOptimizationProperties.m_def, repeats=False)
    molecular_dynamics = SubSection(sub_section=MolecularDynamicsProperties.m_def, repeats=False)
    vibrational = SubSection(sub_section=VibrationalProperties.m_def, repeats=False)
    electronic = SubSection(sub_section=ElectronicProperties.m_def, repeats=False)
    mechanical = SubSection(sub_section=MechanicalProperties.m_def, repeats=False)
    spectroscopy = SubSection(sub_section=SpectroscopyProperties.m_def, repeats=False)

    n_calculations = Quantity(
        type=int,
        description='''
        The number of performed single configuration calculations.'
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type, metrics=dict(n_calculations='sum')),
    )
    available_properties = Quantity(
        type=str,
        shape=['0..*'],
        description='List of all properties which are present in this section_results',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )


class Results(MSection):
    m_def = Section(
        description='''
        Contains a summary of the entry contents.
        '''
    )
    material = SubSection(sub_section=Material.m_def, repeats=False)
    method = SubSection(sub_section=Method.m_def, repeats=False)
    properties = SubSection(sub_section=Properties.m_def, repeats=False)


m_package.__init_metainfo__()
