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

from typing import List
import numpy as np
from elasticsearch_dsl import Text

from ase.data import chemical_symbols

from nomad import config
from nomad.datamodel.metainfo.common import ProvenanceTracker, PropertySection
from nomad.datamodel.metainfo.simulation.system import Atoms, System as SystemRun
from nomad.datamodel.metainfo.simulation.workflow import (
    EquationOfStateResults,
    EOSFit,
    RadialDistributionFunction as RDFWorkflow,
    RadialDistributionFunctionValues,
    MeanSquaredDisplacement as MSDWorkflow,
    MeanSquaredDisplacementValues,
    DiffusionConstantValues
)
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

from nomad.datamodel.optimade import Species as OptimadeSpecies  # noqa
from nomad.datamodel.metainfo.simulation.calculation import (
    RadiusOfGyration as RgCalculation,
    RadiusOfGyrationValues,
    Dos,
    BandStructure as BandStructureCalculation,
    GreensFunctions as GreensFunctionsCalculation,
    BandEnergies,
    DosValues,
    BandGapDeprecated as CalcBandGapDeprecated,
    BandGap as CalcBandGap,
    Calculation,
    ElectronicStructureProvenance
)  # noqa
from nomad.datamodel.metainfo.simulation.method import (
    Scf, Electronic, Smearing, ExcitedStateMethodology as XSMethod,
    GW as GWMethod, BSE as BSEMethod, HubbardKanamoriModel as HubbardKanamori,
    AtomParameters, DMFT as DMFTMethod, BasisSet, BasisSetContainer
)  # noqa
from nomad.datamodel.metainfo.simulation.workflow import (
    GeometryOptimization as MGeometryOptimization,
    GeometryOptimizationMethod,
    GeometryOptimizationResults,
    ThermodynamicsResults, MolecularDynamicsMethod
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
orbitals = ['s', 'p', 'd', 'f']
orbitals += ['{}{}'.format(n, orbital) for n in range(1, 10) for orbital in orbitals]
xc_treatments = {  # Note: respect the ordering (Python3.7+), is used in DFTMethod.xc_functional_type()
    'lda': 'LDA',
    'gga': 'GGA',
    'mgg': 'meta-GGA',
    'hyb_mgg': 'hyper-GGA',
    'hyb': 'hybrid',
}
xc_treatments_extended = {**xc_treatments, 'hf_': 'HF'}
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

structure_name_map = {
    'CaTiO<sub>3</sub> Pnma Perovskite Structure': 'perovskite',
    'Hypothetical Tetrahedrally Bonded Carbon with 4&ndash;Member Rings': '4-member ring',
    'In (A6) Structure': 'fct',
    '$\\alpha$&ndash;Pa (A<sub>a</sub>) Structure': 'bct',
    'Hypothetical BCT5 Si Structure': 'bct5',
    'Wurtzite (ZnS, B4) Structure': 'wurtzite',
    'Hexagonal Close Packed (Mg, A3) Structure': 'hcp',
    'Half&ndash;Heusler (C1<sub>b</sub>) Structure': 'half-Heusler',
    'Zincblende (ZnS, B3) Structure': 'zincblende',
    'Cubic Perovskite (CaTiO<sub>3</sub>, E2<sub>1</sub>) Structure': 'perovskite',
    '$\\alpha$&ndash;Po (A<sub>h</sub>) Structure': 'simple cubic',
    'Si<sub>46</sub> Clathrate Structure': 'clathrate',
    'Cuprite (Cu<sub>2</sub>O, C3) Structure': 'cuprite',
    'Heusler (L2<sub>1</sub>) Structure': 'Heusler',
    'Rock Salt (NaCl, B1) Structure': 'rock salt',
    'Face&ndash;Centered Cubic (Cu, A1) Structure': 'fcc',
    'Diamond (A4) Structure': 'diamond',
    'Body&ndash;Centered Cubic (W, A2) Structure': 'bcc',
}


def get_formula_hill(formula: str) -> str:
    '''
    Converts the given chemical formula into the Hill format.

    Args:
        formula: Original formula.

    Returns:
        Chemical formula in the Hill format.
    '''
    from nomad.atomutils import Formula
    return None if formula is None else Formula(formula).format('hill')


def get_formula_iupac(formula: str) -> str:
    '''
    Converts the given chemical formula into the IUPAC format.

    Args:
        formula: Original formula.

    Returns:
        Chemical formula in the IUPAC format.
    '''
    from nomad.atomutils import Formula
    return None if formula is None else Formula(formula).format('iupac')


def available_properties(root: MSection) -> List[str]:
    '''Returns a list of property names that are available in results.properties.

    Args:
        root: The metainfo section containing the properties

    Returns:
        List of property names that are present
    '''
    from nomad.utils import traverse_reversed
    available_property_names = {
        'electronic.band_structure_electronic.band_gap': 'electronic.band_structure_electronic.band_gap',
        'electronic.band_structure_electronic': 'band_structure_electronic',
        'electronic.dos_electronic': 'dos_electronic',
        'electronic.greens_functions_electronic': 'greens_functions_electronic',
        'vibrational.dos_phonon': 'dos_phonon',
        'vibrational.band_structure_phonon': 'band_structure_phonon',
        'vibrational.energy_free_helmholtz': 'energy_free_helmholtz',
        'vibrational.heat_capacity_constant_volume': 'heat_capacity_constant_volume',
        'thermodynamic.trajectory': 'trajectory',
        'structural.radial_distribution_function': 'radial_distribution_function',
        'dynamical.mean_squared_displacement': 'mean_squared_displacement',
        'structural.radius_of_gyration': 'radius_of_gyration',
        'geometry_optimization': 'geometry_optimization',
        'mechanical.bulk_modulus': 'bulk_modulus',
        'mechanical.shear_modulus': 'shear_modulus',
        'mechanical.energy_volume_curve': 'energy_volume_curve',
        'spectroscopic.spectra': 'spectra',
        'optoelectronic.solar_cell': 'solar_cell',
    }
    available_properties: List[str] = []
    for path, shortcut in available_property_names.items():
        for _ in traverse_reversed(root, path.split('.')):
            available_properties.append(shortcut)
            break
    return sorted(available_properties)


tokenizer_formula = get_tokenizer(r'[A-Z][a-z]?\d*')


class BandGapDeprecated(CalcBandGapDeprecated):
    label = Quantity(
        type=str,
        description='''
        Label to identify the band gap data, e.g. method employed.
        '''
    )
    index = CalcBandGapDeprecated.index.m_copy()
    index.m_annotations['elasticsearch'] = [Elasticsearch(material_entry_type)]
    value = CalcBandGapDeprecated.value.m_copy()
    value.m_annotations['elasticsearch'] = [Elasticsearch(material_entry_type)]
    type = CalcBandGapDeprecated.type.m_copy()
    type.m_annotations['elasticsearch'] = [Elasticsearch(material_entry_type)]

    label = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary label for GUI generating GUI artifacts.
        Can be left unpopulated.
        '''
    )


class BandGap(CalcBandGap):
    label = Quantity(
        type=str,
        description='''
        Label to identify the band gap data, e.g. the method employed.
        ''')
    index = CalcBandGap.index.m_copy()
    index.m_annotations['elasticsearch'] = [Elasticsearch(material_entry_type)]
    value = CalcBandGap.value.m_copy()
    value.m_annotations['elasticsearch'] = [Elasticsearch(material_entry_type)]
    type = CalcBandGap.type.m_copy()
    type.m_annotations['elasticsearch'] = [Elasticsearch(material_entry_type)]


class ElementalComposition(MSection):
    m_def = Section(
        description='''
        Section containing information about the concentration of an element,
        given by its atomic and mass fraction within the system or material.
        ''',
        label_quantity='element',
    )
    element = Quantity(
        type=MEnum(chemical_symbols[1:]),
        description='''
        The symbol of the element, e.g. 'Pb'.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion='simple')
        ],
    )
    mass = Quantity(
        type=np.float64,
        unit='kg',
        description='''
        The (average) mass of the element.
        ''',
    )
    atomic_fraction = Quantity(
        type=np.float64,
        description='''
        The atomic fraction of the element in the system it is contained within.
        Per definition a positive value less than or equal to 1.
        ''',
        a_elasticsearch=Elasticsearch(material_type),
    )
    mass_fraction = Quantity(
        type=np.float64,
        description='''
        The mass fraction of the element in the system it is contained within.
        Per definition a positive value less than or equal to 1.
        ''',
        a_elasticsearch=Elasticsearch(material_type),
    )


class LatticeParameters(MSection):
    m_def = Section(
        description='''
        Lattice parameters of a cell.
        ''',
    )
    a = Quantity(
        type=np.float64,
        unit='m',
        description='''
        Length of the first basis vector.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    b = Quantity(
        type=np.float64,
        unit='m',
        description='''
        Length of the second basis vector.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    c = Quantity(
        type=np.float64,
        unit='m',
        description='''
        Length of the third basis vector.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    alpha = Quantity(
        type=np.float64,
        unit='radian',
        description='''
        Angle between second and third basis vector.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    beta = Quantity(
        type=np.float64,
        unit='radian',
        description='''
        Angle between first and third basis vector.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    gamma = Quantity(
        type=np.float64,
        unit='radian',
        description='''
        Angle between first and second basis vector.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
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
        type=np.float64,
        description='''
        The free parameter x if present.
        '''
    )
    y = Quantity(
        type=np.float64,
        description='''
        The free parameter y if present.
        '''
    )
    z = Quantity(
        type=np.float64,
        description='''
        The free parameter z if present.
        '''
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
        type=np.float64,
        shape=[3, 3],
        unit='m',
        description='''
        The three lattice vectors in Cartesian coordinates.
        '''
    )
    cartesian_site_positions = Quantity(
        type=np.float64,
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
        type=np.float64,
        unit='m ** 3',
        description='''
        Volume of the cell.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    atomic_density = Quantity(
        type=np.float64,
        unit='1 / m ** 3',
        description='''
        Atomic density of the material (atoms/volume).'
        '''
    )
    mass_density = Quantity(
        type=np.float64,
        unit='kg / m ** 3',
        description='''
        Mass density of the material.
        '''
    )
    species = SubSection(sub_section=OptimadeSpecies.m_def, repeats=True)
    lattice_parameters = SubSection(sub_section=LatticeParameters.m_def)
    wyckoff_sets = SubSection(sub_section=WyckoffSet.m_def, repeats=True)


class Structures(MSection):
    m_def = Section(
        description='''
        Contains full atomistic representations of the material in different
        forms.
        ''',
    )
    structure_original = SubSection(
        description='''
        Contains a selected representative structure from the the original
        data.
        ''',
        sub_section=Structure.m_def,
        repeats=False,
    )
    structure_conventional = SubSection(
        description='''
        Contains the conventional structure that is derived from
        structure_original. This conventional stucture has been idealized and
        the conventions employed by spglib are used.
        ''',
        sub_section=Structure.m_def,
        repeats=False,
    )
    structure_primitive = SubSection(
        description='''
        Contains the primitive structure that is derived from
        structure_original. This primitive stucture has been idealized and the
        conventions employed by spglib are used.
        ''',
        sub_section=Structure.m_def,
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
        type=np.int32,
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
        type=np.int32,
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
        type=MEnum(sorted(list(set(structure_name_map.values())))),
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


# =============================================================================
# New topological data
class Cell(MSection):
    m_def = Section(
        description='''
        Properties of a unit cell.
        ''',
    )
    a = Quantity(
        type=np.float64,
        unit='m',
        description='''
        Length of the first basis vector.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    b = Quantity(
        type=np.float64,
        unit='m',
        description='''
        Length of the second basis vector.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    c = Quantity(
        type=np.float64,
        unit='m',
        description='''
        Length of the third basis vector.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    alpha = Quantity(
        type=np.float64,
        unit='radian',
        description='''
        Angle between second and third basis vector.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    beta = Quantity(
        type=np.float64,
        unit='radian',
        description='''
        Angle between first and third basis vector.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    gamma = Quantity(
        type=np.float64,
        unit='radian',
        description='''
        Angle between first and second basis vector.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    volume = Quantity(
        type=np.float64,
        unit='m ** 3',
        description='''
        Volume of the cell.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    atomic_density = Quantity(
        type=np.float64,
        unit='1 / m ** 3',
        description='''
        Atomic density of the material (atoms/volume).'
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    mass_density = Quantity(
        type=np.float64,
        unit='kg / m ** 3',
        description='''
        Mass density of the material.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )


class SymmetryNew(MSection):
    m_def = Section(
        description='''
        Section containing information about the symmetry properties of a
        conventional cell related to a system.
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
        type=np.int32,
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
        type=np.int32,
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
    choice = Quantity(
        type=str,
        shape=[],
        description='''
        String that specifies the centering, origin and basis vector settings of the 3D
        space group that defines the symmetry group of the simulated physical system (see
        section system). Values are as defined by spglib.
        ''')
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
    symmetry_method = Quantity(
        type=str,
        shape=[],
        description='''
        Identifies the source of the symmetry information contained within this
        section. If equal to 'spg_normalized' the information comes from a
        normalization step.
        '''
    )
    origin_shift = Quantity(
        type=np.float64,
        shape=[3],
        description='''
        Vector $\\mathbf{p}$ from the origin of the standardized system to the origin of
        the original system. Together with the matrix $\\mathbf{P}$, found in
        space_group_3D_transformation_matrix, the transformation between the standardized
        coordinates $\\mathbf{x}_s$ and original coordinates $\\mathbf{x}$ is then given
        by $\\mathbf{x}_s = \\mathbf{P} \\mathbf{x} + \\mathbf{p}$.
        '''
    )
    transformation_matrix = Quantity(
        type=np.float64,
        shape=[3, 3],
        description='''
        Matrix $\\mathbf{P}$ that is used to transform the standardized coordinates to the
        original coordinates. Together with the vector $\\mathbf{p}$, found in
        space_group_3D_origin_shift, the transformation between the standardized
        coordinates $\\mathbf{x}_s$ and original coordinates $\\mathbf{x}$ is then given by
        $\\mathbf{x}_s = \\mathbf{P} \\mathbf{x} + \\mathbf{p}$.
        '''
    )
    symmorphic = Quantity(
        type=bool,
        shape=[],
        description='''
        Specifies if the space group is symmorphic. Set to True if all
        translations are zero.
        '''
    )
    prototype_label_aflow = Quantity(
        type=str,
        shape=[],
        description='''
        AFLOW label of the prototype (see
        http://aflowlib.org/CrystalDatabase/prototype_index.html) identified on the basis
        of the space_group and normalized_wyckoff.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion='simple')
        ]
    )
    prototype_name = Quantity(
        type=MEnum(sorted(list(set(structure_name_map.values())))),
        description='''
        A common name for this prototypical structure, e.g. fcc, bcc.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion='default')
        ],
    )
    wyckoff_sets = SubSection(sub_section=WyckoffSet.m_def, repeats=True)


class Relation(MSection):
    '''
    Contains information about the relation between two different systems.
    '''
    m_def = Section(validate=False)
    type = Quantity(
        type=MEnum('root', 'subsystem', 'group', 'primitive_cell', 'conventional_cell'),
        description='''
        The type of relation between a system and it's parent.

        | Value | Description |
        | --------- | ----------------------- |
        | `'root'` | System representing the entire structure, has no parent system. |
        | `'subsystem'` | A single logical entity extracted from the parent system. |
        | `'group'` | A logical group of subsystems within the parent, e.g. a group of molecules in MD. |
        | `'primitive_cell'` | The conventional cell from which the parent is constructed from. |
        | `'conventional_cell'` | The primitive cell from which the parent is constructed from. |
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion='default')
        ]
    )


class Coordination(MSection):
    '''
   Coordination number of an element, which represents the number
   of atoms directly bonded to the element.
    '''
    element = Quantity(
        type=MEnum(chemical_symbols),
        description='''
        Chemical symbol of element, whose coordination number is being determined.
        '''
    )
    coordination_number = Quantity(
        type=int,
        description='''
        The number of neighbours directly connected to an atom
        '''
    )


class System(MSection):
    '''Describes a physical system as identified in an entry. Can be a part of a
    larger structural hierarchy.
    '''
    m_def = Section(
        description='''
        Describes a a structural part that has been identified within the entry.
        May be related to other systems.
        '''
    )
    system_id = Quantity(
        type=str,
        description='''
        That path of this section within the metainfo that is used as a unique
        identifier.
        ''',
        a_elasticsearch=Elasticsearch(material_type)
    )
    label = Quantity(
        type=str,
        description='''
        Descriptive label that identifies this structural part.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion='default')
        ]
    )
    method = Quantity(
        type=MEnum('parser', 'user', 'matid', 'porosity'),
        description='''
        The method used for identifying this system.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion='default')
        ]
    )
    description = Quantity(
        type=str,
        description='''
        A short description about this part of the topology.
        ''',
        a_elasticsearch=Elasticsearch(material_type)
    )
    material_id = Quantity(
        type=str,
        description='''
        A fixed length, unique material identifier in the form of a hash
        digest.
        ''',
        a_elasticsearch=Elasticsearch(material_type)
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
        type=MEnum(structure_classes + ['group', 'molecule', 'monomer']),
        description='''
        Structural class determined from the atomic structure.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion='default')
        ],
    )
    dimensionality = Quantity(
        type=MEnum(['0D', '1D', '2D', '3D']),
        description='''
        Dimensionality of the system. For atomistic systems this is
        automatically evaluated by using the topology-scaling algorithm:
        https://doi.org/10.1103/PhysRevLett.118.106101.

        | Value | Description |
        | --------- | ----------------------- |
        | `'0D'` | Not connected periodically |
        | `'1D'` | Periodically connected in one dimension |
        | `'2D'` | Periodically connected in two dimensions |
        | `'3D'` | Periodically connected in three dimensions |
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion='default')
        ],
    )
    building_block = Quantity(
        type=MEnum(['surface', '2D material', 'molecule', 'monomer']),
        description='''
        More exact classification for this system, i.e. the type of "building
        block" it represents.

        | Value | Description |
        | --------- | ----------------------- |
        | `'surface'` | Structure built from a unit cell that repeats periodically in two directions and at least twice, but not infinitely in a third direction. |
        | `'2D material'` | Structure built from a unit cell that repeats periodically in two directions and only once in a third direction. |
        | `'molecule'` | Molecule defined in the force-field topology |
        | `'monomer'` | Monomer defined in the force-field topology |
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
            Alphabetically sorted chemical formula with reduced integer chemical
            proportion numbers. The proportion number is omitted if it is 1.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion=tokenizer_formula)
        ],
    )
    chemical_formula_hill = Quantity(
        type=str,
        description='''
            The chemical formula for a structure in Hill form with element
            symbols followed by non-reduced integer chemical proportion numbers.
            The proportion number is omitted if it is 1.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type, normalizer=get_formula_hill),
            Elasticsearch(suggestion=tokenizer_formula)
        ],
    )
    chemical_formula_iupac = Quantity(
        type=str,
        description='''
            Formula where the elements are ordered using a formal list loosely
            based on electronegativity as defined in the IUPAC nomenclature of
            inorganic chemistry (2005). Contains reduced integer chemical
            proportion numbers where the proportion number is omitted if it is
            1.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type, normalizer=get_formula_iupac),
            Elasticsearch(suggestion=tokenizer_formula)
        ],
    )
    chemical_formula_anonymous = Quantity(
        type=str,
        description='''
            Formula with the elements ordered by their reduced integer chemical
            proportion number, and the chemical species replaced by
            alphabetically ordered letters. The proportion number is omitted if
            it is 1. E.g.  H2O becomes A2B and H2O2 becomes AB. The letters are
            drawn from the english alphabet that may be extended by increasing
            the number of letters, e.g. A, B, ..., Z, Aa, Ab and so on. This
            definition is in line with the similarly named OPTIMADE definition.
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
            Alphabetically sorted chemical formula with reduced integer chemical
            proportion numbers. The proportion number is omitted if it is 1.
        ''',
        a_elasticsearch=Elasticsearch(material_type, mapping=Text(multi=True)),
    )
    parent_system = Quantity(
        type=str,
        description='''
        Reference to the parent system.
        ''',
        a_elasticsearch=Elasticsearch(material_type)
    )
    child_systems = Quantity(
        type=str,
        shape=['*'],
        description='''
        References to the child systems.
        ''',
        a_elasticsearch=Elasticsearch(material_type)
    )
    atomic_fraction = Quantity(
        type=np.float64,
        description='''
        The atomic fraction of this system in the full structure it is contained in.
        Per definition a positive value less than or equal to 1.
        ''',
        a_elasticsearch=Elasticsearch(material_type),
    )
    mass_fraction = Quantity(
        type=np.float64,
        description='''
        The mass fraction of this system in the full structure it is contained within.
        Per definition a positive value less than or equal to 1.
        ''',
        a_elasticsearch=Elasticsearch(material_type),
    )
    atoms = SubSection(
        description='''
        The atomistic structure that is associated with this
        system.
        ''',
        sub_section=Atoms.m_def,
        repeats=False
    )
    atoms_ref = Quantity(
        type=Atoms,
        description='''
        Reference to an atomistic structure that is associated with this
        system.
        ''',
    )
    n_atoms = Quantity(
        type=int,
        shape=[],
        description='''
        The total number of species (atoms, particles) in the system.
        ''',
        a_elasticsearch=Elasticsearch(material_type)
    )
    indices = Quantity(
        type=np.int64,
        shape=['*', '*'],
        description='''
        Indices of the atoms belonging to this group. These indices refer to the
        system specified in atoms_ref. Each row represents a new instance.
        '''
    )
    elemental_composition = SubSection(
        sub_section=ElementalComposition.m_def,
        repeats=True,
        a_elasticsearch=Elasticsearch(material_type, nested=True)
    )
    system_relation = SubSection(sub_section=Relation.m_def, repeats=False)
    cell = SubSection(sub_section=Cell.m_def, repeats=False)
    symmetry = SubSection(sub_section=SymmetryNew.m_def, repeats=False)

    sbu_type = Quantity(
        type=str,
        shape=[],
        description='''
        The topological representation of the metal secondary building units (sbus).
        The shape of most metal sbus are well defined and form the basis of most
         popular MOFs. The most common example is the paddlewheel, rodlike mofs,
         irmofs, uio66
         ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion='default')
        ],
    )
    largest_cavity_diameter = Quantity(
        type=np.float64,
        alias=['lcd', 'largest_included_sphere'],
        unit='m',
        description='''
        The largest cavity diameter is the largest sphere that can be inserted in a porous
        system without overlapping with any of the atoms in the system.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
        ]
    )
    pore_limiting_diameter = Quantity(
        type=np.float64,
        alias=['pld', 'free_sphere'],
        unit='m',
        description='''
        The pore limiting diameter is the largest sphere that can freely
        diffuse through the porous network without overlapping with any of the
        atoms in the system.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
        ]
    )
    largest_included_sphere_along_free_sphere_path = Quantity(
        type=np.float64,
        alias=['lfpd'],
        unit='m',
        description='''
        The largest included sphere along free sphere path is
        largest sphere that can be inserted in the pore.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
        ]
    )
    accessible_surface_area = Quantity(
        type=np.float64,
        alias=['asa', 'sasa', 'solvent_accessible_surface_area'],
        unit='m ** 2',
        description='''
        The surface area accessible is the area that is accessible to guest molecules
        in a porous system. It is generally considered to be the entire surface area
        that can be spanned by a probe of a specific radius. In NOMAD, by default we use
        a probe that has a radius of 1.86 Angstrom, which correspond to the
        covalent radii of nitrogen gas. For biomolecular system, a radii of
        1.4 Angstrom can be used, which correspond to the covalent radii
        of water.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
        ]
    )
    accessible_volume = Quantity(
        type=np.float64,
        alias=['common_solvent_accessible_volume', 'csav',
               'solvent_accessible_volume', 'sav'],
        unit='m ** 3',
        description='''
        Volume of unoccupied space in a system that can be accessible to
        guest molecules, like solvents.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
        ]
    )

    metal_coordination = SubSection(sub_section=Coordination.m_def, repeats=True)

    void_fraction = Quantity(
        type=np.float64,
        alias=['void_ratio'],
        description='''
        Ratio of the the volume of the unoccupied space in the system
        to the volume of the entire system. It is a good proxy to
        determine how porous a system is. Highly porous systems
        often have a larger void fraction, meanwhile compact or dense
        systems have smaller void fractions.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
        ]
    )
    n_channels = Quantity(
        type=int,
        shape=[],
        description='''
        Number of channels present in the porous system, which correspond to the number of
        pores within the system.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
        ]
    )
    sbu_coordination_number = Quantity(
        type=int,
        alias=['n_point_of_extensions'],
        description='''
        The number of connecting point in the secondary building units(sbu), which corresponds to
        the to the number of point of extension in the secondary building unit. Some common
        terminology include
        1 : monotopic
        2 : ditopic
        3 : tritopic
        4 : tetratopic
        5 : pentatopic
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
        ]
    )
# =============================================================================


class Material(MSection):
    m_def = Section(
        description='''
        Section containing information on the material composition and structure.
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
        Structural class determined from the atomic structure.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion='default')
        ],
    )
    dimensionality = System.dimensionality.m_copy()
    building_block = System.building_block.m_copy()
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
            Alphabetically sorted chemical formula with reduced integer chemical
            proportion numbers. The proportion number is omitted if it is 1.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion=tokenizer_formula)
        ],
    )
    chemical_formula_hill = Quantity(
        type=str,
        description='''
            The chemical formula for a structure in Hill form with element
            symbols followed by non-reduced integer chemical proportion numbers.
            The proportion number is omitted if it is 1.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type, normalizer=get_formula_hill),
            Elasticsearch(suggestion=tokenizer_formula)
        ],
    )
    chemical_formula_iupac = Quantity(
        type=str,
        description='''
            Formula where the elements are ordered using a formal list loosely
            based on electronegativity as defined in the IUPAC nomenclature of
            inorganic chemistry (2005). Contains reduced integer chemical
            proportion numbers where the proportion number is omitted if it is
            1.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type, normalizer=get_formula_iupac),
            Elasticsearch(suggestion=tokenizer_formula)
        ],
    )
    chemical_formula_anonymous = Quantity(
        type=str,
        description='''
            Formula with the elements ordered by their reduced integer chemical
            proportion number, and the chemical species replaced by
            alphabetically ordered letters. The proportion number is omitted if
            it is 1. E.g.  H2O becomes A2B and H2O2 becomes AB. The letters are
            drawn from the english alphabet that may be extended by increasing
            the number of letters, e.g. A, B, ..., Z, Aa, Ab and so on. This
            definition is in line with the similarly named OPTIMADE definition.
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
            Alphabetically sorted chemical formula with reduced integer chemical
            proportion numbers. The proportion number is omitted if it is 1.
        ''',
        a_elasticsearch=Elasticsearch(material_type, mapping=Text(multi=True)),
    )
    elemental_composition = SubSection(
        sub_section=ElementalComposition.m_def,
        repeats=True,
        a_elasticsearch=Elasticsearch(material_type, nested=True)
    )
    symmetry = SubSection(sub_section=Symmetry.m_def, repeats=False)
    topology = SubSection(
        sub_section=System.m_def,
        repeats=True,
        a_elasticsearch=Elasticsearch(material_type, nested=True)
    )


class HubbardKanamoriModel(MSection):
    '''
    Setup of the Hubbard model used in DFT+U
    '''

    m_def = Section(validate=False)

    atom_label = AtomParameters.label.m_copy()
    orbital = HubbardKanamori.orbital.m_copy()
    u_effective = HubbardKanamori.u_effective.m_copy()
    u_effective.m_annotations['elasticsearch'] = [Elasticsearch(material_entry_type)]
    u = HubbardKanamori.u.m_copy()
    u.m_annotations['elasticsearch'] = [Elasticsearch(material_entry_type)]
    j = HubbardKanamori.j.m_copy()
    j.m_annotations['elasticsearch'] = [Elasticsearch(material_entry_type)]
    double_counting_correction = HubbardKanamori.double_counting_correction.m_copy()


class DFT(MSection):
    m_def = Section(
        description='''
        Methodology for a DFT calculation.
        '''
    )
    basis_set_type = Quantity(  # TODO: deprecate after reparsing
        type=MEnum(basis_set_types),
        default=unavailable,
        description='The used basis set functions.',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )
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

    jacobs_ladder = Quantity(
        type=MEnum(list(xc_treatments.values()) + [unavailable, not_processed]),
        default=not_processed,
        description='''Functional classification in line with Jacob\'s Ladder.
        For more information, see https://doi.org/10.1063/1.1390175 (original paper);
        https://doi.org/10.1103/PhysRevLett.91.146401 (meta-GGA);
        and https://doi.org/10.1063/1.1904565 (hyper-GGA).''',
        a_elasticsearch=Elasticsearch(material_entry_type, default_aggregation_size=100),
    )
    xc_functional_type = Quantity(
        type=jacobs_ladder.type,
        default=jacobs_ladder.default,
        description=jacobs_ladder.description,
        a_elasticsearch=Elasticsearch(material_entry_type, default_aggregation_size=100),
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
    exact_exchange_mixing_factor = Quantity(
        type=np.float64,
        description='Amount of exact exchange mixed in with the XC functional (value range = [0,1]).',
        a_elasticsearch=Elasticsearch(material_entry_type)
    )
    hubbard_kanamori_model = SubSection(
        sub_section=HubbardKanamoriModel.m_def, repeats=True,
        a_elasticsearch=[Elasticsearch(material_entry_type, nested=True)])


class Projection(MSection):
    m_def = Section(
        description='''
        Methodology for a Projection calculation.
        '''
    )
    type = Quantity(
        type=MEnum(['wannier', 'slater_koster', 'custom']),
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
        description='''
        Projection type for the virtual orbitals: Wannier or Slater-Koster formalisms, or a
        custom tight-binding model.
        '''
    )
    localization_type = Quantity(
        type=MEnum(['single_shot', 'maximally_localized']),
        description='''
        Localization type of the virtual Wannier orbitals.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )


class ExcitedStateMethodology(MSection):
    m_def = Section(
        description='''
        Methodology for a Excited-State calculation.
        '''
    )
    type = XSMethod.type.m_copy(
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ]
    )
    basis_set_type = Quantity(
        type=MEnum(basis_set_types),
        description='The used basis set functions.',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )
    starting_point_type = Quantity(
        type=MEnum(list(xc_treatments_extended.values()) + [unavailable, not_processed]),
        description='The libXC based xc functional classification used in the starting point DFT simulation.',
        a_elasticsearch=Elasticsearch(material_entry_type)
    )
    starting_point_names = Quantity(
        type=str,
        shape=['*'],
        description='The list of libXC functional names that where used in this entry.',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ]
    )


class GW(ExcitedStateMethodology):
    m_def = Section(
        description='''
        Methodology for a GW calculation.
        '''
    )
    type = GWMethod.type.m_copy(
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ]
    )


class BSE(ExcitedStateMethodology):
    m_def = Section(
        description='''
        Methodology for a BSE calculation.
        '''
    )
    type = BSEMethod.type.m_copy(
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ]
    )
    solver = BSEMethod.solver.m_copy(
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ]
    )
    gw_type = Quantity(
        type=MEnum(GWMethod.type.type._list),  # TODO solve conflict between BSE.gw_type and GW.type when using GWMethod.type.m_copy()
        description=GWMethod.type.description,
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ]
    )


class DMFT(MSection):
    m_def = Section(
        description='''
        Methodology for a DMFT calculation.
        '''
    )
    impurity_solver_type = DMFTMethod.impurity_solver.m_copy(
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ]
    )
    inverse_temperature = DMFTMethod.inverse_temperature.m_copy(
        a_elasticsearch=[
            Elasticsearch(material_entry_type)
        ]
    )
    magnetic_state = DMFTMethod.magnetic_state.m_copy()
    magnetic_state.description = 'Magnetic state in which the DMFT calculation is done.'
    magnetic_state.m_annotations['elasticsearch'] = [
        Elasticsearch(material_entry_type),
        Elasticsearch(suggestion='default')
    ]
    u = HubbardKanamori.u.m_copy(
        a_elasticsearch=[
            Elasticsearch(material_entry_type)
        ]
    )
    jh = HubbardKanamori.jh.m_copy(
        a_elasticsearch=[
            Elasticsearch(material_entry_type)
        ]
    )
    analytical_continuation = Quantity(
        type=MEnum('Pade', 'MaxEnt', 'SVD', ''),
        shape=[],
        description='''
        Analytical continuation used to continuate the imaginary space Green's functions into
        the real frequencies space.

        | Name           | Description         | Reference                        |

        | -------------- | ------------------- | -------------------------------- |

        | `'Pade'` | Pade's approximant  | https://www.sciencedirect.com/science/article/pii/0021999173901277?via%3Dihub |

        | `'MaxEnt'` | Maximum Entropy method | https://journals.aps.org/prb/abstract/10.1103/PhysRevB.41.2380 |

        | `'SVD'` | Singular value decomposition | https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.75.517 |

        | `'Stochastic'` | Stochastic method | https://journals.aps.org/prb/abstract/10.1103/PhysRevB.57.10287 |
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )


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


class Precision(MSection):
    m_def = Section(
        description='''
        Contains parameters for controlling or evaluating the convergence of the electronic structure.
        '''
    )
    k_line_density = Quantity(
        type=np.float64,
        shape=[],
        unit='m',
        description='''
        Amount of sampled k-points per unit reciprocal length along each axis.
        Contains the least precise density out of all axes.
        Should only be compared between calulations of similar dimensionality.
        ''',
        a_elasticsearch=[Elasticsearch(material_entry_type)],
    )
    native_tier = BasisSetContainer.native_tier.m_copy(a_elasticsearch=[
        Elasticsearch(material_entry_type)])
    basis_set = BasisSetContainer.type.m_copy(a_elasticsearch=[
        Elasticsearch(material_entry_type),
        Elasticsearch(suggestion='default'),
    ])
    planewave_cutoff = BasisSet.cutoff.m_copy(a_elasticsearch=[  # TODO: set better names?
        Elasticsearch(material_entry_type)
    ])
    apw_cutoff = BasisSet.cutoff_fractional.m_copy(a_elasticsearch=[  # TODO: set better names?
        Elasticsearch(material_entry_type)
    ])


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
    projection = SubSection(sub_section=Projection.m_def, repeats=False)
    gw = SubSection(sub_section=GW.m_def, repeats=False)
    bse = SubSection(sub_section=BSE.m_def, repeats=False)
    dmft = SubSection(sub_section=DMFT.m_def, repeats=False)
    quantum_cms = SubSection(sub_section=QuantumCMS.m_def, repeats=False)
    precision = SubSection(sub_section=Precision.m_def, repeats=False)


class XRDMethod(MSection):
    m_def = Section(
        description='''
        Methodology for an X-Ray Diffraction measurement.
        '''
    )
    diffraction_method_name = Quantity(
        type=MEnum(
            [
                "Powder X-Ray Diffraction (PXRD)",
                "Single Crystal X-Ray Diffraction (SCXRD)",
                "High-Resolution X-Ray Diffraction (HRXRD)",
                "Small-Angle X-Ray Scattering (SAXS)",
                "X-Ray Reflectivity (XRR)",
                "Grazing Incidence X-Ray Diffraction (GIXRD)",
                config.services.unavailable_value
            ]
        ),
        description='''
        The diffraction method used to obtain the diffraction pattern.
        | X-Ray Diffraction Method                                   | Description                                                                                                                                                                                                 |
        |------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
        | **Powder X-Ray Diffraction (PXRD)**                        | The term "powder" refers more to the random orientation of small crystallites than to the physical form of the sample. Can be used with non-powder samples if they present random crystallite orientations. |
        | **Single Crystal X-Ray Diffraction (SCXRD)**               | Used for determining the atomic structure of a single crystal.                                                                                                                                              |
        | **High-Resolution X-Ray Diffraction (HRXRD)**              | A technique typically used for detailed characterization of epitaxial thin films using precise diffraction measurements.                                                                                    |
        | **Small-Angle X-Ray Scattering (SAXS)**                    | Used for studying nanostructures in the size range of 1-100 nm. Provides information on particle size, shape, and distribution.                                                                             |
        | **X-Ray Reflectivity (XRR)**                               | Used to study thin film layers, interfaces, and multilayers. Provides info on film thickness, density, and roughness.                                                                                       |
        | **Grazing Incidence X-Ray Diffraction (GIXRD)**            | Primarily used for the analysis of thin films with the incident beam at a fixed shallow angle.                                                                                                              |
        ''',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )


class MeasurementMethod(MSection):
    m_def = Section(
        description='''
        Contains method details for a measurement entry.
        '''
    )
    xrd = SubSection(sub_section=XRDMethod.m_def, repeats=False)


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
        type=MEnum(
            [
                'DFT', 'Projection', 'GW', 'DMFT', 'CoreHole', 'BSE', 'EELS', 'XPS',
                'XRD', config.services.unavailable_value
            ]
        ),
        description='''
        Common name for the used method.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )

    # TODO error in registering quantity if Workflow.name.m_copy()
    workflow_name = Quantity(
        type=str
    )
    workflow_name.m_annotations['elasticsearch'] = [
        Elasticsearch(material_entry_type),
        Elasticsearch(suggestion='default')
    ]
    simulation = SubSection(sub_section=Simulation.m_def, repeats=False)
    measurement = SubSection(sub_section=MeasurementMethod.m_def, repeats=False)


class MolecularDynamics(MSection):
    m_def = Section(
        description="""
        Methodology for molecular dynamics.
        """,
    )
    time_step = MolecularDynamicsMethod.integration_timestep.m_copy()
    time_step.m_annotations["elasticsearch"] = Elasticsearch(material_entry_type)
    ensemble_type = MolecularDynamicsMethod.thermodynamic_ensemble.m_copy()
    ensemble_type.m_annotations["elasticsearch"] = Elasticsearch(material_entry_type)


class MDProvenance(ProvenanceTracker):
    m_def = Section(
        description="""
        Contains provenance information for properties derived from molecular
        dynamics simulations.
        """,
    )
    molecular_dynamics = SubSection(sub_section=MolecularDynamics.m_def, repeats=False)


class MDPropertySection(PropertySection):
    m_def = Section(
        description="""
        Base class for referring to molecular dynamics properties.
        """,
    )
    provenance = SubSection(sub_section=MDProvenance.m_def)


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
        Density of states (DOS) values for the entire system and all species.
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
    label = Quantity(
        type=str,
        description='''
        Label to identify the DOS data, e.g. the method employed.
        ''')

    spin_polarized = Quantity(
        type=bool,
        description='''
        Whether the DOS is spin-polarized, i.e. is contains channels for both
        spin values.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    band_gap = SubSection(
        sub_section=BandGapDeprecated.m_def,
        repeats=True,
        a_elasticsearch=Elasticsearch(material_entry_type, nested=True)
    )
    energy_fermi = Quantity(
        type=np.float64,
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
    label = Quantity(
        type=str,
        description='''
        Label to identify the bandstructure data, e.g. the method employed.
        ''')
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
        sub_section=BandGapDeprecated.m_def,
        repeats=True,
        a_elasticsearch=Elasticsearch(material_entry_type, nested=True)
    )
    energy_fermi = Quantity(
        type=np.float64,
        unit='joule',
        shape=[],
        description='''
        Fermi energy.
        '''
    )


class GreensFunctionsElectronic(MSection):
    m_def = Section(
        description='''
        Base class for Green's functions information.
        ''',
    )
    type = GreensFunctionsCalculation.type.m_copy()
    label = Quantity(
        type=str,
        description='''
        Label to identify the Greens functions data, e.g. the method employed.
        ''')
    tau = Quantity(
        type=GreensFunctionsCalculation.tau,
        description='''
        Array containing the set of discrete imaginary times.
        ''',
    )
    matsubara_freq = Quantity(
        type=GreensFunctionsCalculation.matsubara_freq,
        description='''
        Array containing the set of discrete imaginary (Matsubara) frequencies.
        ''',
    )
    frequencies = Quantity(
        type=GreensFunctionsCalculation.frequencies,
        description='''
        Array containing the set of discrete real frequencies.
        ''',
    )
    greens_function_tau = Quantity(
        type=GreensFunctionsCalculation.greens_function_tau,
        description='''
        Green's functions values in imaginary times.
        ''',
    )
    greens_function_iw = Quantity(
        type=GreensFunctionsCalculation.greens_function_iw,
        description='''
        Green's functions values in imaginary (Matsubara) frequencies.
        ''',
    )
    self_energy_iw = Quantity(
        type=GreensFunctionsCalculation.self_energy_iw,
        description='''
        Self-energy values in imaginary (Matsubara) frequencies.
        ''',
    )
    greens_function_freq = Quantity(
        type=GreensFunctionsCalculation.greens_function_freq,
        description='''
        Green's function values in real frequencies.
        ''',
    )
    self_energy_freq = Quantity(
        type=GreensFunctionsCalculation.self_energy_freq,
        description='''
        Self-energy values in real frequencies.
        ''',
    )
    hybridization_function_freq = Quantity(
        type=GreensFunctionsCalculation.hybridization_function_freq,
        description='''
        Hybridization function values in real frequencies.
        ''',
    )
    orbital_occupations = Quantity(
        type=GreensFunctionsCalculation.orbital_occupations,
        description='''
        Orbital occupation per correlated atom in the unit cell and per spin.
        ''',
    )
    quasiparticle_weights = Quantity(
        type=GreensFunctionsCalculation.quasiparticle_weights,
        description='''
        Quasiparticle weights of each orbital per site and spin. Calculated from:
            Z = inv(1.0 - d [Re Sigma] / dw at w=0)
        it takes values ∈ [0.0, 1.0], being Z=1 non-correlated, and Z=0 in a Mott state.
        ''',
    )
    chemical_potential = Quantity(
        type=GreensFunctionsCalculation.chemical_potential,
        description='''
        Chemical potential.
        ''',
    )


class HeatCapacityConstantVolume(MSection):
    m_def = Section(
        description='''
        Contains the values of the specific (per mass) and isochoric (constant
        volume) heat capacity at different temperatures.
        '''
    )
    heat_capacities = Quantity(
        type=ThermodynamicsResults.heat_capacity_c_v,
        shape=[],
        description='''
        Specific heat capacity values at constant volume.
        ''',
    )
    temperatures = Quantity(
        type=ThermodynamicsResults.temperature,
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
        type=ThermodynamicsResults.vibrational_free_energy_at_constant_volume,
        shape=[],
        description='''
        The Helmholtz free energies per atom at constant volume.
        ''',
    )
    temperatures = Quantity(
        type=ThermodynamicsResults.temperature,
        description='''
        The temperatures at which Helmholtz free energies are calculated.
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
            'birch_murnaghan',
            'murnaghan',
        ),
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )
    volumes = Quantity(type=EquationOfStateResults.volumes)
    energies_raw = Quantity(type=EquationOfStateResults.energies)
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
            'birch_murnaghan',
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
        type=np.float64,
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
        type=np.float64,
        description='Shear modulus value.',
        unit='pascal',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )


class GeometryOptimization(MSection):
    m_def = Section(
        description='''
        Geometry optimization results and settings.
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
        type=GeometryOptimizationResults.energies,
        description='''
        List of energy_total values gathered from the single configuration
        calculations that are a part of the optimization trajectory.
        ''',
    )
    system_optimized = Quantity(
        type=SystemRun,
        description='''
        Contains the optimized geometry that is the result of a geometry optimization.
        ''',
    )
    type = MGeometryOptimization.name.m_copy()
    convergence_tolerance_energy_difference = GeometryOptimizationMethod.convergence_tolerance_energy_difference.m_copy()
    convergence_tolerance_energy_difference.m_annotations['elasticsearch'] = Elasticsearch(material_entry_type)
    convergence_tolerance_force_maximum = GeometryOptimizationMethod.convergence_tolerance_force_maximum.m_copy()
    convergence_tolerance_force_maximum.m_annotations['elasticsearch'] = Elasticsearch(material_entry_type)
    final_force_maximum = GeometryOptimizationResults.final_force_maximum.m_copy()
    final_force_maximum.m_annotations['elasticsearch'] = Elasticsearch(material_entry_type)
    final_energy_difference = GeometryOptimizationResults.final_energy_difference.m_copy()
    final_energy_difference.m_annotations['elasticsearch'] = Elasticsearch(material_entry_type)
    final_displacement_maximum = GeometryOptimizationResults.final_displacement_maximum.m_copy()
    final_displacement_maximum.m_annotations['elasticsearch'] = Elasticsearch(material_entry_type)


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
    band_gap = SubSection(
        sub_section=BandGap.m_def,
        repeats=True,
        a_elasticsearch=Elasticsearch(material_entry_type, nested=True)
    )
    dos_electronic = SubSection(sub_section=DOSElectronic.m_def, repeats=True)
    band_structure_electronic = SubSection(sub_section=BandStructureElectronic.m_def, repeats=True)
    greens_functions_electronic = SubSection(sub_section=GreensFunctionsElectronic.m_def, repeats=True)


class QuantityDynamic(MSection):
    m_def = Section(
        description="""
        Contains the values for a quantity at different times.
        """
    )
    time = Quantity(
        type=np.float64,
        shape=[],
        unit='second',
        description="""
        The explicit times at which the values are evaluated. Provide either
        this or time_step and time_start.
        """,
    )
    time_step = Quantity(
        type=np.float64,
        unit='second',
        description="""
        The time step between successive evaluations. Provide either
        this and time_start or the explicit times.
        """,
    )
    time_start = Quantity(
        type=np.float64,
        unit='second',
        description="""
        The time at which the evaluation started. Provide either this and
        time_step or the explicit times.
        """,
    )


class VolumeDynamic(QuantityDynamic):
    m_def = Section(
        description="""
        Contains volume values evaluated at different times.
        """
    )
    value = Quantity(
        type=np.float64,
        shape=[],
        unit="m ** 3",
        description="""
        The volume values.
        """,
    )


class PressureDynamic(QuantityDynamic):
    m_def = Section(
        description="""
        Contains pressure values evaluated at different times.
        """
    )
    value = Quantity(
        type=np.float64,
        shape=[],
        unit="pascal",
        description="""
        The pressure values.
        """,
    )


class TemperatureDynamic(QuantityDynamic):
    m_def = Section(
        description="""
        Contains temperature values evaluated at different times.
        """
    )
    value = Quantity(
        type=np.float64,
        shape=[],
        unit="kelvin",
        description="""
        The temperature value.
        """,
    )


class EnergyDynamic(QuantityDynamic):
    m_def = Section(
        description="""
        Contains energy values evaluated at different times.
        """
    )
    value = Quantity(
        type=np.float64,
        shape=[],
        unit="joule",
        description="""
        The energy values.
        """,
    )


class Trajectory(MDPropertySection):
    m_def = Section(
        description='''
        Thermodynamic properties reported for an ensemble evolving in time.
        ''',
    )
    temperature = SubSection(sub_section=TemperatureDynamic.m_def, repeats=False)
    pressure = SubSection(sub_section=PressureDynamic.m_def, repeats=False)
    volume = SubSection(sub_section=VolumeDynamic.m_def, repeats=False)
    energy_potential = SubSection(sub_section=EnergyDynamic.m_def, repeats=False)
    available_properties = Quantity(
        type=MEnum(
            'temperature',
            'pressure',
            'volume',
            'energy_potential'
        ),
        shape=['0..*'],
        description='Subset of the property names that are present in this trajectory.',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )


class ThermodynamicProperties(MSection):
    m_def = Section(
        description='''
        Thermodynamic properties.
        ''',
    )
    trajectory = SubSection(
        sub_section=Trajectory.m_def,
        repeats=True,
        a_elasticsearch=Elasticsearch(material_entry_type, nested=True)
    )


class RadiusOfGyration(QuantityDynamic, MDPropertySection):
    m_def = Section(
        description='''
        Contains Radius of Gyration values as a trajectory.
        ''',
    )
    kind = RgCalculation.kind.m_copy()
    kind.m_annotations['elasticsearch'] = [
        Elasticsearch(material_entry_type),
        Elasticsearch(suggestion='default')
    ]
    label = RadiusOfGyrationValues.label.m_copy()
    label.m_annotations['elasticsearch'] = [
        Elasticsearch(material_entry_type),
        Elasticsearch(suggestion='default')
    ]
    atomsgroup_ref = RadiusOfGyrationValues.atomsgroup_ref.m_copy()
    value = RadiusOfGyrationValues.value.m_copy()


class RadialDistributionFunction(MDPropertySection):
    m_def = Section(
        description='''
        Radial distribution function.
        ''',
    )
    type = RDFWorkflow.type.m_copy()
    type.m_annotations['elasticsearch'] = [
        Elasticsearch(material_entry_type),
        Elasticsearch(suggestion='default')
    ]
    label = RadialDistributionFunctionValues.label.m_copy()
    label.m_annotations['elasticsearch'] = [
        Elasticsearch(material_entry_type),
        Elasticsearch(suggestion='default')
    ]
    bins = RadialDistributionFunctionValues.bins.m_copy()
    n_bins = RadialDistributionFunctionValues.n_bins.m_copy()
    value = RadialDistributionFunctionValues.value.m_copy()
    frame_start = RadialDistributionFunctionValues.frame_start.m_copy()
    frame_end = RadialDistributionFunctionValues.frame_end.m_copy()


class DiffractionPattern(MSection):
    m_def = Section(
        description='''
        Diffraction pattern.
        ''',
    )
    incident_beam_wavelength = Quantity(
        links=['https://manual.nexusformat.org/classes/base_classes/NXbeam.html#nxbeam-incident-wavelength-field'],
        description='''
        The wavelength of the incident beam.
        ''',
        type=np.float64,
        unit='m',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    two_theta_angles = Quantity(
        # type=Dos.energies,
        type=np.float64,  # TODO convert to reference when schema is in place
        unit='degree',
        shape=['*'],
        description='''
        Array containing the set of 2-theta angles.
        ''',
    )
    intensity = Quantity(
        # type=DosValues,
        type=np.float64,  # TODO convert to reference when schema is in place
        shape=['*'],
        description='''
        Array containing the set of intensities.
        ''',
    )
    q_vector = Quantity(
        type=np.dtype(np.float64),
        shape=['*'],
        unit='meter**(-1)',
        description='The scattering vector *Q*.')


class StructuralProperties(MSection):
    m_def = Section(
        description='''
        Structural properties.
        ''',
    )
    radial_distribution_function = SubSection(
        sub_section=RadialDistributionFunction.m_def,
        repeats=True,
        a_elasticsearch=Elasticsearch(material_entry_type, nested=True)
    )
    radius_of_gyration = SubSection(
        sub_section=RadiusOfGyration.m_def,
        repeats=True,
        a_elasticsearch=Elasticsearch(material_entry_type, nested=True)
    )
    diffraction_pattern = SubSection(
        sub_section=DiffractionPattern.m_def,
        repeats=True,
        a_elasticsearch=Elasticsearch(material_entry_type, nested=True)
    )


class MeanSquaredDisplacement(MDPropertySection):
    m_def = Section(
        description='''
        Mean Squared Displacements.
        ''',
    )
    type = MSDWorkflow.type.m_copy()
    type.m_annotations['elasticsearch'] = [
        Elasticsearch(material_entry_type),
        Elasticsearch(suggestion='default')
    ]
    direction = MSDWorkflow.direction.m_copy()
    error_type = MSDWorkflow.error_type.m_copy()
    label = MeanSquaredDisplacementValues.label.m_copy()
    label.m_annotations['elasticsearch'] = [
        Elasticsearch(material_entry_type),
        Elasticsearch(suggestion='default')
    ]
    n_times = MeanSquaredDisplacementValues.n_times.m_copy()
    times = MeanSquaredDisplacementValues.times.m_copy()
    value = MeanSquaredDisplacementValues.value.m_copy()
    errors = MeanSquaredDisplacementValues.errors.m_copy()

    diffusion_constant_value = DiffusionConstantValues.value.m_copy()
    diffusion_constant_error_type = DiffusionConstantValues.error_type.m_copy()
    diffusion_constant_errors = DiffusionConstantValues.errors.m_copy()


class DynamicalProperties(MSection):
    m_def = Section(
        description='''
        Dynamical properties.
        ''',
    )
    mean_squared_displacement = SubSection(
        sub_section=MeanSquaredDisplacement.m_def,
        repeats=True,
        a_elasticsearch=Elasticsearch(material_entry_type, nested=True)
    )


class SolarCell(MSection):
    m_def = Section(
        description='''
        Properties of solar cells.
        '''
    )
    efficiency = Quantity(
        type=np.float64,
        shape=[],
        description='''
        Power conversion effciency of a solar cell in percentage %.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type)
    )
    fill_factor = Quantity(
        type=np.float64,
        shape=[],
        description='''
        Fill factor of a solar cell in absolute values (from 0 to 1).
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type)
    )
    open_circuit_voltage = Quantity(
        type=np.float64,
        unit='V',
        shape=[],
        description='''
        Open circuit voltage of a solar cell.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type)
    )
    short_circuit_current_density = Quantity(
        type=np.float64,
        unit='A / m**2',
        shape=[],
        description='''
        Short circuit current density of a solar cell.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type)
    )
    illumination_intensity = Quantity(
        type=np.float64,
        unit=('W/m**2'),
        shape=[],
        description='''
        The light intensity during the IV measurement.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type)
    )
    device_area = Quantity(
        type=np.float64,
        unit=('m**2'),
        shape=[],
        description='''
        The total area of the cell during IV and stability measurements under illumination.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type)
    )
    device_architecture = Quantity(
        type=str,
        description='''
        Device architecture of the solar cell. Examples are:
        `pn-Heterojunction`, `pin`, `nip`, ...
        ''',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )
    device_stack = Quantity(
        type=str,
        shape=['0..*'],
        description='''
        Layers of the entire device.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )
    absorber = Quantity(
        type=str,
        shape=['0..*'],
        description='''
        Absorber layers used in the solar cell.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )
    absorber_fabrication = Quantity(
        type=str,
        shape=['0..*'],
        description='''
        Technique describing the fabrication of the absorber layer. Examples are:
        `Spin-coating`, `Evaporation`, `Doctor blading`, ...
        ''',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )
    electron_transport_layer = Quantity(
        type=str,
        shape=['0..*'],
        description='''
        Electron selective contact layers used in the solar cell.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )
    hole_transport_layer = Quantity(
        type=str,
        shape=['0..*'],
        description='''
        Hole selective contact layers used in the solar cell.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )
    substrate = Quantity(
        type=str,
        shape=['0..*'],
        description='''
        Substrate layers used in the solar cell.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )
    back_contact = Quantity(
        type=str,
        shape=['0..*'],
        description='''
        Back contact layers used in the solar cell.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )


class OptoelectronicProperties(MSection):
    m_def = Section(
        description='''
        Optoelectronic properties.
        '''
    )
    solar_cell = SubSection(
        sub_section=SolarCell.m_def,
        repeats=False
    )


class EELSInstrument(MSection):
    m_def = Section(
        description='''
        Base class for an EELS instrument.
        ''',
    )
    detector_type = Quantity(
        type=str,
        description='''
        Detector type.
        ''')
    max_energy = Quantity(
        type=np.float64,
        unit='joule',
        description='''
        Maximum energy of the detector.
        ''')
    min_energy = Quantity(
        type=np.float64,
        unit='joule',
        description='''
        Minimum energy of the detector.
        ''')
    guntype = Quantity(
        type=str,
        description='''
        Gun type of the detector.
        ''')
    beam_energy = Quantity(
        type=np.float64,
        unit='volt',
        description='''
        Incoming beam energy.
        ''')
    beam_current = Quantity(
        type=str,
        description='''
        Incoming beam current.
        ''')
    resolution = Quantity(
        type=np.float64,
        unit='joule',
        description='''
        Energy resolution of the detector.
        ''')
    step_size = Quantity(
        type=str,
        description='''
        Step size for axes in units of energy / pixel.
        ''')
    acquisition_mode = Quantity(
        type=str,
        description='''
        Acquisition mode for the counts in the detector.
        ''')
    dark_current = Quantity(
        type=bool,
        description='''
        Is dark current or noise to be substract included in the output?
        ''')


class EELSMethodology(MSection):
    m_def = Section(
        description='''
        Base class for the EELS methodology.
        ''',
    )
    detector_type = EELSInstrument.detector_type.m_copy(
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion="default")])
    resolution = EELSInstrument.resolution.m_copy(
        a_elasticsearch=[
            Elasticsearch(material_entry_type)])
    max_energy = EELSInstrument.max_energy.m_copy(
        a_elasticsearch=[
            Elasticsearch(material_entry_type)])
    min_energy = EELSInstrument.min_energy.m_copy(
        a_elasticsearch=[
            Elasticsearch(material_entry_type)])


class SpectraProvenance(ProvenanceTracker):
    m_def = Section(
        description="""
        Contains provenance information (mainly the methodology section) for spectra properties
        derived from an experiment or a calculation.
        """,
    )
    eels = SubSection(
        sub_section=EELSMethodology.m_def,
    )
    electronic_structure = SubSection(
        sub_section=ElectronicStructureProvenance.m_def,
        repeats=True
    )


class Spectra(MSection):
    m_def = Section(
        description='''
        Base class for Spectra calculation information as obtained from an experiment or a computation.
        ''',
    )
    type = Quantity(
        type=MEnum('EELS', 'XAS', 'XANES', 'EXAFS', 'XES', 'XPS', 'RXIS', config.services.unavailable_value),
        description='''
        Identifier for the methodology done to obtain the spectra data: EELS, XAS, XPS, etc.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )
    label = Quantity(
        type=MEnum('computation', 'experiment'),
        description='''
        Identifier for the source of the spectra data, either 'computation' or 'experiment'.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )
    n_energies = Quantity(
        type=np.int32,
        description='''
        Number of excitation energies.
        ''')
    energies = Quantity(
        type=np.float64,
        shape=['n_energies'],
        unit='joule',
        description='''
        Excitation energies for which the spectra is obtained.
        ''')
    intensities = Quantity(
        type=np.float64,
        shape=['n_energies'],
        description='''
        Intensitites obtained at each excitation energy. This can be computationally calculated,
        or electron counts coming from an experiment. In arbitrary units.
        ''')
    intensities_units = Quantity(
        type=str,
        description='''
        Units in which the intensities of the spectra are returned. It can be `F/m` as for
        the dielectric constant, or `counts` for the data of a CCD device.
        ''')
    provenance = SubSection(sub_section=SpectraProvenance.m_def)


class SpectroscopicProperties(MSection):
    m_def = Section(
        description='''
        Spectroscopic properties.
        ''',
    )
    spectra = SubSection(
        sub_section=Spectra.m_def,
        repeats=True,
        a_elasticsearch=Elasticsearch(material_type, nested=True)
    )


class Properties(MSection):
    m_def = Section(
        description='''
        Contains the physical properties that have been calculated or used in
        this entry.
        '''
    )
    structural = SubSection(sub_section=StructuralProperties.m_def, repeats=False)
    dynamical = SubSection(sub_section=DynamicalProperties.m_def, repeats=False)
    structures = SubSection(sub_section=Structures.m_def, repeats=False)
    vibrational = SubSection(sub_section=VibrationalProperties.m_def, repeats=False)
    electronic = SubSection(sub_section=ElectronicProperties.m_def, repeats=False)
    optoelectronic = SubSection(sub_section=OptoelectronicProperties.m_def, repeats=False)
    mechanical = SubSection(sub_section=MechanicalProperties.m_def, repeats=False)
    thermodynamic = SubSection(sub_section=ThermodynamicProperties.m_def, repeats=False)
    spectroscopic = SubSection(sub_section=SpectroscopicProperties.m_def, repeats=False)
    geometry_optimization = SubSection(sub_section=GeometryOptimization.m_def, repeats=False)

    n_calculations = Quantity(
        type=int,
        description='''
        The number of performed single configuration calculations.'
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type, metrics=dict(n_calculations='sum')),
    )
    available_properties = Quantity(
        type=str,
        default=[],
        derived=lambda a: available_properties(a),
        shape=['0..*'],
        description='Subset of the property names that are present in this entry.',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )


class ELN(MSection):
    sections = Quantity(
        type=str, shape=['*'],
        description='''
            The type of sections used in entries to search for. By default these are the names
            of the used section definitions.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type))

    tags = Quantity(
        type=str, shape=['*'],
        description='''
            Short tags that are useful to quickly search based on various
            user defined criteria.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type))

    names = Quantity(
        type=str, shape=['*'],
        description='''
            Short human readable and descriptive names that appear in
            ELN entries.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type, mapping='text'))

    descriptions = Quantity(
        type=str, shape=['*'],
        description='''
            'Human descriptions that appear in ELN entries.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type, mapping='text'))

    instruments = Quantity(
        type=str, shape=['*'],
        description='''
            The name or type of instrument used in an activity, e.g. process or
            measurement.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type))

    methods = Quantity(
        type=str, shape=['*'],
        description='''
            The name or the applied method in an activity, e.g. process or measurement
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type))

    lab_ids = Quantity(
        type=str, shape=['*'],
        description='''
            The laboratory specific id for any item, e.g. sample, chemical, instrument.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')])


class Results(MSection):
    m_def = Section(
        description='''
        Contains a summary of the entry contents.
        '''
    )
    material = SubSection(sub_section=Material.m_def, repeats=False)
    method = SubSection(sub_section=Method.m_def, repeats=False)
    properties = SubSection(sub_section=Properties.m_def, repeats=False)
    eln = SubSection(sub_section=ELN.m_def, repeats=False)


m_package.__init_metainfo__()
