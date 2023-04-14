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
from nomad.utils import traverse_reversed
from nomad.atomutils import Formula
from nomad.datamodel.metainfo.measurements import Spectrum
from nomad.datamodel.metainfo.simulation.system import Atoms
from nomad.datamodel.metainfo.workflow import (
    EquationOfState,
    EOSFit,
    RadialDistributionFunction as RDFWorkflow,
    RadialDistributionFunctionValues,
    MeanSquaredDisplacement as MSDWorkflow,
    MeanSquaredDisplacementValues,
    DiffusionConstantValues,
    Workflow
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
    Calculation
)  # noqa
from nomad.datamodel.metainfo.simulation.method import (
    BasisSet, Scf, Electronic, Smearing,
    GW as GWMethod, HubbardKanamoriModel as HubbardKanamori, AtomParameters,
    DMFT as DMFTMethod
)  # noqa
from nomad.datamodel.metainfo.workflow import (
    GeometryOptimization as MGeometryOptimization,
    MolecularDynamics as MMolecularDynamics,
    Thermodynamics, IntegrationParameters
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
    return None if formula is None else Formula(formula).format('hill')


def get_formula_iupac(formula: str) -> str:
    '''
    Converts the given chemical formula into the IUPAC format.

    Args:
        formula: Original formula.

    Returns:
        Chemical formula in the IUPAC format.
    '''
    return None if formula is None else Formula(formula).format('iupac')


def available_properties(root: MSection) -> List[str]:
    '''Returns a list of property names that are available in results.properties.

    Args:
        root: The metainfo section containing the properties

    Returns:
        List of property names that are present
    '''
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
        'spectroscopy.eels': 'eels',
        'optoelectronic.solar_cell': 'solar_cell',
    }
    available_properties: List[str] = []
    for path, shortcut in available_property_names.items():
        for _ in traverse_reversed(root, path.split('.')):
            available_properties.append(shortcut)
            break
    return sorted(available_properties)


tokenizer_formula = get_tokenizer(r'[A-Z][a-z]?\d*')


class BandGap(MSection):
    m_def = Section(
        description='''
        Band gap information.
        '''
    )
    label = Quantity(
        type=str,
        description='''
        Label to identify the band gap data, e.g. method employed.
        '''
    )
    index = Quantity(
        type=np.int32,
        description='''
        Index of the data, e.g. spin channel index.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    value = Quantity(
        type=np.float64,
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


class BandGapElectronic(BandGap):
    m_def = Section(
        description='''
        Band gap information for electronic structure.
        '''
    )

    energy_highest_occupied = Quantity(
        type=np.float64,
        unit='joule',
        shape=[],
        description='''
        The highest occupied energy.
        ''',
    )
    energy_lowest_unoccupied = Quantity(
        type=np.float64,
        unit='joule',
        shape=[],
        description='''
        The lowest unoccupied energy.
        ''',
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
    pbc = Quantity(
        type=bool,
        shape=[3],
        description='''
        Periodic boundary conditions of the cell lattice vectors, given in order a, b, c.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )


class Prototype(MSection):
    '''
    Information on the prototype corresponding to the current section.
    '''
    m_def = Section(validate=False)

    aflow_id = Quantity(
        type=str,
        shape=[],
        description='''
        AFLOW id of the prototype (see
        http://aflowlib.org/CrystalDatabase/prototype_index.html) identified on the basis
        of the space_group and normalized_wyckoff.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion='simple')
        ]
    )
    assignment_method = Quantity(
        type=str,
        shape=[],
        description='''
        Method used to identify the prototype.
        '''
    )
    label = Quantity(
        type=str,
        shape=[],
        description='''
        Label of the prototype identified on the basis of the space_group and
        normalized_wyckoff. The label is in the same format as in the read_prototypes
        function: <space_group_number>-<prototype_name>-<Pearson's symbol>).
        '''
    )
    name = Quantity(
        type=str,
        description='''
        A common name for this prototypical structure, e.g. fcc, bcc.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_type),
            Elasticsearch(suggestion='default')
        ],
    )
    formula = Quantity(
        type=str,
        description='''
        The formula of the prototypical material for this structure.
        ''',
        a_elasticsearch=Elasticsearch(material_type),
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


class Relation(MSection):
    '''
    Contains information about the relation between two different systems.
    '''
    m_def = Section(validate=False)
    type = Quantity(
        type=MEnum('subsystem', 'idealization'),
        description='''
        The type of relation.
        ''',
        a_elasticsearch=Elasticsearch(material_type),
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
        a_elasticsearch=Elasticsearch(material_type)
    )
    method = Quantity(
        type=MEnum('parser', 'user', 'matid'),
        description='''
        The method used for identifying this system.
        ''',
        a_elasticsearch=Elasticsearch(material_type)
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
    atoms = SubSection(
        description='''
        The atomistic structure that is associated with this
        system'.
        ''',
        sub_section=Atoms.m_def,
        repeats=False
    )
    atoms_ref = Quantity(
        type=Atoms,
        description='''
        Reference to an atomistic structure that is associated with this
        system'.
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
        Indices of the atoms belonging to this group. These indices refer to
        the original system. Each row represents a new instance.
        '''
    )
    system_relation = SubSection(sub_section=Relation.m_def, repeats=False)
    cell = SubSection(sub_section=Cell.m_def, repeats=False)
    symmetry = SubSection(sub_section=SymmetryNew.m_def, repeats=False)
    prototype = SubSection(sub_section=Prototype.m_def, repeats=False)


# =============================================================================


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
        Structural class determined from the atomic structure.
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
    basis_set_type = Quantity(
        type=MEnum(basis_set_types),
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
    starting_point_type = Quantity(
        type=MEnum(list(xc_treatments.values()) + [unavailable, not_processed]),
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


class DMFT(MSection):
    m_def = Section(
        description='''
        Methodology for a DMFT calculation.
        '''
    )
    impurity_solver_type = DMFTMethod.impurity_solver.m_copy()
    impurity_solver_type.m_annotations['elasticsearch'] = [
        Elasticsearch(material_entry_type),
        Elasticsearch(suggestion='default')
    ]
    total_filling = Quantity(
        type=np.float64,
        description='''
        Total filling of the correlated atoms in the unit cell per spin ∈[0.0, 1.0]. E.g., half-filling
        is defined as 0.5.
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    inverse_temperature = DMFTMethod.inverse_temperature.m_copy()
    inverse_temperature.m_annotations['elasticsearch'] = [Elasticsearch(material_entry_type)]
    magnetic_state = DMFTMethod.magnetic_state.m_copy()
    magnetic_state.description = 'Magnetic state in which the DMFT calculation is done.'
    magnetic_state.m_annotations['elasticsearch'] = [
        Elasticsearch(material_entry_type),
        Elasticsearch(suggestion='default')
    ]
    u = HubbardKanamori.u.m_copy()
    u.m_annotations['elasticsearch'] = [Elasticsearch(material_entry_type)]
    hunds_hubbard_ratio = Quantity(
        type=np.float64,
        description='''
        Ratio JH/U, with JH being the Hunds coupling and U being the Hubbard local interaction.
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
        '''
    )


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
    dmft = SubSection(sub_section=DMFT.m_def, repeats=False)
    quantum_cms = SubSection(sub_section=QuantumCMS.m_def, repeats=False)
    precision = SubSection(sub_section=Precision.m_def, repeats=False)


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
        type=MEnum(['DFT', 'Projection', 'GW', 'DMFT', 'CoreHole', 'BSE', 'EELS', 'XPS', config.services.unavailable_value]),
        description='''
        Common name for the used method.
        ''',
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion='default')
        ],
    )

    workflow_name = Workflow.type.m_copy()
    workflow_name.m_annotations['elasticsearch'] = [
        Elasticsearch(material_entry_type),
        Elasticsearch(suggestion='default')
    ]
    simulation = SubSection(sub_section=Simulation.m_def, repeats=False)


class MolecularDynamics(MSection):
    m_def = Section(
        description="""
        Methodology for molecular dynamics.
        """,
    )
    time_step = IntegrationParameters.integration_timestep.m_copy()
    time_step.m_annotations["elasticsearch"] = Elasticsearch(material_entry_type)
    ensemble_type = MMolecularDynamics.thermodynamic_ensemble.m_copy()
    ensemble_type.m_annotations["elasticsearch"] = Elasticsearch(material_entry_type)


class Methodology(MSection):
    m_def = Section(
        description="""
        Contains methodological information and can be attached to any physical
        property.
        """,
    )
    molecular_dynamics = SubSection(sub_section=MolecularDynamics.m_def, repeats=False)


class PropertySection(MSection):
    m_def = Section(
        description="""
        Base class for that can be used to attach a specific methodology to a
        physical property.
        """,
    )
    methodology = SubSection(sub_section=Methodology.m_def, repeats=False)


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
        sub_section=BandGapElectronic.m_def,
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
        sub_section=BandGapElectronic.m_def,
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
    tau = GreensFunctionsCalculation.tau.m_copy()
    real_greens_function_tau = Quantity(
        type=np.float64,
        shape=['n_atoms_per_unit_cell', 2, 'n_correlated_orbitals', 'n_tau'],
        description='''
        Real part (extraction done in normalizer) of the Green's function in tau (imaginary time).
        '''
    )
    matsubara_freq = GreensFunctionsCalculation.matsubara_freq.m_copy()
    imag_self_energy_iw = Quantity(
        type=np.float64,
        shape=['n_atoms_per_unit_cell', 2, 'n_correlated_orbitals', '2 * n_matsubara_freq'],
        description='''
        Imaginary part (extraction done in normalizer) of the Self energy in Matsubara (imaginary frequency).
        '''
    )
    orbital_occupations = GreensFunctionsCalculation.orbital_occupations.m_copy()
    quasiparticle_weights = GreensFunctionsCalculation.quasiparticle_weights.m_copy()
    chemical_potential = GreensFunctionsCalculation.chemical_potential.m_copy()


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
        type=MGeometryOptimization.energies,
        description='''
        List of energy_total values gathered from the single configuration
        calculations that are a part of the optimization trajectory.
        ''',
    )
    structure_optimized = SubSection(
        description='''
        Contains a structure that is the result of a geometry optimization.
        ''',
        sub_section=Structure.m_def,
        repeats=False,
    )
    type = MGeometryOptimization.type.m_copy()
    convergence_tolerance_energy_difference = MGeometryOptimization.convergence_tolerance_energy_difference.m_copy()
    convergence_tolerance_energy_difference.m_annotations['elasticsearch'] = Elasticsearch(material_entry_type)
    convergence_tolerance_force_maximum = MGeometryOptimization.convergence_tolerance_force_maximum.m_copy()
    convergence_tolerance_force_maximum.m_annotations['elasticsearch'] = Elasticsearch(material_entry_type)
    final_force_maximum = MGeometryOptimization.final_force_maximum.m_copy()
    final_force_maximum.m_annotations['elasticsearch'] = Elasticsearch(material_entry_type)
    final_energy_difference = MGeometryOptimization.final_energy_difference.m_copy()
    final_energy_difference.m_annotations['elasticsearch'] = Elasticsearch(material_entry_type)
    final_displacement_maximum = MGeometryOptimization.final_displacement_maximum.m_copy()
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
    band_structure_electronic = SubSection(sub_section=BandStructureElectronic.m_def, repeats=True)
    dos_electronic = SubSection(sub_section=DOSElectronic.m_def, repeats=True)
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


class Trajectory(PropertySection):
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


class RadiusOfGyration(QuantityDynamic, PropertySection):
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


class RadialDistributionFunction(PropertySection):
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


class MeanSquaredDisplacement(PropertySection):
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
    structural = SubSection(sub_section=StructuralProperties.m_def, repeats=False)
    dynamical = SubSection(sub_section=DynamicalProperties.m_def, repeats=False)
    structures = SubSection(sub_section=Structures.m_def, repeats=False)
    vibrational = SubSection(sub_section=VibrationalProperties.m_def, repeats=False)
    electronic = SubSection(sub_section=ElectronicProperties.m_def, repeats=False)
    optoelectronic = SubSection(sub_section=OptoelectronicProperties.m_def, repeats=False)
    mechanical = SubSection(sub_section=MechanicalProperties.m_def, repeats=False)
    thermodynamic = SubSection(sub_section=ThermodynamicProperties.m_def, repeats=False)
    spectroscopy = SubSection(sub_section=SpectroscopyProperties.m_def, repeats=False)
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
