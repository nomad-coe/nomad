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
from nomad.datamodel.metainfo.workflow import EquationOfState, EOSFit, Workflow
from nomad.metainfo.elasticsearch_extension import (
    Elasticsearch,
    material_type,
    material_entry_type,
    entry_type as search_entry_type,
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
    Dos,
    BandStructure as BandStructureCalculation,
    BandEnergies,
    DosValues,
    Calculation)  # noqa
from nomad.datamodel.metainfo.simulation.method import (
    BasisSet, Scf, Electronic, Smearing, GW as GWMethod
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


tokenizer_formula = get_tokenizer(r'[A-Z][a-z]?\d*')


class BandGap(MSection):
    m_def = Section(
        description='''
        Band gap information for each spin channel.
        '''
    )
    index = Quantity(
        type=np.dtype(np.int32),
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
    species = SubSection(sub_section=OptimadeSpecies.m_def, repeats=True)
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


# =============================================================================
# New topological data
class Cell(MSection):
    m_def = Section(
        description='''
        Properties of a unit cell.
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
    volume = Quantity(
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
        ''',
        a_elasticsearch=Elasticsearch(material_entry_type),
    )
    mass_density = Quantity(
        type=np.dtype(np.float64),
        unit='kg / m ** 3',
        description='''
        Mass density of the material.
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
        type=np.dtype(np.float64),
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
        type=np.dtype(np.float64),
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


class Species(MSection):
    '''
    Contains information about a particle species. Note that the particle can
    also be something else than atoms, e.g. coarse-grained particle, isotopes,
    etc.
    '''
    m_def = Section(validate=False)
    name = Quantity(
        type=str,
        description='''
        Name that uniquely identifies this species within a system.
        ''',
        a_elasticsearch=Elasticsearch(material_type)
    )
    mass = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='kilogram',
        description='''
        Mass of the species.
        '''
    )
    atomic_number = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        The atomic number of the species if available.
        '''
    )


class Atoms(MSection):
    '''
    Describes the atomic structure of the physical system. This includes the atom
    positions, lattice vectors, etc.
    '''
    m_def = Section(validate=False)

    concentrations = Quantity(
        type=np.dtype(np.float64),
        shape=['n_atoms'],
        description='''
        Concentrations of the species defined by labels which can be assigned for systems
        with variable compositions.
        '''
    )
    labels = Quantity(
        type=str,
        shape=['n_atoms'],
        description='''
        List containing the labels of the atoms. In the usual case, these correspond to
        the chemical symbols of the atoms. One can also append an index if there is a
        need to distinguish between species with the same symbol, e.g., atoms of the
        same species assigned to different atom-centered basis sets or pseudo-potentials,
        or simply atoms in different locations in the structure such as those in the bulk
        and on the surface. In the case where a species is not an atom, and therefore
        cannot be representated by a chemical symbol, the label can simply be the name of
        the particles.
        '''
    )
    positions = Quantity(
        type=np.dtype(np.float64),
        shape=['n_atoms', 3],
        unit='meter',
        description='''
        Positions of all the species, in cartesian coordinates. This metadata defines a
        configuration and is therefore required. For alloys where concentrations of
        species are given for each site in the unit cell, it stores the position of the
        sites.
        '''
    )
    velocities = Quantity(
        type=np.dtype(np.float64),
        shape=['n_atoms', 3],
        unit='meter / second',
        description='''
        Velocities of the nuclei, defined as the change in cartesian coordinates of the
        nuclei with respect to time.
        '''
    )
    lattice_vectors = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        unit='meter',
        description='''
        Lattice vectors in cartesian coordinates of the simulation cell. The
        last (fastest) index runs over the $x,y,z$ Cartesian coordinates, and the first
        index runs over the 3 lattice vectors.
        '''
    )
    lattice_vectors_reciprocal = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        unit='1/meter',
        description='''
        Reciprocal lattice vectors in cartesian coordinates of the simulation cell. The
        first index runs over the $x,y,z$ Cartesian coordinates, and the second index runs
        over the 3 lattice vectors.
        '''
    )
    local_rotations = Quantity(
        type=np.dtype(np.float64),
        shape=['n_atoms', 3, 3],
        description='''
        A rotation matrix defining the orientation of each atom. If the rotation matrix
        cannot be specified for an atom, the remaining atoms should set it to
        the zero matrix (not the identity!)
        '''
    )
    periodic = Quantity(
        type=bool,
        shape=[3],
        description='''
        Denotes if periodic boundary condition is applied to each of the lattice vectors.'
        '''
    )
    supercell_matrix = Quantity(
        type=np.dtype(np.int32),
        shape=[3, 3],
        description='''
        Specifies the matrix that transforms the unit-cell into the super-cell in which
        the actual calculation is performed.
        '''
    )
    species = SubSection(sub_section=Species.m_def, repeats=False)
    wyckoff_sets = SubSection(sub_section=WyckoffSet.m_def, repeats=True)


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
    '''Describes a physical structure as identified in an entry. Can be a part
    of a larger structural hierarchy, i.e. the topology.
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
        The structural classification for this system.
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
    formula_hill = Quantity(
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
    formula_reduced = Quantity(
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
    formula_anonymous = Quantity(
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
    formula_reduced_fragments = Quantity(
        type=str,
        shape=['*'],
        description='''
        The reduced formula separated into individual terms containing both the atom
        type and count. Used for searching parts of a formula.
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
    atoms = Quantity(
        type=Structure,
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
        type=np.dtype(np.int64),
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
    topology = SubSection(
        sub_section=System.m_def,
        repeats=True,
        a_elasticsearch=Elasticsearch(material_type, nested=True)
    )


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

    workflow_name = Workflow.type.m_copy()
    workflow_name.shape = ['*']
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
        sub_section=StructureOptimized.m_def,
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
    band_structure_electronic = SubSection(sub_section=BandStructureElectronic.m_def, repeats=False)
    dos_electronic = SubSection(sub_section=DOSElectronic.m_def, repeats=False)


class SpectroscopyProperties(MSection):
    m_def = Section(
        description='''
        Spectroscopic properties.
        ''',
    )
    spectrum = Quantity(type=Spectrum)


class QuantityDynamic(MSection):
    m_def = Section(
        description="""
        Contains the values for a quantity at different times.
        """
    )
    time = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description="""
        The explicit times at which the values are evaluated. Provide either
        this or time_step and time_start.
        """,
    )
    time_step = Quantity(
        type=np.dtype(np.float64),
        unit='second',
        description="""
        The time step between successive evaluations. Provide either
        this and time_start or the explicit times.
        """,
    )
    time_start = Quantity(
        type=np.dtype(np.float64),
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
        type=np.dtype(np.float64),
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
        type=np.dtype(np.float64),
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
        type=np.dtype(np.float64),
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
        type=np.dtype(np.float64),
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


class SolarCell(MSection):

    m_def = Section(
        description='''
        Properties of solar cells.
        ''',
    )

    open_circuit_voltage = Quantity(
        type=np.dtype(np.float64),
        unit='V',
        shape=[],
        description="""
        Open circuit voltage of a solar cell.
                    """)

    short_circuit_current_density = Quantity(
        type=np.dtype(np.float64),
        unit='A / m**2',
        shape=[],
        description="""
        Short circuit current density of a solar cell.
                    """)

    fill_factor = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description="""
        Fill factor of a solar cell in absolute values (from 0 to 1).
                    """)

    efficiency = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description="""
        Power conversion effciency of a solar cell in percentage %.
                    """)

    illumination_intensity = Quantity(
        type=np.dtype(np.float64),
        unit=('W/m**2'),
        shape=[],
        description="""
    The light intensity during the IV measurement.
                    """)

    substrate = Quantity(
        type=str,
        description="""
        Substrate used in the solar cell. Might be a stack in which seperate layers are
        separated by a vertical bar " | ".
                    """)

    back_contact = Quantity(
        type=str,
        description="""
        Back contact used in the solar cell. Might be a stack in which seperate layers are
        separated by a vertical bar " | ".
                    """)

    electron_transport_layer = Quantity(
        type=str,
        description="""
        Electron selective contact used in the solar cell. Might be a stack in which
        seperate layers are separated by a vertical bar " | ".
                    """)

    hole_transport_layer = Quantity(
        type=str,
        description="""
        Hole selective contact used in the solar cell. Might be a stack in which seperate layers are
        separated by a vertical bar " | ".
                    """)

    absorber = Quantity(
        type=str,
        description="""
        Absorber layer used in the solar cell. Might be an stack in which seperate layers are
        separated by a vertical bar " | ".
                    """)

    device_stack = Quantity(
        type=str,
        description="""
        Substrate used in the solar cell. Might be an stack in which seperate layers are
        separated by a vertical bar " | ".
                    """)

    device_architecture = Quantity(
        type=str,
        description="""
        Device architecture of the solar cell. Examples are:
        `pn-Heterojunction`, `pin`, `nip`, ...
                    """)

    device_area = Quantity(
        type=np.dtype(np.float64),
        unit=('m**2'),
        shape=[],
        description="""
        The total area of the cell during IV and stability measurements under illumination.
                    """)

    absorber_fabrication = Quantity(
        type=str,
        description="""
        Technique describing the fabrication of the absorber layer. Examples are:
        `Spin-coating`, `Evaporation`, `Doctor blading`, ...
                    """)


class Properties(MSection):
    m_def = Section(
        description='''
        Contains the physical properties that have been calculated or used in
        this entry.
        '''
    )
    structures = SubSection(sub_section=Structures.m_def, repeats=False)
    vibrational = SubSection(sub_section=VibrationalProperties.m_def, repeats=False)
    electronic = SubSection(sub_section=ElectronicProperties.m_def, repeats=False)
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
        a_elasticsearch=Elasticsearch(search_entry_type))

    tags = Quantity(
        type=str, shape=['*'],
        description='''
            Short tags that are useful to quickly search based on various
            user defined criteria.
        ''',
        a_elasticsearch=Elasticsearch(search_entry_type))

    names = Quantity(
        type=str, shape=['*'],
        description='''
            Short human readable and descriptive names that appear in
            ELN entries.
        ''',
        a_elasticsearch=Elasticsearch(search_entry_type, mapping='text'))

    descriptions = Quantity(
        type=str, shape=['*'],
        description='''
            'Human descriptions that appear in ELN entries.
        ''',
        a_elasticsearch=Elasticsearch(search_entry_type, mapping='text'))

    instruments = Quantity(
        type=str, shape=['*'],
        description='''
            The name or type of instrument used in an activity, e.g. process or
            measurement.
        ''',
        a_elasticsearch=Elasticsearch(search_entry_type))

    methods = Quantity(
        type=str, shape=['*'],
        description='''
            The name or the applied method in an activity, e.g. process or measurement
        ''',
        a_elasticsearch=Elasticsearch(search_entry_type))

    lab_ids = Quantity(
        type=str, shape=['*'],
        description='''
            The laboratory specific id for any item, e.g. sample, chemical, instrument.
        ''',
        a_elasticsearch=Elasticsearch(search_entry_type))


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
    solarcell = SubSection(sub_section=SolarCell.m_def, repeats=False)


m_package.__init_metainfo__()
