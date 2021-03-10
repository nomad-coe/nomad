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
from elasticsearch_dsl import Text, Keyword, Integer, Boolean, Double

from ase.data import chemical_symbols

from nomad import config
from nomad.metainfo.search_extension import Search
from nomad.metainfo import (
    MSection,
    Section,
    SubSection,
    Quantity,
    MEnum,
    Reference,
    SectionProxy,
    QuantityReference,
)
from nomad.datamodel.optimade import Species  # noqa
from nomad.datamodel.metainfo.common_dft import (  # noqa
    Method as section_method,
    Dos,
    KBand,
    KBandSegment,
    ThermodynamicalProperties,
)

structure_classes = [
    "1D",
    "2D",
    "atom",
    "bulk",
    "molecule / cluster",
    "surface",
    "unavailable",
    "not processed",
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
    'unavailable',
    'not processed',
]
core_electron_treatments = [
    "full all electron",
    "all electron frozen core",
    "pseudopotential",
    "unavailable",
]


class WyckoffSet(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Section for storing Wyckoff set information. Only available for
        conventional cells that have undergone symmetry analysis.
        """
    )
    wyckoff_letter = Quantity(
        type=str,
        description="""
        The Wyckoff letter for this set.
        """
    )
    indices = Quantity(
        type=np.dtype("i4"),
        shape=["1..*"],
        description="""
        Indices of the atoms belonging to this group.
        """
    )
    element = Quantity(
        type=str,
        description="""
        Chemical element at this Wyckoff position.
        """
    )
    x = Quantity(
        type=np.dtype(np.float64),
        description="""
        The free parameter x if present.
        """
    )
    y = Quantity(
        type=np.dtype(np.float64),
        description="""
        The free parameter y if present.
        """
    )
    z = Quantity(
        type=np.dtype(np.float64),
        description="""
        The free parameter z if present.
        """
    )


class LatticeParameters(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Lattice parameters of a cell.
        """,
        a_search="lattice_parameters"
    )
    a = Quantity(
        type=float,
        description="""
        Length of the first basis vector.
        """,
        a_search=Search()
    )
    b = Quantity(
        type=float,
        description="""
        Length of the second basis vector.
        """,
        a_search=Search()
    )
    c = Quantity(
        type=float,
        description="""
        Length of the third basis vector.
        """,
        a_search=Search()
    )
    alpha = Quantity(
        type=float,
        description="""
        Angle between second and third basis vector.
        """,
        a_search=Search()
    )
    beta = Quantity(
        type=float,
        description="""
        Angle between first and third basis vector.
        """,
        a_search=Search()
    )
    gamma = Quantity(
        type=float,
        description="""
        Angle between first and second basis vector.
        """,
        a_search=Search()
    )


class Structure(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Describes an atomistic structure.
        """
    )
    dimension_types = Quantity(
        type=int,
        shape=[3],
        default=[0, 0, 0],
        description="""
        List of three integers. For each of the three directions indicated by
        the three lattice vectors (see property lattice_vectors). This list
        indicates if the direction is periodic (value 1) or non-periodic (value
        0). Note: the elements in this list each refer to the direction of the
        corresponding entry in lattice_vectors and not the Cartesian x, y, z
        directions.
        """
    )
    nperiodic_dimensions = Quantity(
        type=int,
        derived=lambda a: sum(a.dimension_types),
        a_search=Search(mapping=Integer()),
        description="""
        An integer specifying the number of periodic dimensions in the
        structure, equivalent to the number of non-zero entries in
        dimension_types.
        """
    )
    lattice_vectors = Quantity(
        type=np.dtype("float64"),
        shape=[3, 3],
        unit="m",
        description="""
        The three lattice vectors in Cartesian coordinates.
        """
    )
    cartesian_site_positions = Quantity(
        type=np.dtype("float64"),
        shape=["nsites", 3],
        unit="m",
        description="""
        Cartesian positions of each site. A site is an atom, a site potentially
        occupied by an atom, or a placeholder for a virtual mixture of atoms
        (e.g., in a virtual crystal approximation).
        """
    )
    nsites = Quantity(
        type=int,
        default=0,
        derived=lambda a: len(a.cartesian_site_positions),
        a_search=Search(mapping=Integer()),
        description="""
        An integer specifying the length of the cartesian_site_positions property.
        """
    )
    species_at_sites = Quantity(
        type=str,
        shape=["nsites"],
        description="""
        Name of the species at each site (where values for sites are specified with the same
        order of the cartesian_site_positions property). The properties of the species are
        found in the species property.
        """
    )
    cell_volume = Quantity(
        type=np.dtype(np.float64),
        unit="m ** 3",
        description="""
        Volume of the cell. The cell volume can only be reported consistently
        after idealization and may not perfectly correspond to the original
        simulation cell.
        """,
        a_search=Search()
    )
    species = SubSection(sub_section=Species.m_def, repeats=True)
    lattice_parameters = SubSection(sub_section=LatticeParameters.m_def)


class StructureOriginal(Structure):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Contains a selected representative structure from the the original
        data.
        """
    )


class StructurePrimitive(Structure):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Contains the primitive structure that is derived from
        structure_original. This primitive stucture has been idealized and the
        conventions employed by spglib are used.
        """
    )


class StructureConventional(Structure):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Contains the conventional structure that is derived from
        structure_original. This conventional stucture has been idealized and
        the conventions employed by spglib are used.
        """
    )
    wyckoff_sets = SubSection(sub_section=WyckoffSet.m_def, repeats=True)


class StructureOptimized(Structure):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Contains a structure that is the result of a geometry optimization.
        """
    )


class Symmetry(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Section containing information about the symmetry of the material. All
        of these properties are derived by running a symmetry analysis on a
        representative geometry from the original data. This original geometry
        is stored in results.properties together with the primitive and
        conventional structures.
        """
    )
    bravais_lattice = Quantity(
        type=str,
        shape=[],
        description="""
        Identifier for the Bravais lattice in Pearson notation. The first lowercase letter
        identifies the crystal family and can be one of the following: a (triclinic), b
        (monoclinic), o (orthorhombic), t (tetragonal), h (hexagonal) or c (cubic). The
        second uppercase letter identifies the centring and can be one of the following: P
        (primitive), S (face centred), I (body centred), R (rhombohedral centring) or F
        (all faces centred).
        """,
    )
    crystal_system = Quantity(
        type=str,
        shape=[],
        description="""
        Name of the crystal system. Can be one of the following: triclinic, monoclinic,
        orthorhombic, tetragonal, trigonal, hexagonal or cubic.
        """,
    )
    hall_number = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description="""
        The Hall number for this system.
        """,
    )
    hall_symbol = Quantity(
        type=str,
        shape=[],
        description="""
        The Hall symbol for this system.
        """,
    )
    point_group = Quantity(
        type=str,
        shape=[],
        description="""
        Symbol of the crystallographic point group in the Hermann-Mauguin notation.
        """,
    )
    space_group_number = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description="""
        Specifies the International Union of Crystallography (IUC) number of the 3D space
        group of this system.
        """,
    )
    space_group_symbol = Quantity(
        type=str,
        shape=[],
        description="""
        The International Union of Crystallography (IUC) short symbol of the 3D
        space group of this system.
        """,
    )
    prototype_formula = Quantity(
        type=str,
        description="""
        The formula of the prototypical material for this structure.
        """,
        a_search=Search()
    )
    prototype_aflow_id = Quantity(
        type=str,
        description="""
        The identifier of this structure in the AFLOW encyclopedia of
        crystallographic prototypes:
        http://www.aflowlib.org/prototype-encyclopedia/index.html
        """,
        a_search=Search()
    )
    structure_name = Quantity(
        type=str,
        description="""
        A common name for this structure, e.g. fcc, bcc.
        """,
        a_search=Search()
    )
    strukturbericht_designation = Quantity(
        type=str,
        description="""
        Classification of the material according to the historically grown
        "strukturbericht".
        """,
        a_search=Search()
    )


class Material(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Contains information that is specific to bulk crystalline materials.
        """
    )
    material_id = Quantity(
        type=str,
        description="""
        A fixed length, unique material identifier in the form of a hash
        digest.
        """,
        a_search=Search()
    )
    material_name = Quantity(
        type=str,
        description="""
        Meaningful names for this a material if any can be assigned.
        """,
        a_search=Search()
    )
    type_structural = Quantity(
        type=MEnum(structure_classes), default="not processed",
        description="""
        Classification based on structural features.
        """,
        a_search=Search(statistic_values=structure_classes)
    )
    type_functional = Quantity(
        type=str,
        shape=["*"],
        description="""
        Classification based on the functional properties.
        """,
        a_search=Search()
    )
    type_compound = Quantity(
        type=str,
        shape=["*"],
        description="""
        Classification based on the chemical formula.
        """,
        a_search=Search()
    )
    elements = Quantity(
        type=MEnum(chemical_symbols), shape=['1..*'],
        description="""
        Names of the different elements present in the structure.
        """,
        a_search=Search(mapping=Text(multi=True, fields={"keyword": Keyword()}))
    )
    nelements = Quantity(
        type=int,
        default=0,
        a_search=Search(mapping=Integer()),
        description="""
        Number of different elements in the structure as an integer.
        """
    )
    chemical_formula_descriptive = Quantity(
        type=str,
        a_search=Search(),
        description="""
            The chemical formula for a structure as a string in a form chosen by the API
            implementation.
        """
    )
    chemical_formula_reduced = Quantity(
        type=str,
        a_search=Search(),
        description="""
            The reduced chemical formula for a structure as a string with element symbols and
            integer chemical proportion numbers. The proportion number MUST be omitted if it is 1.
        """
    )
    chemical_formula_hill = Quantity(
        type=str,
        a_search=Search(),
        description="""
            The chemical formula for a structure in Hill form with element symbols followed by
            integer chemical proportion numbers. The proportion number MUST be omitted if it is 1.
        """
    )
    chemical_formula_anonymous = Quantity(
        type=str,
        a_search=Search(),
        description="""
            The anonymous formula is the chemical_formula_reduced, but where the elements are
            instead first ordered by their chemical proportion number, and then, in order left to
            right, replaced by anonymous symbols A, B, C, ..., Z, Aa, Ba, ..., Za, Ab, Bb, ... and
            so on.
        """
    )
    chemical_formula_reduced_fragments = Quantity(
        type=str,
        shape=["*"],
        description="""
        The reduced formula separated into individual terms containing both the atom
        type and count. Used for searching parts of a formula.
        """,
        a_search=Search(mapping=Text(multi=True))
    )
    symmetry = SubSection(sub_section=Symmetry.m_def, repeats=False, a_search="symmetry")


class DFT(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Methodology for a DFT calculation.
        """
    )
    basis_set_type = Quantity(
        type=MEnum(basis_set_types), default='not processed',
        description='The used basis set functions.',
        a_search=Search(statistic_values=basis_set_types)
    )
    basis_set_name = section_method.basis_set.m_copy()
    core_electron_treatment = Quantity(
        type=MEnum(core_electron_treatments),
        description="""
        How the core electrons are described.
        """,
        a_search=Search()
    )
    spin_polarized = Quantity(
        type=bool,
        description="""
        Whether the calculation is spin-polarized.
        """,
        a_search=Search(mapping=Boolean()),
    )
    scf_threshold_energy_change = section_method.scf_threshold_energy_change.m_copy()
    van_der_Waals_method = section_method.van_der_Waals_method.m_copy()
    relativity_method = section_method.relativity_method.m_copy()
    smearing_kind = section_method.smearing_kind.m_copy()
    smearing_width = section_method.smearing_width.m_copy()
    xc_functional_type = Quantity(
        type=str, default='not processed',
        description='The libXC based xc functional classification used in the simulation.',
        a_search=Search(
            statistic_values=list(xc_treatments.values()) + ['unavailable', 'not processed'],
            statistic_size=100
        )
    )
    xc_functional_names = Quantity(
        type=str, default=[], shape=['*'],
        description='The list of libXC functional names that where used in this entry.',
        a_search=Search(many_and='append')
    )


class GW(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Methodology for a GW calculation.
        """
    )
    gw_type = section_method.gw_type.m_copy()
    starting_point = Quantity(
        type=str, default=[], shape=['*'],
        description='The list of libXC functional names that were used for the ground state calculation.',
        a_search=Search(many_and='append')
    )


class Simulation(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Contains method details for a simulation entry.
        """
    )
    program_name = Quantity(
        type=str, default='not processed',
        description='The name of the used program.',
        a_search=Search()
    )
    program_version = Quantity(
        type=str, default='not processed',
        description='The version of the used program.',
        a_search=Search()
    )
    dft = SubSection(sub_section=DFT.m_def, repeats=False, a_search="dft")
    gw = SubSection(sub_section=GW.m_def, repeats=False, a_search="gw")
    phonon = Quantity(
        type=Reference(SectionProxy('Phonon')),
        shape=[],
        description='''
        Reference to phonon calculation methodology.
        ''',
        a_search="phonon",
    )
    geometry_optimization = Quantity(
        type=Reference(SectionProxy('GeometryOptimization')),
        shape=[],
        description='''
        Reference to geometry optimization methodology.
        ''',
        a_search="geometry_optimization",
    )
    molecular_dynamics = Quantity(
        type=Reference(SectionProxy('MolecularDynamics')),
        shape=[],
        description='''
        Reference to molecular dynamics methodology.
        ''',
        a_search="molecular_dynamics",
    )


class Method(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Contains a summary of the methodology that has been used in this entry.
        This methodology applies to all of the reported properties and
        determines the result of a single energy evalution. The individual
        properties may be further methodological details affect e.g. the
        sampling.
        """
    )
    method_id = Quantity(
        type=str,
        description="""
        Identifier for the used method. Only available if the method could be
        identified precisely.
        """,
        a_search=Search()
    )
    method_name = Quantity(
        type=MEnum(["DFT", "GW", config.services.unavailable_value]),
        description="""
        Common name for the used method.
        """,
        a_search=Search()
    )
    simulation = SubSection(sub_section=Simulation.m_def, repeats=False, a_search="simulation")


class HeatCapacityConstantVolume(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Contains the values of the specific (per mass) and isochoric (constant
        volume) heat capacity at different temperatures.
        """
    )
    # values = Quantity(
    #     type=np.dtype("float64"),
    #     shape=["*"],
    #     unit='joule / (kelvin * kilogram)',
    #     description="""
    #     Heat capacity values.
    #     """,
    # )
    temperature = Quantity(
        type=QuantityReference(ThermodynamicalProperties.thermodynamical_property_temperature),
        description='''
        Reference to the temperature values.
        ''',
    )


class HelmholtzFreeEnergy(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Contains the values of the specific (per mass) Helmholtz free energy at
        different temperatures.
        """
    )
    # values = Quantity(
    #     type=QuantityReference(SectionProxy('section_dos')),
    #     shape=[],
    #     description='''
    #     Reference to the values of Helmholtz free energy.
    #     ''',
    # )
    temperature = Quantity(
        type=QuantityReference(ThermodynamicalProperties.thermodynamical_property_temperature),
        description='''
        Reference to the temperature values.
        ''',
    )


class DOS(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Base class for density of states information.
        """,
    )
    energies = Quantity(
        type=Dos.dos_energies,
        description="""
        Array containing the set of discrete energy values for the density of
        states (DOS).
        """,
    )
    densities = Quantity(
        type=Dos.dos_values_normalized,
        description="""
        Density of states (DOS) values normalized with unit cell volume and
        number of atoms.
        """,
    )


class DOSPhonon(DOS):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Contains the total phonon density of states.
        """,
    )


class DOSElectronic(DOS):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Contains the total electronic density of states.
        """,
    )
    spin_polarized = Quantity(
        type=bool,
        description="""
        Whether the DOS is spin-polarized, i.e. is contains channels for both
        spin values.
        """,
        a_search=Search(mapping=Boolean()),
    )
    energy_fermi = Quantity(
        type=np.dtype(np.float64),
        unit='joule',
        shape=["n_spin_channels"],
        description="""
        Fermi energy for each spin channel.
        """,
    )
    energy_highest_occupied = Quantity(
        type=np.dtype(np.float64),
        unit='joule',
        shape=["n_spin_channels"],
        description="""
        The highest occupied energy for each spin channel.
        """,
    )
    energy_lowest_unoccupied = Quantity(
        type=np.dtype(np.float64),
        unit='joule',
        shape=["n_spin_channels"],
        description="""
        The lowest unoccupied energy for each spin channel.
        """,
    )


class BandStructure(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Base class for band structure information.
        """,
    )
    reciprocal_cell = Quantity(
        type=KBand.reciprocal_cell,
        description="""
        The reciprocal cell within which the band structure is calculated.
        """,
    )
    segments = Quantity(
        type=KBandSegment,
        shape=["*"],
        description="""
        Collection of linear path segments in the reciprocal space. The
        segments are represented as third-order tensors: one dimension for the
        spin channels, one for the sequence of reciprocal space points for the
        segment, and one for the sequence of eigenvalues at a given point.
        """,
    )
    path_standard = Quantity(
        type=str,
        shape=[],
        description="""
        String that identifies the possible standard used in sampling the
        reciprocal space.
        """,
    )


class BandStructurePhonon(BandStructure):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        This section stores information on a vibrational band structure
        evaluation along one-dimensional pathways in the reciprocal space.
        """
    )


class BandStructureElectronic(BandStructure):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        This section stores information on a electonic band structure
        evaluation along one-dimensional pathways in the reciprocal space.
        """
    )
    spin_polarized = Quantity(
        type=bool,
        description="""
        Whether the band structure is spin-polarized, i.e. is contains channels
        for both spin values.
        """,
        a_search=Search(mapping=Boolean()),
    )
    energy_fermi = Quantity(
        type=np.dtype(np.float64),
        unit='joule',
        shape=["n_spin_channels"],
        description="""
        Fermi energy for each spin channel.
        """,
    )
    energy_highest_occupied = Quantity(
        type=np.dtype(np.float64),
        unit='joule',
        shape=["n_spin_channels"],
        description="""
        The highest occupied energy for each spin channel.
        """,
    )
    energy_lowest_unoccupied = Quantity(
        type=np.dtype(np.float64),
        unit='joule',
        shape=["n_spin_channels"],
        description="""
        The lowest unoccupied energy for each spin channel.
        """,
    )
    band_gap = Quantity(
        type=np.dtype(np.float64),
        unit='joule',
        shape=["n_spin_channels"],
        description="""
        Band gap value for each spin channel. If no gap is found, the band gap
        is reported as zero.
        """,
        a_search=Search(mapping=Double),
    )
    band_gap_type = Quantity(
        type=MEnum("direct", "indirect", "no_gap"),
        shape=["n_spin_channels"],
        description="""
        Band gap type for each spin channel.
        """,
        a_search=Search(),
    )


class Properties(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Contains the physical properties that have been calculated or used in
        this entry.
        """
    )
    structure_original = SubSection(
        sub_section=StructureOriginal.m_def,
        repeats=False,
        a_search="structure_original"
    )
    structure_conventional = SubSection(
        sub_section=StructureConventional.m_def,
        repeats=False,
        a_search="structure_conventional"
    )
    structure_primitive = SubSection(
        sub_section=StructurePrimitive.m_def,
        repeats=False,
        a_search="structure_primitive"
    )
    structure_optimized = SubSection(
        sub_section=StructureOptimized.m_def,
        repeats=False,
        a_search="structure_optimized"
    )
    band_structure_electronic = SubSection(sub_section=BandStructureElectronic.m_def, repeats=False, a_search="band_structure_electronic")
    dos_electronic = SubSection(sub_section=DOSElectronic.m_def, repeats=False, a_search="dos_electronic")
    band_structure_phonon = SubSection(sub_section=BandStructurePhonon.m_def, repeats=False, a_search="band_structure_phonon")
    dos_phonon = SubSection(sub_section=DOSPhonon.m_def, repeats=False, a_search="dos_phonon")
    # heat_capacity_constant_volume = SubSection(sub_section=HeatCapacityConstantVolume.m_def, repeats=False, a_search="heat_capacity_constant_volume")
    # helmholtz_free_energy = SubSection(sub_section=HelmholtzFreeEnergy.m_def, repeats=False, a_search="helmholtz_free_energy")


class Results(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        a_search=Search(),
        description="""
        Contains a summary of the entry contents.
        """
    )
    material = SubSection(sub_section=Material.m_def, repeats=False, a_search="material")
    method = SubSection(sub_section=Method.m_def, repeats=False, a_search="method")
    properties = SubSection(sub_section=Properties.m_def, repeats=False, a_search="properties")
