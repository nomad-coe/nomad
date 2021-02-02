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
from elasticsearch_dsl import Text, Keyword

from nomad import config
from nomad.metainfo import MSection, Section, SubSection, Quantity, MEnum
from nomad.metainfo.mongoengine_extension import Mongo, MongoDocument
from nomad.metainfo.search_extension import Search
from nomad.metainfo.elastic_extension import ElasticDocument
from nomad.datamodel.datamodel import MongoMetadata, EditableUserMetadata


class DOSSimilarity(MSection):
    m_def = Section(
        a_mongo=MongoDocument(),
        description="""
        Contains information about the similarity of the density of states of
        this material compared to other materials.
        """,
    )
    material_ids = Quantity(
        type=str,
        shape=["n_similar_materials"],
        a_mongo=Mongo(),
    )
    values = Quantity(
        type=float,
        shape=["n_similar_materials"],
        a_mongo=Mongo(),
    )


class Similarity(MSection):
    m_def = Section(
        a_mongo=MongoDocument(),
        description="""
        Contains information about the similarity of this material to other
        materials.
        """,
    )
    electronic_dos = SubSection(sub_section=DOSSimilarity.m_def, repeats=False)


class Energies(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        a_search="energies",
        description="""
        Contains different types of energies extracted from this entry. The
        energies are extracted from a representative calculation: for geometry
        optimization it is the last optimization step.
        """
    )
    energy_total = Quantity(
        type=np.dtype(np.float64),
        unit="eV",
        description="""
        Total energy.
        """,
        a_search=Search()
    )
    energy_total_T0 = Quantity(
        type=np.dtype(np.float64),
        unit="eV",
        description="""
        Total energy projected to T=0.
        """,
        a_search=Search()
    )
    energy_free = Quantity(
        type=np.dtype(np.float64),
        unit="eV",
        description="""
        Free energy.
        """,
        a_search=Search()
    )


class Properties(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        a_search="properties",
        description="""
        Contains derived physical properties that are specific to the NOMAD
        Encyclopedia.
        """
    )
    atomic_density = Quantity(
        type=np.dtype(np.float64),
        unit="1 / m ** 3",
        description="""
        Atomic density of the material (atoms/volume)."
        """,
        a_search=Search()
    )
    mass_density = Quantity(
        type=np.dtype(np.float64),
        unit="kg / m ** 3",
        description="""
        Mass density of the material.
        """,
        a_search=Search()
    )
    band_gap = Quantity(
        type=np.dtype(np.float64),
        unit="eV",
        description="""
        Band gap value. If multiple spin channels are present, this value is
        taken from the channel with smallest band gap value.
        """,
        a_search=Search()
    )
    band_gap_direct = Quantity(
        type=bool,
        description="""
        Whether band gap is direct or not. If multiple spin channels are
        present, this value is taken from the channel with smallest band gap
        value.
        """,
        a_search=Search()
    )
    energies = SubSection(sub_section=Energies.m_def, repeats=False, a_search=Search())
    has_electronic_band_structure = Quantity(
        type=bool,
        shape=[],
        description="""
        True if the calculation contains an electronic band structure.
        """,
        a_search=Search()
    )
    has_electronic_dos = Quantity(
        type=bool,
        shape=[],
        description="""
        True if the calculation contains an electronic density of states.
        """,
        a_search=Search()
    )
    has_thermodynamical_properties = Quantity(
        type=bool,
        shape=[],
        description="""
        True if the calculation contains thermodynamical properties.
        """,
        a_search=Search()
    )


class Method(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        a_search="method",
        description="""
        Contains an overview of the methodology that was detected in this
        entry.
        """
    )
    program_name = Quantity(
        type=MEnum("ABINIT", "Amber", "ASAP", "ATK", "BAND", "BigDFT", "CASTEP", "Charmm", "CP2K", "CPMD", "Crystal", "DFTB+", "DL_POLY", "DMol3", "elastic", "elk", "exciting", "FHI-aims", "fleur", "fplo", "GAMESS", "Gaussian", "GPAW", "Gromacs", "Gromos", "gulp", "LAMMPS", "libAtoms", "MOLCAS", "MOPAC", "Namd", "NWChem", "Octopus", "ONETEP", "OpenKIM", "ORCA", "Phonopy", "qbox", "Quantum Espresso", "Siesta", "TINKER", "turbomole", "VASP", "WIEN2k"),
        a_search=Search(),
        description="""
        Name of the program used for this calculation.
        """,
    )
    program_version = Quantity(
        type=str,
        a_search=Search(),
        description="""
        Version of the program used for this calculation.
        """,
    )
    basis_set = Quantity(
        type=MEnum('numeric AOs', 'LCAO', 'gaussians', '(L)APW+lo', 'plane waves', 'psinc functions', 'real-space grid'),
        description='The used basis set functions.',
        a_search=Search(),
    )
    method_type = Quantity(
        type=MEnum("DFT", "GW", "unavailable", DFTU="DFT+U"),
        description="""
        Generic name for the used methodology.
        """,
        a_search=Search()
    )
    core_electron_treatment = Quantity(
        type=MEnum("full all electron", "all electron frozen core", "pseudopotential", "unavailable"),
        description="""
        How the core electrons are described.
        """,
        a_search=Search()
    )
    functional_long_name = Quantity(
        type=str,
        description="""
        Full identified for the used exchange-correlation functional.
        """,
        a_search=Search()
    )
    functional_type = Quantity(
        type=MEnum("GGA", "LDA", "hybrid-GGA", "hybrid-meta-GGA", "HF", "GW", "meta-GGA"),
        description="""
        Basic type of the used exchange-correlation functional.
        """,
        a_search=Search()
    )
    method_id = Quantity(
        type=str,
        description="""
        A fixed length, unique method identifier in the form of a hash digest.
        The hash is created by using several method settings as seed. This hash
        is only defined if a set of well-defined method settings is available
        for the used program.
        """
    )
    group_eos_id = Quantity(
        type=str,
        description="""
        A fixed length, unique identifier for equation-of-state calculations.
        Only calculations within the same upload and with a method hash
        available will be grouped under the same hash.
        """,
        a_search=Search()
    )
    group_parametervariation_id = Quantity(
        type=str,
        description="""
        A fixed length, unique identifier for calculations where structure is
        identical but the used computational parameters are varied. Only
        calculations within the same upload and with a method hash available
        will be grouped under the same hash.
        """,
        a_search=Search()
    )
    gw_starting_point = Quantity(
        type=str,
        description="""
        The exchange-correlation functional that was used as a starting point
        for this GW calculation.
        """
    )
    gw_type = Quantity(
        type=MEnum("G0W0", "scGW"),
        description="""
        Basic type of GW calculation.
        """
    )
    smearing_kind = Quantity(
        type=str,
        description="""
        Smearing function used for the electronic structure calculation.
        """
    )
    smearing_parameter = Quantity(
        type=np.dtype(np.float64),
        description="""
        Parameter for smearing, usually the width.
        """
    )


class Workflow(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        a_search="calculation",
        description="""
        Contains an overview of the type of workflow that was detected in
        this entry.
        """
    )
    workflow_type = Quantity(
        type=MEnum(
            single_point="single point",
            geometry_optimization="geometry optimization",
            molecular_dynamics="molecular dynamics",
            phonon_calculation="phonon calculation",
            elastic_constants="elastic constants",
            qha_calculation="QHA calculation",
            gw_calculation="GW calculation",
            equation_of_state="equation of state",
            parameter_variation="parameter variation",
            unavailable="unavailable"),
        description="""
        Defines the type of workflow that was detected for this entry.
        """,
        a_search=Search()
    )


class WyckoffVariables(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Contains the variables associated with a Wyckoff set.
        """
    )
    x = Quantity(
        type=np.dtype(np.float64),
        description="""
        The x variable if present.
        """
    )
    y = Quantity(
        type=np.dtype(np.float64),
        description="""
        The y variable if present.
        """
    )
    z = Quantity(
        type=np.dtype(np.float64),
        description="""
        The z variable if present.
        """
    )


class WyckoffSet(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Section for storing Wyckoff set information.
        """
    )
    wyckoff_letter = Quantity(
        type=str,
        description="""
        The Wyckoff letter for this set.
        """
    )
    indices = Quantity(
        type=np.dtype('i4'),
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
    variables = SubSection(sub_section=WyckoffVariables.m_def, repeats=False)


class LatticeParameters(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Lattice parameters of the idealized cell. The lattice parameters can
        only be reported consistently after idealization and may not perfectly
        correspond to the original simulation cell.
        """,
    )
    a = Quantity(
        type=float,
        description="""
        Length of the first basis vector.
        """,
    )
    b = Quantity(
        type=float,
        description="""
        Length of the second basis vector.
        """,
    )
    c = Quantity(
        type=float,
        description="""
        Length of the third basis vector.
        """,
    )
    alpha = Quantity(
        type=float,
        description="""
        Angle between second and third basis vector.
        """,
    )
    beta = Quantity(
        type=float,
        description="""
        Angle between first and third basis vector.
        """,
    )
    gamma = Quantity(
        type=float,
        description="""
        Angle between first and second basis vector.
        """,
    )


class IdealizedStructure(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Contains structural information for an idealized representation of the
        material used in the calculation. This idealization is used for
        visualizing the material and for calculating the structural properties.
        The properties of the idealized structure may slightly vary from the
        original structure used in the calculation.
        """,
        a_search="idealized_structure",
    )
    atom_labels = Quantity(
        type=str,
        shape=['1..*'],
        description="""
        Type (element, species) of each atom.
        """
    )
    atom_positions = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        description="""
        Atom positions given in coordinates that are relative to the idealized
        cell.
        """
    )
    lattice_vectors = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        description="""
        Lattice vectors of the idealized structure. For bulk materials it is
        the Bravais cell. This cell is representative and is idealized to match
        the detected symmetry properties.
        """
    )
    lattice_vectors_primitive = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        description="""
        Lattice vectors of the the primitive unit cell in a form to be visualized
        within the idealized cell. This cell is representative and is
        idealized to match the detected symmemtry properties.
        """
    )
    periodicity = Quantity(
        type=np.bool_,
        shape=[3],
        description="""
        Automatically detected true periodicity of each lattice direction. May
        not correspond to the periodicity used in the calculation.
        """
    )
    number_of_atoms = Quantity(
        type=int,
        description="""
        Number of atoms in the idealized structure."
        """
    )
    cell_volume = Quantity(
        type=np.dtype(np.float64),
        unit="m ** 3",
        description="""
        Volume of the idealized cell. The cell volume can only be reported
        consistently after idealization and may not perfectly correspond to the
        original simulation cell.
        """,
        a_search=Search()
    )
    wyckoff_sets = SubSection(sub_section=WyckoffSet.m_def, repeats=True)
    lattice_parameters = SubSection(sub_section=LatticeParameters.m_def)


class Calculation(MSection):
    m_def = Section(
        description="""
        List of all calculations related to this material.
        """
    )
    calc_id = Quantity(
        type=str,
        a_search=Search()
    )
    upload_id = Quantity(
        type=str,
        a_search=Search()
    )
    published = Quantity(
        type=bool, default=False,
        description='Indicates if the entry is published',
        categories=[MongoMetadata],
        a_search=Search()
    )
    with_embargo = Quantity(
        type=bool, default=False, categories=[MongoMetadata, EditableUserMetadata],
        description='Indicated if this entry is under an embargo',
        a_search=Search()
    )
    owners = Quantity(
        type=str, shape=['0..*'],
        description='The ids of the users that have access (upload and shared with users) to this calculation',
        a_search=Search(many_or='append')
    )
    method = SubSection(sub_section=Method.m_def, repeats=False, a_search="method")
    workflow = SubSection(sub_section=Workflow.m_def, repeats=False, a_search="workflow")
    properties = SubSection(sub_section=Properties.m_def, repeats=False, a_search="properties")
    idealized_structure = SubSection(sub_section=IdealizedStructure.m_def, repeats=False, a_search="idealized_structure")


class Bulk(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Contains information that is specific to bulk crystalline materials.
        """
    )
    bravais_lattice = Quantity(
        type=MEnum("aP", "mP", "mS", "oP", "oS", "oI", "oF", "tP", "tI", "hR", "hP", "cP", "cI", "cF"),
        description="""
        The Bravais lattice type in the Pearson notation, where the first
        lowercase letter indicates the crystal system, and the second uppercase
        letter indicates the lattice type. The value can only be one of the 14
        different Bravais lattices in three dimensions.

        Crystal system letters:

        a = Triclinic
        m = Monoclinic
        o = Orthorhombic
        t = Tetragonal
        h = Hexagonal and Trigonal
        c = Cubic

        Lattice type letters:

        P = Primitive
        S (A, B, C) = One side/face centred
        I = Body centered
        R = Rhombohedral centring
        F = All faces centred
        """,
        a_search=Search()
    )
    crystal_system = Quantity(
        type=MEnum("triclinic", "monoclinic", "orthorhombic", "tetragonal", "trigonal", "hexagonal", "cubic"),
        description="""
        The detected crystal system. One of seven possibilities in three dimensions.
        """,
        a_search=Search()
    )
    has_free_wyckoff_parameters = Quantity(
        type=bool,
        description="""
        Whether the material has any Wyckoff sites with free parameters. If a
        materials has free Wyckoff parameters, at least some of the atoms are
        not bound to a particular location in the structure but are allowed to
        move with possible restrictions set by the symmetry.
        """,
        a_search=Search()
    )
    point_group = Quantity(
        type=MEnum("1", "-1", "2", "m", "2/m", "222", "mm2", "mmm", "4", "-4", "4/m", "422", "4mm", "-42m", "4/mmm", "3", "-3", "32", "3m", "-3m", "6", "-6", "6/m", "622", "6mm", "-6m2", "6/mmm", "23", "m-3", "432", "-43m", "m-3m"),
        description="""
        Point group in Hermann-Mauguin notation, part of crystal structure
        classification. There are 32 point groups in three dimensional space.
        """,
        a_search=Search()
    )
    space_group_number = Quantity(
        type=int,
        description="""
        Integer representation of the space group, part of crystal structure
        classification, part of material definition.
        """,
        a_search=Search()
    )
    space_group_international_short_symbol = Quantity(
        type=str,
        description="""
        International short symbol notation of the space group.
        """,
        a_search=Search()
    )
    structure_prototype = Quantity(
        type=str,
        description="""
        The prototypical material for this crystal structure.
        """,
        a_search=Search()
    )
    structure_type = Quantity(
        type=str,
        description="""
        Classification according to known structure type, considering the point
        group of the crystal and the occupations with different atom types.
        """,
        a_search=Search()
    )
    strukturbericht_designation = Quantity(
        type=str,
        description="""
        Classification of the material according to the historically grown "strukturbericht".
        """,
        a_search=Search()
    )
    program_name = Quantity(
        type=str,
        description="""
        The prototypical material for this crystal structure.
        """,
        a_search=Search()
    )
    structure_type = Quantity(
        type=str,
        description="""
        Classification according to known structure type, considering the point
        group of the crystal and the occupations with different atom types.
        """,
        a_search=Search()
    )
    strukturbericht_designation = Quantity(
        type=str,
        description="""
        Classification of the material according to the historically grown "strukturbericht".
        """,
        a_search=Search()
    )


class Material(MSection):
    m_def = Section(
        a_elastic=ElasticDocument(index_name=config.elastic.materials_index_name, id=lambda x: x.material_id),
        a_mongo=MongoDocument()
    )
    material_id = Quantity(
        type=str,
        description="""
        A fixed length, unique material identifier in the form of a hash
        digest.
        """,
        a_mongo=Mongo(primary_key=True),
        a_search=Search()
    )
    material_type = Quantity(
        type=MEnum(bulk="bulk", two_d="2D", one_d="1D"),
        description="""
        Broad structural classification for the material.
        """,
        a_search=Search()
    )
    material_name = Quantity(
        type=str,
        description="""
        Most meaningful name for a material if one could be assigned
        """,
        a_search=Search()
    )
    material_classification = Quantity(
        type=str,
        description="""
        Contains the compound class and classification of the material
        according to springer materials in JSON format.
        """
    )
    formula = Quantity(
        type=str,
        description="""
        Formula giving the composition and occurrences of the elements in the
        Hill notation. For periodic materials the formula is calculated fom the
        primitive unit cell.
        """,
        a_search=Search()
    )
    formula_reduced = Quantity(
        type=str,
        description="""
        Formula giving the composition and occurrences of the elements in the
        Hill notation where the number of occurences have been divided by the
        greatest common divisor.
        """,
        a_search=Search()
    )
    species_and_counts = Quantity(
        type=str,
        description="""
        The formula separated into individual terms containing both the atom
        type and count. Used for searching parts of a formula.
        """,
        a_search=Search(mapping=Text(multi=True, fields={'keyword': Keyword()}))
    )
    species = Quantity(
        type=str,
        description="""
        The formula separated into individual terms containing only unique atom
        species. Used for searching materials containing specific elements.
        """,
        a_search=Search(mapping=Text(multi=True, fields={'keyword': Keyword()}))
    )
    similarity = SubSection(sub_section=Similarity.m_def, repeats=False)
    calculations = SubSection(
        sub_section=Calculation.m_def,
        repeats=True,
        a_search=Search("calculations", nested=True),
        description="""
        List of all calculations related to this material.
        """
    )
    bulk = SubSection(sub_section=Bulk.m_def, repeats=False, a_search="bulk")
