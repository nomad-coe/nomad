import numpy as np
from elasticsearch_dsl import InnerDoc
from nomad.metainfo import MSection, Section, SectionProxy, SubSection, Quantity, Reference, MEnum, units
from nomad.datamodel.metainfo.public import section_dos, section_k_band


class WyckoffVariables(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        a_elastic=dict(type=InnerDoc),
        description="""
        Contains the variables associated with a Wyckoff set.
        """
    )
    x = Quantity(
        type=float,
        description="""
        The x variable if present.
        """
    )
    y = Quantity(
        type=float,
        description="""
        The y variable if present.
        """
    )
    z = Quantity(
        type=float,
        description="""
        The z variable if present.
        """
    )


class WyckoffSet(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        a_elastic=dict(type=InnerDoc),
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


class IdealizedStructure(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        a_elastic=dict(type=InnerDoc),
        description="""
        Contains structural information for an idealized representation of the
        material used in the calculation. This idealization is used for
        visualizing the material and for calculating the structural properties.
        The properties of the idealized structure may slightly vary from the
        original structure used in the calculation.
        """
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
    lattice_parameters = Quantity(
        type=np.dtype(np.float64),
        shape=[6],
        description="""
        Lattice parameters of the idealized cell. The lattice parameters can
        only be reported consistently after idealization and may not perfectly
        correspond to the original simulation cell.
        """
    )
    periodicity = Quantity(
        type=np.bool,
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
        type=float,
        description="""
        Volume of the idealized cell. The cell volume can only be reported
        consistently after idealization and may not perfectly correspond to the
        original simulation cell.
        """
    )


class Bulk(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        a_elastic=dict(type=InnerDoc),
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
        """
    )
    crystal_system = Quantity(
        type=MEnum("triclinic", "monoclinic", "orthorhombic", "tetragonal", "trigonal", "hexagonal", "cubic"),
        description="""
        The detected crystal system. One of seven possibilities in three dimensions.
        """
    )
    has_free_wyckoff_parameters = Quantity(
        type=bool,
        description="""
        Whether the material has any Wyckoff sites with free parameters. If a
        materials has free Wyckoff parameters, at least some of the atoms are
        not bound to a particular location in the structure but are allowed to
        move with possible restrictions set by the symmetry.
        """
    )
    point_group = Quantity(
        type=MEnum("1", "-1", "2", "m", "2/m", "222", "mm2", "mmm", "4", "-4", "4/m", "422", "4mm", "-42m", "4/mmm", "3", "-3", "32", "3m", "-3m", "6", "-6", "6/m", "622", "6mm", "-6m2", "6/mmm", "23", "m-3", "432", "-43m", "m-3m"),
        description="""
        Point group in Hermann-Mauguin notation, part of crystal structure
        classification. There are 32 point groups in three dimensional space.
        """
    )
    space_group_number = Quantity(
        type=int,
        description="""
        Integer representation of the space group, part of crystal structure
        classification, part of material definition.
        """
    )
    space_group_international_short_symbol = Quantity(
        type=str,
        description="""
        International short symbol notation of the space group.
        """
    )
    structure_prototype = Quantity(
        type=str,
        description="""
        The prototypical material for this crystal structure.
        """
    )
    structure_type = Quantity(
        type=str,
        description="""
        Classification according to known structure type, considering the point
        group of the crystal and the occupations with different atom types.
        """
    )
    strukturbericht_designation = Quantity(
        type=str,
        description="""
        Classification of the material according to the historically grown "strukturbericht".
        """
    )
    wyckoff_sets = SubSection(sub_section=WyckoffSet.m_def, repeats=True)


class Material(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        a_elastic=dict(type=InnerDoc),
        description="""
        Contains an overview of the type of material that was detected in this
        entry.
        """
    )
    material_type = Quantity(
        type=MEnum(bulk="bulk", two_d="2D", one_d="1D", unavailable="unavailable"),
        description="""
        "Broad structural classification for the material, e.g. bulk, 2D, 1D... ",
        """
    )
    material_hash = Quantity(
        type=str,
        description="""
        A fixed length, unique material identifier in the form of a hash
        digest.
        """
    )
    material_name = Quantity(
        type=str,
        description="""
        Most meaningful name for a material.
        """
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
        """
    )
    formula_reduced = Quantity(
        type=str,
        description="""
        Formula giving the composition and occurrences of the elements in the
        Hill notation whre the number of occurences have been divided by the
        greatest common divisor.
        """
    )

    # The idealized structure for this material
    idealized_structure = SubSection(sub_section=IdealizedStructure.m_def, repeats=False)

    # Bulk-specific properties
    bulk = SubSection(sub_section=Bulk.m_def, repeats=False)


class Method(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        a_elastic=dict(type=InnerDoc),
        description="""
        Contains an overview of the methodology that was detected in this
        entry.
        """
    )
    method_type = Quantity(
        type=MEnum("DFT", "GW", "unavailable", DFTU="DFT+U"),
        description="""
        Generic name for the used methodology.
        """
    )
    basis_set_type = Quantity(
        type=MEnum("Numeric AOs", "Gaussians", "(L)APW+lo", "FLAPW (full-potential linearized augmented planewave)", "Plane waves", "Real-space grid", "Local-orbital minimum-basis"),
        description="""
        Basic type of the used basis set.
        """
    )
    core_electron_treatment = Quantity(
        type=MEnum("full all electron", "all electron frozen core", "pseudopotential", "unavailable"),
        description="""
        How the core electrons are described.
        """
    )
    functional_long_name = Quantity(
        type=str,
        description="""
        Full identified for the used exchange-correlation functional.
        """
    )
    functional_type = Quantity(
        type=str,
        description="""
        Basic type of the used exchange-correlation functional.
        """
    )
    method_hash = Quantity(
        type=str,
        description="""
        A fixed length, unique method identifier in the form of a hash digest.
        The hash is created by using several method settings as seed. This hash
        is only defined if a set of well-defined method settings is available
        for the used program.
        """
    )
    group_eos_hash = Quantity(
        type=str,
        description="""
        A fixed length, unique identifier for equation-of-state calculations.
        Only calculations within the same upload and with a method hash
        available will be grouped under the same hash.
        """
    )
    group_parametervariation_hash = Quantity(
        type=str,
        description="""
        A fixed length, unique identifier for calculations where structure is
        identical but the used computational parameters are varied. Only
        calculations within the same upload and with a method hash available
        will be grouped under the same hash.
        """
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
        type=float,
        description="""
        Parameter for smearing, usually the width.
        """
    )


class Calculation(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        a_elastic=dict(type=InnerDoc),
        description="""
        Contains an overview of the type of calculation that was detected in
        this entry.
        """
    )
    calculation_type = Quantity(
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
        Defines the type of calculation that was detected for this entry.
        """
    )


class Properties(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        a_elastic=dict(type=InnerDoc),
        description="""
        Contains a summary of the physical properties that have been calculated
        in this entry.
        """
    )
    atomic_density = Quantity(
        type=float,
        unit=units.m**(-3),
        description="""
        Atomic density of the material (atoms/volume)."
        """
    )
    mass_density = Quantity(
        type=float,
        unit=units.kg / units.m**3,
        description="""
        Mass density of the material.
        """
    )
    energies = Quantity(
        type=str,
        description="""
        Code dependent energy values, corrected to be per formula unit.
        """
    )
    electronic_band_structure = Quantity(
        type=Reference(section_k_band.m_def),
        shape=[],
        description="""
        Reference to an electronic band structure.
        """
    )
    electronic_dos = Quantity(
        type=Reference(section_dos.m_def),
        shape=[],
        description="""
        Reference to an electronic density of states.
        """
    )
    phonon_band_structure = Quantity(
        type=Reference(section_k_band.m_def),
        shape=[],
        description="""
        Reference to a phonon band structure.
        """
    )
    phonon_dos = Quantity(
        type=Reference(section_dos.m_def),
        shape=[],
        description="""
        Reference to a phonon density of states.
        """
    )


class section_encyclopedia(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        a_elastic=dict(type=InnerDoc),
        description="""
        Section which stores information for the NOMAD Encyclopedia.
        """
    )
    material = SubSection(sub_section=Material.m_def, repeats=False)
    method = SubSection(sub_section=Method.m_def, repeats=False)
    properties = SubSection(sub_section=Properties.m_def, repeats=False)
    calculation = SubSection(sub_section=Calculation.m_def, repeats=False)
