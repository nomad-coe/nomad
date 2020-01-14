import numpy as np
from elasticsearch_dsl import InnerDoc
from nomad.metainfo import MSection, Section, SubSection, Quantity, MEnum, units


class Material(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        a_elastic=dict(type=InnerDoc),
        description="""
        Section for storing the data that links this entry into a specific material.
        """
    )
    material_hash = Quantity(
        type=str,
        description="""
        A unique material identifier. For crystals the hash
        identifier is constructed from formula, space group and
        wyckoff_position_population.
        """
    )
    system_type = Quantity(
        type=MEnum(
            bulk="bulk",
            two_d="2D",
            one_d="1D",
            unavailable="unavailable"
        ),
        description="""
        "Character of physical system's geometry, e.g. bulk, surface... ",
        """
    )
    number_of_atoms = Quantity(
        type=int,
        description="""
        Number of atoms in the bravais cell."
        """
    )
    atom_labels = Quantity(
        type=str,
        shape=['1..*'],
        description="""
        Type (element, species) of each atom,
        """
    )
    atom_positions = Quantity(
        type=np.dtype('f8'),
        shape=['number_of_atoms', 3],
        description="""
        Position of each atom, given in relative coordinates.
        """
    )
    bravais_lattice = Quantity(
        type=str,
        description="""
        The Bravais lattice type in the Pearson notation, where the first
        lowercase letter indicates the crystal system, and the second uppercase
        letter indicates the lattice type.

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


class Calculation(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        a_elastic=dict(type=InnerDoc),
        description="""
        Section for storing data related to a calculation that is identified
        from this entry.
        """
    )
    run_type = Quantity(
        type=MEnum(
            single_point="single point",
            geometry_optimization="geometry optimization",
            molecular_dynamics="molecular dynamics",
            phonon_calculation="phonon calculation",
            elastic_constants="elastic constants",
            qha_calculation="QHA calculation",
            qw_calculation="GW calculation",
            equation_of_state="equation of state",
            parameter_variation="parameter variation",
            unavailable="unavailable"),
        description="""
        Defines the type of run identified for this entry.
        """
    )
    atomic_density = Quantity(
        type=float,
        unit=units.m**3,
        description="""
        Atomic density of the material (atoms/volume)."
        """
    )


class Encyclopedia(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        a_elastic=dict(type=InnerDoc)
    )
    material = SubSection(sub_section=Material.m_def, repeats=False)
    calculation = SubSection(sub_section=Calculation.m_def, repeats=False)
