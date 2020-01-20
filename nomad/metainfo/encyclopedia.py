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
        type=MEnum(bulk="bulk", two_d="2D", one_d="1D", unavailable="unavailable"),
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
    cell_normalized = Quantity(
        type=np.dtype('f8'),
        shape=[3, 3],
        description="""
        Unit cell in normalized form, meaning the bravais cell. This cell is
        representative and is idealized to match the detected symmetry
        properties.
        """
    )
    cell_primitive = Quantity(
        type=np.dtype('f8'),
        shape=[3, 3],
        description="""
        Definition of the primitive unit cell in a form to be visualized well
        within the normalized cell. This cell is representative and is
        idealized to match the detected symmemtry properties.
        """
    )
    crystal_system = Quantity(
        type=MEnum("triclinic", "monoclinic", "orthorhombic", "tetragonal", "trigonal", "hexagonal", "cubic"),
        description="""
        The detected crystal system. One of seven possibilities in three dimensions.
        """
    )
    formula = Quantity(
        type=str,
        description="""
        Formula giving the composition and occurrences of the elements in the
        Hill notation for the irreducible unit cell.
        """
    )
    formula_reduced = Quantity(
        type=str,
        description="""
        Formula giving the composition and occurrences of the elements in the
        Hill notation for the irreducible unit cell. In this reduced form the
        number of occurences have been divided by the greatest common divisor.
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
    material_name = Quantity(
        type=str,
        description="""
        Most meaningful name for a material.
        """
    )
    periodicity = Quantity(
        type=np.dtype('i1'),
        shape=["1..*"],
        description="""
        The indices of the periodic dimensions.
        """
    )
    point_group = Quantity(
        type=MEnum("1", "-1", "2", "m", "2/m", "222", "mm2", "mmm", "4", "-4", "4/m", "422", "4mm", "-42m", "4/mmm", "3", "-3", "32", "3m", "-3m", "6", "-6", "6/m", "622", "6mm", "-6m2", "6/mmm", "23", "m-3", "432", "-43m", "m-3m"),
        description="""
        Point group in Hermann-Mauguin notation, part of crystal structure
        classification. There are 32 point groups in three dimensional space.
        """
    )
    wyckoff_groups = Quantity(
        type=str,
        description="""
        Returns a list of information about the Wyckoff groups in the JSON format.

        An example of the output:
            [
                {
                    'wyckoff_letter': 'a',
                    'variables': {'z': 0.0},
                    'indices': [0, 6, 12],
                    'element': 'Bi'
                },
                {
                    'wyckoff_letter': 'b',
                    'variables': {'x': 0.50155295, 'z': 0.87461175999999996},
                    'indices': [1, 3, 4, 7, 9, 10, 13, 15, 16],
                    'element': 'Ga'
                }, ...
            ]
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
        unit=units.m**(-3),
        description="""
        Atomic density of the material (atoms/volume)."
        """
    )
    cell_angles_string = Quantity(
        type=str,
        description="""
        A summary of the cell angles, part of material definition.
        """
    )
    cell_volume = Quantity(
        type=float,
        description="""
        Cell volume for a specific calculation. The cell volume can only be
        reported consistently after normalization. Thus it corresponds to the
        normalized cell that is idealized to fit the detected symmetry and may
        not perfectly correspond to the original simulation cell.
        """
    )
    lattice_parameters = Quantity(
        type=np.dtype('f8'),
        shape=[6],
        description="""
        Lattice parameters of a specific calculation. The lattice parameters
        can only be reported consistently after normalization. Thus they
        correspond to the normalized cell that is idealized to fit the detected
        symmetry and may not perfectly correspond to the original simulation
        cell.
        """
    )
    mass_density = Quantity(
        type=float,
        unit=units.kg / units.m**3,
        description="""
        Mass density of the material based on the structural information.
        """
    )


class Encyclopedia(MSection):
    m_def = Section(
        name="encyclopedia",
        a_flask=dict(skip_none=True),
        a_elastic=dict(type=InnerDoc)
    )
    material = SubSection(sub_section=Material.m_def, repeats=False)
    calculation = SubSection(sub_section=Calculation.m_def, repeats=False)
