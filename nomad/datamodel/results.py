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
from elasticsearch_dsl import Text, Keyword, Integer

from ase.data import chemical_symbols

from nomad.metainfo import MSection, Section, SubSection, Quantity, MEnum
from nomad.metainfo.search_extension import Search
from nomad import atomutils

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
        Section containing information about the symmetry properties of the
        system.
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
        derived=lambda a: atomutils.get_volume(a.lattice_vectors.magnitude),
        unit="m ** 3",
        description="""
        Volume of the cell. The cell volume can only be reported consistently
        after idealization and may not perfectly correspond to the original
        simulation cell.
        """,
        a_search=Search()
    )
    lattice_parameters = SubSection(sub_section=LatticeParameters.m_def)
    wyckoff_sets = SubSection(sub_section=WyckoffSet.m_def, repeats=True)


class Symmetry(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        description="""
        Section containing information about the symmetry/crystallography
        related properties of the system.
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
    prototype = Quantity(
        type=str,
        description="""
        The prototypical material for this crystal structure.
        """,
        a_search=Search()
    )
    prototype_aflow_label = Quantity(
        type=str,
        description="""
        The label of this structure in the AFLOW encyclopedia of
        crystallographic prototypes:
        http://www.aflowlib.org/prototype-encyclopedia/index.html
        """,
        a_search=Search()
    )
    structure_name = Quantity(
        type=str,
        description="""
        A common name for this structure.
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
    structural_type = Quantity(
        type=MEnum(structure_classes), default="not processed",
        description="""
        Classification based on structural features.
        """,
        a_search=Search(statistic_values=structure_classes)
    )
    functional_type = Quantity(
        type=str,
        shape=["*"],
        description="""
        Classification based on the functional properties.
        """,
        a_search=Search()
    )
    compound_type = Quantity(
        type=str,
        shape=["*"],
        description="""
        Classification based on the chemical formula.
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
    structure_original = SubSection(sub_section=Structure.m_def, repeats=False, a_search="structure_original")
    structure_conventional = SubSection(sub_section=Structure.m_def, repeats=False, a_search="structure_conventional")
    structure_primitive = SubSection(sub_section=Structure.m_def, repeats=False, a_search="structure_primitive")


class Results(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        a_search=Search(),
        description="""
        Contains a summary of the entry contents.
        """
    )
    material = SubSection(sub_section=Material.m_def, repeats=False, a_search="material")
