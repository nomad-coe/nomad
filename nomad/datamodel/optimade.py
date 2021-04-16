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

from ase.data import chemical_symbols
from elasticsearch_dsl import Keyword, Float, InnerDoc, Nested, Integer
import numpy as np

from nomad.units import ureg
from nomad.metainfo import MSection, Section, Quantity, SubSection, MEnum, DefinitionAnnotation
from nomad.metainfo.search_extension import Search


# TODO move the module

def optimade_links(section: str):
    return [
        'https://github.com/Materials-Consortia/OPTiMaDe/blob/develop/optimade.md#%s' %
        section]


class ElementRatio(InnerDoc):
    elements = Keyword()
    elements_ratios = Float()

    @staticmethod
    def from_structure_entry(entry: 'OptimadeEntry'):
        return [
            ElementRatio(elements=entry.elements[i], elements_ratios=entry.elements_ratios[i])
            for i in range(0, entry.nelements)]


class Optimade(DefinitionAnnotation):
    def __init__(self, query: bool = False, entry: bool = False, sortable: bool = False, type: str = None):
        self.sortable = sortable
        self.type = type


class Species(MSection):
    '''
    Used to describe the species of the sites of this structure. Species can be pure
    chemical elements, or virtual-crystal atoms representing a statistical occupation of a
    given site by multiple chemical elements.
    '''

    m_def = Section(links=optimade_links('h.6.2.13'))

    name = Quantity(
        type=str, a_optimade=Optimade(entry=True, sortable=False, type='string'), description='''
            The name of the species; the name value MUST be unique in the species list.
        ''')

    chemical_symbols = Quantity(
        type=MEnum(chemical_symbols + ['x', 'vacancy']), shape=['1..*'],
        a_optimade=Optimade(entry=True, sortable=False, type='list'), description='''
            A list of strings of all chemical elements composing this species.

            It MUST be one of the following:

            - a valid chemical-element name, or
            - the special value "X" to represent a non-chemical element, or
            - the special value "vacancy" to represent that this site has a non-zero probability

            of having a vacancy (the respective probability is indicated in the concentration
            list, see below).

            If any one entry in the species list has a chemical_symbols list that is longer than 1
            element, the correct flag MUST be set in the list structure_features (see
            structure_features)
        ''')

    concentration = Quantity(
        type=float, shape=['1..*'],
        a_optimade=Optimade(entry=True, sortable=False, type='list'), description='''
            A list of floats, with same length as chemical_symbols. The numbers represent the
            relative concentration of the corresponding chemical symbol in this species. The
            numbers SHOULD sum to one. Cases in which the numbers do not sum to one typically fall
            only in the following two categories:

            - Numerical errors when representing float numbers in fixed precision, e.g. for two
            chemical symbols with concentrations 1/3 and 2/3, the concentration might look
            something like [0.33333333333, 0.66666666666]. If the client is aware that the sum
            is not one because of numerical precision, it can renormalize the values so that the
            sum is exactly one.
            - Experimental errors in the data present in the database. In this case, it is the
            responsibility of the client to decide how to process the data.

            Note that concentrations are uncorrelated between different sites (even of the same
            species).
        ''')

    mass = Quantity(type=float, unit=ureg.amu, a_optimade=dict(entry='optional', sortable=False, type='float'))

    original_name = Quantity(type=str, a_optimade=dict(entry='optional', sortable=False, type='string'), description='''
        Can be any valid Unicode string, and SHOULD contain (if specified) the name of the
        species that is used internally in the source database.

        Note: With regards to "source database", we refer to the immediate source being
        queried via the OPTiMaDe API implementation. The main use of this field is for source
        databases that use species names, containing characters that are not allowed (see
        description of the species_at_sites list).
        ''')


class OptimadeEntry(MSection):
    m_def = Section(
        links=optimade_links('h.6.2'),
        a_flask=dict(skip_none=True))

    elements = Quantity(
        type=MEnum(chemical_symbols), shape=['1..*'],
        links=optimade_links('h.6.2.1'),
        a_search=Search(),
        a_optimade=Optimade(query=True, entry=True, sortable=False, type='list'),
        description='''
            Names of the different elements present in the structure.
        ''')

    nelements = Quantity(
        type=int, default=0,
        links=optimade_links('h.6.2.2'),
        a_search=Search(),
        a_optimade=Optimade(query=True, entry=True, sortable=True, type='integer'),
        description='''
            Number of different elements in the structure as an integer.
        ''')

    elements_ratios = Quantity(
        type=float, shape=['nelements'],
        links=optimade_links('h.6.2.3'),
        a_search=Search(mapping=Nested(ElementRatio), value=ElementRatio.from_structure_entry),
        a_optimade=Optimade(query=True, entry=True, sortable=False, type='list'),
        description='''
            Relative proportions of different elements in the structure.
        ''')

    chemical_formula_descriptive = Quantity(
        type=str,
        links=optimade_links('h.6.2.4'),
        a_search=Search(),
        a_optimade=Optimade(query=True, entry=True, sortable=True, type='string'),
        description='''
            The chemical formula for a structure as a string in a form chosen by the API
            implementation.
        ''')

    chemical_formula_reduced = Quantity(
        type=str,
        links=optimade_links('h.6.2.5'),
        a_search=Search(),
        a_optimade=Optimade(query=True, entry=True, sortable=True, type='string'),
        description='''
            The reduced chemical formula for a structure as a string with element symbols and
            integer chemical proportion numbers. The proportion number MUST be omitted if it is 1.
        ''')

    chemical_formula_hill = Quantity(
        type=str,
        links=optimade_links('h.6.2.6'),
        a_search=Search(),
        a_optimade=Optimade(query=True, entry=False, sortable=True, type='string'),
        description='''
            The chemical formula for a structure in Hill form with element symbols followed by
            integer chemical proportion numbers. The proportion number MUST be omitted if it is 1.
        ''')

    chemical_formula_anonymous = Quantity(
        type=str,
        links=optimade_links('h.6.2.7'),
        a_search=Search(),
        a_optimade=Optimade(query=True, entry=True, sortable=True, type='string'),
        description='''
            The anonymous formula is the chemical_formula_reduced, but where the elements are
            instead first ordered by their chemical proportion number, and then, in order left to
            right, replaced by anonymous symbols A, B, C, ..., Z, Aa, Ba, ..., Za, Ab, Bb, ... and
            so on.
        ''')

    dimension_types = Quantity(
        type=int, shape=[3], default=[0, 0, 0],
        links=optimade_links('h.6.2.8'),
        a_optimade=Optimade(query=False, entry=True, sortable=False, type='list'),
        description='''
            List of three integers. For each of the three directions indicated by the three lattice
            vectors (see property lattice_vectors). This list indicates if the direction is
            periodic (value 1) or non-periodic (value 0). Note: the elements in this list each
            refer to the direction of the corresponding entry in lattice_vectors and not
            the Cartesian x, y, z directions.
        ''')

    nperiodic_dimensions = Quantity(
        type=int, derived=lambda a: sum(a.dimension_types),
        links=optimade_links('h.6.2.8'),
        a_search=Search(mapping=Integer()),
        a_optimade=Optimade(query=True, entry=True, sortable=True, type='integer'),
        description='''
            An integer specifying the number of periodic dimensions in the structure, equivalent
            to the number of non-zero entries in dimension_types.
        ''')

    lattice_vectors = Quantity(
        type=np.dtype('float64'), shape=[3, 3], unit=ureg.angstrom,
        links=optimade_links('h.6.2.9'),
        a_optimade=Optimade(query=False, entry=True, sortable=False, type='list'),
        description='''
            The three lattice vectors in Cartesian coordinates, in ångström (Å).
        ''')

    cartesian_site_positions = Quantity(
        type=np.dtype('float64'), shape=['nsites', 3], unit=ureg.angstrom,
        links=optimade_links('h.6.2.10'),
        a_optimade=Optimade(query=False, entry=True, sortable=False, type='list'), description='''
            Cartesian positions of each site. A site is an atom, a site potentially occupied by
            an atom, or a placeholder for a virtual mixture of atoms (e.g., in a virtual crystal
            approximation).
        ''')

    nsites = Quantity(
        type=int, default=0,
        links=optimade_links('h.6.2.11'),
        a_search=Search(),
        a_optimade=Optimade(query=True, entry=True, sortable=True, type='integer'), description='''
            An integer specifying the length of the cartesian_site_positions property.
        ''')

    species_at_sites = Quantity(
        type=str, shape=['nsites'],
        links=optimade_links('h.6.2.12'),
        a_search=Search(),
        a_optimade=Optimade(query=False, entry=True, sortable=False, type='list'), description='''
            Name of the species at each site (where values for sites are specified with the same
            order of the cartesian_site_positions property). The properties of the species are
            found in the species property.
        ''')

    # TODO assemblies

    structure_features = Quantity(
        type=MEnum(['disorder', 'unknown_positions', 'assemblies']), shape=['1..*'], default=[],
        links=optimade_links('h.6.2.15'),
        a_search=Search(),
        a_optimade=Optimade(query=True, entry=True, sortable=False, type='list'), description='''
            A list of strings that flag which special features are used by the structure.

            - disorder: This flag MUST be present if any one entry in the species list has a
            chemical_symbols list that is longer than 1 element.
            - unknown_positions: This flag MUST be present if at least one component of the
            cartesian_site_positions list of lists has value null.
            - assemblies: This flag MUST be present if the assemblies list is present.

        ''')

    species = SubSection(sub_section=Species.m_def, repeats=True)
