""" An example metainfo package. """

import numpy as np

from nomad.metainfo import MSection, MCategory, Section, Quantity, Enum, Package, SubSection, units

m_package = Package(links=['http://metainfo.nomad-coe.eu'])


class SystemHash(MCategory):
    """ All quantities that contribute to what makes a system unique. """


class Parsing(MSection):
    """ All data that describes the NOMAD parsing of this run. """

    parser_name = Quantity(type=str)
    parser_version = Quantity(type=str)
    nomad_version = Quantity(type=str)
    warnings = Quantity(type=str, shape=['0..*'])


class System(MSection):
    """ All data that describes a simulated system. """

    n_atoms = Quantity(
        type=int, default=0,
        description='Number of atoms in the simulated system.')

    atom_labels = Quantity(
        type=str, shape=['n_atoms'], categories=[SystemHash.m_def],
        description='The atoms in the simulated systems.')

    atom_positions = Quantity(
        type=np.dtype('f'), shape=['n_atoms', 3], unit=units.m, categories=[SystemHash.m_def],
        description='The atom positions in the simulated system.')

    lattice_vectors = Quantity(
        type=np.dtype('f'), shape=[3, 3], unit=units.m, categories=[SystemHash.m_def],
        description='The lattice vectors of the simulated unit cell.')

    unit_cell = Quantity(synonym_for='lattice_vectors')

    periodic_dimensions = Quantity(
        type=bool, shape=[3], categories=[SystemHash.m_def],
        description='A vector of booleans indicating in which dimensions the unit cell is repeated.')


class Run(MSection):
    """ All data that belongs to a single code run. """

    code_name = Quantity(type=str, description='The name of the code that was run.')
    code_version = Quantity(type=str, description='The version of the code that was run.')

    systems = SubSection(sub_section=System.m_def, repeats=True)
    parsing = SubSection(sub_section=Parsing.m_def)


# class VaspRun(MSection):
#     """ All VASP specific quantities for section Run. """
#     m_def = Section(extends=Run.m_def)

#     x_vasp_raw_format = Quantity(
#         type=Enum(['xml', 'outcar']),
#         description='The file format of the parsed VASP mainfile.')


if __name__ == '__main__':
    # Demonstration of how to reflect on the definitions

    # All definitions are metainfo data themselves, and they can be accessed like any other
    # metainfo data. E.g. all section definitions are sections themselves.

    # To get quantities of a given section
    print(Run.m_def.m_sub_sections(Quantity))

    # Or all Sections in the package
    print(m_package.m_sub_sections(Section))  # type: ignore, pylint: disable=undefined-variable

    # There are also some definition specific helper methods.
    # For example to get all attributes (Quantities and possible sub-sections) of a section.
    print(Run.m_def.all_properties)

    # Demonstration on how to use the definitions, e.g. to create a run with system:
    run = Run()
    run.code_name = 'VASP'
    run.code_version = '1.0.0'

    system = run.m_create(System)
    system.n_atoms = 3
    system.atom_labels = ['H', 'H', 'O']

    # Or to read data from existing metainfo data:
    print(system.atom_labels)

    # To serialize the data:
    print(run.m_to_json(indent=2))

    print(m_package.m_to_json(indent=2))  # type: ignore, pylint: disable=undefined-variable
