""" An example metainfo package. """

import numpy as np

from nomad.metainfo import MSection, MCategory, Section, Quantity, Package, SubSection, Enum, units

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
        type=int, derived=lambda system: len(system.atom_labels),
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
        type=bool, shape=[3], default=[False, False, False], categories=[SystemHash.m_def],
        description='A vector of booleans indicating in which dimensions the unit cell is repeated.')


class SCC(MSection):

    energy_total = Quantity(type=float, default=0.0, unit=units.J)

    system = Quantity(type=System.m_def, description='The system that this calculation is based on.')


class Run(MSection):
    """ All data that belongs to a single code run. """

    code_name = Quantity(type=str, description='The name of the code that was run.')
    code_version = Quantity(type=str, description='The version of the code that was run.')

    parsing = SubSection(sub_section=Parsing.m_def)
    systems = SubSection(sub_section=System.m_def, repeats=True)
    sccs = SubSection(sub_section=SCC.m_def, repeats=True)


class VaspRun(Run):
    """ All VASP specific quantities for section Run. """
    m_def = Section(extends_base_section=True)

    x_vasp_raw_format = Quantity(
        type=Enum(['xml', 'outcar']),
        description='The file format of the parsed VASP mainfile.')


if __name__ == '__main__':
    # Demonstration of how to reflect on the definitions

    # All definitions are metainfo data themselves, and they can be accessed like any other
    # metainfo data. E.g. all section definitions are sections themselves.

    # To get quantities of a given section
    print(Run.m_def.m_get_sub_sections(Section.quantities))

    # Or all Sections in the package
    print(m_package.m_get_sub_sections(Package.section_definitions))  # type: ignore, pylint: disable=undefined-variable

    # There are also some definition specific helper methods.
    # For example to get all attributes (Quantities and possible sub-sections) of a section.
    print(Run.m_def.all_properties)

    # Demonstration on how to use the definitions, e.g. to create a run with system:
    run = Run()
    run.code_name = 'VASP'
    run.code_version = '1.0.0'
    run.m_as(VaspRun).x_vasp_raw_format = 'outcar'
    # The same as
    run.x_vasp_raw_format = 'outcar'  # type: ignore

    system = run.m_create(System)
    system.atom_labels = ['H', 'H', 'O']

    calc = run.m_create(SCC)
    calc.energy_total = 1.23e-10
    calc.system = system

    # Or to read data from existing metainfo data:
    print(system.atom_labels)
    print(system.n_atoms)

    # To serialize the data:
    serializable = run.m_to_dict()
    # or
    print(run.m_to_json(indent=2))

    # To deserialize data
    run = Run.m_from_dict(serializable)
    print(run.sccs[0].system)

    # print(m_package.m_to_json(indent=2))  # type: ignore, pylint: disable=undefined-variable
