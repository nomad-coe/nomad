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

''' An example metainfo package. '''

import numpy as np
from datetime import datetime

from nomad.units import ureg
from nomad.metainfo import (
    MSection, MCategory, Section, Quantity, Package, SubSection, MEnum,
    Datetime, constraint)

m_package = Package(links=['https://nomad-lab.eu/prod/rae/docs/metainfo.html'])


class SystemHash(MCategory):
    ''' All quantities that contribute to what makes a system unique. '''


class Parsing(MSection):
    ''' All data that describes the NOMAD parsing of this run.

    Quantities can also be documented like this:

    Args:
        parser_name: 'Name of the used parser'
        parser_version: 'Version of the used parser'
    '''

    parser_name = Quantity(type=str)
    parser_version = Quantity(type=str)
    nomad_version = Quantity(type=str, default='latest')
    warnings = Quantity(type=str, shape=['0..*'])
    parse_time = Quantity(type=Datetime)


class System(MSection):
    ''' All data that describes a simulated system. '''

    n_atoms = Quantity(
        type=int, derived=lambda system: len(system.atom_labels),
        description='Number of atoms in the simulated system.')

    atom_labels = Quantity(
        type=str, shape=['n_atoms'], categories=[SystemHash],
        description='The atoms in the simulated systems.')

    atom_positions = Quantity(
        type=np.dtype('f'), shape=['n_atoms', 3], unit=ureg.m, categories=[SystemHash],
        description='The atom positions in the simulated system.')

    lattice_vectors = Quantity(
        type=np.dtype('f'), shape=[3, 3], unit=ureg.m, categories=[SystemHash],
        aliases=['unit_cell'],
        description='The lattice vectors of the simulated unit cell.')

    periodic_dimensions = Quantity(
        type=bool, shape=[3], default=[False, False, False], categories=[SystemHash],
        description='A vector of booleans indicating in which dimensions the unit cell is repeated.')

    system_type = Quantity(type=str)


class SCC(MSection):

    energy_total = Quantity(type=float, default=0.0, unit=ureg.J)
    energy_total_0 = Quantity(type=np.dtype(np.float32), default=0.0, unit=ureg.J)
    an_int = Quantity(type=np.dtype(np.int32))

    system = Quantity(type=System, description='The system that this calculation is based on.')


class Run(MSection):
    ''' All data that belongs to a single code run. '''

    code_name = Quantity(type=str, description='The name of the code that was run.')
    code_version = Quantity(type=str, description='The version of the code that was run.')

    parsing = SubSection(sub_section=Parsing)
    systems = SubSection(sub_section=System, repeats=True)
    sccs = SubSection(sub_section=SCC, repeats=True)

    @constraint
    def one_scc_per_system(self):
        assert self.m_sub_section_count(Run.systems) == self.m_sub_section_count(Run.sccs),\
            'Numbers of system does not match numbers of calculations.'


class VaspRun(Run):
    ''' All VASP specific quantities for section Run. '''
    m_def = Section(extends_base_section=True)

    x_vasp_raw_format = Quantity(
        type=MEnum(['xml', 'outcar']),
        description='The file format of the parsed VASP mainfile.')


if __name__ == '__main__':
    # Demonstration of how to reflect on the definitions

    # All definitions are metainfo data themselves, and they can be accessed like any other
    # metainfo data. E.g. all section definitions are sections themselves.

    # To get quantities of a given section
    print(Run.m_def.m_get_sub_sections(Section.quantities))

    # Or all Sections in the package
    print(m_package.m_get_sub_sections(Package.section_definitions))

    # There are also some definition specific helper methods.
    # For example to get all attributes (Quantities and possible sub-sections) of a section.
    print(Run.m_def.all_properties)

    # Demonstration on how to use the definitions, e.g. to create a run with system:
    run = Run()
    run.code_name = 'VASP'
    run.code_version = '1.0.0'

    parsing = run.m_create(Parsing)
    parsing.parse_time = datetime.now()

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

    # To validate dimensions and custom constraints
    print('errors: %s' % run.m_all_validate())

    # To serialize the data:
    serializable = run.m_to_dict()
    # or
    print(run.m_to_json(indent=2))

    # To deserialize data
    run = Run.m_from_dict(serializable)
    print(run.sccs[0].system)

    # print(m_package.m_to_json(indent=2))  # type: ignore, pylint: disable=undefined-variable
