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

from nomad.metainfo import MSection, Quantity, Package, Datetime, SubSection

m_package = Package()


class QuantumCircuit(MSection):
    processors = Quantity(type=str, shape=['0..*'])
    number_of_registers = Quantity(type=int)
    simulated = Quantity(type=bool)


class QuantumCMS(MSection):
    '''
    The root section for all (meta)data that belongs to a single calculation.
    '''
    chemical_formula = Quantity(
        type=str,
        description=''' The chemical formula that describes the simulated material ''')

    chemical_name = Quantity(
        type=str,
        description=''' The chemical name that describes the simulated material ''')

    atom_labels = Quantity(
        type=str, shape=['1..*'],
        description=''' Labels for the atoms/elements that comprise the simulated material ''')

    transformation = Quantity(type=str)
    quantum_computer_system = Quantity(type=str)
    quantum_computing_libraries = Quantity(type=str, shape=['0..*'])
    computation_datetime = Quantity(type=Datetime)

    number_of_shots = Quantity(type=int)
    quantum_volume = Quantity(type=int)

    section_quantum_circuit = SubSection(sub_section=QuantumCircuit)


m_package.__init_metainfo__()
