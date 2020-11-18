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

'''
Quantum computational materials science metadata
'''

from nomad import config
from nomad.metainfo import Quantity, MSection, Section, Datetime
from nomad.metainfo.search_extension import Search


class QCMSMetadata(MSection):
    m_def = Section(a_domain='qcms')

    # sample quantities
    chemical = Quantity(type=str, default='not processed', a_search=Search())
    quantum_computer_system = Quantity(type=str, a_search=Search())
    quantum_computing_libraries = Quantity(type=str, shape=['0..*'], a_search=Search())
    computation_datetime = Quantity(type=Datetime, a_search=Search())

    # TODO move
    quantities = Quantity(type=str, shape=['0..*'], default=[], a_search=Search())

    def apply_domain_metadata(self, entry_archive):
        if entry_archive is None:
            return

        entry = self.m_parent

        root_section = entry_archive.section_quantum_cms
        entry.formula = root_section.chemical_formula
        if not entry.formula:
            entry.formula = config.services.unavailable_value
        atoms = root_section.atom_labels
        if hasattr(atoms, 'tolist'):
            atoms = atoms.tolist()
        entry.n_atoms = len(atoms)

        atoms = list(set(atoms))
        atoms.sort()
        entry.atoms = atoms

        self.chemical = root_section.chemical_name
        if not self.chemical:
            self.chemical = config.services.unavailable_value

        self.quantum_computer_system = root_section.quantum_computer_system
        if root_section.quantum_computing_libraries is not None:
            self.quantum_computing_libraries = root_section.quantum_computing_libraries
        self.computation_datetime = root_section.computation_datetime

        quantities = set()

        quantities.add(root_section.m_def.name)
        for _, property_def, _ in root_section.m_traverse():
            quantities.add(property_def.name)

        self.quantities = list(quantities)
