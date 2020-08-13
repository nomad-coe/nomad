# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

'''
Experimental material science specific metadata
'''

from nomad import config
from nomad.metainfo import Quantity, MSection, Section, Datetime
from nomad.metainfo.search_extension import Search


class EMSMetadata(MSection):
    m_def = Section(a_domain='ems')

    # sample quantities
    chemical = Quantity(type=str, default='not processed', a_search=Search())
    sample_constituents = Quantity(type=str, default='not processed', a_search=Search())
    sample_microstructure = Quantity(type=str, default='not processed', a_search=Search())

    # general metadata
    experiment_summary = Quantity(type=str, default='not processed', a_search=Search())
    experiment_location = Quantity(type=str, default='not processed', a_search=Search())
    experiment_time = Quantity(type=Datetime, a_search=Search())

    # method
    method = Quantity(type=str, default='not processed', a_search=Search())
    probing_method = Quantity(type=str, default='not processed', a_search=Search())

    # data metadata
    repository_name = Quantity(type=str, default='not processed', a_search=Search())
    repository_url = Quantity(type=str, default='not processed', a_search=Search())
    entry_repository_url = Quantity(type=str, default='not processed', a_search=Search())
    preview_url = Quantity(type=str, default='not processed', a_search=Search())

    # TODO move
    quantities = Quantity(type=str, shape=['0..*'], default=[], a_search=Search())
    group_hash = Quantity(type=str, a_search=Search())

    def apply_domain_metadata(self, backend):
        from nomad import utils

        if backend is None:
            return

        entry = self.m_parent

        root_section = backend.entry_archive.section_experiment
        entry.formula = root_section.section_sample[0].sample_chemical_formula
        atoms = root_section.section_sample[0].sample_atom_labels
        if hasattr(atoms, 'tolist'):
            atoms = atoms.tolist()
        entry.n_atoms = len(atoms)

        atoms = list(set(atoms))
        atoms.sort()
        entry.atoms = atoms

        self.chemical = root_section.section_sample[0].sample_chemical_name
        self.sample_microstructure = root_section.section_sample[0].sample_microstructure
        self.sample_constituents = root_section.section_sample[0].sample_constituents

        self.experiment_summary = root_section.experiment_summary
        self.experiment_location = root_section.experiment_location
        experiment_time = root_section.experiment_time
        if experiment_time != config.services.unavailable_value:
            self.experiment_time = experiment_time

        self.method = root_section.section_method[0].experiment_method_name
        self.probing_method = root_section.section_method[0].probing_method

        self.repository_name = root_section.section_data[0].data_repository_name
        self.repository_url = root_section.section_data[0].data_repository_url
        self.preview_url = root_section.section_data[0].data_preview_url
        self.entry_repository_url = root_section.section_data[0].entry_repository_url

        self.group_hash = utils.hash(
            entry.formula,
            self.method,
            self.experiment_location,
            entry.with_embargo,
            entry.uploader)

        quantities = set()

        quantities.add(root_section.m_def.name)
        for _, property_def, _ in root_section.m_traverse():
            quantities.add(property_def.name)

        self.quantities = list(quantities)
