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
Experimental material science specific metadata
'''
from nomad import config
from nomad.metainfo import Quantity, MSection, Section, Datetime
from nomad.metainfo.search_extension import Search


def _unavailable(value):
    if value is None:
        return config.services.unavailable_value

    return value


class EMSMetadata(MSection):
    m_def = Section(a_domain='ems')

    # sample quantities
    chemical = Quantity(type=str, a_search=Search())
    sample_constituents = Quantity(type=str, a_search=Search())
    sample_microstructure = Quantity(type=str, a_search=Search())

    # general metadata
    experiment_summary = Quantity(type=str, a_search=Search())
    origin_time = Quantity(type=Datetime, a_search=Search())
    experiment_location = Quantity(type=str, a_search=Search())

    # method
    method = Quantity(type=str, a_search=Search())
    data_type = Quantity(type=str, a_search=Search())
    probing_method = Quantity(type=str, a_search=Search())

    # data metadata
    repository_name = Quantity(type=str, a_search=Search())
    repository_url = Quantity(type=str, a_search=Search())
    entry_repository_url = Quantity(type=str, a_search=Search())
    preview_url = Quantity(type=str, a_search=Search())

    # TODO move
    quantities = Quantity(type=str, shape=['0..*'], default=[], a_search=Search())
    group_hash = Quantity(type=str, a_search=Search())

    def apply_domain_metadata(self, entry_archive):
        from nomad import utils

        if entry_archive is None:
            return

        entry = self.m_parent

        root_section = entry_archive.section_experiment
        entry.formula = root_section.section_sample.section_material.chemical_formula
        atoms = root_section.section_sample.section_material.atom_labels

        if atoms is None:
            entry.atoms = []
        else:
            if hasattr(atoms, 'tolist'):
                atoms = atoms.tolist()
            entry.n_atoms = len(atoms)

            atoms = list(set(atoms))
            atoms.sort()
            entry.atoms = atoms

        self.chemical = _unavailable(root_section.section_sample.section_material.chemical_name)
        self.sample_microstructure = _unavailable(root_section.section_sample.sample_microstructure)
        self.sample_constituents = _unavailable(root_section.section_sample.sample_constituents)

        self.experiment_summary = root_section.experiment_summary
        location = root_section.experiment_location
        if location is not None:
            location_str = ', '.join([
                getattr(location, prop)
                for prop in ['facility', 'institution', 'address']
                if getattr(location, prop) is not None])
            self.experiment_location = location_str

        if root_section.experiment_time:
            self.origin_time = root_section.experiment_time
        elif root_section.experiment_publish_time:
            self.origin_time = root_section.experiment_publish_time
        else:
            self.origin_time = self.m_parent.upload_time

        self.data_type = _unavailable(root_section.section_method.data_type)
        self.method = _unavailable(root_section.section_method.method_name)
        self.probing_method = _unavailable(root_section.section_method.probing_method)

        self.repository_name = _unavailable(root_section.section_data.repository_name)
        self.repository_url = root_section.section_data.repository_url
        self.preview_url = root_section.section_data.preview_url
        self.entry_repository_url = root_section.section_data.entry_repository_url

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

        if self.m_parent.references is None:
            self.m_parent.references = [self.entry_repository_url]
