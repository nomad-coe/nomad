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

from .common import get_optional_backend_value


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
        logger = utils.get_logger(__name__).bind(
            upload_id=entry.upload_id, calc_id=entry.calc_id, mainfile=entry.mainfile)

        entry.formula = get_optional_backend_value(
            backend, 'sample_chemical_formula', 'section_sample', logger=logger)
        atoms = get_optional_backend_value(
            backend, 'sample_atom_labels', 'section_sample', logger=logger)
        if hasattr(atoms, 'tolist'):
            atoms = atoms.tolist()
        entry.n_atoms = len(atoms)

        atoms = list(set(atoms))
        atoms.sort()
        entry.atoms = atoms

        self.chemical = get_optional_backend_value(
            backend, 'sample_chemical_name', 'section_sample', logger=logger)
        self.sample_microstructure = get_optional_backend_value(
            backend, 'sample_microstructure', 'section_sample', logger=logger)
        self.sample_constituents = get_optional_backend_value(
            backend, 'sample_constituents', 'section_sample', logger=logger)

        self.experiment_summary = get_optional_backend_value(
            backend, 'experiment_summary', 'section_experiment', logger=logger)
        self.experiment_location = get_optional_backend_value(
            backend, 'experiment_location', 'section_experiment', logger=logger)
        experiment_time = get_optional_backend_value(
            backend, 'experiment_time', 'section_experiment', None, logger=logger)
        if experiment_time != config.services.unavailable_value:
            self.experiment_time = experiment_time

        self.method = get_optional_backend_value(
            backend, 'experiment_method_name', 'section_method', logger=logger)
        self.probing_method = get_optional_backend_value(
            backend, 'probing_method', 'section_method', logger=logger)

        self.repository_name = get_optional_backend_value(
            backend, 'data_repository_name', 'section_data', logger=logger)
        self.repository_url = get_optional_backend_value(
            backend, 'data_repository_url', 'section_data', logger=logger)
        self.preview_url = get_optional_backend_value(
            backend, 'data_preview_url', 'section_data', logger=logger)
        self.entry_repository_url = get_optional_backend_value(
            backend, 'entry_repository_url', 'section_data', logger=logger)

        self.group_hash = utils.hash(
            entry.formula,
            self.method,
            self.experiment_location,
            entry.with_embargo,
            entry.uploader)

        quantities = set()

        root_section = backend.entry_archive.section_experiment
        quantities.add(root_section.m_def.name)
        for _, property_def, _ in root_section.m_traverse():
            quantities.add(property_def.name)

        self.quantities = list(quantities)
