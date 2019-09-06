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

"""
Experimental material science specific metadata
"""

from typing import List
import ase.data

from nomad import utils

from .base import CalcWithMetadata, DomainQuantity, Domain, get_optional_backend_value


class EMSEntryWithMetadata(CalcWithMetadata):

    def __init__(self, **kwargs):
        # sample quantities
        self.formula: str = None
        self.atoms: List[str] = []
        self.n_atoms: int = 0
        self.chemical: str = None
        self.sample_constituents: str = None
        self.sample_microstructure: str = None

        # general metadata
        self.experiment_summary: str = None
        self.experiment_location: str = None
        self.experiment_time: str = None

        # method
        self.method: str = None
        self.probing_method: str = None

        # data metadata
        self.repository_name: str = None
        self.repository_url: str = None
        self.preview_url: str = None

        self.quantities = []
        self.group_hash: str = None

        super().__init__(**kwargs)

    def apply_domain_metadata(self, backend):
        logger = utils.get_logger(__name__).bind(
            upload_id=self.upload_id, calc_id=self.calc_id, mainfile=self.mainfile)

        self.formula = get_optional_backend_value(
            backend, 'sample_chemical_formula', 'section_sample', logger=logger)
        self.atoms = get_optional_backend_value(
            backend, 'sample_atom_labels', 'section_sample', logger=logger)
        if hasattr(self.atoms, 'tolist'):
            self.atoms = self.atoms.tolist()
        self.n_atoms = len(self.atoms)
        self.atoms = list(set(self.atoms))
        self.atoms.sort()
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
        self.experiment_time = get_optional_backend_value(
            backend, 'experiment_time', 'section_experiment', logger=logger)

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

        self.group_hash = utils.hash(
            self.formula,
            self.method,
            self.experiment_location,
            self.with_embargo,
            self.comment,
            self.references,
            self.uploader,
            self.coauthors)

        quantities = set()

        for meta_info, _, _ in backend._delegate.results.traverse(root_section='section_experiment'):
            quantities.add(meta_info)

        self.quantities = list(quantities)


Domain(
    'EMS', EMSEntryWithMetadata,
    root_sections=['section_experiment', 'section_entry_info'],
    metainfo_all_package='all.experimental.nomadmetainfo.json',
    quantities=dict(
        formula=DomainQuantity(
            'The chemical (hill) formula of the simulated system.',
            order_default=True),
        atoms=DomainQuantity(
            'The atom labels of all atoms in the simulated system.',
            aggregations=len(ase.data.chemical_symbols)),
        method=DomainQuantity(
            'The experimental method used.', aggregations=20),
        probing_method=DomainQuantity(
            'The used probing method.', aggregations=10),
        sample_microstructure=DomainQuantity(
            'The sample micro structure.', aggregations=10),
        sample_constituents=DomainQuantity(
            'The sample constituents.', aggregations=10),
        quantities=DomainQuantity(
            'All quantities that are used by this calculation')),
    metrics=dict(
        quantities=('quantities', 'value_count')),
    default_statistics=[
        'method', 'probing_method', 'sample_microstructure', 'sample_constituents'])
