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
DFT specific metadata
"""

from typing import List
import re
from elasticsearch_dsl import Integer, Object
import ase.data

from nomadcore.local_backend import ParserEvent

from nomad import utils, config
from nomad.metainfo import optimade

from .base import CalcWithMetadata, DomainQuantity, Domain, get_optional_backend_value


xc_treatments = {
    'gga': 'GGA',
    'hf_': 'HF',
    'oep': 'OEP',
    'hyb': 'hybrid',
    'mgg': 'meta-GGA',
    'vdw': 'vdW',
    'lda': 'LDA',
}
""" https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/XC-functional """

basis_sets = {
    'gaussians': 'gaussians',
    'realspacegrid': 'real-space grid',
    'planewaves': 'plane waves'
}

version_re = re.compile(r'(\d+(\.\d+(\.\d+)?)?)')


def map_functional_name_to_xc_treatment(name):
    if name == config.services.unavailable_value:
        return name

    return xc_treatments.get(name[:3].lower(), name)


def map_basis_set_to_basis_set_label(name):
    key = name.replace('_', '').replace('-', '').replace(' ', '').lower()
    return basis_sets.get(key, name)


def simplify_version(version):
    match = version_re.search(version)
    if match is None:
        return version
    else:
        return match.group(0)


class DFTCalcWithMetadata(CalcWithMetadata):

    def __init__(self, **kwargs):
        self.formula: str = None
        self.atoms: List[str] = []
        self.n_atoms: int = 0
        self.basis_set: str = None
        self.xc_functional: str = None
        self.system: str = None
        self.crystal_system: str = None
        self.spacegroup: str = None
        self.spacegroup_symbol: str = None
        self.code_name: str = None
        self.code_version: str = None

        self.n_geometries = 0
        self.n_calculations = 0
        self.n_total_energies = 0
        self.n_quantities = 0
        self.quantities = []
        self.geometries = []
        self.group_hash: str = None

        self.optimade: optimade.OptimadeStructureEntry = None

        super().__init__(**kwargs)

    def update(self, **kwargs):
        super().update(**kwargs)

        if self.optimade is not None and isinstance(self.optimade, dict):
            self.optimade = optimade.OptimadeStructureEntry.m_from_dict(self.optimade)

    def apply_domain_metadata(self, backend):
        from nomad.normalizing.system import normalized_atom_labels

        logger = utils.get_logger(__name__).bind(
            upload_id=self.upload_id, calc_id=self.calc_id, mainfile=self.mainfile)

        self.code_name = backend.get_value('program_name', 0)
        try:
            self.code_version = simplify_version(backend.get_value('program_version', 0))
        except KeyError:
            self.code_version = config.services.unavailable_value

        self.atoms = get_optional_backend_value(backend, 'atom_labels', 'section_system', [], logger=logger)
        if hasattr(self.atoms, 'tolist'):
            self.atoms = self.atoms.tolist()
        self.n_atoms = len(self.atoms)
        self.atoms = list(set(normalized_atom_labels(set(self.atoms))))
        self.atoms.sort()

        self.crystal_system = get_optional_backend_value(
            backend, 'crystal_system', 'section_symmetry', logger=logger)
        self.spacegroup = get_optional_backend_value(
            backend, 'space_group_number', 'section_symmetry', 0, logger=logger)
        self.spacegroup_symbol = get_optional_backend_value(
            backend, 'international_short_symbol', 'section_symmetry', 0, logger=logger)
        self.basis_set = map_basis_set_to_basis_set_label(
            get_optional_backend_value(backend, 'program_basis_set_type', 'section_run', logger=logger))
        self.system = get_optional_backend_value(
            backend, 'system_type', 'section_system', logger=logger)
        self.formula = get_optional_backend_value(
            backend, 'chemical_composition_bulk_reduced', 'section_system', logger=logger)
        self.xc_functional = map_functional_name_to_xc_treatment(
            get_optional_backend_value(backend, 'XC_functional_name', 'section_method', logger=logger))

        self.group_hash = utils.hash(
            self.formula,
            self.spacegroup,
            self.basis_set,
            self.xc_functional,
            self.code_name,
            self.code_version,
            self.with_embargo,
            self.comment,
            self.references,
            self.uploader,
            self.coauthors)

        quantities = set()
        geometries = set()
        n_quantities = 0
        n_calculations = 0
        n_total_energies = 0
        n_geometries = 0

        for meta_info, event, value in backend._delegate.results.traverse():
            quantities.add(meta_info)

            if event == ParserEvent.add_value or event == ParserEvent.add_array_value:
                n_quantities += 1

                if meta_info == 'energy_total':
                    n_total_energies += 1

                if meta_info == 'configuration_raw_gid':
                    geometries.add(value)

            elif event == ParserEvent.open_section:
                if meta_info == 'section_single_configuration_calculation':
                    n_calculations += 1

                if meta_info == 'section_system':
                    n_geometries += 1

        self.quantities = list(quantities)
        self.geometries = list(geometries)
        self.n_quantities = n_quantities
        self.n_calculations = n_calculations
        self.n_total_energies = n_total_energies
        self.n_geometries = n_geometries

        self.optimade = backend.get_mi2_section(optimade.OptimadeStructureEntry.m_def)


def only_atoms(atoms):
    numbers = [ase.data.atomic_numbers[atom] for atom in atoms]
    only_atoms = [ase.data.chemical_symbols[number] for number in sorted(numbers)]
    return ''.join(only_atoms)


Domain(
    'DFT', DFTCalcWithMetadata,
    quantities=dict(
        formula=DomainQuantity(
            'The chemical (hill) formula of the simulated system.',
            order_default=True),
        atoms=DomainQuantity(
            'The atom labels of all atoms in the simulated system.',
            aggregations=len(ase.data.chemical_symbols), multi=True, zero_aggs=False),
        only_atoms=DomainQuantity(
            'The atom labels concatenated in species-number order. Used with keyword search '
            'to facilitate exclusive searches.',
            elastic_value=only_atoms, metadata_field='atoms', multi=True),
        basis_set=DomainQuantity(
            'The used basis set functions.', aggregations=20),
        xc_functional=DomainQuantity(
            'The xc functional type used for the simulation.', aggregations=20),
        system=DomainQuantity(
            'The system type of the simulated system.', aggregations=10),
        crystal_system=DomainQuantity(
            'The crystal system type of the simulated system.', aggregations=10),
        code_name=DomainQuantity(
            'The code name.', aggregations=40),
        spacegroup=DomainQuantity('The spacegroup of the simulated system as number'),
        spacegroup_symbol=DomainQuantity('The spacegroup as international short symbol'),
        geometries=DomainQuantity(
            'Hashes that describe unique geometries simulated by this code run.', multi=True),
        quantities=DomainQuantity(
            'All quantities that are used by this calculation',
            metric=('quantities', 'value_count'), multi=True),
        n_total_energies=DomainQuantity(
            'Number of total energy calculations',
            elastic_mapping=Integer()),
        n_calculations=DomainQuantity(
            'Number of single configuration calculation sections',
            elastic_mapping=Integer()),
        n_quantities=DomainQuantity(
            'Number of overall parsed quantities',
            elastic_mapping=Integer()),
        n_geometries=DomainQuantity(
            'Number of unique geometries',
            elastic_mapping=Integer()),
        n_atoms=DomainQuantity(
            'Number of atoms in the simulated system',
            elastic_mapping=Integer()),
        optimade=DomainQuantity(
            'Data for the optimade API',
            elastic_mapping=Object(optimade.ESOptimadeEntry),
            elastic_value=lambda entry: optimade.elastic_obj(entry, optimade.ESOptimadeEntry)
        )),
    metrics=dict(
        total_energies=('n_total_energies', 'sum'),
        calculations=('n_calculations', 'sum'),
        quantities=('n_quantities', 'sum'),
        geometries=('n_geometries', 'sum'),
        unique_geometries=('geometries', 'cardinality')
    ),
    default_statistics=[
        'atoms', 'basis_set', 'xc_functional', 'system', 'crystal_system', 'code_name'])
