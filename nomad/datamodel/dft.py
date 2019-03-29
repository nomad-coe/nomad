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
from elasticsearch_dsl import Integer
import ase.data

from nomad import utils, config

from .base import CalcWithMetadata, DomainQuantity, Domain


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

        self.n_total_energies = 0
        self.n_geometries = 0
        self.quantities = []
        self.geometries = []
        self.group_hash: str = None

        super().__init__(**kwargs)

    def apply_domain_metadata(self, backend):
        logger = utils.get_logger(__name__).bind(
            upload_id=self.upload_id, calc_id=self.calc_id, mainfile=self.mainfile)

        def get_optional_value(key, section, unavailable_value=None):
            # Section is section_system, section_symmetry, etc...
            val = None  # Initialize to None, so we can compare section values.
            # Loop over the sections with the name section in the backend.
            for section_index in backend.get_sections(section):
                try:
                    new_val = backend.get_value(key, section_index)
                except KeyError:
                    new_val = None

                # Compare values from iterations.
                if val is not None and new_val is not None:
                    if val.__repr__() != new_val.__repr__():
                        logger.warning(
                            'The values for %s differ between different %s: %s vs %s' %
                            (key, section, str(val), str(new_val)))

                val = new_val if new_val is not None else val

            if val is None:
                logger.warning(
                    'The values for %s where not available in any %s' % (key, section))
                return unavailable_value if unavailable_value is not None else config.services.unavailable_value
            else:
                return val

        if self.calc_id is None:
            self.calc_id = backend.get_value('calc_id')
        if self.upload_id is None:
            self.upload_id = backend.get_value('upload_id')
        if self.mainfile is None:
            self.mainfile = backend.get_value('main_file')

        self.code_name = backend.get_value('program_name', 0)
        self.code_version = simplify_version(backend.get_value('program_version', 0))

        self.atoms = get_optional_value('atom_labels', 'section_system')
        if hasattr(self.atoms, 'tolist'):
            self.atoms = self.atoms.tolist()
        self.n_atoms = len(self.atoms)
        self.atoms = list(set(self.atoms))
        self.atoms.sort()

        self.crystal_system = get_optional_value('crystal_system', 'section_symmetry')
        self.spacegroup = get_optional_value('space_group_number', 'section_symmetry', 0)
        self.spacegroup_symbol = get_optional_value('international_short_symbol', 'section_symmetry', 0)
        self.basis_set = map_basis_set_to_basis_set_label(
            get_optional_value('program_basis_set_type', 'section_run'))
        self.system = get_optional_value('system_type', 'section_system')
        self.formula = get_optional_value('chemical_composition_bulk_reduced', 'section_system')
        self.xc_functional = map_functional_name_to_xc_treatment(
            get_optional_value('XC_functional_name', 'section_method'))

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
        n_total_energies = 0
        n_geometries = 0

        for meta_info, _, value in backend._delegate.results.traverse():
            quantities.add(meta_info)
            if meta_info == 'energy_total':
                n_total_energies += 1
            if meta_info == 'section_system':
                n_geometries += 1
            if meta_info == 'configuration_raw_gid':
                geometries.add(value)

        self.quantities = list(quantities)
        self.geometries = list(geometries)
        self.n_total_energies = n_total_energies
        self.n_geometries = n_geometries


Domain.register_domain(DFTCalcWithMetadata, 'DFT', quantities=dict(
    formula=DomainQuantity(
        'The chemical (hill) formula of the simulated system.',
        order_default=True),
    atoms=DomainQuantity(
        'The atom labels of all atoms in the simulated system.',
        aggregations=len(ase.data.chemical_symbols)),
    basis_set=DomainQuantity(
        'The used basis set functions.', aggregations=10),
    xc_functional=DomainQuantity(
        'The xc functional type used for the simulation.', aggregations=10),
    system=DomainQuantity(
        'The system type of the simulated system.', aggregations=10),
    crystal_system=DomainQuantity(
        'The crystal system type of the simulated system.', aggregations=10),
    code_name=DomainQuantity(
        'The code name.', aggregations=10),
    spacegroup=DomainQuantity('The spacegroup of the simulated system as number'),
    spacegroup_symbol=DomainQuantity('The spacegroup as international short symbol'),
    n_total_energies=DomainQuantity(
        'Number of total energy calculations',
        metric=('total_energies', 'sum'),
        elastic_mapping=Integer()),
    n_geometries=DomainQuantity(
        'Number of unique geometries',
        metric=('geometries', 'cardinality'),
        elastic_mapping=Integer()),
    n_atoms=DomainQuantity('Number of atoms in the simulated system', elastic_mapping=Integer())))
