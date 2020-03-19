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
DFT specific metadata
'''

import re

from nomad import utils, config
from nomad.metainfo import MSection, Section, Quantity, MEnum, SubSection
from nomad.metainfo.search_extension import Search

from .common import get_optional_backend_value
from .optimade import OptimadeEntry
from .metainfo.public import section_run


xc_treatments = {
    'gga': 'GGA',
    'hf_': 'HF',
    'oep': 'OEP',
    'hyb': 'hybrid',
    'mgg': 'meta-GGA',
    'vdw': 'vdW',
    'lda': 'LDA',
}
''' https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/XC-functional '''

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


class Label(MSection):
    '''
    Label that further classify a structure.

    Attributes:
        label: The label as a string
        type: The type of the label
        source: The source that this label was taken from.

    '''
    label = Quantity(type=str, a_search=Search())

    type = Quantity(type=MEnum(
        'compound_class', 'classification', 'prototype', 'prototype_id'),
        a_search=Search())

    source = Quantity(
        type=MEnum('springer', 'aflow_prototype_library'),
        a_search=Search())


class DFTMetadata(MSection):
    m_def = Section(a_domain='dft')

    basis_set = Quantity(
        type=str, default='not processed',
        description='The used basis set functions.',
        a_search=Search(statistic_size=20, default_statistic=True))

    xc_functional = Quantity(
        type=str, default='not processed',
        description='The libXC based xc functional classification used in the simulation.',
        a_search=Search(statistic_size=20, default_statistic=True))

    system = Quantity(
        type=str, default='not processed',
        description='The system type of the simulated system.',
        a_search=Search(default_statistic=True))

    crystal_system = Quantity(
        type=str, default='not processed',
        description='The crystal system type of the simulated system.',
        a_search=Search(default_statistic=True))

    spacegroup = Quantity(
        type=int, default='not processed',
        description='The spacegroup of the simulated system as number.',
        a_search=Search())

    spacegroup_symbol = Quantity(
        type=str, default='not processed',
        description='The spacegroup as international short symbol.',
        a_search=Search())

    code_name = Quantity(
        type=str, default='not processed',
        description='The name of the used code.',
        a_search=Search(statistic_size=40, default_statistic=True))

    code_version = Quantity(
        type=str, default='not processed',
        description='The version of the used code.',
        a_search=Search())

    n_geometries = Quantity(
        type=int, description='Number of unique geometries.',
        a_search=Search(metric_name='geometries', metric='sum'))

    n_calculations = Quantity(
        type=int,
        description='Number of single configuration calculation sections',
        a_search=Search(metric_name='calculations', metric='sum'))

    n_total_energies = Quantity(
        type=int, description='Number of total energy calculations',
        a_search=Search(metric_name='total_energies', metric='sum'))

    n_quantities = Quantity(
        type=int, description='Number of metainfo quantities parsed from the entry.',
        a_search=Search(metric='sum', metric_name='quantities'))

    quantities = Quantity(
        type=str, shape=['0..*'],
        description='All quantities that are used by this entry.',
        a_search=Search(
            metric_name='distinct_quantities', metric='cardinality', many_and='append'))

    geometries = Quantity(
        type=str, shape=['0..*'],
        description='Hashes for each simulated geometry',
        a_search=Search(metric_name='unique_geometries', metric='cardinality'))

    group_hash = Quantity(
        type=str,
        description='Hashes that describe unique geometries simulated by this code run.',
        a_search=Search(many_or='append', group='groups_grouped', metric_name='groups', metric='cardinality'))

    labels = SubSection(
        sub_section=Label, repeats=True,
        description='The labels taken from AFLOW prototypes and springer.',
        a_search='labels')

    optimade = SubSection(
        sub_section=OptimadeEntry,
        description='Metadata used for the optimade API.',
        a_search='optimade')

    def apply_domain_metadata(self, backend):
        from nomad.normalizing.system import normalized_atom_labels
        entry = self.m_parent

        logger = utils.get_logger(__name__).bind(
            upload_id=entry.upload_id, calc_id=entry.calc_id, mainfile=entry.mainfile)

        # code and code specific ids
        self.code_name = backend.get_value('program_name', 0)
        try:
            self.code_version = simplify_version(backend.get_value('program_version', 0))
        except KeyError:
            self.code_version = config.services.unavailable_value

        raw_id = get_optional_backend_value(backend, 'raw_id', 'section_run', None)
        if raw_id is not None:
            entry.raw_id = raw_id

        # metadata (system, method, chemistry)
        atoms = get_optional_backend_value(backend, 'atom_labels', 'section_system', [], logger=logger)
        if hasattr(atoms, 'tolist'):
            atoms = atoms.tolist()
        entry.n_atoms = len(atoms)
        atoms = list(set(normalized_atom_labels(set(atoms))))
        atoms.sort()
        entry.atoms = atoms

        self.crystal_system = get_optional_backend_value(
            backend, 'crystal_system', 'section_symmetry', logger=logger)
        self.spacegroup = get_optional_backend_value(
            backend, 'space_group_number', 'section_symmetry', 0, logger=logger)
        self.spacegroup_symbol = get_optional_backend_value(
            backend, 'international_short_symbol', 'section_symmetry', logger=logger)
        self.basis_set = map_basis_set_to_basis_set_label(
            get_optional_backend_value(backend, 'program_basis_set_type', 'section_run', logger=logger))
        self.system = get_optional_backend_value(
            backend, 'system_type', 'section_system', logger=logger)
        entry.formula = get_optional_backend_value(
            backend, 'chemical_composition_bulk_reduced', 'section_system', logger=logger)
        self.xc_functional = map_functional_name_to_xc_treatment(
            get_optional_backend_value(backend, 'XC_functional_name', 'section_method', logger=logger))

        # grouping
        self.group_hash = utils.hash(
            entry.formula,
            self.spacegroup,
            self.basis_set,
            self.xc_functional,
            self.code_name,
            self.code_version,
            entry.with_embargo,
            entry.uploader)

        # metrics and quantities
        quantities = set()
        geometries = set()
        n_quantities = 0
        n_calculations = 0
        n_total_energies = 0
        n_geometries = 0

        for root_section in backend.resource.contents:
            if not root_section.m_follows(section_run.m_def):
                continue

            quantities.add(root_section.m_def.name)
            n_quantities += 1

            for section, property_def, _ in root_section.m_traverse():
                property_name = property_def.name
                quantities.add(property_name)
                n_quantities += 1

                if property_name == 'energy_total':
                    n_total_energies += 1

                if property_name == 'configuration_raw_gid':
                    geometries.add(section.m_get(property_def))

                if property_name == 'section_single_configuration_calculation':
                    n_calculations += 1

                if property_name == 'section_system':
                    n_geometries += 1

        self.quantities = list(quantities)
        self.geometries = list(geometries)
        self.n_quantities = n_quantities
        self.n_calculations = n_calculations
        self.n_total_energies = n_total_energies
        self.n_geometries = n_geometries

        # labels
        compounds = set()
        classifications = set()
        for index in backend.get_sections('section_springer_material'):
            compounds.update(backend.get_value('springer_compound_class', index))
            classifications.update(backend.get_value('springer_classification', index))

        for compound in compounds:
            self.labels.append(Label(label=compound, type='compound_class', source='springer'))
        for classification in classifications:
            self.labels.append(Label(label=classification, type='classification', source='springer'))

        aflow_id = get_optional_backend_value(backend, 'prototype_aflow_id', 'section_prototype')
        aflow_label = get_optional_backend_value(backend, 'prototype_label', 'section_prototype')

        if aflow_id is not None and aflow_label is not None:
            self.labels.append(Label(label=aflow_label, type='prototype', source='aflow_prototype_library'))
            self.labels.append(Label(label=aflow_id, type='prototype_id', source='aflow_prototype_library'))

        # optimade
        optimade = backend.get_mi2_section(OptimadeEntry.m_def)
        if optimade is not None:
            self.optimade = optimade.m_copy()
