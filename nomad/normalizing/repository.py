# Copyright 2018 Fawzi Mohamed, Danio Brambila, Markus Scheidgen
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

import re

from nomad.parsing import BadContextURI

from .normalizer import Normalizer

unavailable_label = 'unavailable'


class RepositoryNormalizer(Normalizer):
    """
    The normalizer that turnes normalized parse results into a set of metadata
    quantities for the repository.
    """
    xc_treatments = {
        'gga': 'GGA',
        'hf_': 'HF',
        'oep': 'OEP',
        'hyb': 'hybrid',
        'mgga': 'meta-GGA',
        'vdw': 'vdW',
        'lda': 'LDA',
    }
    """ https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/XC-functional """

    version_re = re.compile(r'(\d+(\.\d+(\.\d+)?)?)')

    def map_functional_name_to_xc_treatment(self, name):
        if name == unavailable_label:
            return name

        return RepositoryNormalizer.xc_treatments.get(name[:3].lower(), name)

    def simplify_version(self, version):
        match = RepositoryNormalizer.version_re.search(version)
        if match is None:
            return version
        else:
            return match.group(0)

    def get_optional_value(self, key):
        try:
            return self._backend.get_value(key, 0)
        except KeyError:
            return unavailable_label

    def normalize(self, logger=None) -> None:
        super().normalize(logger)
        b = self._backend
        repository_info_context = '/section_repository_info/0'
        try:
            b.openContext(repository_info_context)
        except BadContextURI:
            b.openNonOverlappingSection('section_repository_info')
            repository_info_context = None
        b.openNonOverlappingSection('section_repository_parserdata')

        b.addValue('repository_checksum', b.get_value('calc_hash', 0))
        b.addValue('repository_program_name', b.get_value('program_name', 0))
        b.addValue(
            'repository_code_version',
            self.simplify_version(b.get_value('program_version', 0)))
        b.addValue('repository_parser_id', b.get_value('parser_name', 0))
        b.addValue('repository_chemical_formula', b.get_value('chemical_composition_bulk_reduced', 0))
        atom_labels = b.get_value('atom_labels', 0)
        b.addValue('repository_atomic_elements', list(set(atom_labels)))
        b.addValue('repository_atomic_elements_count', len(atom_labels))
        b.addValue('repository_system_type', b.get_value('system_type', 0))

        b.addValue('repository_crystal_system', self.get_optional_value('crystal_system'))
        b.addValue('repository_spacegroup_nr', self.get_optional_value('space_group_number'))

        b.addValue('repository_basis_set_type', self.get_optional_value('program_basis_set_type'))
        b.addValue(
            'repository_xc_treatment',
            self.map_functional_name_to_xc_treatment(self.get_optional_value(('XC_functional_name'))))

        b.closeNonOverlappingSection('section_repository_parserdata')
        if repository_info_context is None:
            b.closeNonOverlappingSection('section_repository_info')
        else:
            b.closeContext(repository_info_context)
        b.finishedParsingSession("ParseSuccess", None)
