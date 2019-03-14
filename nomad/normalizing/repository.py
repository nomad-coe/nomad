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

    def map_functional_name_to_xc_treatment(self, name):
        if name == unavailable_label:
            return name

        return RepositoryNormalizer.xc_treatments.get(name[:3].lower(), name)

    def map_basis_set_to_basis_set_label(self, name):
        key = name.replace('_', '').replace('-', '').replace(' ', '').lower()
        return RepositoryNormalizer.basis_sets.get(key, name)

    def simplify_version(self, version):
        match = RepositoryNormalizer.version_re.search(version)
        if match is None:
            return version
        else:
            return match.group(0)

    def get_optional_value(self, key, section, unavailable_value=None):
        # Section is section_system, section_symmetry, etc...
        val = None  # Initialize to None, so we can compare section values.
        diff_flag = False  # Flag to store whether vals differ between sections.
        # Loop over the sections with the name section in the backend.
        for section_index in self._backend.get_sections(section):
            try:
                new_val = self._backend.get_value(key, section_index)
            except KeyError:
                continue

            # Compare values from iterations.
            diff_bool = new_val != val

            if type(diff_bool) is bool:
                if diff_bool and val is not None:
                    diff_flag = True
            elif diff_bool.all() and (val is not None):
                # Then we have an array, and diff bool has multiple values since
                # each item in array has been compared item for item.
                diff_flag = True

            val = new_val

        if diff_flag is True:
            self.logger.warning(
                'The values for %s differ between different %s' % (key, section))

        if val is None:
            self.logger.warning(
                'The values for %s where not available in any %s' % (key, section))
            return unavailable_value if unavailable_value is not None else unavailable_label
        else:
            return val

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
        # We add these values as optional as some parser developers may create parsers
        # where the first section contains only properties that are constant over
        # several simulations and use the other sections to define simulation specific
        # parameters. Ex. in first section_system CASTEP parsers defines # of elections
        # and in subsequent sections it defines atom labels, positions, etc...
        atom_labels = self.get_optional_value('atom_labels', 'section_system')
        b.addValue('repository_atomic_elements', list(set(atom_labels)))
        b.addValue('repository_atomic_elements_count', len(atom_labels))
        b.addValue(
            'repository_crystal_system',
            self.get_optional_value('crystal_system', 'section_symmetry'))
        b.addValue(
            'repository_spacegroup_nr',
            self.get_optional_value('space_group_number', 'section_symmetry', 0))
        b.addValue(
            'repository_spacegroup_symbol',
            self.get_optional_value('international_short_symbol', 'section_symmetry', 0))
        b.addValue(
            'repository_basis_set_type',
            self.map_basis_set_to_basis_set_label(
                self.get_optional_value('program_basis_set_type', 'section_run')))
        b.addValue(
            'repository_system_type',
            self.get_optional_value('system_type', 'section_system'))
        b.addValue(
            'repository_chemical_formula',
            self.get_optional_value('chemical_composition_bulk_reduced', 'section_system'))
        b.addValue(
            'repository_xc_treatment',
            self.map_functional_name_to_xc_treatment(
                self.get_optional_value('XC_functional_name', 'section_method')))

        b.closeNonOverlappingSection('section_repository_parserdata')
        if repository_info_context is None:
            b.closeNonOverlappingSection('section_repository_info')
        else:
            b.closeContext(repository_info_context)
        b.finishedParsingSession("ParseSuccess", None)
