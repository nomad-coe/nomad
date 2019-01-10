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

from nomad import utils

from .normalizer import Normalizer


class RepositoryNormalizer(Normalizer):
    """
    The normalizer that turnes normalized parse results into a set of metadata
    quantities for the repository.
    """
    def normalize(self, logger=None) -> None:
        super().normalize(logger)
        b = self._backend

        b.openNonOverlappingSection('section_repository_info')
        b.openNonOverlappingSection('section_repository_parserdata')

        b.addValue('repository_checksum', utils.archive.calc_hash(b.get_value('archive_id', 0)))
        b.addValue('repository_chemical_formula', b.get_value('chemical_composition_bulk_reduced', 0))
        b.addValue('repository_parser_id', b.get_value('parser_name', 0))
        atoms = b.get_value('atom_labels', 0)
        b.addValue('repository_atomic_elements', atoms)
        b.addValue('repository_atomic_elements_count', len(atoms))
        b.addValue('repository_basis_set_type', b.get_value('program_basis_set_type', 0))
        b.addValue('repository_crystal_system', b.get_value('crystal_system', 0))
        b.addValue('repository_program_name', b.get_value('program_name', 0))
        # TODO shorten and normalize the code version
        b.addValue('repository_code_version', b.get_value('program_version', 0))
        b.addValue('repository_spacegroup_nr', b.get_value('space_group_number', 0))
        b.addValue('repository_system_type', b.get_value('system_type', 0))
        # TODO shorten and normalize to functional type
        b.addValue('repository_xc_treatment', b.get_value('XC_functional_name', 0))

        b.closeNonOverlappingSection('section_repository_parserdata')
        b.closeNonOverlappingSection('section_repository_info')
        b.finishedParsingSession("ParseSuccess", None)
