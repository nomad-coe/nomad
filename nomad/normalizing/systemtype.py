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

from nomad.normalizing.normalizer import SystemBasedNormalizer
from systemtypenormalizer import classify_structure


class SystemTypeNormalizer(SystemBasedNormalizer):
    def __init__(self, backend):
        super().__init__(backend, all_sections=True)

    def normalize_system(self, section_system) -> None:
        classify_structure.logger = self.logger

        structure = classify_structure.ClassifyStructure(section_system)
        structure.classify()
        structure_type = structure.classification
        print("Structure type form the system type normalizer is:")
        print(structure_type)
        self._backend.addValue('system_type', structure_type)
