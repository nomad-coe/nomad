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
from systemtypenormalizer.classify_structure import ClassifyStructure


class SystemTypeNormalizer(SystemBasedNormalizer):

    def normalize_system(self, section_system) -> None:
        structure = ClassifyStructure(section_system)
        structure.classify()
        structure_type = structure.classification
        self._backend.addValue('system_type', structure_type)
