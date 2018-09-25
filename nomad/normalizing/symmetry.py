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
from symmetrynormalizer import symmetry_analysis


class SymmetryNormalizer(SystemBasedNormalizer):
    """
    This is basically a copy of the legace NOMAD-coe symmetry normalizer.
    """
    def __init__(self, backend):
        super().__init__(backend, all_sections=True)

    def normalize_system(self, section_system) -> None:
        symmetry_analysis.logging = self.logger
        symmetry_analysis.normalize(self._backend, section_system)
