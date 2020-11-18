#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

class Context():
    """A simple class for holding the context related to an Encylopedia entry.
    """
    def __init__(
        self,
        material_type: str,
        method_type: str,
        calc_type: str,
        representative_system,
        representative_method,
        representative_scc,
        representative_scc_idx,
    ):
        self.material_type = material_type
        self.method_type = method_type
        self.calc_type = calc_type
        self.representative_system = representative_system
        self.representative_method = representative_method
        self.representative_scc = representative_scc
        self.representative_scc_idx = representative_scc_idx
        self.greatest_common_divisor: int = None
