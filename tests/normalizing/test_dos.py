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

import numpy as np

from tests.test_parsing import parse_file
from tests.normalizing.conftest import run_normalize


vasp_parser_dos = (
    'parsers/vasp', 'tests/data/parsers/vasp/vasp_dos.xml')


def test_dos_normalizer():
    """
    Ensure the DOS normalizer acted on the DOS values. We take a VASP example.
    """
    backend = parse_file(vasp_parser_dos)
    backend = run_normalize(backend)

    # Check if 'dos_values' were indeed normalized
    # 'dvn' stands for 'dos_values_normalized'
    backend_dvn = backend.get_value('dos_values_normalized', 0)
    last_value = backend_dvn[0, -1]
    expected = 1.7362195274239454e+47
    # Compare floats properly with numpy (delta tolerance involved)
    assert np.allclose(last_value, expected)
