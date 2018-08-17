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

from nomad.parsing import JSONStreamWriter, parser_dict
from nomad.normalizing import normalizers
import sys


def test_normalizer():
    vasp_parser = parser_dict['parsers/vasp']
    example_mainfile = '.dependencies/parsers/vasp/test/examples/xml/perovskite.xml'
    parser_backend = vasp_parser.run(example_mainfile)
    status, _ = parser_backend.status

    assert status == 'ParseSuccess'

    for normalizer_class in normalizers:
        normalizer = normalizer_class(parser_backend)
        normalizer.normalize()

    print(parser_backend)

    assert parser_backend.get_value('atom_species', 0) is not None
    assert parser_backend.get_value('system_type', 0) is not None
    assert parser_backend.get_value('crystal_system', 0) is not None
    assert parser_backend.get_value('space_group_number', 0) is not None
    assert parser_backend.get_value('XC_functional_name', 0) is not None
