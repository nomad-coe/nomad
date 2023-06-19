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

from tests.normalizing.conftest import run_normalize
from tests.normalizing.conftest import run_processing


def test_processing(raw_files, no_warn):
    directory = 'tests/data/datamodel/metainfo/eln/structure_file'
    mainfile = 'structure_file.archive.json'
    mainfile_schema = 'eln_with_structure.schema.archive.yaml'

    test_archive_schema = run_processing(directory, mainfile_schema)
    test_archive = run_processing(directory, mainfile)
    # assert archive for schema and solar cell entry
    assert len(test_archive_schema.definitions.section_definitions) == 1
    assert len(test_archive.results.material.chemical_formula_reduced) > 0
