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
    directory = 'tests/data/datamodel/metainfo/eln/solar_cells'
    mainfile = 'solar_cell.archive.json'
    mainfile_schema = 'solar_cell_eln.schema.archive.yaml'

    test_archive_schema = run_processing(directory, mainfile_schema)
    run_normalize(test_archive_schema)

    test_archive = run_processing(directory, mainfile)
    run_normalize(test_archive)
    # assert archive for schema and solar cell entry
    assert len(test_archive_schema.definitions.section_definitions) == 1
    assert test_archive_schema.metadata.entry_type == 'Schema'
    assert test_archive.metadata.entry_type == 'SolarCell'
    assert test_archive.results.properties.optoelectronic.solar_cell.efficiency >= 0
    assert len(test_archive.results.material.chemical_formula_reduced) > 0
    assert len(test_archive.results.properties.optoelectronic.solar_cell.device_stack) > 0
    assert len(test_archive.results.properties.optoelectronic.solar_cell.absorber) > 0
    assert test_archive.results.properties.optoelectronic.band_gap[0].value.magnitude >= 0
