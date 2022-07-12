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

from tests.normalizing.conftest import run_processing


def test_processing(raw_files, no_warn):
    directory = 'tests/data/datamodel/metainfo/eln/perovskite_database'
    mainfile = 'example.archive.json'
    test_archive = run_processing(directory, mainfile)

    # assert archive
    assert test_archive.metadata.entry_type == 'PerovskiteSolarCell'
    assert test_archive.results.properties.optoelectronic.solar_cell.efficiency > -1
    assert len(test_archive.results.material.chemical_formula_reduced) > 0
    assert len(test_archive.data.jv.jv_curve[0].current_density) > -1
    # assert min(0.0, 10) < test_archive.data.eqe.bandgap_eqe < max(0.0, 10)
