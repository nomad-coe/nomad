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

import pytest

from tests.normalizing.conftest import run_processing


@pytest.mark.parametrize('mainfile, assert_xpaths', [
    pytest.param('schema.archive.yaml', [], id='schema'),
    pytest.param('sample.archive.json', ['data.processes.pvd_evaporation.time'], id='sample'),
    pytest.param('PVD-P.archive.json', [], id='instrument'),
    pytest.param('Zinc_Selenide.archive.json', [], id='chemical')
])
def test_eln(mainfile, assert_xpaths, raw_files, no_warn):
    mainfile_directory = 'examples/data/eln'
    archive = run_processing(mainfile_directory, mainfile)

    for xpath in assert_xpaths:
        assert archive.m_xpath(xpath) is not None


@pytest.mark.parametrize('mainfile, assert_xpaths', [
    pytest.param('tabular-parser-col-mode.archive.yaml', ['data.My_Quantity'], id='col_mode'),
    pytest.param('tabular-parser-row-mode.archive.yaml', ['data.My_Subsection.My_Section[4].My_Quantity'],
                 id='row_mode'),
    pytest.param('tabular-parser-entry-mode.archive.yaml', [], id='entry_mode'),
])
def test_sample_tabular(mainfile, assert_xpaths, raw_files, no_warn):
    mainfile_directory = 'examples/data/docs'
    archive = run_processing(mainfile_directory, mainfile)

    for xpath in assert_xpaths:
        assert archive.m_xpath(xpath) is not None
