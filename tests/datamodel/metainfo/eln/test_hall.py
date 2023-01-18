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

from numpy.testing import assert_array_almost_equal
from tests.normalizing.conftest import run_processing


def test_reading_of_measurement_file(raw_files, with_warn):
    directory = 'tests/data/datamodel/metainfo/eln/hall'
    mainfile = 'hall_python.data.archive.yaml'
    test_archive = run_processing(directory, mainfile)

    assert test_archive.metadata.entry_type == "Hall_experiment"
    assert len(test_archive.data.tasks) == 1
    assert len(test_archive.data.tasks[0].measurements) == 7
    assert_array_almost_equal(
        test_archive.data.tasks[0].measurements[0].data[0].field.magnitude,
        [3000., 2999.7, 3000.1, 3000.2, 3000.1]
    )


def test_reading_of_instrument_file(raw_files, with_warn):
    directory = 'tests/data/datamodel/metainfo/eln/hall'
    mainfile = 'hall_python.data.archive.yaml'
    test_archive = run_processing(directory, mainfile)

    assert test_archive.data.instrument
    assert test_archive.data.instrument.data_file
    assert isinstance(test_archive.data.instrument.instrument.use_instruments, bool)
    assert isinstance(test_archive.data.instrument.instrument.numberofsamples, int)
    assert isinstance(test_archive.data.instrument.instrument.ac_hall_type, str)
    assert len(test_archive.data.instrument.instrument.temperature_domain) == 8
