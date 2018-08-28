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

import pytest
from typing import Generator
from datetime import datetime
import time

from nomad.parsing import LocalBackend
from nomad.search import Calc, key_mappings

from tests.test_normalizing import normalized_vasp_example  # pylint: disable=unused-import
from tests.test_parsing import parsed_vasp_example  # pylint: disable=unused-import


@pytest.fixture(scope='function')
def example_entry(normalized_vasp_example: LocalBackend) -> Generator[Calc, None, None]:
    entry = Calc.add_from_backend(
        normalized_vasp_example,
        upload_hash='test_upload_hash',
        calc_hash='test_calc_hash',
        upload_id='test_upload_id',
        mainfile='/test/mainfile',
        upload_time=datetime.now())
    time.sleep(1)  # eventually consistent?
    yield entry
    entry.delete()


def test_add(example_entry: Calc):
    result = Calc.get(id='%s/%s' % (example_entry.upload_hash, example_entry.calc_hash))

    assert result is not None
    for property in Calc._doc_type.mapping:
        property = key_mappings.get(property, property)
        assert getattr(result, property) is not None
