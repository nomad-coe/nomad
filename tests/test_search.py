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
from elasticsearch_dsl import Index
from datetime import datetime

from nomad.parsing import LocalBackend
from nomad.search import Calc

from tests.test_normalizing import normalized_vasp_example  # pylint: disable=unused-import
from tests.test_parsing import parsed_vasp_example  # pylint: disable=unused-import


@pytest.fixture(scope='function', autouse=True)
def index():
    yield
    Index('calcs').delete()


def test_add(normalized_vasp_example: LocalBackend):
    Calc.add_from_backend(
        normalized_vasp_example,
        upload_hash='test_upload_hash',
        calc_hash='test_calc_hash',
        mainfile='/test/mainfile',
        upload_time=datetime.now())

    result = Calc.get(id='%s/%s' % ('test_upload_hash', 'test_calc_hash'))

    assert result is not None
    for property in Calc._doc_type.mapping:
        assert getattr(result, property) is not None
