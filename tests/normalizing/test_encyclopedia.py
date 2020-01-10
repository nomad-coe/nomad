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

from nomad.parsing import LocalBackend
from nomad.normalizing import normalizers

from tests.test_parsing import parsed_vasp_example  # pylint: disable=unused-import
from tests.test_parsing import parsed_template_example  # pylint: disable=unused-import
from tests.test_parsing import parsed_example  # pylint: disable=unused-import
from tests.test_parsing import parse_file
from tests.utils import assert_log


def run_normalize(backend: LocalBackend) -> LocalBackend:
    status, _ = backend.status

    assert status == 'ParseSuccess'

    for normalizer_class in normalizers:
        normalizer = normalizer_class(backend)
        normalizer.normalize()
    return backend


@pytest.fixture
def single_point(parsed_template_example) -> LocalBackend:
    return run_normalize(parsed_template_example)


def test_single_point(single_point: LocalBackend):
    """Tests that single point calculations are correctly processed."
    """
    pass
