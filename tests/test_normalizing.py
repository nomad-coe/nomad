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
from tests.test_parsing import parsed_faulty_unknown_matid_example  # pylint: disable=unused-import
from tests.utils import assert_log


def run_normalize(backend: LocalBackend) -> LocalBackend:
    status, _ = backend.status

    assert status == 'ParseSuccess'

    for normalizer_class in normalizers:
        normalizer = normalizer_class(backend)
        normalizer.normalize()

    return backend


@pytest.fixture
def normalized_vasp_example(parsed_vasp_example: LocalBackend) -> LocalBackend:
    return run_normalize(parsed_vasp_example)


@pytest.fixture
def normalized_example(parsed_example: LocalBackend) -> LocalBackend:
    return run_normalize(parsed_example)


@pytest.fixture
def normalized_template_example(parsed_template_example) -> LocalBackend:
    return run_normalize(parsed_template_example)


def test_template_example_normalizer(parsed_template_example, no_warn, caplog):
    run_normalize(parsed_template_example)
    print(str(caplog.records))


def assert_normalized(backend):
    with open("test_file_name.json", "wt") as file:
        backend.write_json(file)
    metadata = backend.metadata()['section_repository_info']['section_repository_parserdata']
    count = 0
    for metainfo in backend.metaInfoEnv().infoKindEls():
        if 'section_repository_parserdata' in metainfo.superNames:
            count += 1
            assert backend.get_value(metainfo.name, 0) is not None
            assert metadata.get(metainfo.name, None) is not None
    assert count > 0


def test_normalizer(normalized_example: LocalBackend, no_warn):
    assert_normalized(normalized_example)


def test_normalizer_faulty_matid(
        parsed_faulty_unknown_matid_example: LocalBackend, caplog):
    """ Runs normalizer on an example w/ bools for atom pos. Should force matid error."""
    run_normalize(parsed_faulty_unknown_matid_example)
    unknown_class_error = (
        'Matid classfication has given us an unexpected type')

    wrong_class_for_no_sim_cell = (
        'Matid classified more than 1D despite having no simulation_cell')

    assert_log(caplog, 'ERROR', unknown_class_error)
    assert_log(caplog, 'ERROR', wrong_class_for_no_sim_cell)
