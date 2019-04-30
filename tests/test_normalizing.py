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

from nomad import datamodel, config
from nomad.parsing import LocalBackend
from nomad.normalizing import normalizers

from tests.test_parsing import parsed_vasp_example  # pylint: disable=unused-import
from tests.test_parsing import parsed_template_example  # pylint: disable=unused-import
from tests.test_parsing import parsed_example  # pylint: disable=unused-import
from tests.test_parsing import parsed_faulty_unknown_matid_example  # pylint: disable=unused-import
from tests.test_parsing import parsed_single_string_atom_labels_test  # pylint: disable=unused-import
from tests.test_parsing import parsed_unknown_atom_label_test  # pylint: disable=unused-import
from tests.utils import assert_log


symmetry_keys = ['spacegroup', 'spacegroup_symbol', 'crystal_system']
calc_metadata_keys = [
    'code_name', 'code_version', 'basis_set', 'xc_functional', 'system', 'formula'] + symmetry_keys

parser_exceptions = {
    'parsers/wien2k': ['xc_functional'],
    'parsers/nwchem': symmetry_keys,
    'parsers/bigdft': symmetry_keys,
    'parsers/gaussian': symmetry_keys,
    'parsers/abinit': ['formula', 'system'] + symmetry_keys,
    'parsers/dl-poly': ['formula', 'basis_set', 'xc_functional', 'system'] + symmetry_keys,
    'parsers/lib-atoms': ['basis_set', 'xc_functional'],
    'parsers/orca': symmetry_keys,
    'parsers/octopus': symmetry_keys,
    'parsers/phonopy': ['basis_set', 'xc_functional'],
    'parsers/gpaw2': symmetry_keys,
    'parsers/gamess': ['formula', 'system'] + symmetry_keys,
    'parsers/gulp': ['formula', 'xc_functional', 'system'] + symmetry_keys,
    'parsers/turbomole': symmetry_keys,
    'parsers/elastic': ['basis_set', 'xc_functional', 'system'] + symmetry_keys
}
"""
Keys that the normalizer for certain parsers might not produce. In an ideal world this
map would be empty.
"""


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


def assert_normalized(backend: LocalBackend):
    metadata = datamodel.DFTCalcWithMetadata()
    metadata.apply_domain_metadata(backend)
    assert metadata.formula is not None
    assert metadata.code_name is not None
    assert metadata.code_version is not None
    assert metadata.basis_set is not None
    assert metadata.xc_functional is not None
    assert metadata.system is not None
    assert metadata.crystal_system is not None
    assert len(metadata.atoms) > 0
    assert metadata.spacegroup is not None

    exceptions = parser_exceptions.get(backend.get_value('parser_name'), [])

    for key in calc_metadata_keys:
        if key not in exceptions:
            assert getattr(metadata, key) != config.services.unavailable_value


def test_normalizer(normalized_example: LocalBackend):
    assert_normalized(normalized_example)


def test_normalizer_faulty_matid(
        parsed_faulty_unknown_matid_example: LocalBackend, caplog):
    """ Runs normalizer on an example w/ bools for atom pos. Should force matid error."""
    run_normalize(parsed_faulty_unknown_matid_example)

    assert_log(caplog, 'ERROR', 'matid project system classification failed')
    assert_log(caplog, 'ERROR', 'no lattice vectors but periodicity')


def test_normalizer_single_string_atom_labels(
        parsed_single_string_atom_labels_test: LocalBackend, caplog):
    """ Runs normalizer on ['Br1SiSiK'] expects error that it is formatted wrong."""
    run_normalize(parsed_single_string_atom_labels_test)
    assert_log(caplog, 'ERROR', 'Atom labels cannot be recognized.')


def test_normalizer_unknown_atom_label(
        parsed_unknown_atom_label_test: LocalBackend, caplog):
    """ Runs normalizer on ['Br','Si','Si','Za'], expects Za throws an error"""
    run_normalize(parsed_unknown_atom_label_test)
    assert_log(caplog, 'ERROR', 'Atom labels cannot be recognized.')
