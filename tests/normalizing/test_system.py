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

from nomad import datamodel, config
from nomad.parsing import LocalBackend

from tests.test_parsing import parse_file
from tests.normalizing.conftest import run_normalize, two_d   # pylint: disable=unused-import
from tests.utils import assert_log


boolean_positions = (
    'parsers/template', 'tests/data/normalizers/no_sim_cell_boolean_positions.json')

single_string_atom_labels = (
    'parsers/template', 'tests/data/normalizers/single_string_atom_labels.json')

unknown_atom_label = (
    'parsers/template', 'tests/data/normalizers/unknown_atom_label_test.json')

fcc_symmetry = (
    'parsers/template', 'tests/data/normalizers/fcc_crystal_structure.json')

two_d = (
    'parsers/template', 'tests/data/normalizers/fcc_crystal_structure.json')

vasp_parser = (
    'parsers/vasp', 'tests/data/parsers/vasp/vasp.xml')

glucose_atom_labels = (
    'parsers/template', 'tests/data/normalizers/glucose_atom_labels.json')

calc_metadata_keys = [
    'code_name', 'code_version', 'basis_set', 'xc_functional', 'system', 'formula']

parser_exceptions = {
    'parsers/wien2k': ['xc_functional'],
    'parsers/abinit': ['formula', 'system'],
    'parsers/dl-poly': ['formula', 'basis_set', 'xc_functional', 'system'],
    'parsers/lib-atoms': ['basis_set', 'xc_functional'],
    'parsers/phonopy': ['basis_set', 'xc_functional'],
    'parsers/gamess': ['formula', 'system'],
    'parsers/gulp': ['formula', 'xc_functional', 'system'],
    'parsers/elastic': ['basis_set', 'xc_functional', 'system'],
    'parsers/dmol': ['system'],
    'parsers/band': ['system'],
    'parsers/qbox': ['xc_functional'],
    'parser/onetep': ['formula', 'basis_set', 'xc_functional', 'system']
}
"""
Keys that the normalizer for certain parsers might not produce. In an ideal world this
map would be empty.
"""


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
    assert len(metadata.atoms) is not None

    exceptions = parser_exceptions.get(backend.get_value('parser_name'), [])

    if metadata.formula != config.services.unavailable_value:
        assert len(metadata.atoms) > 0

    for key in calc_metadata_keys:
        if key not in exceptions:
            assert getattr(metadata, key) != config.services.unavailable_value


def test_normalizer(normalized_example: LocalBackend):
    assert_normalized(normalized_example)


def test_normalizer_faulty_matid(caplog):
    """Runs normalizer on an example w/ bools for atom pos. Should force matid error."""
    # assert isinstance(backend, LocalBackend)
    backend = parse_file(boolean_positions)
    run_normalize(backend)
    assert_log(caplog, 'ERROR', 'matid project system classification failed')
    assert_log(caplog, 'ERROR', 'no lattice vectors but periodicity')


def test_normalizer_single_string_atom_labels(caplog):
    """
    Runs normalizer on ['Br1SiSiK'] expects error. Should replace the label with 'X' and
    the numbers of postitions should not match the labels.
    """
    backend = parse_file(single_string_atom_labels)
    run_normalize(backend)
    assert_log(caplog, 'ERROR', 'len of atom position does not match number of atoms')


def test_normalizer_unknown_atom_label(caplog, no_warn):
    """Runs normalizer on ['Br','Si','Si','Za'], for normalization Za will be
    replaced, but stays int the labels.
    """
    backend = parse_file(unknown_atom_label)
    run_normalize(backend)
    assert backend.get_value('atom_labels')[3] == 'Za'


def test_symmetry_classification_fcc():
    """Runs normalizer where lattice vectors should give fcc symmetry."""
    backend = parse_file(fcc_symmetry)
    backend = run_normalize(backend)
    expected_crystal_system = 'cubic'
    expected_bravais_lattice = 'cF'
    expected_point_group = 'm-3m'
    expected_origin_shift = [0, 0, 0]
    cyrstal_system = backend.get_value('crystal_system')
    assert cyrstal_system == expected_crystal_system
    bravais_lattice = backend.get_value('bravais_lattice')
    assert bravais_lattice == expected_bravais_lattice
    point_group = backend.get_value('point_group')
    assert point_group == expected_point_group
    origin_shift = backend.get_value('origin_shift')
    assert all(origin_shift == expected_origin_shift)


def test_system_classification(bulk, two_d):
    """Tests that the system classification is correct for different kind of systems
    """
    # Bulk system
    assert bulk.get_value('system_type') == "bulk"
    # 2D system
    assert two_d.get_value('system_type') == "2D"


def test_reduced_chemical_formula():
    "Ensure we get the right reduced chemical formula for glucose atom labels"
    backend = parse_file(glucose_atom_labels)
    backend = run_normalize(backend)
    expected_red_chem_formula = 'C6H12O6'
    reduced_chemical_formula = backend.get_value('chemical_composition_bulk_reduced')
    assert expected_red_chem_formula == reduced_chemical_formula


def test_vasp_incar_system():
    """
    Ensure we can test an incar value in the VASP example
    """
    backend = parse_file(vasp_parser)
    backend = run_normalize(backend)
    expected_value = 'SrTiO3'  # material's formula in vasp.xml

    # backend_value = backend.get_value('x_vasp_unknown_incars')  # OK
    # backend_value = backend.get_value('x_vasp_atom_kind_refs')  # OK
    backend_value = backend.get_value('x_vasp_incar_SYSTEM')  # OK

    print("backend_value: ", backend_value)
    assert expected_value == backend_value


def test_springer_normalizer():
    """
    Ensure the Springer normalizer works well with the VASP example.
    """
    backend = parse_file(vasp_parser)
    backend = run_normalize(backend)

    backend_value = backend.get_value('springer_id', 89)
    expected_value = 'sd_1932539'
    assert expected_value == backend_value

    backend_value = backend.get_value('springer_alphabetical_formula', 89)
    expected_value = 'O3SrTi'
    assert expected_value == backend_value

    backend_value = backend.get_value('springer_url', 89)
    expected_value = 'http://materials.springer.com/isp/crystallographic/docs/sd_1932539'
    assert expected_value == backend_value
