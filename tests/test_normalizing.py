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
import numpy as np
from ase import Atoms
import ase.build

from nomad import datamodel, config
from nomad.parsing import LocalBackend
from nomad.normalizing import normalizers

from tests.test_parsing import parsed_vasp_example  # pylint: disable=unused-import
from tests.test_parsing import parsed_template_example  # pylint: disable=unused-import
from tests.test_parsing import parsed_example  # pylint: disable=unused-import
from tests.test_parsing import parsed_template_no_system  # pylint: disable=unused-import
from tests.test_parsing import parse_file
from tests.utils import assert_log


boolean_positions = (
    'parsers/template', 'tests/data/normalizers/no_sim_cell_boolean_positions.json')

single_string_atom_labels = (
    'parsers/template', 'tests/data/normalizers/single_string_atom_labels.json')

unknown_atom_label = (
    'parsers/template', 'tests/data/normalizers/unknown_atom_label_test.json')

fcc_symmetry = (
    'parsers/template', 'tests/data/normalizers/fcc_crystal_structure.json')

vasp_parser = (
    'parsers/vasp', 'tests/data/parsers/vasp/vasp.xml')

vasp_parser_dos = (
    'parsers/vasp', 'tests/data/parsers/vasp/vasp_dos.xml')

glucose_atom_labels = (
    'parsers/template', 'tests/data/normalizers/glucose_atom_labels.json')

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
    'parsers/gulp': ['formula', 'xc_functional', 'system', 'basis_set'] + symmetry_keys,
    'parsers/turbomole': symmetry_keys,
    'parsers/elastic': ['basis_set', 'xc_functional', 'system'] + symmetry_keys,
    'parsers/dmol': ['system'] + symmetry_keys,
    'parser/molcas': symmetry_keys,
    'parsers/band': ['system'] + symmetry_keys,
    'parsers/qbox': ['xc_functional'],
    'parser/onetep': ['formula', 'basis_set', 'xc_functional', 'system'] + symmetry_keys
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
    if parsed_example.domain != 'ems':
        return run_normalize(parsed_example)

    return None


@pytest.fixture
def normalized_template_example(parsed_template_example) -> LocalBackend:
    return run_normalize(parsed_template_example)


def run_normalize_for_structure(atoms: Atoms) -> LocalBackend:
    template = parsed_template_no_system()

    # Fill structural information
    gid = template.openSection("section_system")
    template.addArrayValues("atom_positions", atoms.get_positions() * 1E-10)
    template.addArrayValues("atom_labels", atoms.get_chemical_symbols())
    template.addArrayValues("simulation_cell", atoms.get_cell() * 1E-10)
    template.addArrayValues("configuration_periodic_dimensions", atoms.get_pbc())
    template.closeSection("section_system", gid)

    return run_normalize(template)


@pytest.fixture(scope='session')
def single_point(two_d) -> LocalBackend:
    return two_d


@pytest.fixture(scope='session')
def geometry_optimization() -> LocalBackend:
    parser_name = "parsers/template"
    filepath = "tests/data/normalizers/fcc_crystal_structure.json"
    backend = parse_file((parser_name, filepath))
    backend = run_normalize(backend)
    return backend


@pytest.fixture(scope='session')
def molecular_dynamics(bulk) -> LocalBackend:
    return bulk


@pytest.fixture(scope='session')
def phonon() -> LocalBackend:
    parser_name = "parsers/phonopy"
    filepath = "tests/data/parsers/phonopy/phonopy-FHI-aims-displacement-01/control.in"
    backend = parse_file((parser_name, filepath))
    backend = run_normalize(backend)
    return backend


@pytest.fixture(scope='session')
def bulk() -> LocalBackend:
    parser_name = "parsers/cp2k"
    filepath = "tests/data/normalizers/cp2k_bulk_md/si_md.out"
    backend = parse_file((parser_name, filepath))
    backend = run_normalize(backend)
    return backend


@pytest.fixture(scope='session')
def two_d() -> LocalBackend:
    parser_name = "parsers/fhi-aims"
    filepath = "tests/data/normalizers/fhiaims_2d_singlepoint/aims.out"
    backend = parse_file((parser_name, filepath))
    backend = run_normalize(backend)
    return backend


@pytest.fixture(scope='session')
def surface() -> LocalBackend:
    parser_name = "parsers/fhi-aims"
    filepath = "tests/data/normalizers/fhiaims_surface_singlepoint/PBE-light+tight-rho2.out"
    backend = parse_file((parser_name, filepath))
    backend = run_normalize(backend)
    return backend


@pytest.fixture(scope='session')
def molecule() -> LocalBackend:
    parser_name = "parsers/fhi-aims"
    filepath = "tests/data/normalizers/fhiaims_molecule_singlepoint/aims.out"
    backend = parse_file((parser_name, filepath))
    backend = run_normalize(backend)
    return backend


@pytest.fixture(scope='session')
def atom() -> LocalBackend:
    parser_name = "parsers/gaussian"
    filepath = "tests/data/normalizers/gaussian_atom_singlepoint/m9b7.out"
    backend = parse_file((parser_name, filepath))
    backend = run_normalize(backend)
    return backend


@pytest.fixture(scope='session')
def one_d() -> LocalBackend:
    parser_name = "parsers/exciting"
    filepath = "tests/data/normalizers/exciting_1d_singlepoint/INFO.OUT"
    backend = parse_file((parser_name, filepath))
    backend = run_normalize(backend)
    return backend


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
    assert len(metadata.atoms) is not None
    assert metadata.spacegroup is not None

    exceptions = parser_exceptions.get(backend.get_value('parser_name'), [])

    if metadata.formula != config.services.unavailable_value:
        assert len(metadata.atoms) > 0

    for key in calc_metadata_keys:
        if key not in exceptions:
            assert getattr(metadata, key) != config.services.unavailable_value


def test_normalizer(normalized_example: LocalBackend):
    if normalized_example is not None:  # non dft cacls are not normalized
        assert_normalized(normalized_example)


def test_normalizer_faulty_matid(caplog):
    """ Runs normalizer on an example w/ bools for atom pos. Should force matid error."""
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
    """ Runs normalizer on ['Br','Si','Si','Za'], for normalizeation Za will be replaced,
        but stays int the labels.
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


def test_system_classification(atom, molecule, one_d, two_d, surface, bulk):
    # Atom
    assert atom['system_type'] == "atom"
    # Molecule / cluster
    assert molecule['system_type'] == "molecule / cluster"
    # 1D
    assert one_d['system_type'] == "1D"
    # 2D
    assert two_d['system_type'] == "2D"
    # Surface
    assert surface['system_type'] == "surface"
    # Bulk
    assert bulk['system_type'] == "bulk"


def test_representative_systems(single_point, molecular_dynamics, geometry_optimization, phonon):
    """Checks that the representative systems are correctly identified and
    processed by SystemNormalizer.
    """
    def check_representative_frames(backend):
        # For systems with multiple frames the first and two last should be processed.
        try:
            frames = backend["frame_sequence_local_frames_ref"]
        except KeyError:
            sccs = backend["section_single_configuration_calculation"]
            scc = sccs[-1]
            repr_system_idx = scc["single_configuration_calculation_to_system_ref"]
        else:
            sampling_method = backend["sampling_method"]
            if sampling_method == "molecular_dynamics":
                idx = 0
            else:
                idx = -1
            scc_idx = frames[idx]
            scc = backend["section_single_configuration_calculation"][scc_idx]
            repr_system_idx = scc["single_configuration_calculation_to_system_ref"]

        # Check that only the representative system has been labels with
        # "is_representative"
        for i, system in enumerate(backend["section_system"]):
            if i == repr_system_idx:
                assert system["is_representative"] is True
            else:
                with pytest.raises(KeyError):
                    system["is_representative"]

    check_representative_frames(single_point)
    check_representative_frames(molecular_dynamics)
    check_representative_frames(geometry_optimization)
    check_representative_frames(phonon)


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


def test_aflow_prototypes():
    """Tests that some basis structures are matched with the correct AFLOW prototypes
    """
    # No prototype info for non-bulk structures
    backend = run_normalize_for_structure(ase.build.molecule("H2O"))
    assert len(backend["section_prototype"]) == 0

    # No prototype info for bulk structure without match
    rattled = ase.build.bulk("C", crystalstructure="diamond", a=3.57, cubic=True)
    rattled.rattle(stdev=2, seed=42)
    rattled.wrap()
    backend = run_normalize_for_structure(rattled)
    assert len(backend["section_prototype"]) == 0

    # Diamond
    diamond = ase.build.bulk("C", crystalstructure="diamond", a=3.57)
    backend = run_normalize_for_structure(diamond)
    prototype_aflow_id = backend.get_value("prototype_aflow_id")
    prototype_label = backend.get_value("prototype_label")
    assert prototype_aflow_id == "A_cF8_227_a"
    assert prototype_label == "227-C-cF8"

    # BCC
    bcc = ase.build.bulk("Fe", crystalstructure="bcc", a=2.856)
    backend = run_normalize_for_structure(bcc)
    prototype_aflow_id = backend.get_value("prototype_aflow_id")
    prototype_label = backend.get_value("prototype_label")
    assert prototype_aflow_id == "A_cI2_229_a"
    assert prototype_label == "229-W-cI2"

    # FCC
    fcc = ase.build.bulk("Ge", crystalstructure="fcc", a=5.658)
    backend = run_normalize_for_structure(fcc)
    prototype_aflow_id = backend.get_value("prototype_aflow_id")
    prototype_label = backend.get_value("prototype_label")
    assert prototype_aflow_id == "A_cF4_225_a"
    assert prototype_label == "225-Cu-cF4"

    # Rocksalt
    rocksalt = ase.build.bulk("NaCl", crystalstructure="rocksalt", a=5.64)
    backend = run_normalize_for_structure(rocksalt)
    prototype_aflow_id = backend.get_value("prototype_aflow_id")
    prototype_label = backend.get_value("prototype_label")
    assert prototype_aflow_id == "AB_cF8_225_a_b"
    assert prototype_label == "225-ClNa-cF8"

    # Zincblende
    zincblende = ase.build.bulk("ZnS", crystalstructure="zincblende", a=5.42, cubic=True)
    backend = run_normalize_for_structure(zincblende)
    prototype_aflow_id = backend.get_value("prototype_aflow_id")
    prototype_label = backend.get_value("prototype_label")
    assert prototype_aflow_id == "AB_cF8_216_c_a"
    assert prototype_label == "216-SZn-cF8"

    # Wurtzite
    wurtzite = ase.build.bulk("SiC", crystalstructure="wurtzite", a=3.086, c=10.053)
    backend = run_normalize_for_structure(wurtzite)
    prototype_aflow_id = backend.get_value("prototype_aflow_id")
    prototype_label = backend.get_value("prototype_label")
    assert prototype_aflow_id == "AB_hP4_186_b_b"
    assert prototype_label == "186-SZn-hP4"


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


def test_dos_normalizer():
    """
    Ensure the DOS normalizer acted on the DOS values. We take a VASP example.
    """
    backend = parse_file(vasp_parser_dos)
    backend = run_normalize(backend)

    # Check if 'dos_values' were indeed normalized
    # 'dvn' stands for 'dos_values_normalized'
    backend_dvn = backend.get_value('dos_values_normalized', 0)
    last_value = backend_dvn[0, -1]
    expected = 1.7362195274239454e+47
    # Compare floats properly with numpy (delta tolerance involved)
    assert np.allclose(last_value, expected)
