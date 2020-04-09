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

import ase.build

from nomad import datamodel, config
from nomad.parsing import Backend
from nomad.app import dump_json

from tests.test_parsing import parsed_vasp_example  # pylint: disable=unused-import
from tests.test_parsing import parsed_template_example  # pylint: disable=unused-import
from tests.test_parsing import parsed_example  # pylint: disable=unused-import
from tests.test_parsing import parsed_template_no_system  # pylint: disable=unused-import
from tests.test_parsing import parse_file
from tests.normalizing.conftest import run_normalize, run_normalize_for_structure   # pylint: disable=unused-import
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
    'dft.code_name', 'dft.code_version', 'dft.basis_set', 'dft.xc_functional', 'dft.system', 'formula']

parser_exceptions = {
    'parsers/wien2k': ['dft.xc_functional'],
    'parsers/abinit': ['formula', 'dft.system'],
    'parsers/dl-poly': ['formula', 'dft.basis_set', 'dft.xc_functional', 'dft.system'],
    'parsers/lib-atoms': ['dft.basis_set', 'dft.xc_functional'],
    'parsers/phonopy': ['dft.basis_set', 'dft.xc_functional'],
    'parsers/gamess': ['formula', 'dft.system', 'dft.xc_functional'],
    'parsers/gulp': ['formula', 'dft.xc_functional', 'dft.system', 'dft.basis_set'],
    'parsers/elastic': ['dft.basis_set', 'dft.xc_functional', 'dft.system'],
    'parsers/dmol': ['dft.system'],
    'parsers/band': ['dft.system'],
    'parsers/qbox': ['dft.xc_functional'],
    'parser/onetep': ['formula', 'dft.basis_set', 'dft.xc_functional', 'dft.system']
}
'''
Keys that the normalizer for certain parsers might not produce. In an ideal world this
map would be empty.
'''


def test_template_example_normalizer(parsed_template_example, no_warn, caplog):
    run_normalize(parsed_template_example)


def assert_normalized(backend: Backend):

    metadata = datamodel.EntryMetadata(domain=backend.domain)
    metadata.apply_domain_metadata(backend)
    assert metadata.formula is not None
    assert len(metadata.atoms) is not None

    if backend.domain == 'dft':
        assert metadata.dft.code_name is not None
        assert metadata.dft.code_version is not None
        assert metadata.dft.basis_set is not None
        assert metadata.dft.xc_functional is not None
        assert metadata.dft.system is not None

    parser_name = backend.entry_archive.section_metadata.parser_name
    exceptions = parser_exceptions.get(parser_name, [])

    if metadata.formula != config.services.unavailable_value:
        assert len(metadata.atoms) > 0

    for key in calc_metadata_keys:
        if key in exceptions:
            continue

        if '.' in key and not key.startswith(backend.domain):
            continue

        assert metadata[key] != config.services.unavailable_value, '%s must not be unavailable' % key

    # check if the result can be dumped
    dump_json(backend.entry_archive.m_to_dict())


def test_normalizer(normalized_example: Backend):
    assert_normalized(normalized_example)


def test_normalizer_faulty_matid(caplog):
    '''Runs normalizer on an example w/ bools for atom pos. Should force matid error.
    '''
    backend = parse_file(boolean_positions)
    run_normalize(backend)
    assert_log(caplog, 'ERROR', 'matid project system classification failed')
    assert_log(caplog, 'ERROR', 'no lattice vectors but periodicity')


def test_normalizer_single_string_atom_labels(caplog):
    '''
    Runs normalizer on ['Br1SiSiK'] expects error. Should replace the label with 'X' and
    the numbers of postitions should not match the labels.
    '''
    backend = parse_file(single_string_atom_labels)
    run_normalize(backend)
    assert_log(caplog, 'ERROR', 'len of atom position does not match number of atoms')


def test_normalizer_unknown_atom_label(caplog, no_warn):
    '''Runs normalizer on ['Br','Si','Si','Za'], for normalization Za will be replaced,
    but stays in the labels.
    '''
    backend = parse_file(unknown_atom_label)
    run_normalize(backend)
    assert backend['atom_labels'][3] == 'Za'


def test_symmetry_classification_fcc():
    '''Runs normalizer where lattice vectors should give fcc symmetry.'''
    backend = parse_file(fcc_symmetry)
    backend = run_normalize(backend)
    expected_crystal_system = 'cubic'
    expected_bravais_lattice = 'cF'
    expected_point_group = 'm-3m'
    expected_origin_shift = [0, 0, 0]
    sec_symmetry = backend.entry_archive.section_run[0].section_system[0].section_symmetry[0]
    crystal_system = sec_symmetry.crystal_system
    assert crystal_system == expected_crystal_system
    bravais_lattice = sec_symmetry.bravais_lattice
    assert bravais_lattice == expected_bravais_lattice
    point_group = sec_symmetry.point_group
    assert point_group == expected_point_group
    origin_shift = sec_symmetry.origin_shift
    assert all(origin_shift == expected_origin_shift)


def test_system_classification(atom, molecule, one_d, two_d, surface, bulk):
    """Tests that the system classification is correct for different kind of systems
    """
    # Atom
    assert atom.get_value("system_type") == "atom"
    # Molecule
    assert molecule.get_value("system_type") == "molecule / cluster"
    # 1D system
    assert one_d.get_value("system_type") == "1D"
    # 2D system
    assert two_d.get_value("system_type") == "2D"
    # Surface
    assert surface.get_value("system_type") == "surface"
    # Bulk system
    assert bulk.get_value("system_type") == "bulk"


def test_representative_systems(single_point, molecular_dynamics, geometry_optimization, phonon):
    '''Checks that the representative systems are correctly identified and
    processed by SystemNormalizer.
    '''
    def check_representative_frames(backend):
        # For systems with multiple frames the first and two last should be processed.
        try:
            frames = backend.entry_archive.section_run[0].section_frame_sequence[0].frame_sequence_local_frames_ref
        except Exception:
            scc = backend.entry_archive.section_run[0].section_single_configuration_calculation[-1]
            repr_system = scc.single_configuration_calculation_to_system_ref
        else:
            sampling_method = backend.entry_archive.section_run[0].section_sampling_method[0].sampling_method
            if sampling_method == "molecular_dynamics":
                scc = frames[0]
            else:
                scc = frames[-1]
            repr_system = scc.single_configuration_calculation_to_system_ref

        # Check that only the representative system has been labels with
        # "is_representative"
        for system in backend.entry_archive.section_run[0].section_system:
            if system.m_parent_index == repr_system.m_parent_index:
                assert system.is_representative is True
            else:
                assert system.is_representative is None

    check_representative_frames(single_point)
    check_representative_frames(molecular_dynamics)
    check_representative_frames(geometry_optimization)
    check_representative_frames(phonon)


def test_reduced_chemical_formula():
    "Ensure we get the right reduced chemical formula for glucose atom labels"
    backend = parse_file(glucose_atom_labels)
    backend = run_normalize(backend)
    expected_red_chem_formula = 'C6H12O6'
    sec_system = backend.entry_archive.section_run[0].section_system[0]
    reduced_chemical_formula = sec_system.chemical_composition_bulk_reduced
    assert expected_red_chem_formula == reduced_chemical_formula


def test_vasp_incar_system():
    '''
    Ensure we can test an incar value in the VASP example
    '''
    backend = parse_file(vasp_parser)
    backend = run_normalize(backend)
    expected_value = 'SrTiO3'  # material's formula in vasp.xml

    # backend_value = backend.get_value('x_vasp_unknown_incars')  # OK
    # backend_value = backend.get_value('x_vasp_atom_kind_refs')  # OK
    backend_value = backend.get_value('x_vasp_incar_SYSTEM')  # OK

    assert expected_value == backend_value


def test_aflow_prototypes():
    '''Tests that some basis structures are matched with the correct AFLOW prototypes
    '''
    def get_proto(atoms):
        backend = run_normalize_for_structure(atoms)
        try:
            sec_proto = backend.entry_archive.section_run[0].section_system[0].section_prototype[0]
            prototype_aflow_id = sec_proto.prototype_aflow_id
            prototype_label = sec_proto.prototype_label
        except Exception:
            prototype_aflow_id = None
            prototype_label = None

        return prototype_aflow_id, prototype_label

    # No prototype info for non-bulk structures
    water = ase.build.molecule("H2O")
    prototype_aflow_id, prototype_label = get_proto(water)
    assert prototype_aflow_id is None
    assert prototype_label is None

    # No prototype info for bulk structure without match
    rattled = ase.build.bulk("C", crystalstructure="diamond", a=3.57, cubic=True)
    rattled.rattle(stdev=2, seed=42)
    rattled.wrap()
    aflow_id, prototype_label = get_proto(rattled)
    assert aflow_id is None
    assert prototype_label is None

    # Diamond
    diamond = ase.build.bulk("C", crystalstructure="diamond", a=3.57)
    prototype_aflow_id, prototype_label = get_proto(diamond)
    assert prototype_aflow_id == "A_cF8_227_a"
    assert prototype_label == "227-C-cF8"

    # BCC
    bcc = ase.build.bulk("Fe", crystalstructure="bcc", a=2.856)
    prototype_aflow_id, prototype_label = get_proto(bcc)
    assert prototype_aflow_id == "A_cI2_229_a"
    assert prototype_label == "229-W-cI2"

    # FCC
    fcc = ase.build.bulk("Ge", crystalstructure="fcc", a=5.658)
    prototype_aflow_id, prototype_label = get_proto(fcc)
    assert prototype_aflow_id == "A_cF4_225_a"
    assert prototype_label == "225-Cu-cF4"

    # Rocksalt
    rocksalt = ase.build.bulk("NaCl", crystalstructure="rocksalt", a=5.64)
    prototype_aflow_id, prototype_label = get_proto(rocksalt)
    assert prototype_aflow_id == "AB_cF8_225_a_b"
    assert prototype_label == "225-ClNa-cF8"

    # Zincblende
    zincblende = ase.build.bulk("ZnS", crystalstructure="zincblende", a=5.42, cubic=True)
    prototype_aflow_id, prototype_label = get_proto(zincblende)
    assert prototype_aflow_id == "AB_cF8_216_c_a"
    assert prototype_label == "216-SZn-cF8"

    # Wurtzite
    wurtzite = ase.build.bulk("SiC", crystalstructure="wurtzite", a=3.086, c=10.053)
    prototype_aflow_id, prototype_label = get_proto(wurtzite)
    assert prototype_aflow_id == "AB_hP4_186_b_b"
    assert prototype_label == "186-SZn-hP4"


def test_springer_normalizer():
    '''
    Ensure the Springer normalizer works well with the VASP example.
    '''
    backend = parse_file(vasp_parser)
    backend = run_normalize(backend)

    gindex = 0

    backend_value = backend.get_value('springer_id', gindex)
    expected_value = 'sd_0305232'
    assert expected_value == backend_value

    backend_value = backend.get_value('springer_alphabetical_formula', gindex)
    expected_value = 'O3SrTi'
    assert expected_value == backend_value

    backend_value = backend.get_value('springer_url', gindex)
    expected_value = 'http://materials.springer.com/isp/crystallographic/docs/sd_0305232'
    assert expected_value == backend_value
