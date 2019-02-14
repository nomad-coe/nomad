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

import os
from io import StringIO
import json
import pytest

from nomadcore.local_meta_info import loadJsonFile
import nomad_meta_info

from nomad import utils, files
from nomad.parsing import JSONStreamWriter, parser_dict, match_parser
from nomad.parsing import LocalBackend, BadContextURI

parser_examples = [
    ('parsers/random', 'test/data/parsers/random_0'),
    ('parsers/template', 'tests/data/parsers/template.json'),
    ('parsers/exciting', 'tests/data/parsers/exciting/Ag/INFO.OUT'),
    ('parsers/exciting', 'tests/data/parsers/exciting/GW/INFO.OUT'),
    ('parsers/vasp', 'tests/data/parsers/vasp/vasp.xml'),
    ('parsers/vasp', 'tests/data/parsers/vasp_compressed/vasp.xml.gz'),
    ('parsers/vaspoutcar', 'tests/data/parsers/vasp_outcar/OUTCAR'),
    ('parsers/fhi-aims', 'tests/data/parsers/fhi-aims/aims.out'),
    ('parsers/cp2k', 'tests/data/parsers/cp2k/si_bulk8.out'),
    ('parsers/crystal', 'tests/data/parsers/crystal/si.out'),
    ('parsers/cpmd', 'tests/data/parsers/cpmd/geo_output.out'),
    ('parsers/nwchem', 'tests/data/parsers/nwchem/single_point/output.out'),
    ('parsers/bigdft', 'tests/data/parsers/bigdft/n2_output.out'),
    ('parsers/wien2k', 'tests/data/parsers/wien2k/AlN/AlN_ZB.scf'),
    ('parsers/band', 'tests/data/parsers/band_adf.out'),
    ('parsers/gaussian', 'tests/data/parsers/gaussian/aniline.out'),
    ('parsers/quantumespresso', 'tests/data/parsers/quantum-espresso/benchmark.out')
]

faulty_unknown_one_d_matid_example = [
    ('parsers/template', 'tests/data/normalizers/no_sim_cell_boolean_positions.json')
]

correct_num_output_files = 19


class TestLocalBackend(object):

    @pytest.fixture(scope='session')
    def meta_info(self):
        file_dir = os.path.dirname(os.path.abspath(nomad_meta_info.__file__))
        path = os.path.join(file_dir, 'all.nomadmetainfo.json')
        meta_info, _ = loadJsonFile(path)
        return meta_info

    @pytest.fixture(scope='function')
    def backend(self, meta_info):
        return LocalBackend(meta_info, debug=True)

    def test_meta_info(self, meta_info, no_warn):
        assert 'section_topology' in meta_info

    def test_metadata(self, backend, no_warn):
        g_index = backend.openSection('section_calculation_info')
        assert g_index == 0
        backend.addValue('calc_id', 't0')
        backend.closeSection('section_calculation_info', 0)
        g_index = backend.openSection('section_repository_info')
        backend.addValue('repository_calc_id', 1)
        backend.closeSection('section_repository_info', 0)
        assert json.dumps(backend.metadata()) is not None

    def test_section(self, backend, no_warn):
        g_index = backend.openSection('section_run')
        assert g_index == 0
        backend.addValue('program_name', 't0')
        backend.closeSection('section_run', 0)

        g_index = backend.openSection('section_run')
        assert g_index == 1

        g_index = backend.openSection('section_run')
        assert g_index == 2

        backend.addValue('program_name', 't1', 1)
        backend.addValue('program_name', 't2', 2)

        backend.closeSection('section_run', 1)
        backend.closeSection('section_run', 2)

        assert backend.get_sections('section_run') == [0, 1, 2]
        for i in range(0, 3):
            assert backend.get_value('program_name', i) == 't%d' % i

    def test_two_sections(self, backend, no_warn):
        g_index = backend.openSection('section_run')
        assert g_index == 0
        backend.addValue('program_name', 't0')
        backend.closeSection('section_run', 0)

        g_index = backend.openSection('section_calculation_info')
        assert g_index == 0
        backend.addValue('parser_name', 'p0')
        backend.closeSection('section_calculation_info', 0)

        assert backend.get_sections('section_run') == [0]
        assert backend.get_sections('section_calculation_info') == [0]

        output = StringIO()
        backend.write_json(output)
        archive = json.loads(output.getvalue())
        assert 'section_run' in archive
        assert 'section_calculation_info' in archive

    def test_subsection(self, backend: LocalBackend, no_warn):
        backend.openSection('section_run')
        backend.openSection('section_method')
        backend.closeSection('section_method', -1)

        backend.openSection('section_method')
        backend.closeSection('section_method', -1)

        backend.openSection('section_run')
        backend.closeSection('section_run', 0)
        backend.closeSection('section_run', 1)

        backend.openSection('section_method')
        backend.closeSection('section_method', -1)

        runs = backend.data['section_run']
        assert len(runs) == 2
        assert len(runs[0]['section_method']) == 2
        assert len(runs[1]['section_method']) == 1

    def test_context(self, backend: LocalBackend, no_warn):
        backend.openSection('section_run')
        backend.openSection('section_method')
        backend.closeSection('section_method', -1)
        backend.closeSection('section_run', -1)

        backend.openSection('section_run')
        backend.closeSection('section_run', -1)

        backend.openContext('/section_run/0')
        backend.addValue('program_name', 't1')
        backend.closeContext('/section_run/0')

        backend.openContext('/section_run/1')
        backend.addValue('program_name', 't2')
        backend.closeContext('/section_run/1')

        backend.openContext('/section_run/0/section_method/0')
        backend.closeContext('/section_run/0/section_method/0')

        runs = backend.data['section_run']
        assert runs[0]['program_name'] == 't1'
        assert runs[1]['program_name'] == 't2'

    def test_multi_context(self, backend: LocalBackend, no_warn):
        backend.openSection('section_run')
        backend.closeSection('section_run', -1)

        backend.openContext('/section_run/0')
        backend.openSection('section_method')
        backend.closeSection('section_method', -1)
        backend.closeContext('/section_run/0')

        backend.openContext('/section_run/0')
        backend.openSection('section_method')
        backend.closeSection('section_method', -1)
        backend.closeContext('/section_run/0')

        assert len(backend.data['section_method']) == 1

    def test_bad_context(self, backend: LocalBackend, no_warn):
        try:
            backend.openContext('section_run/0')
            assert False
        except BadContextURI:
            pass

        try:
            backend.openContext('dsfds')
            assert False
        except BadContextURI:
            pass


def create_reference(data, pretty):
    if (pretty):
        return json.dumps(data, indent=2)
    else:
        return json.dumps(data, separators=(',', ':'))


@pytest.mark.parametrize("pretty", [False, True])
def test_stream_generator(pretty, no_warn):
    example_data = [
        {
            'key1': 'value',
            'key2': 1
        },
        {
            'key': {
                'key': 'value'
            }
        }
    ]

    out = StringIO()
    writer = JSONStreamWriter(out, pretty=pretty)
    writer.open_array()
    writer.open_object()
    writer.key('key1')
    writer.value('value')
    writer.key('key2')
    writer.value(1)
    writer.close_object()
    writer.open_object()
    writer.key('key')
    writer.open_object()
    writer.key('key')
    writer.value('value')
    writer.close_object()
    writer.close_object()
    writer.close_array()
    writer.close()

    assert create_reference(example_data, pretty) == out.getvalue()


def assert_parser_result(backend):
    status, errors = backend.status
    assert status == 'ParseSuccess'
    assert errors is None or len(errors) == 0


def run_parser(parser_name, mainfile):
    parser = parser_dict[parser_name]
    result = parser.run(mainfile, logger=utils.get_logger(__name__))
    return add_calculation_info(result)


@pytest.fixture
def parsed_vasp_example() -> LocalBackend:
    return run_parser(
        'parsers/vasp', 'dependencies/parsers/vasp/test/examples/xml/perovskite.xml')


@pytest.fixture
def parsed_template_example() -> LocalBackend:
    return run_parser(
        'parsers/template', 'tests/data/parsers/template.json')


@pytest.fixture(
    params=faulty_unknown_one_d_matid_example, ids=lambda spec: '%s-%s' % spec)
def parsed_faulty_unknown_matid_example(caplog, request) -> LocalBackend:
    parser_name, mainfile = request.param
    return run_parser(parser_name, mainfile)


@pytest.fixture(params=parser_examples, ids=lambda spec: '%s-%s' % spec)
def parsed_example(request) -> LocalBackend:
    parser_name, mainfile = request.param
    result = run_parser(parser_name, mainfile)
    return result


def add_calculation_info(backend: LocalBackend) -> LocalBackend:
    backend.openNonOverlappingSection('section_calculation_info')
    backend.addValue('upload_id', 'test_upload_id')
    backend.addValue('calc_id', 'test_calc_id')
    backend.addValue('calc_hash', 'test_calc_hash')
    backend.addValue('main_file', 'test/mainfile.txt')
    backend.addValue('parser_name', 'testParser')
    backend.closeNonOverlappingSection('section_calculation_info')
    return backend


@pytest.mark.parametrize('parser_name, mainfile', parser_examples)
def test_parser(parser_name, mainfile):
    parsed_example = run_parser(parser_name, mainfile)
    assert_parser_result(parsed_example)


def test_match(raw_files, no_warn):
    example_upload_id = 'example_upload_id'
    upload_files = files.StagingUploadFiles(example_upload_id, create=True, is_authorized=lambda: True)
    upload_files.add_rawfiles('tests/data/parsers')

    count = 0
    for mainfile in upload_files.raw_file_manifest():
        parser = match_parser(mainfile, upload_files)
        if parser is not None:
            count += 1

    assert count == correct_num_output_files
