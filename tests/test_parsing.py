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

from nomad.parsing import JSONStreamWriter, parser_dict
from nomad.parsing import LocalBackend, BadContextURI


class TestLocalBackend(object):

    @pytest.fixture(scope='session')
    def meta_info(self):
        path = '.dependencies/nomad-meta-info/meta_info/nomad_meta_info/all.nomadmetainfo.json'
        meta_info, _ = loadJsonFile(path)
        return meta_info

    @pytest.fixture(scope='function')
    def backend(self, meta_info):
        return LocalBackend(meta_info, debug=True)

    def test_meta_info(self, meta_info):
        assert 'section_topology' in meta_info

    def test_section(self, backend):
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

    def test_subsection(self, backend: LocalBackend):
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

    def test_context(self, backend: LocalBackend):
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

    def test_multi_context(self, backend: LocalBackend):
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

    def test_bad_context(self, backend: LocalBackend):
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
def test_stream_generator(pretty):
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


@pytest.fixture
def parsed_vasp_example() -> LocalBackend:
    vasp_parser = parser_dict['parsers/vasp']
    example_mainfile = '.dependencies/parsers/vasp/test/examples/xml/perovskite.xml'
    return vasp_parser.run(example_mainfile)


def test_vasp_parser(parsed_vasp_example: LocalBackend):
    status, errors = parsed_vasp_example.status

    assert status == 'ParseSuccess'
    assert errors is None or len(errors) == 0


exciting_examples = [
    '.dependencies/parsers/exciting/test/examples/Ag/INFO.OUT',
    '.dependencies/parsers/exciting/test/examples/lithium/INFO.OUT',
    '.dependencies/parsers/exciting/test/examples/GW/INFO.OUT',
]


@pytest.mark.parametrize('example_mainfile', exciting_examples)
def test_exciting_parser(example_mainfile):
    exciting_parser = parser_dict['parsers/exciting']
    parsed_exciting_example = exciting_parser.run(example_mainfile)

    status, errors = parsed_exciting_example.status

    assert status == 'ParseSuccess'
    assert errors is None or len(errors) == 0


def test_match():
    vasp_parser = parser_dict['parsers/vasp']
    directory = './data/examples_vasp_6'

    count = 0
    for dirpath, _, filenames in os.walk(directory):
        for filename in filenames:
            fullname = os.path.join(dirpath, filename)
            if vasp_parser.is_mainfile(fullname, lambda fn: open(fn)):
                count += 1

    assert count == 6
