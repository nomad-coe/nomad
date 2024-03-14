#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import pytest
import json
import os
import os.path

from nomad import config
from nomad.parsing.parser import ArchiveParser
from nomad.datamodel import EntryArchive, Context


def test_archive_parser(raw_files_function):
    archive_data = {
        'definitions': {
            'section_definitions': [
                {
                    'name': 'TestSection',
                    'base_sections': ['nomad.datamodel.data.EntryData'],
                    'quantities': [{'name': 'test_quantity', 'type': 'str'}],
                }
            ]
        },
        'data': {
            'm_def': '#/definitions/section_definitions/0',
            'test_quantity': 'test_value',
        },
    }

    mainfile = os.path.join(config.fs.tmp, 'test_mainfile.archive.json')
    with open(mainfile, 'wt') as f:
        json.dump(archive_data, f)

    archive = EntryArchive()
    archive.m_context = Context()
    ArchiveParser().parse(mainfile, archive)

    assert archive.data.m_to_dict() == archive_data['data']


def get_file_parameter():
    example_files = [
        'schema.archive.yaml',
        'schema.archive.json',
        'intra-entry.archive.json',
    ]
    path = os.walk(os.path.join(os.path.dirname(__file__), '../../examples/data'))
    for root, _, files in path:
        for file in files:
            if os.path.basename(file) in example_files:
                yield pytest.param(os.path.join(root, file), id=file)


@pytest.mark.parametrize('mainfile', get_file_parameter())
def test_example_data(mainfile, no_warn):
    archive = EntryArchive()
    archive.m_context = Context()
    ArchiveParser().parse(mainfile, archive)
    archive.m_to_dict()
