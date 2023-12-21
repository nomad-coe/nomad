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

import os.path
import importlib
import sys

from nomad.utils.exampledata import ExampleData

from tests.parsing.test_parsing import run_singular_parser
from tests.normalizing.conftest import run_normalize
from tests.test_client import async_api_v1


def test_archive_query(async_api_v1, elastic, raw_files, mongo, test_user, capsys):
    mainfile = os.path.join(
        __file__, '..', '..', 'data', 'examples', 'archive_query_vasprun.xml.gz'
    )
    archive = run_singular_parser('parsers/vasp', mainfile)
    run_normalize(archive)
    archive.metadata.apply_archive_metadata(archive)

    data = ExampleData(main_author=test_user)
    data.create_upload('test_upload_id', published=True)
    data.create_entry(upload_id='test_upload_id', entry_archive=archive)
    data.save()

    sys.path.append('examples/archive')
    importlib.import_module('archive_query', 'nomad')

    captured = capsys.readouterr()
    assert captured.err == ''
    assert 'CaO4Ti2' in captured.out
