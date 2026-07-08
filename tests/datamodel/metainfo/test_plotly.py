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

import json

from tests.normalizing.conftest import run_processing


def test_plotly_snapshot(raw_files_function):
    directory = 'tests/data/datamodel/metainfo/plotly'
    mainfile = 'plotly.schema.archive.yaml'
    plotly_archive = run_processing(directory, mainfile)

    f = open('tests/data/datamodel/metainfo/plotly/snapshot.archive.json')
    snapshot = json.load(f)
    f.close()

    figures = plotly_archive.data['figures']
    snapshot_figures = snapshot['figures']
    for i in range(0, 4):
        assert json.dumps(figures[i].figure, sort_keys=True) == json.dumps(
            snapshot_figures[i], sort_keys=True
        )
