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
import numpy as np

from nomad.units import ureg
from nomad.datamodel.datamodel import EntryArchive
from nomad.datamodel.results import Results, Material, System
from nomad.datamodel.metainfo.simulation.run import Run
from nomad.datamodel.metainfo.simulation.system import System as SystemRun, Atoms
from nomad.utils.exampledata import ExampleData

from .common import assert_response


@pytest.fixture(scope="module")
def example_data_systems(elastic_module, mongo_module, test_user):
    data = ExampleData(main_author=test_user)
    upload_id = "systems_upload"

    data.create_upload(
        upload_id=upload_id,
        published=True
    )
    atoms = Atoms(
        n_atoms=2,
        labels=["C", "H"],
        species=[6, 1],
        positions=np.array([[0, 0, 0], [1, 1, 1]]) * ureg.angstrom,
        lattice_vectors=np.array([[5, 0, 0], [0, 5, 0], [0, 0, 5]]) * ureg.angstrom,
        periodic=[True, True, True],
    )
    data.create_entry(
        upload_id=upload_id,
        entry_id="systems_entry_1",
        mainfile="test_content/test_entry/main-file.json",
        entry_archive=EntryArchive(
            run=[Run(system=[SystemRun(atoms=atoms)])],
            results=Results(
                material=Material(
                    topology=[System(atoms=atoms)]
                )
            )
        )
    )

    data.save()

    yield

    data.delete()
    from nomad.search import search
    assert search(query=dict(upload_id=upload_id)).pagination.total == 0


def run_query(entry_id, path, format, client):
    response = client.get(f'systems/{entry_id}/?path={path}&format={format}', headers={})
    return response


@pytest.mark.parametrize("path", [
    pytest.param('run/0/system/0', id='explicit indexing'),
    pytest.param('run/0/system/-1', id='negative indexing'),
    pytest.param('/run/0/system/0', id='start with slash'),
    pytest.param('/run/0/system/0', id='end with slash'),
    pytest.param('results/material/topology/0', id='saved in topology'),
])
def test_paths(path, client, example_data_systems):
    response = run_query('systems_entry_1', path, 'pdb', client)
    assert_response(response, 200)
