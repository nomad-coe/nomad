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

from nomad.app.v1.models import WithQuery
from nomad.search import update_by_query
from nomad.search.v1 import quantity_values, search
from nomad.metainfo.elasticsearch_extension import entry_type, entry_index, material_index

from tests.utils import ExampleData


@pytest.fixture()
def example_data(elastic, test_user, other_test_user):
    data = ExampleData(uploader=test_user)

    for i in range(0, 4):
        data.create_entry(
            upload_id='test_upload_id',
            calc_id=f'test_entry_id_{i}',
            mainfile='test_content/test_embargo_entry/mainfile.json',
            shared_with=[],
            with_embargo=True)

    data.save(with_files=False, with_mongo=False)


def test_index(indices, example_data):
    assert material_index.get(id='test_material_id') is not None
    assert entry_index.get(id='test_entry_id_0') is not None


@pytest.fixture()
def indices(elastic):
    pass


def test_indices(indices):
    assert entry_type.quantities.get('entry_id') is not None
    assert entry_type.quantities.get('upload_id') is not None


@pytest.mark.parametrize('api_query, total', [
    pytest.param('{}', 4, id='empty'),
    pytest.param('{"results.method.simulation.program_name": "VASP"}', 4, id="match"),
    pytest.param('{"results.method.simulation.program_name": "VASP", "results.method.simulation.dft.xc_functional_type": "dne"}', 0, id="match_all"),
    pytest.param('{"and": [{"results.method.simulation.program_name": "VASP"}, {"results.method.simulation.dft.xc_functional_type": "dne"}]}', 0, id="and"),
    pytest.param('{"or":[{"results.method.simulation.program_name": "VASP"}, {"results.method.simulation.dft.xc_functional_type": "dne"}]}', 4, id="or"),
    pytest.param('{"not":{"results.method.simulation.program_name": "VASP"}}', 0, id="not"),
    pytest.param('{"results.method.simulation.program_name": {"all": ["VASP", "dne"]}}', 0, id="all"),
    pytest.param('{"results.method.simulation.program_name": {"any": ["VASP", "dne"]}}', 4, id="any"),
    pytest.param('{"results.method.simulation.program_name": {"none": ["VASP", "dne"]}}', 0, id="none"),
    pytest.param('{"results.method.simulation.program_name": {"gte": "VASP"}}', 4, id="gte"),
    pytest.param('{"results.method.simulation.program_name": {"gt": "A"}}', 4, id="gt"),
    pytest.param('{"results.method.simulation.program_name": {"lte": "VASP"}}', 4, id="lte"),
    pytest.param('{"results.method.simulation.program_name": {"lt": "A"}}', 0, id="lt"),
])
def test_search_query(indices, example_data, api_query, total):
    api_query = json.loads(api_query)
    results = search(owner='all', query=WithQuery(query=api_query).query)
    assert results.pagination.total == total  # pylint: disable=no-member


def test_update_by_query(indices, example_data):
    update_by_query(
        update_script='''
            ctx._source.entry_id = "other test id";
        ''',
        owner='all', query={}, index='v1')

    entry_index.refresh()

    results = search(owner='all', query=dict(entry_id='other test id'))
    assert results.pagination.total == 4


def test_quantity_values(indices, example_data):
    results = list(quantity_values('entry_id', page_size=1, owner='all'))
    assert results == ['test_entry_id_0', 'test_entry_id_1', 'test_entry_id_2', 'test_entry_id_3']
