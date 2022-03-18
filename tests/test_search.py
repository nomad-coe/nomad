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

from typing import List, Dict, Any, Union, Iterable
import pytest
import json

from nomad import config, utils, infrastructure
from nomad.app.v1.models import WithQuery
from nomad.search import quantity_values, search, update_by_query, refresh
from nomad.metainfo.elasticsearch_extension import entry_type, entry_index, material_index
from nomad.utils.exampledata import ExampleData


def assert_search_upload(
        entries: Union[int, Iterable] = -1,
        additional_keys: List[str] = [],
        upload_id: str = None,
        **kwargs):

    if isinstance(entries, list):
        size = len(entries)
    elif isinstance(entries, int):
        size = entries
    else:
        assert False

    keys = ['entry_id', 'upload_id', 'mainfile']
    refresh()
    body: Dict[str, Any] = {}
    body.update(size=10)
    if upload_id is not None:
        body['query'] = dict(match=dict(upload_id=upload_id))

    search_results = infrastructure.elastic_client.search(
        index=config.elastic.entries_index, body=body)['hits']

    if size != -1:
        assert search_results['total']['value'] == size

    if search_results['total']['value'] > 0:
        for hit in search_results['hits']:
            hit = utils.flat(hit['_source'])
            for key, value in kwargs.items():
                assert hit.get(key, None) == value, key

            if 'pid' in hit:
                assert int(hit.get('pid')) > 0

            for key in keys:
                assert key in hit, f'{key} is missing'

            for key in additional_keys:
                assert key in hit, f'{key} is missing'
                assert hit[key] != config.services.unavailable_value

            for coauthor in hit.get('entry_coauthors', []):
                assert coauthor.get('name', None) is not None


def test_mapping_compatibility(elastic_infra):
    from nomad.infrastructure import elastic_client

    v0 = elastic_client.indices.get(config.elastic.entries_index)
    v1 = elastic_client.indices.get(config.elastic.entries_index)

    def get_mapping(index):
        assert len(index) == 1
        index = index[next(iter(index))]
        assert len(index['mappings']) == 1
        return index['mappings'][next(iter(index['mappings']))]

    v0, v1 = get_mapping(v0), get_mapping(v1)

    def compare(a, b, path='', results=None):
        if results is None:
            results = []
        if path != '':
            path += '.'
        for key in set(list(a.keys()) + list(b.keys())):
            if key in a and key in b:
                next_a, next_b = a[key], b[key]
                if isinstance(next_a, dict) and isinstance(next_b, dict):
                    compare(next_a, next_b, f'{path}{key}', results=results)
                    continue

                if next_a == next_b:
                    continue

            results.append(f"{'v0' if key in a else 'v1'}:{path}{key}")

        return results

    for diff in compare(v0, v1):
        # assert that there are only top-level differences and mapping types and fields are
        # the same
        assert len([c for c in diff if c == '.']) == 1, diff


@pytest.fixture()
def example_data(elastic, test_user):
    data = ExampleData(main_author=test_user)
    data.create_upload(upload_id='test_upload_id', published=True, embargo_length=12)
    for i in range(0, 4):
        data.create_entry(
            upload_id='test_upload_id',
            entry_id=f'test_entry_id_{i}',
            mainfile='test_content/test_embargo_entry/mainfile.json')

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
