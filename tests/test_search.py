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
from datetime import datetime

from nomad.utils import deep_get
from nomad import config, utils, infrastructure
from nomad.app.v1.models import WithQuery, MetadataRequired, MetadataPagination, Direction
from nomad.datamodel.datamodel import EntryArchive, EntryData, EntryMetadata
from nomad.metainfo.metainfo import Datetime, Quantity
from nomad.metainfo.util import MEnum
from nomad.search import quantity_values, search, update_by_query, refresh
from nomad.metainfo.elasticsearch_extension import entry_type, entry_index, material_index, schema_separator, dtype_separator
from nomad.utils.exampledata import ExampleData
from tests.config import yaml_schema_name, python_schema_name


def split(path):
    return [int(x) if x.isdigit() else x for x in path.split('.')]


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
            hit = utils.flatten_dict(hit['_source'])
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


def get_schema_fixture(type, request):
    return request.getfixturevalue(f'example_data_schema_{type}')


def get_schema_quantity(type, quantity):
    if quantity is None:
        return None
    if not quantity.startswith('data.'):
        return quantity
    dtype = ''
    if type == 'python':
        name = python_schema_name
    elif type == 'yaml':
        name = yaml_schema_name
        dtype = {
            'data.name': 'str',
            'data.count': 'int',
            'data.frequency': 'float'
        }[quantity]
        dtype = f'{dtype_separator}{dtype}'

    return f'{quantity}{schema_separator}{name}{dtype}'


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


@pytest.fixture()
def example_eln_data(elastic, test_user):
    class DataSection(EntryData):
        text = Quantity(type=str)
        keyword = Quantity(type=MEnum('one', 'two'))
        long = Quantity(type=int)
        double = Quantity(type=float)
        date = Quantity(type=Datetime)
        boolean = Quantity(type=bool)

    data = ExampleData(main_author=test_user)
    data.create_upload(upload_id='test_upload_id', published=True, embargo_length=12)

    parameters = [
        ('text', 'test value'),
        ('keyword', 'one'),
        ('long', 1),
        ('double', 1.2),
        ('date', datetime.fromtimestamp(0)),
        ('boolean', False)]

    for index, item in enumerate(parameters):
        quantity, value = item
        archive = EntryArchive(metadata=EntryMetadata(), data=DataSection())
        archive.data.m_set(DataSection.m_def.all_quantities[quantity], value)
        archive.metadata.apply_archive_metadata(archive)

        data.create_entry(
            entry_archive=archive,
            upload_id='test_upload_id',
            entry_id=f'test_entry_id_{index}',
            mainfile=f'test_content/test_embargo_entry/mainfile_{index}.archive.json')

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


@pytest.mark.parametrize('api_query, total', [
    pytest.param({}, 6, id='empty'),
    pytest.param(
        {"search_quantities": {"id": f"data.text{schema_separator}tests.test_search.DataSection", "str_value": "test"}},
        1, id='simple-positive'),
    pytest.param(
        {"search_quantities": {"id": f"data.text{schema_separator}tests.test_search.DataSection", "str_value": "wrong"}},
        0, id='simple-negative'),
    pytest.param(
        {"search_quantities": {"id": f"data.keyword{schema_separator}tests.test_search.DataSection", "str_value": "one"}}, 1, id='keyword-as-text'),
    pytest.param(
        {"search_quantities": {"id": f"data.long{schema_separator}tests.test_search.DataSection", "int_value": {"lte": 1}}}, 1, id='int'),
    pytest.param(
        {"search_quantities": {"id": f"data.double{schema_separator}tests.test_search.DataSection", "float_value": {"lt": 1.3}}}, 1, id='float'),
    pytest.param(
        {"search_quantities": {"id": f"data.date{schema_separator}tests.test_search.DataSection", "datetime_value": {"lt": "1971-01-01"}}}, 1, id='datetime'),
    pytest.param(
        {"search_quantities": {"id": f"data.boolean{schema_separator}tests.test_search.DataSection", "bool_value": False}}, 1, id='bool'),
    pytest.param(
        {"search_quantities": {"id": f"data.text{schema_separator}tests.test_search.DataSection", "str_value": "one"}}, 0, id='uses-nested'),
])
def test_search_query_dynamic(indices, example_eln_data, api_query, total):
    '''Tests that search queries targeting dynamic quantities work for different data types.'''
    results = search(owner='all', query=WithQuery(query=api_query).query)
    assert results.pagination.total == total  # pylint: disable=no-member


@pytest.mark.parametrize('include, exclude, include_response, exclude_response', [
    pytest.param(
        None,
        None,
        ['search_quantities.0.str_value', f'data.name'],
        [],
        id='no-required'
    ),
    pytest.param(
        [],
        None,
        [],
        ['search_quantities.0.str_value', f'data.name'],
        id='include-nothing'
    ),
    pytest.param(
        None,
        ['search_quantities*'],
        [f'data.name'],
        ['search_quantities.0.str_value'],
        id='exclude-searchable-quantities'
    ),
    pytest.param(
        None,
        ['search_quantities.str_value'],
        [f'data.name', 'search_quantities.0.id'],
        ['search_quantities/0/str_value'],
        id='exclude-searchable-quantities-subset'
    ),
    pytest.param(
        None,
        [f'data.name'],
        ['search_quantities.0.str_value'],
        [f'data.name'],
        id='exclude-dynamic'
    ),
    pytest.param(
        ['search_quantities.str_value'],
        None,
        ['search_quantities.0.str_value'],
        [f'data.name', 'search_quantities.0.id'],
        id='include-searchable-quantities'
    ),
    pytest.param(
        [f'data.name'],
        None,
        [f'data.name'],
        ['search_quantities/0/id'],
        id='include-dynamic'
    ),
    pytest.param(
        [f'data.count', f'data.frequency'],
        None,
        [f'data.count', f'data.frequency'],
        ['search_quantities.0.int_value', 'search_quantities.0.float_value'],
        id='include-dynamic-multiple'
    ),
    pytest.param(
        None,
        [f'data.name', 'search_quantities*'],
        [],
        [f'data.name', 'search_quantities.0.id'],
        id='exclude-both'
    ),
    pytest.param(
        [f'data.name', 'search_quantities*'],
        None,
        [f'data.name', 'search_quantities.0.id'],
        [],
        id='include-both'
    ),
])
@pytest.mark.parametrize('schema_type', ['python', 'yaml'])
def test_search_query_dynamic_required(indices, plugin_schema, example_data_schema_python, schema_type, include, exclude, include_response, exclude_response):
    '''Tests that the requiring works correctly for dynamic fields.'''
    if include:
        include = [get_schema_quantity(schema_type, x) for x in include]
    if exclude:
        exclude = [get_schema_quantity(schema_type, x) for x in exclude]
    results = search(owner='all', required=MetadataRequired(include=include, exclude=exclude))

    for path in include_response:
        assert deep_get(results.data, 0, *split(path)) is not None
    for path in exclude_response:
        with pytest.raises(ValueError):
            deep_get(results.data, 0, *split(path))


@pytest.mark.parametrize('include', [
    pytest.param([(f'data.name', 'test0')], id='root-section'),
    pytest.param([(f'data.child.name', 'test_child0')], id='nested-section'),
    pytest.param([(f'data.child_repeating.0.name', 'test_child_repeating0')], id='nested-repeating-section'),
])
@pytest.mark.parametrize('schema_type', ['python', 'yaml'])
def test_search_hits_dynamic(indices, plugin_schema, schema_type, include, request):
    '''Tests that the response hit structure is properly reconstructed for
    dynamic quantities that target a schema.
    '''
    get_schema_fixture(schema_type, request)
    results = search(owner='all')
    for path, value in include:
        assert deep_get(results.data, 0, *split(path)) == value


@pytest.mark.parametrize('order_by, order, expected', [
    pytest.param(None, None, None, id='default'),
    pytest.param(
        f'data.name',
        Direction.asc,
        [f'test{i}' for i in [0, 1, 10, 11, 12, 13, 14, 2, 3, 4]],
        id='sort-string-single-asc'
    ),
    pytest.param(
        f'data.name',
        Direction.desc,
        [f'test{i}' for i in [9, 8, 7, 6, 5, 4, 3, 2, 14, 13]],
        id='sort-string-single-desc'
    ),
    pytest.param(
        f'data.count',
        Direction.asc,
        range(10),
        id='sort-int-single-asc'
    ),
    pytest.param(
        f'data.count',
        Direction.desc,
        range(14, 4, -1),
        id='sort-int-single-desc'
    ),
    pytest.param(
        f'data.frequency',
        Direction.asc,
        [x + 0.5 for x in range(10)],
        id='sort-float-single-asc'
    ),
    pytest.param(
        f'data.frequency',
        Direction.desc,
        [x + 0.5 for x in range(14, 4, -1)],
        id='sort-float-single-desc'
    ),
])
@pytest.mark.parametrize('schema_type', ['python', 'yaml'])
def test_pagination_dynamic(indices, plugin_schema, schema_type, order_by, order, expected, request):
    '''Tests that sorting by a dynamic field works as expected.
    '''
    get_schema_fixture(schema_type, request)
    order_by = get_schema_quantity(schema_type, order_by)
    pagination = MetadataPagination(order_by=order_by, page_size=10)
    if order is not None:
        pagination.order = order
    results = search(owner='all', pagination=pagination)
    if expected:
        assert len(results.data) == len(expected)
        path = order_by.split(schema_separator)[0]
        for i, data in enumerate(results.data):
            assert deep_get(data, *split(path)) == expected[i]
