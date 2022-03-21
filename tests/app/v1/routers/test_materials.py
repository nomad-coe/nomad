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
from urllib.parse import urlencode

from nomad.metainfo.elasticsearch_extension import material_entry_type

from tests.test_files import example_mainfile_contents  # pylint: disable=unused-import

from .common import (
    aggregation_exclude_from_search_test_parameters, assert_pagination, assert_metadata_response, assert_required, assert_aggregations,
    perform_metadata_test, perform_owner_test, owner_test_parameters,
    post_query_test_parameters, get_query_test_parameters, pagination_test_parameters,
    aggregation_test_parameters)
from tests.conftest import example_data  # pylint: disable=unused-import

'''
These are the tests for all API operations below ``entries``. The tests are organized
using the following type of methods: fixtures, ``perform*_test``, ``assert_*``, and
``test_*``. While some ``test_*`` methods test individual API operations, some
test methods will test multiple API operations that use common aspects like
supporting queries, pagination, or the owner parameter. The test methods will use
``perform_*_test`` methods as an parameter. Similarly, the ``assert_*`` methods allow
to assert for certain aspects in the responses.
'''


def perform_materials_metadata_test(*args, **kwargs):
    kwargs.update(endpoint='materials')
    return perform_metadata_test(*args, **kwargs)


program_name = 'entries.results.method.simulation.program_name'


@pytest.mark.parametrize(
    'aggregation, total, size, status_code, user',
    aggregation_test_parameters(
        entity_id='material_id', material_prefix='', entry_prefix='entries.', total=6))
def test_materials_aggregations(client, example_data, test_user_auth, aggregation, total, size, status_code, user):
    headers = {}
    if user == 'test_user':
        headers = test_user_auth

    aggregations = {'test_agg_name': aggregation}

    response_json = perform_materials_metadata_test(
        client, headers=headers, owner='visible', aggregations=aggregations,
        pagination=dict(page_size=0),
        status_code=status_code, http_method='post')

    if response_json is None:
        return

    for aggregation_obj in aggregation.values():
        assert_aggregations(
            response_json, 'test_agg_name', aggregation_obj, total=total, size=size,
            default_key='material_id')

        # make sure that (nested) terms aggregation produce sensible counts
        if 'terms' in aggregation:
            for bucket in response_json['aggregations']['test_agg_name']['terms']['data']:
                assert bucket['count'] <= 6  # the total number of materials, counting entries we would surpass this


@pytest.mark.parametrize(
    'query,agg_data,total,status_code',
    aggregation_exclude_from_search_test_parameters(resource='materials', total_per_entity=3, total=6))
def test_materials_aggregations_exclude_from_search(client, example_data, query, agg_data, total, status_code):
    aggs, types, lengths = agg_data
    response_json = perform_materials_metadata_test(
        client, owner='visible',
        query=query, aggregations=aggs,
        pagination=dict(page_size=0),
        status_code=status_code, http_method='post')

    if response_json is None:
        return

    assert response_json['pagination']['total'] == total
    for i, (type, length) in enumerate(zip(types, lengths)):
        response_agg = response_json['aggregations'][f'agg_{i}'][type]
        assert len(response_agg['data']) == length


@pytest.mark.parametrize('required, status_code', [
    pytest.param({'include': ['material_id', program_name]}, 200, id='include'),
    pytest.param({'include': ['entries.*', program_name]}, 200, id='include-section'),
    pytest.param({'exclude': [program_name]}, 200, id='exclude'),
    pytest.param({'exclude': ['missspelled', program_name]}, 422, id='bad-quantitiy'),
    pytest.param({'exclude': ['material_id']}, 200, id='exclude-id'),
    pytest.param({'exclude': ['entries.results.*']}, 200, id='exclude-sub-section'),
    pytest.param({'exclude': [program_name, 'entries.results.method.*']}, 200, id='exclude-multiple'),
    pytest.param({'include': [program_name]}, 200, id='include-id')
])
@pytest.mark.parametrize('http_method', ['post', 'get'])
def test_materials_required(client, example_data, required, status_code, http_method):
    response_json = perform_materials_metadata_test(
        client, required=required, pagination={'page_size': 1}, status_code=status_code, http_method=http_method)

    if response_json is None:
        return

    assert_required(response_json['data'][0], required, default_key='material_id')


@pytest.mark.parametrize('material_id, required, status_code', [
    pytest.param('id_01', {}, 200, id='id'),
    pytest.param('doesnotexist', {}, 404, id='404'),
    pytest.param('id_01', {'include': ['material_id', 'n_elements']}, 200, id='include'),
    pytest.param('id_01', {'exclude': ['n_elements']}, 200, id='exclude'),
    pytest.param('id_01', {'exclude': ['material_id', 'n_elements']}, 200, id='exclude-id')
])
def test_material_metadata(client, example_data, material_id, required, status_code):
    response = client.get('materials/%s?%s' % (material_id, urlencode(required, doseq=True)))
    response_json = assert_metadata_response(response, status_code=status_code)

    if response_json is None:
        return

    assert_required(response_json['data'], required, default_key='material_id')


@pytest.mark.parametrize(
    'query, status_code, total',
    post_query_test_parameters('material_id', total=6, material_prefix='', entry_prefix='entries.') + [
        pytest.param({'entries.entry_id': 'id_01'}, 200, 1, id='entries-single'),
        pytest.param({'entries.entry_id:any': ['id_01', 'id_02']}, 200, 1, id='any-entry-same-material'),
        pytest.param({'entries.entry_id:any': ['id_01', 'id_05']}, 200, 2, id='any-entry-diff-material'),
        pytest.param({'entries.entry_id:all': ['id_01', 'id_02']}, 200, 0, id='all-entry-same-material'),
        pytest.param({'entries.entry_id:all': ['id_01', 'id_05']}, 200, 0, id='all-entry-diff-material'),
        pytest.param({
            'and': [
                {'entries.entry_id': 'id_01'},
                {'entries.entry_id': 'id_02'}
            ]
        }, 200, 1, id='per-entry-same-material'),
        pytest.param({
            'and': [
                {'entries.entry_id': 'id_01'},
                {'entries.entry_id': 'id_05'}
            ]
        }, 200, 0, id='per-entry-diff-material'),
        pytest.param({'entries': {'entry_id:any': ['id_01', 'id_02']}}, 200, 1, id='alt-any-entry-same-material'),
        pytest.param({'entries': {'entry_id:any': ['id_01', 'id_05']}}, 200, 2, id='alt-any-entry-diff-material'),
        pytest.param({'entries': {'entry_id:all': ['id_01', 'id_02']}}, 200, 0, id='alt-all-entry-same-material'),
        pytest.param({'entries': {'entry_id:all': ['id_01', 'id_05']}}, 200, 0, id='alt-all-entry-diff-material'),
        pytest.param({'entry_id': 'id_01'}, 422, 0, id='not-material-quantity'),
        pytest.param({'entries.material_id': 'id_01'}, 422, 0, id='not-entry-quantity')
    ])
def test_materials_post_query(client, example_data, query, status_code, total):
    response_json = perform_materials_metadata_test(
        client, query=query, status_code=status_code, total=total,
        http_method='post')

    response = client.post('materials/query', json={'query': query})
    response_json = assert_metadata_response(response, status_code=status_code)

    if response_json is None:
        return

    if 'pagination' not in response_json:
        return

    pagination = response_json['pagination']
    assert pagination['total'] == total
    assert pagination['page_size'] == 10
    assert pagination['order_by'] == 'material_id'
    assert pagination['order'] == 'asc'
    assert ('next_page_after_value' in pagination) == (total > 10)


@pytest.mark.parametrize('query, status_code, total', get_query_test_parameters(
    'material_id', total=6, material_prefix='', entry_prefix='entries.'))
def test_materials_get_query(client, example_data, query, status_code, total):
    assert 'entries.upload_create_time' in material_entry_type.quantities

    response_json = perform_materials_metadata_test(
        client, query=query, status_code=status_code, total=total, http_method='get')

    if response_json is None:
        return

    if 'pagination' not in response_json:
        return

    response = client.get('materials?%s' % urlencode(query, doseq=True))

    response_json = assert_metadata_response(response, status_code=status_code)

    if response_json is None:
        return

    pagination = response_json['pagination']
    assert pagination['total'] == total
    assert pagination['page_size'] == 10
    assert pagination['order_by'] == 'material_id'
    assert pagination['order'] == 'asc'
    assert ('next_page_after_value' in pagination) == (total > 10)


@pytest.mark.parametrize(
    'owner, user, status_code, total_entries, total_mainfiles, total_materials',
    owner_test_parameters())
@pytest.mark.parametrize('http_method', ['post', 'get'])
@pytest.mark.parametrize('test_method', [
    pytest.param(perform_materials_metadata_test, id='metadata')])
def test_materials_owner(
        client, example_data, test_user_auth, other_test_user_auth, admin_user_auth,
        owner, user, status_code, total_entries, total_mainfiles, total_materials,
        http_method, test_method):

    perform_owner_test(
        client, test_user_auth, other_test_user_auth, admin_user_auth,
        owner, user, status_code, total_materials, http_method, test_method)


@pytest.mark.parametrize('pagination, response_pagination, status_code', pagination_test_parameters(
    elements='elements', n_elements='n_elements', crystal_system='symmetry.crystal_system',
    total=6))
@pytest.mark.parametrize('http_method', ['post', 'get'])
def test_materials_pagination(client, example_data, pagination, response_pagination, status_code, http_method):
    response_json = perform_materials_metadata_test(
        client, pagination=pagination, status_code=status_code, http_method=http_method)

    if response_json is None:
        return

    assert_pagination(pagination, response_json['pagination'], response_json['data'], is_get=(http_method == 'get'))
