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

from nomad.datamodel import results
from nomad.metainfo.elasticsearch_extension import material_entry_type

from tests.test_files import example_mainfile_contents  # pylint: disable=unused-import

from .common import (
    assert_pagination, assert_metadata_response, assert_required,
    perform_metadata_test, perform_owner_test, owner_test_parameters,
    post_query_test_parameters, get_query_test_parameters, pagination_test_parameters)
from ..conftest import example_data as data  # pylint: disable=unused-import

'''
These are the tests for all API operations below ``entries``. The tests are organized
using the following type of methods: fixtures, ``perfrom_*_test``, ``assert_*``, and
``test_*``. While some ``test_*`` methods test individual API operations, some
test methods will test multiple API operations that use common aspects like
supporting queries, pagination, or the owner parameter. The test methods will use
``perform_*_test`` methods as an parameter. Similarely, the ``assert_*`` methods allow
to assert for certain aspects in the responses.
'''


def perform_materials_metadata_test(*args, **kwargs):
    kwargs.update(endpoint='materials')
    return perform_metadata_test(*args, **kwargs)


n_code_names = results.Simulation.program_name.a_elasticsearch.statistics_size
program_name = 'entries.results.method.simulation.program_name'


# @pytest.mark.parametrize('statistic, size, status_code, user', [
#     pytest.param({'quantity': program_name}, n_code_names, 200, None, id='fixed-values'),
#     pytest.param({'quantity': program_name, 'metrics': ['uploads']}, n_code_names, 200, None, id='metrics'),
#     pytest.param({'quantity': program_name, 'metrics': ['does not exist']}, -1, 422, None, id='bad-metric'),
#     pytest.param({'quantity': 'entry_id', 'size': 1000}, 23, 200, None, id='size-to-large'),
#     pytest.param({'quantity': 'entry_id', 'size': 10}, 10, 200, None, id='size'),
#     pytest.param({'quantity': 'entry_id', 'size': -1}, -1, 422, None, id='bad-size-1'),
#     pytest.param({'quantity': 'entry_id', 'size': 0}, -1, 422, None, id='bad-size-2'),
#     pytest.param({'quantity': 'entry_id'}, 10, 200, None, id='size-default'),
#     pytest.param({'quantity': 'entry_id', 'value_filter': '_0'}, 9, 200, None, id='filter'),
#     pytest.param({'quantity': 'entry_id', 'value_filter': '.*_0.*'}, -1, 422, None, id='bad-filter'),
#     pytest.param({'quantity': 'upload_id', 'order': {'type': 'values'}}, 3, 200, 'test_user', id='order-type'),
#     pytest.param({'quantity': 'upload_id', 'order': {'direction': 'asc'}}, 3, 200, 'test_user', id='order-direction'),
#     pytest.param({'quantity': 'does not exist'}, -1, 422, None, id='bad-quantity')])
# def test_entries_statistics(client, data, test_user_auth, statistic, size, status_code, user):
#     statistics = {'test_statistic': statistic}
#     headers = {}
#     if user == 'test_user':
#         headers = test_user_auth

#     response_json = perform_materials_metadata_test(
#         client, headers=headers, owner='visible', statistics=statistics,
#         status_code=status_code, http_method='post')

#     if response_json is None:
#         return

#     assert_statistic(response_json, 'test_statistic', statistic, size=size)


# # TODO is this really the desired behavior
# def test_entries_statistics_ignore_size(client, data):
#     statistic = {'quantity': program_name, 'size': 10}
#     statistics = {'test_statistic': statistic}
#     response_json = perform_materials_metadata_test(
#         client, statistics=statistics, status_code=200, http_method='post')
#     statistic.update(size=n_code_names)
#     assert_statistic(response_json, 'test_statistic', statistic, size=n_code_names)


# def test_entries_all_statistics(client, data):
#     statistics = {
#         quantity: {'quantity': quantity, 'metrics': [metric for metric in entry_type.metrics]}
#         for quantity in entry_type.quantities if entry_type.quantities[quantity].aggregateable}
#     response_json = perform_materials_metadata_test(
#         client, statistics=statistics, status_code=200, http_method='post')
#     for name, statistic in statistics.items():
#         assert_statistic(response_json, name, statistic)


# @pytest.mark.parametrize('aggregation, total, size, status_code', [
#     pytest.param({'quantity': 'upload_id', 'pagination': {'order_by': 'uploader.user_id'}}, 3, 3, 200, id='order-str'),
#     pytest.param({'quantity': 'upload_id', 'pagination': {'order_by': 'upload_time'}}, 3, 3, 200, id='order-date'),
#     pytest.param({'quantity': 'upload_id', 'pagination': {'order_by': 'results.properties.n_calculations'}}, 3, 3, 200, id='order-int'),
#     pytest.param({'quantity': 'results.material.symmetry.structure_name'}, 0, 0, 200, id='no-results'),
#     pytest.param({'quantity': 'upload_id', 'pagination': {'page_after_value': 'id_published'}}, 3, 1, 200, id='after'),
#     pytest.param({'quantity': 'upload_id', 'pagination': {'order_by': 'uploader.name', 'page_after_value': 'Sheldon Cooper:id_published'}}, 3, 1, 200, id='after-order'),
#     pytest.param({'quantity': 'upload_id', 'entries': {'size': 10}}, 3, 3, 200, id='entries'),
#     pytest.param({'quantity': 'upload_id', 'entries': {'size': 1}}, 3, 3, 200, id='entries-size'),
#     pytest.param({'quantity': 'upload_id', 'entries': {'size': 0}}, -1, -1, 422, id='bad-entries'),
#     pytest.param({'quantity': 'upload_id', 'entries': {'size': 10, 'required': {'include': ['entry_id', 'uploader.*']}}}, 3, 3, 200, id='entries-include'),
#     pytest.param({'quantity': 'upload_id', 'entries': {'size': 10, 'required': {'exclude': ['files', 'mainfile']}}}, 3, 3, 200, id='entries-exclude')
# ])
# def test_entries_aggregations(client, data, test_user_auth, aggregation, total, size, status_code):
#     headers = test_user_auth
#     aggregations = {'test_agg_name': aggregation}
#     response_json = perform_materials_metadata_test(
#         client, headers=headers, owner='visible', aggregations=aggregations,
#         pagination=dict(page_size=0),
#         status_code=status_code, http_method='post')

#     if response_json is None:
#         return

#     assert_aggregations(response_json, 'test_agg_name', aggregation, total=total, size=size)


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
def test_entries_required(client, data, required, status_code, http_method):
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
def test_entry_metadata(client, data, material_id, required, status_code):
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
def test_materials_post_query(client, data, query, status_code, total):
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
def test_materials_get_query(client, data, query, status_code, total):
    assert 'entries.upload_time' in material_entry_type.quantities

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


@pytest.mark.parametrize('owner, user, status_code, total', owner_test_parameters(total=6))
@pytest.mark.parametrize('http_method', ['post', 'get'])
@pytest.mark.parametrize('test_method', [
    pytest.param(perform_materials_metadata_test, id='metadata')])
def test_materials_owner(
        client, data, test_user_auth, other_test_user_auth, admin_user_auth,
        owner, user, status_code, total, http_method, test_method):

    perform_owner_test(
        client, test_user_auth, other_test_user_auth, admin_user_auth,
        owner, user, status_code, total, http_method, test_method)


@pytest.mark.parametrize('pagination, response_pagination, status_code', pagination_test_parameters(
    elements='elements', n_elements='n_elements', crystal_system='symmetry.crystal_system',
    total=6))
@pytest.mark.parametrize('http_method', ['post', 'get'])
def test_materials_pagination(client, data, pagination, response_pagination, status_code, http_method):
    response_json = perform_materials_metadata_test(
        client, pagination=pagination, status_code=status_code, http_method=http_method)

    if response_json is None:
        return

    assert_pagination(pagination, response_json['pagination'], response_json['data'], is_get=(http_method == 'get'))
