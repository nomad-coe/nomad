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

from devtools import debug
from urllib.parse import urlencode

from nomad.metainfo.elasticsearch_extension import DocumentType

from tests.utils import assert_at_least, assert_url_query_args


def assert_response(response, status_code=None):
    ''' General assertions for status_code and error messages '''
    if status_code and response.status_code != status_code:
        try:
            debug(response.json())
        except Exception:
            pass

    if status_code is not None:
        if response.status_code != status_code and response.status_code == 422:
            print(response.json()['detail'])
        assert response.status_code == status_code

    if status_code == 422:
        response_json = response.json()
        details = response_json['detail']
        assert len(details) > 0
        for detail in details:
            assert 'loc' in detail
            assert 'msg' in detail
        return

    if 400 <= status_code < 500:
        response_json = response.json()
        assert 'detail' in response_json


def assert_base_metadata_response(response, status_code=None):
    assert_response(response, status_code)

    if status_code != 200 or response.status_code != 200:
        return None

    response_json = response.json()
    assert 'es_query' not in response_json
    assert 'data' in response_json
    return response_json


def assert_metadata(response_json):
    if isinstance(response_json['data'], list):
        metadatas = response_json['data']
    else:
        metadatas = [response_json['data']]

    for metadata in metadatas:
        if 'required' not in response_json:
            assert 'license' in metadata


def assert_metadata_response(response, status_code=None):
    response_json = assert_base_metadata_response(response, status_code=status_code)
    if response_json is not None:
        assert_metadata(response_json)
    return response_json


def assert_statistic(response_json, name, statistic, doc_type: DocumentType, size=-1):
    assert 'statistics' in response_json
    assert name in response_json['statistics']
    statistic_response = response_json['statistics'][name]
    for key in ['data', 'size', 'order', 'quantity']:
        assert key in statistic_response

    assert_at_least(statistic, statistic_response)

    default_size = doc_type.quantities[statistic['quantity']].statistics_size
    assert statistic.get('size', default_size) >= len(statistic_response['data'])

    if size != -1:
        assert len(statistic_response['data']) == size

    values = list(statistic_response['data'].keys())
    for index, value in enumerate(values):
        data = statistic_response['data'][value]
        assert 'entries' in data
        for metric in statistic.get('metrics', []):
            assert metric in data

        if index < len(values) - 1:

            def order_value(value, data):
                if statistic_response['order']['type'] == 'entries':
                    return data['entries']
                else:
                    return value

            if statistic_response['order']['direction'] == 'asc':
                assert order_value(value, data) <= order_value(values[index + 1], statistic_response['data'][values[index + 1]])
            else:
                assert order_value(value, data) >= order_value(values[index + 1], statistic_response['data'][values[index + 1]])

    if 'order' in statistic:
        assert statistic_response['order']['type'] == statistic['order'].get('type', 'entries')
        assert statistic_response['order']['direction'] == statistic['order'].get('direction', 'desc')


def assert_required(data, required, default_key: str):
    if 'include' in required:
        for key in data:
            assert key in required['include'] or '%s.*' % key in required['include'] or key == default_key
    if 'exclude' in required:
        for key in required['exclude']:
            assert key not in data or key == default_key


def assert_aggregations(response_json, name, agg, total: int, size: int, default_key: str):
    assert 'aggregations' in response_json
    assert name in response_json['aggregations']
    agg_response = response_json['aggregations'][name]

    for key in ['data', 'pagination', 'quantity']:
        assert key in agg_response

    assert_at_least(agg, agg_response)

    n_data = len(agg_response['data'])
    assert agg.get('pagination', {}).get('page_size', 10) >= n_data
    assert agg_response['pagination']['total'] >= n_data
    for item in agg_response['data'].values():
        for key in ['size']:
            assert key in item
            assert item['size'] > 0
    if size >= 0:
        assert n_data == size
    if total >= 0:
        assert agg_response['pagination']['total'] == total

    if 'entries' in agg:
        agg_data = [item['data'][0] for item in agg_response['data'].values()]
    else:
        agg_data = [{agg['quantity']: value} for value in agg_response['data']]

    if 'pagination' in agg:
        assert_pagination(agg['pagination'], agg_response['pagination'], agg_data, is_get=False)
    else:
        assert_pagination({}, agg_response['pagination'], agg_data, order_by=agg['quantity'], is_get=False)

    if 'entries' in agg:
        for item in agg_response['data'].values():
            assert 'data' in item
            assert agg['entries'].get(size, 10) >= len(item['data']) > 0
            if 'required' in agg['entries']:
                for entry in item['data']:
                    assert_required(entry, agg['entries']['required'], default_key=default_key)


def assert_pagination(pagination, pagination_response, data, order_by=None, order=None, is_get=True):
    assert_at_least(pagination, pagination_response)
    assert len(data) <= pagination_response['page_size']
    assert len(data) <= pagination_response['total']

    if order is None:
        order = pagination_response.get('order', 'asc')
    if order_by is None:
        order_by = pagination_response.get('order_by')

    if order_by is not None:
        for index, item in enumerate(data):
            if index < len(data) - 1 and order_by in item:
                if order == 'desc':
                    assert item[order_by] >= data[index + 1][order_by]
                else:
                    assert item[order_by] <= data[index + 1][order_by]

    if is_get:
        page_size = pagination_response['page_size']
        page = pagination_response.get('page')
        page_url = pagination_response.get('page_url')
        first_page_url = pagination_response.get('first_page_url')
        prev_page_url = pagination_response.get('prev_page_url')
        next_page_url = pagination_response.get('next_page_url')
        next_page_after_value = pagination_response.get('next_page_after_value')

        assert page_url
        if page_size:
            assert first_page_url
            assert_url_query_args(first_page_url, page_after_value=None, page=None)
        if next_page_after_value:
            assert next_page_url
            assert_url_query_args(next_page_url, page_after_value=next_page_after_value, page=None)
        if page and page > 1:
            assert prev_page_url
            assert_url_query_args(prev_page_url, page=page - 1, page_after_value=None)


def perform_metadata_test(
        client, endpoint: str, owner=None, headers={}, status_code=200,
        total=None, http_method='get', **kwargs):

    if http_method == 'get':
        params = {}
        if owner is not None:
            params['owner'] = owner
        for value in kwargs.values():
            params.update(**value)
        response = client.get(
            f'{endpoint}?{urlencode(params, doseq=True)}', headers=headers)

    elif http_method == 'post':
        body = dict(**kwargs)
        if owner is not None:
            body['owner'] = owner
        response = client.post(f'{endpoint}/query', headers=headers, json=body)

    else:
        assert False

    response_json = assert_metadata_response(response, status_code=status_code)

    if response_json is None:
        return

    assert 'pagination' in response_json
    if total is not None and total >= 0:
        assert response_json['pagination']['total'] == total

    return response_json
