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
import zipfile
import io
import json

from nomad.metainfo.search_extension import search_quantities
from nomad.app.v1.models import AggregateableQuantity, Metric

from tests.utils import assert_at_least, assert_url_query_args
from tests.test_files import example_mainfile_contents, append_raw_files  # pylint: disable=unused-import

from .common import assert_response
from tests.app.conftest import example_data as data  # pylint: disable=unused-import

'''
These are the tests for all API operations below ``entries``. The tests are organized
using the following type of methods: fixtures, ``perfrom_*_test``, ``assert_*``, and
``test_*``. While some ``test_*`` methods test individual API operations, some
test methods will test multiple API operations that use common aspects like
supporting queries, pagination, or the owner parameter. The test methods will use
``perform_*_test`` methods as an parameter. Similarely, the ``assert_*`` methods allow
to assert for certain aspects in the responses.
'''


def perform_entries_metadata_test(
        client, owner=None, headers={}, status_code=200,
        entries=None, http_method='get', **kwargs):

    if http_method == 'get':
        params = {}
        if owner is not None:
            params['owner'] = owner
        for value in kwargs.values():
            params.update(**value)
        response = client.get(
            'entries?%s' % urlencode(params, doseq=True), headers=headers)

    elif http_method == 'post':
        body = dict(**kwargs)
        if owner is not None:
            body['owner'] = owner
        response = client.post('entries/query', headers=headers, json=body)

    else:
        assert False

    response_json = assert_entries_metadata_response(response, status_code=status_code)

    if response_json is None:
        return

    assert 'pagination' in response_json
    if entries is not None and entries >= 0:
        assert response_json['pagination']['total'] == entries

    return response_json


def perform_entries_raw_download_test(
        client, headers={}, query={}, owner=None, files={}, entries=-1, files_per_entry=5,
        status_code=200, http_method='get'):

    if owner == 'all':
        # This operation is not allow for owner 'all'
        status_code = 401

    if http_method == 'post':
        body = {'query': query, 'files': files}
        if owner is not None:
            body['owner'] = owner
        response = client.post('entries/raw/download/query', headers=headers, json=body)

    elif http_method == 'get':
        params = dict(**query)
        params.update(**files)
        if owner is not None:
            params['owner'] = owner
        response = client.get('entries/raw/download?%s' % urlencode(params, doseq=True), headers=headers)

    else:
        assert False

    assert_response(response, status_code)
    if status_code == 200:
        assert_raw_zip_file(
            response, files=entries * files_per_entry + 1, manifest_entries=entries,
            compressed=files.get('compress', False))


def perform_entries_raw_test(
        client, owner=None, headers={}, status_code=200,
        entries=None, http_method='get', files_per_entry=-1, **kwargs):

    if owner == 'all':
        # This operation is not allow for owner 'all'
        status_code = 401

    if http_method == 'get':
        params = {}
        if owner is not None:
            params['owner'] = owner
        for value in kwargs.values():
            params.update(**value)
        response = client.get(
            'entries/raw?%s' % urlencode(params, doseq=True), headers=headers)

    elif http_method == 'post':
        body = dict(**kwargs)
        if owner is not None:
            body['owner'] = owner
        response = client.post('entries/raw/query', headers=headers, json=body)

    else:
        assert False

    response_json = assert_entries_raw_metadata_response(response, status_code=status_code)

    if response_json is None:
        return None

    assert 'pagination' in response_json
    if entries is not None:
        assert response_json['pagination']['total'] == entries

    assert_entries_raw_response(response_json, files_per_entry=files_per_entry)

    return response_json


def perform_entries_archive_download_test(
        client, headers={}, query={}, owner=None, files={},
        entries=-1, status_code=200, http_method='get'):

    if owner == 'all':
        # This operation is not allow for owner 'all'
        status_code = 401

    if http_method == 'post':
        body = {'query': query, 'files': files}
        if owner is not None:
            body['owner'] = owner
        response = client.post('entries/archive/download/query', headers=headers, json=body)

    elif http_method == 'get':
        params = dict(**query)
        params.update(**files)
        if owner is not None:
            params['owner'] = owner
        response = client.get('entries/archive/download?%s' % urlencode(params, doseq=True), headers=headers)

    else:
        assert False

    assert_response(response, status_code)
    if status_code == 200:
        assert_archive_zip_file(
            response, entries=entries,
            compressed=files.get('compress', False))


def perform_entries_archive_test(
        client, headers={}, entries=-1, status_code=200, http_method='get', **kwargs):

    if kwargs.get('owner') == 'all':
        # This operation is not allow for owner 'all'
        status_code = 401

    if http_method == 'get':
        assert 'required' not in kwargs
        params = {}
        if 'owner' in kwargs: params.update(owner=kwargs['owner'])
        if 'query' in kwargs: params.update(**kwargs['query'])
        if 'pagination' in kwargs: params.update(**kwargs['pagination'])
        response = client.get('entries/archive?%s' % urlencode(params, doseq=True), headers=headers)

    else:
        body = dict(**kwargs)
        response = client.post('entries/archive/query', headers=headers, json=body)

    assert_response(response, status_code)
    if status_code != 200:
        return None

    json_response = response.json()
    if entries >= 0:
        assert json_response['pagination']['total'] == entries
    for archive_data in json_response['data']:
        required = kwargs.get('required', '*')
        archive = archive_data['archive']
        if required == '*':
            for key in ['section_metadata', 'section_run']:
                assert key in archive
        else:
            for key in required: assert key in archive
            for key in archive: assert key in required

    return json_response


def assert_entry_metadata(response_json):
    if isinstance(response_json['data'], list):
        entries = response_json['data']
    else:
        entries = [response_json['data']]

    for entry in entries:
        if 'required' not in response_json:
            assert 'license' in entry

        if 'uploader' in entry:
            assert 'email' not in entry['uploader']


def assert_entries_metadata_response(response, status_code=None):
    response_json = assert_entries_raw_metadata_response(response, status_code=status_code)
    if response_json is not None:
        assert_entry_metadata(response_json)
    return response_json


def assert_entries_raw_metadata_response(response, status_code=None):
    assert_response(response, status_code)

    if status_code != 200 or response.status_code != 200:
        return None

    response_json = response.json()
    assert 'es_query' not in response_json
    assert 'data' in response_json
    return response_json


def assert_statistic(response_json, name, statistic, size=-1):
    assert 'statistics' in response_json
    assert name in response_json['statistics']
    statistic_response = response_json['statistics'][name]
    for key in ['data', 'size', 'order', 'quantity']:
        assert key in statistic_response

    assert_at_least(statistic, statistic_response)

    default_size = search_quantities[statistic['quantity']].statistic_size
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


def assert_required(data, required):
    if 'include' in required:
        for key in data:
            assert key in required['include'] or '%s.*' % key in required['include'] or key == 'calc_id'
    if 'exclude' in required:
        for key in required['exclude']:
            assert key not in data or key == 'calc_id'


def assert_aggregations(response_json, name, agg, total: int, size: int):
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
                    assert_required(entry, agg['entries']['required'])


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


def assert_raw_zip_file(
        response, files: int = -1, manifest_entries: int = -1, compressed: bool = False):

    manifest_keys = ['calc_id', 'upload_id', 'mainfile']

    assert len(response.content) > 0
    with zipfile.ZipFile(io.BytesIO(response.content)) as zip_file:
        with zip_file.open('manifest.json', 'r') as f:
            manifest = json.load(f)

        with_missing_files = any(entry['calc_id'] == 'id_02' for entry in manifest)
        with_overlapping_files = any(entry['calc_id'] == 'id_11' for entry in manifest)

        assert zip_file.testzip() is None
        zip_files = set(zip_file.namelist())
        if files >= 0:
            if with_missing_files or with_overlapping_files:
                assert files - (5 if with_missing_files else 0) - (4 if with_overlapping_files else 0) <= len(zip_files) < files
            else:
                assert len(zip_files) == files
            assert (zip_file.getinfo(zip_file.namelist()[0]).compress_type > 0) == compressed

        for path in zip_files:
            assert path == 'manifest.json' or path.startswith('id_')

        if manifest_entries >= 0:
            assert len(manifest) == manifest_entries

            for entry in manifest:
                if 'mainfile' in manifest:
                    manifest['mainfile'] in zip_files
                assert all(key in entry for key in manifest_keys)
                assert all(key in manifest_keys for key in entry)


def assert_entries_raw_response(response_json, files_per_entry: int = -1):
    assert 'data' in response_json
    for entry in response_json['data']:
        assert_entry_raw(entry, files_per_entry)


def assert_entry_raw_response(response_json, files_per_entry: int = -1):
    for key in ['entry_id', 'data']:
        assert key in response_json
    assert_entry_raw(response_json['data'], files_per_entry=files_per_entry)


def assert_entry_raw(data, files_per_entry: int = -1):
    for key in ['upload_id', 'calc_id', 'files']:
        assert key in data
    files = data['files']
    if files_per_entry >= 0:
        if data['calc_id'] == 'id_02':
            # missing files
            assert len(files) == 0
        elif data['calc_id'] in ['id_10', 'id_11']:
            # overlapping files
            assert len(files) == files_per_entry + 1
        else:
            assert len(files) == files_per_entry
    for file_ in files:
        assert 'size' in file_
        assert 'path' in file_


def assert_archive_zip_file(response, entries: int = -1, compressed: bool = False):
    manifest_keys = ['calc_id', 'upload_id', 'path', 'parser_name']

    assert len(response.content) > 0
    with zipfile.ZipFile(io.BytesIO(response.content)) as zip_file:
        assert zip_file.testzip() is None
        with zip_file.open('manifest.json', 'r') as f:
            manifest = json.load(f)

        with_missing_files = any(entry['calc_id'] == 'id_02' for entry in manifest)

        zip_files = set(zip_file.namelist())
        if entries >= 0:
            assert len(zip_files) == entries + 1 - (1 if with_missing_files else 0)
            assert (zip_file.getinfo(zip_file.namelist()[0]).compress_type > 0) == compressed

        for path in zip_files:
            assert path.endswith('.json')
            with zip_file.open(path, 'r') as f:
                data = json.load(f)
                if path != 'manifest.json':
                    for key in ['calc_id', 'archive']:
                        assert key in data
                    assert_archive(data['archive'])

        if entries >= 0:
            assert len(manifest) == entries

            for entry in manifest:
                if 'mainfile' in manifest:
                    manifest['path'] in zip_files
                assert all(key in entry for key in manifest_keys)
                assert all(key in manifest_keys for key in entry)


def assert_archive_response(response_json, required=None):
    for key in ['entry_id', 'required', 'data']:
        assert key in response_json
    if required is not None:
        assert required == response_json['required']
    for key in ['calc_id', 'upload_id', 'parser_name', 'archive']:
        assert key in response_json['data']
    assert_archive(response_json['data']['archive'], required=required)


def assert_archive(archive, required=None):
    for key in ['section_metadata']:
        assert key in archive


n_code_names = search_quantities['dft.code_name'].statistic_size


@pytest.mark.parametrize('statistic, size, status_code, user', [
    pytest.param({'quantity': 'dft.code_name'}, n_code_names, 200, None, id='fixed-values'),
    pytest.param({'quantity': 'dft.code_name', 'metrics': ['uploads']}, n_code_names, 200, None, id='metrics'),
    pytest.param({'quantity': 'dft.code_name', 'metrics': ['does not exist']}, -1, 422, None, id='bad-metric'),
    pytest.param({'quantity': 'calc_id', 'size': 1000}, 23, 200, None, id='size-to-large'),
    pytest.param({'quantity': 'calc_id', 'size': 10}, 10, 200, None, id='size'),
    pytest.param({'quantity': 'calc_id', 'size': -1}, -1, 422, None, id='bad-size-1'),
    pytest.param({'quantity': 'calc_id', 'size': 0}, -1, 422, None, id='bad-size-2'),
    pytest.param({'quantity': 'calc_id'}, 20, 200, None, id='size-default'),
    pytest.param({'quantity': 'calc_id', 'value_filter': '_0'}, 9, 200, None, id='filter'),
    pytest.param({'quantity': 'calc_id', 'value_filter': '.*_0.*'}, -1, 422, None, id='bad-filter'),
    pytest.param({'quantity': 'upload_id', 'order': {'type': 'values'}}, 3, 200, 'test_user', id='order-type'),
    pytest.param({'quantity': 'upload_id', 'order': {'direction': 'asc'}}, 3, 200, 'test_user', id='order-direction'),
    pytest.param({'quantity': 'does not exist'}, -1, 422, None, id='bad-quantity')])
def test_entries_statistics(client, data, test_user_auth, statistic, size, status_code, user):
    statistics = {'test_statistic': statistic}
    headers = {}
    if user == 'test_user':
        headers = test_user_auth

    response_json = perform_entries_metadata_test(
        client, headers=headers, owner='visible', statistics=statistics,
        status_code=status_code, http_method='post')

    if response_json is None:
        return

    assert_statistic(response_json, 'test_statistic', statistic, size=size)


def test_entries_statistics_ignore_size(client, data):
    statistic = {'quantity': 'dft.code_name', 'size': 10}
    statistics = {'test_statistic': statistic}
    response_json = perform_entries_metadata_test(
        client, statistics=statistics, status_code=200, http_method='post')
    statistic.update(size=n_code_names)
    assert_statistic(response_json, 'test_statistic', statistic, size=n_code_names)


def test_entries_all_statistics(client, data):
    statistics = {
        quantity.value: {'quantity': quantity.value, 'metrics': [metric.value for metric in Metric]}
        for quantity in AggregateableQuantity}
    response_json = perform_entries_metadata_test(
        client, statistics=statistics, status_code=200, http_method='post')
    for name, statistic in statistics.items():
        assert_statistic(response_json, name, statistic)


@pytest.mark.parametrize('aggregation, total, size, status_code', [
    pytest.param({'quantity': 'upload_id', 'pagination': {'order_by': 'uploader'}}, 3, 3, 200, id='order-str'),
    pytest.param({'quantity': 'upload_id', 'pagination': {'order_by': 'upload_time'}}, 3, 3, 200, id='order-date'),
    pytest.param({'quantity': 'upload_id', 'pagination': {'order_by': 'dft.n_calculations'}}, 3, 3, 200, id='order-int'),
    pytest.param({'quantity': 'dft.labels_springer_classification'}, 0, 0, 200, id='no-results'),
    pytest.param({'quantity': 'upload_id', 'pagination': {'page_after_value': 'id_published'}}, 3, 1, 200, id='after'),
    pytest.param({'quantity': 'upload_id', 'pagination': {'order_by': 'uploader', 'page_after_value': 'Sheldon Cooper:id_published'}}, 3, 1, 200, id='after-order'),
    pytest.param({'quantity': 'upload_id', 'entries': {'size': 10}}, 3, 3, 200, id='entries'),
    pytest.param({'quantity': 'upload_id', 'entries': {'size': 1}}, 3, 3, 200, id='entries-size'),
    pytest.param({'quantity': 'upload_id', 'entries': {'size': 0}}, -1, -1, 422, id='bad-entries'),
    pytest.param({'quantity': 'upload_id', 'entries': {'size': 10, 'required': {'include': ['calc_id', 'uploader']}}}, 3, 3, 200, id='entries-include'),
    pytest.param({'quantity': 'upload_id', 'entries': {'size': 10, 'required': {'exclude': ['files', 'mainfile']}}}, 3, 3, 200, id='entries-exclude')
])
def test_entries_aggregations(client, data, test_user_auth, aggregation, total, size, status_code):
    headers = test_user_auth
    aggregations = {'test_agg_name': aggregation}
    response_json = perform_entries_metadata_test(
        client, headers=headers, owner='visible', aggregations=aggregations,
        pagination=dict(page_size=0),
        status_code=status_code, http_method='post')

    if response_json is None:
        return

    assert_aggregations(response_json, 'test_agg_name', aggregation, total=total, size=size)


@pytest.mark.parametrize('required, status_code', [
    pytest.param({'include': ['calc_id', 'upload_id']}, 200, id='include'),
    pytest.param({'include': ['dft.*', 'upload_id']}, 200, id='include-section'),
    pytest.param({'exclude': ['upload_id']}, 200, id='exclude'),
    pytest.param({'exclude': ['missspelled', 'upload_id']}, 422, id='bad-quantitiy'),
    pytest.param({'exclude': ['calc_id']}, 200, id='exclude-id'),
    pytest.param({'exclude': ['dft.optimade']}, 200, id='exclude-sub-section'),
    pytest.param({'exclude': ['files', 'dft.optimade', 'dft.quantities']}, 200, id='exclude-multiple'),
    pytest.param({'include': ['upload_id']}, 200, id='include-id')
])
@pytest.mark.parametrize('http_method', ['post', 'get'])
def test_entries_required(client, data, required, status_code, http_method):
    response_json = perform_entries_metadata_test(
        client, required=required, pagination={'page_size': 1}, status_code=status_code, http_method=http_method)

    if response_json is None:
        return

    assert_required(response_json['data'][0], required)


@pytest.mark.parametrize('entry_id, required, status_code', [
    pytest.param('id_01', {}, 200, id='id'),
    pytest.param('doesnotexist', {}, 404, id='404'),
    pytest.param('id_01', {'include': ['calc_id', 'upload_id']}, 200, id='include'),
    pytest.param('id_01', {'exclude': ['upload_id']}, 200, id='exclude'),
    pytest.param('id_01', {'exclude': ['calc_id', 'upload_id']}, 200, id='exclude-calc-id')
])
def test_entry_metadata(client, data, entry_id, required, status_code):
    response = client.get('entries/%s?%s' % (entry_id, urlencode(required, doseq=True)))
    response_json = assert_entries_metadata_response(response, status_code=status_code)

    if response_json is None:
        return

    assert_required(response_json['data'], required)


@pytest.mark.parametrize('query, files, entries, files_per_entry, status_code', [
    pytest.param({}, {}, 23, 5, 200, id='all'),
    pytest.param({'calc_id': 'id_01'}, {}, 1, 5, 200, id='all'),
    pytest.param({'dft.code_name': 'DOESNOTEXIST'}, {}, 0, 5, 200, id='empty')
])
@pytest.mark.parametrize('http_method', ['post', 'get'])
def test_entries_raw(client, data, query, files, entries, files_per_entry, status_code, http_method):
    perform_entries_raw_test(
        client, status_code=status_code, query=query, files=files, entries=entries,
        files_per_entry=files_per_entry, http_method=http_method)


@pytest.mark.parametrize('query, files, entries, files_per_entry, status_code', [
    pytest.param({}, {}, 23, 5, 200, id='all'),
    pytest.param({'dft.code_name': 'DOESNOTEXIST'}, {}, 0, 5, 200, id='empty'),
    pytest.param({}, {'glob_pattern': '*.json'}, 23, 1, 200, id='glob'),
    pytest.param({}, {'re_pattern': '[a-z]*\\.aux'}, 23, 4, 200, id='re'),
    pytest.param({}, {'re_pattern': 'test_entry_02'}, 1, 5, 200, id='re-filter-entries'),
    pytest.param({}, {'re_pattern': 'test_entry_02/.*\\.json'}, 1, 1, 200, id='re-filter-entries-and-files'),
    pytest.param({}, {'glob_pattern': '*.json', 're_pattern': '.*\\.aux'}, 23, 4, 200, id='re-overwrites-glob'),
    pytest.param({}, {'re_pattern': '**'}, -1, -1, 422, id='bad-re-pattern'),
    pytest.param({}, {'compress': True}, 23, 5, 200, id='compress')
])
@pytest.mark.parametrize('http_method', ['post', 'get'])
def test_entries_download_raw(client, data, query, files, entries, files_per_entry, status_code, http_method):
    perform_entries_raw_download_test(
        client, status_code=status_code, query=query, files=files, entries=entries,
        files_per_entry=files_per_entry, http_method=http_method)


@pytest.mark.parametrize('http_method', ['post', 'get'])
@pytest.mark.parametrize('test_method', [
    pytest.param(perform_entries_raw_download_test, id='raw-download'),
    pytest.param(perform_entries_archive_download_test, id='archive-download')])
def test_entries_download_max(monkeypatch, client, data, test_method, http_method):
    monkeypatch.setattr('nomad.config.max_entry_download', 20)

    test_method(client, status_code=400, http_method=http_method)


@pytest.mark.parametrize('entry_id, files_per_entry, status_code', [
    pytest.param('id_01', 5, 200, id='id'),
    pytest.param('id_embargo', -1, 404, id='404'),
    pytest.param('doesnotexist', -1, 404, id='404')])
def test_entry_raw(client, data, entry_id, files_per_entry, status_code):
    response = client.get('entries/%s/raw' % entry_id)
    assert_response(response, status_code)
    if status_code == 200:
        assert_entry_raw_response(response.json(), files_per_entry=files_per_entry)


@pytest.mark.parametrize('entry_id, files, files_per_entry, status_code', [
    pytest.param('id_01', {}, 5, 200, id='id'),
    pytest.param('doesnotexist', {}, -1, 404, id='404'),
    pytest.param('id_01', {'glob_pattern': '*.json'}, 1, 200, id='glob'),
    pytest.param('id_01', {'re_pattern': '[a-z]*\\.aux'}, 4, 200, id='re'),
    pytest.param('id_01', {'re_pattern': '**'}, -1, 422, id='bad-re-pattern'),
    pytest.param('id_01', {'compress': True}, 5, 200, id='compress')])
def test_entry_raw_download(client, data, entry_id, files, files_per_entry, status_code):
    response = client.get('entries/%s/raw/download?%s' % (entry_id, urlencode(files, doseq=True)))
    assert_response(response, status_code)
    if status_code == 200:
        assert_raw_zip_file(
            response, files=files_per_entry + 1, manifest_entries=1,
            compressed=files.get('compress', False))


@pytest.fixture(scope='module')
def data_with_compressed_files(data):
    append_raw_files(
        'id_published', 'tests/data/api/mainfile.xz', 'test_content/subdir/test_entry_02/mainfile.xz')
    append_raw_files(
        'id_published', 'tests/data/api/mainfile.gz', 'test_content/subdir/test_entry_02/mainfile.gz')


@pytest.mark.parametrize('entry_id, path, params, status_code', [
    pytest.param('id_01', 'mainfile.json', {}, 200, id='id'),
    pytest.param('doesnotexist', 'mainfile.json', {}, 404, id='404-entry'),
    pytest.param('id_01', 'doesnot.exist', {}, 404, id='404-file'),
    pytest.param('id_01', 'mainfile.json', {'offset': 10, 'length': 10}, 200, id='offset-length'),
    pytest.param('id_01', 'mainfile.json', {'length': 1000000}, 200, id='length-too-large'),
    pytest.param('id_01', 'mainfile.json', {'offset': 1000000}, 200, id='offset-too-large'),
    pytest.param('id_01', 'mainfile.json', {'offset': -1}, 422, id='bad-offset'),
    pytest.param('id_01', 'mainfile.json', {'length': -1}, 422, id='bad-length'),
    pytest.param('id_01', 'mainfile.json', {'decompress': True}, 200, id='decompress-json'),
    pytest.param('id_02', 'mainfile.xz', {'decompress': True}, 200, id='decompress-xz'),
    pytest.param('id_02', 'mainfile.gz', {'decompress': True}, 200, id='decompress-gz'),
    pytest.param('id_unpublished', 'mainfile.json', {}, 404, id='404-unpublished'),
    pytest.param('id_embargo', 'mainfile.json', {}, 404, id='404-embargo'),
    pytest.param('id_embargo', 'mainfile.json', {'user': 'test-user'}, 200, id='embargo'),
    pytest.param('id_embargo', 'mainfile.json', {'user': 'other-test-user'}, 404, id='404-embargo-shared'),
    pytest.param('id_embargo_shared', 'mainfile.json', {'user': 'other-test-user'}, 200, id='embargo-shared')
])
def test_entry_raw_download_file(
        client, data_with_compressed_files, example_mainfile_contents, test_user_auth, other_test_user_auth,
        entry_id, path, params, status_code):

    user = params.get('user')
    if user:
        del(params['user'])
        if user == 'test-user':
            headers = test_user_auth
        elif user == 'other-test-user':
            headers = other_test_user_auth
    else:
        headers = {}

    response = client.get(
        f'entries/{entry_id}/raw/download/{path}?{urlencode(params, doseq=True)}',
        headers=headers)

    assert_response(response, status_code)
    if status_code == 200:
        content = response.text
        if path.endswith('.json'):
            offset = params.get('offset', 0)
            length = params.get('length', len(example_mainfile_contents) - offset)
            assert content == example_mainfile_contents[offset:offset + length]
        else:
            assert content == 'test content\n'


@pytest.mark.parametrize('query, files, entries, status_code', [
    pytest.param({}, {}, 23, 200, id='all'),
    pytest.param({'dft.code_name': 'DOESNOTEXIST'}, {}, -1, 200, id='empty'),
    pytest.param({}, {'compress': True}, 23, 200, id='compress')
])
@pytest.mark.parametrize('http_method', ['post', 'get'])
def test_entries_archive_download(client, data, query, files, entries, status_code, http_method):
    perform_entries_archive_download_test(
        client, status_code=status_code, query=query, files=files, entries=entries,
        http_method=http_method)


@pytest.mark.parametrize('required, status_code', [
    pytest.param('*', 200, id='full'),
    pytest.param({'section_metadata': '*'}, 200, id='partial'),
    pytest.param({'section_run': {'section_system[NOTANINT]': '*'}}, 422, id='bad-required-1'),
    pytest.param({'section_metadata': {'owners[NOTANINT]': '*'}}, 422, id='bad-required-2'),
    pytest.param({'DOESNOTEXIST': '*'}, 422, id='bad-required-3')
])
def test_entries_archive(client, data, required, status_code):
    perform_entries_archive_test(
        client, status_code=status_code, required=required, http_method='post')


@pytest.mark.parametrize('entry_id, status_code', [
    pytest.param('id_01', 200, id='id'),
    pytest.param('id_02', 404, id='404-not-visible'),
    pytest.param('doesnotexist', 404, id='404-does-not-exist')])
def test_entry_archive(client, data, entry_id, status_code):
    response = client.get('entries/%s/archive' % entry_id)
    assert_response(response, status_code)
    if status_code == 200:
        assert_archive_response(response.json())


@pytest.mark.parametrize('entry_id, required, status_code', [
    pytest.param('id_01', '*', 200, id='full'),
    pytest.param('id_02', '*', 404, id='404'),
    pytest.param('id_01', {'section_metadata': '*'}, 200, id='partial'),
    pytest.param('id_01', {'section_run': {'section_system[NOTANINT]': '*'}}, 422, id='bad-required-1'),
    pytest.param('id_01', {'section_metadata': {'owners[NOTANINT]': '*'}}, 422, id='bad-required-2'),
    pytest.param('id_01', {'DOESNOTEXIST': '*'}, 422, id='bad-required-3'),
    pytest.param('id_01', {'resolve-inplace': 'NotBool', 'section_workflow': '*'}, 422, id='bad-required-4'),
    pytest.param('id_01', {'resolve-inplace': True, 'section_metadata': 'include-resolved'}, 200, id='resolve-inplace')
])
def test_entry_archive_query(client, data, entry_id, required, status_code):
    response = client.post('entries/%s/archive/query' % entry_id, json={
        'required': required
    })
    assert_response(response, status_code)
    if status_code == 200:
        assert_archive_response(response.json(), required=required)


def perform_entries_owner_test(
        client, test_user_auth, other_test_user_auth, admin_user_auth,
        owner, user, status_code, total, http_method, test_method):

    headers = None
    if user == 'test_user':
        headers = test_user_auth
    elif user == 'other_test_user':
        headers = other_test_user_auth
    elif user == 'admin_user':
        headers = admin_user_auth
    elif user == 'bad_user':
        headers = {'Authorization': 'Bearer NOTATOKEN'}

    test_method(
        client, headers=headers, owner=owner, status_code=status_code, entries=total,
        http_method=http_method)


@pytest.mark.parametrize('query, status_code, total', [
    pytest.param({}, 200, 23, id='empty'),
    pytest.param('str', 422, -1, id='not-dict'),
    pytest.param({'calc_id': 'id_01'}, 200, 1, id='match'),
    pytest.param({'mispelled': 'id_01'}, 422, -1, id='not-quantity'),
    pytest.param({'calc_id': ['id_01', 'id_02']}, 200, 0, id='match-list-0'),
    pytest.param({'calc_id': 'id_01', 'atoms': ['H', 'O']}, 200, 1, id='match-list-1'),
    pytest.param({'calc_id:any': ['id_01', 'id_02']}, 200, 2, id='any-short'),
    pytest.param({'calc_id': {'any': ['id_01', 'id_02']}}, 200, 2, id='any'),
    pytest.param({'calc_id': {'any': 'id_01'}}, 422, -1, id='any-not-list'),
    pytest.param({'calc_id:any': 'id_01'}, 422, -1, id='any-short-not-list'),
    pytest.param({'calc_id:gt': 'id_01'}, 200, 22, id='gt-short'),
    pytest.param({'calc_id': {'gt': 'id_01'}}, 200, 22, id='gt'),
    pytest.param({'calc_id': {'gt': ['id_01']}}, 422, 22, id='gt-list'),
    pytest.param({'calc_id': {'missspelled': 'id_01'}}, 422, -1, id='not-op'),
    pytest.param({'calc_id:lt': ['id_01']}, 422, -1, id='gt-shortlist'),
    pytest.param({'calc_id:misspelled': 'id_01'}, 422, -1, id='not-op-short'),
    pytest.param({'or': [{'calc_id': 'id_01'}, {'calc_id': 'id_02'}]}, 200, 2, id='or'),
    pytest.param({'or': {'calc_id': 'id_01', 'dft.code_name': 'VASP'}}, 422, -1, id='or-not-list'),
    pytest.param({'and': [{'calc_id': 'id_01'}, {'calc_id': 'id_02'}]}, 200, 0, id='and'),
    pytest.param({'not': {'calc_id': 'id_01'}}, 200, 22, id='not'),
    pytest.param({'not': [{'calc_id': 'id_01'}]}, 422, -1, id='not-list'),
    pytest.param({'not': {'not': {'calc_id': 'id_01'}}}, 200, 1, id='not-nested-not'),
    pytest.param({'not': {'calc_id:any': ['id_01', 'id_02']}}, 200, 21, id='not-nested-any'),
    pytest.param({'and': [{'calc_id:any': ['id_01', 'id_02']}, {'calc_id:any': ['id_02', 'id_03']}]}, 200, 1, id='and-nested-any'),
    pytest.param({'and': [{'not': {'calc_id': 'id_01'}}, {'not': {'calc_id': 'id_02'}}]}, 200, 21, id='not-nested-not')
])
@pytest.mark.parametrize('test_method', [
    pytest.param(perform_entries_metadata_test, id='metadata'),
    pytest.param(perform_entries_raw_download_test, id='raw-download'),
    pytest.param(perform_entries_raw_test, id='raw'),
    pytest.param(perform_entries_archive_test, id='archive'),
    pytest.param(perform_entries_archive_download_test, id='archive-download')])
def test_entries_post_query(client, data, query, status_code, total, test_method):
    response_json = test_method(client, query=query, status_code=status_code, entries=total, http_method='post')

    response = client.post('entries/query', json={'query': query})
    response_json = assert_entries_metadata_response(response, status_code=status_code)

    if response_json is None:
        return

    if 'pagination' not in response_json:
        return

    pagination = response_json['pagination']
    assert pagination['total'] == total
    assert pagination['page_size'] == 10
    assert pagination['order_by'] == 'calc_id'
    assert pagination['order'] == 'asc'
    assert ('next_page_after_value' in pagination) == (total > 10)


@pytest.mark.parametrize('query, status_code, total', [
    pytest.param({}, 200, 23, id='empty'),
    pytest.param({'calc_id': 'id_01'}, 200, 1, id='match'),
    pytest.param({'mispelled': 'id_01'}, 200, 23, id='not-quantity'),
    pytest.param({'calc_id': ['id_01', 'id_02']}, 200, 2, id='match-many-or'),
    pytest.param({'atoms': ['H', 'O']}, 200, 23, id='match-list-many-and-1'),
    pytest.param({'atoms': ['H', 'O', 'Zn']}, 200, 0, id='match-list-many-and-2'),
    pytest.param({'n_atoms': 2}, 200, 23, id='match-int'),
    pytest.param({'n_atoms__gt': 2}, 200, 0, id='gt-int'),
    pytest.param({'calc_id__any': ['id_01', 'id_02']}, 200, 2, id='any'),
    pytest.param({'calc_id__any': 'id_01'}, 200, 1, id='any-not-list'),
    pytest.param({'domain': ['dft', 'ems']}, 422, -1, id='list-not-op'),
    pytest.param({'calc_id__gt': 'id_01'}, 200, 22, id='gt'),
    pytest.param({'calc_id__gt': ['id_01', 'id_02']}, 422, -1, id='gt-list'),
    pytest.param({'calc_id__missspelled': 'id_01'}, 422, -1, id='not-op'),
    pytest.param({'q': 'calc_id__id_01'}, 200, 1, id='q-match'),
    pytest.param({'q': 'missspelled__id_01'}, 422, -1, id='q-bad-quantity'),
    pytest.param({'q': 'bad_encoded'}, 422, -1, id='q-bad-encode'),
    pytest.param({'q': 'n_atoms__2'}, 200, 23, id='q-match-int'),
    pytest.param({'q': 'n_atoms__gt__2'}, 200, 0, id='q-gt'),
    pytest.param({'q': 'dft.workflow.section_geometry_optimization.final_energy_difference__1e-24'}, 200, 0, id='foat'),
    pytest.param({'q': 'domain__dft'}, 200, 23, id='enum'),
    pytest.param({'q': 'upload_time__gt__2014-01-01'}, 200, 23, id='datetime'),
    pytest.param({'q': ['atoms__all__H', 'atoms__all__O']}, 200, 23, id='q-all'),
    pytest.param({'q': ['atoms__all__H', 'atoms__all__X']}, 200, 0, id='q-all')
])
@pytest.mark.parametrize('test_method', [
    pytest.param(perform_entries_metadata_test, id='metadata'),
    pytest.param(perform_entries_raw_download_test, id='raw-download'),
    pytest.param(perform_entries_raw_test, id='raw'),
    pytest.param(perform_entries_archive_test, id='archive'),
    pytest.param(perform_entries_archive_download_test, id='archive-download')])
def test_entries_get_query(client, data, query, status_code, total, test_method):
    response_json = test_method(
        client, query=query, status_code=status_code, entries=total, http_method='get')

    if response_json is None:
        return

    if 'pagination' not in response_json:
        return

    response = client.get('entries?%s' % urlencode(query, doseq=True))

    response_json = assert_entries_metadata_response(response, status_code=status_code)

    if response_json is None:
        return

    pagination = response_json['pagination']
    assert pagination['total'] == total
    assert pagination['page_size'] == 10
    assert pagination['order_by'] == 'calc_id'
    assert pagination['order'] == 'asc'
    assert ('next_page_after_value' in pagination) == (total > 10)


@pytest.mark.parametrize('owner, user, status_code, total', [
    pytest.param('user', None, 401, -1, id='user-wo-auth'),
    pytest.param('staging', None, 401, -1, id='staging-wo-auth'),
    pytest.param('visible', None, 200, 23, id='visible-wo-auth'),
    pytest.param('admin', None, 401, -1, id='admin-wo-auth'),
    pytest.param('shared', None, 401, -1, id='shared-wo-auth'),
    pytest.param('public', None, 200, 23, id='public-wo-auth'),

    pytest.param('user', 'test_user', 200, 27, id='user-test-user'),
    pytest.param('staging', 'test_user', 200, 2, id='staging-test-user'),
    pytest.param('visible', 'test_user', 200, 27, id='visible-test-user'),
    pytest.param('admin', 'test_user', 401, -1, id='admin-test-user'),
    pytest.param('shared', 'test_user', 200, 27, id='shared-test-user'),
    pytest.param('public', 'test_user', 200, 23, id='public-test-user'),

    pytest.param('user', 'other_test_user', 200, 0, id='user-other-test-user'),
    pytest.param('staging', 'other_test_user', 200, 1, id='staging-other-test-user'),
    pytest.param('visible', 'other_test_user', 200, 25, id='visible-other-test-user'),
    pytest.param('shared', 'other_test_user', 200, 2, id='shared-other-test-user'),
    pytest.param('public', 'other_test_user', 200, 23, id='public-other-test-user'),

    pytest.param('all', None, 200, 25, id='metadata-all-wo-auth'),
    pytest.param('all', 'test_user', 200, 27, id='metadata-all-test-user'),
    pytest.param('all', 'other_test_user', 200, 26, id='metadata-all-other-test-user'),

    pytest.param('admin', 'admin_user', 200, 27, id='admin-admin-user'),
    pytest.param('all', 'bad_user', 401, -1, id='bad-credentials')
])
@pytest.mark.parametrize('http_method', ['post', 'get'])
@pytest.mark.parametrize('test_method', [
    pytest.param(perform_entries_metadata_test, id='metadata'),
    pytest.param(perform_entries_raw_download_test, id='raw-download'),
    pytest.param(perform_entries_raw_test, id='raw'),
    pytest.param(perform_entries_archive_test, id='archive'),
    pytest.param(perform_entries_archive_download_test, id='archive-download')])
def test_entries_owner(
        client, data, test_user_auth, other_test_user_auth, admin_user_auth,
        owner, user, status_code, total, http_method, test_method):

    perform_entries_owner_test(
        client, test_user_auth, other_test_user_auth, admin_user_auth,
        owner, user, status_code, total, http_method, test_method)


@pytest.mark.parametrize('pagination, response_pagination, status_code', [
    pytest.param({}, {'total': 23, 'page_size': 10, 'next_page_after_value': 'id_10'}, 200, id='empty'),
    pytest.param({'page_size': 1}, {'total': 23, 'page_size': 1, 'next_page_after_value': 'id_01'}, 200, id='size'),
    pytest.param({'page_size': 0}, {'total': 23, 'page_size': 0}, 200, id='size-0'),
    pytest.param({'page_size': 1, 'page_after_value': 'id_01'}, {'page_after_value': 'id_01', 'next_page_after_value': 'id_02'}, 200, id='after'),
    pytest.param({'page_size': 1, 'page_after_value': 'id_02', 'order': 'desc'}, {'next_page_after_value': 'id_01'}, 200, id='after-desc'),
    pytest.param({'page_size': 1, 'order_by': 'n_atoms'}, {'next_page_after_value': '2:id_01'}, 200, id='order-by-after-int'),
    pytest.param({'page_size': 1, 'order_by': 'dft.code_name'}, {'next_page_after_value': 'VASP:id_01'}, 200, id='order-by-after-nested'),
    pytest.param({'page_size': -1}, None, 422, id='bad-size'),
    pytest.param({'order': 'misspelled'}, None, 422, id='bad-order'),
    pytest.param({'order_by': 'misspelled'}, None, 422, id='bad-order-by'),
    pytest.param({'order_by': 'atoms', 'page_after_value': 'H:id_01'}, None, 422, id='order-by-list'),
    pytest.param({'order_by': 'n_atoms', 'page_after_value': 'some'}, None, 400, id='order-by-bad-after'),
    pytest.param({'page': 1, 'page_size': 1}, {'total': 23, 'page_size': 1, 'next_page_after_value': 'id_02', 'page': 1}, 200, id='page-1'),
    pytest.param({'page': 2, 'page_size': 1}, {'total': 23, 'page_size': 1, 'next_page_after_value': 'id_03', 'page': 2}, 200, id='page-2'),
    pytest.param({'page': 1000, 'page_size': 10}, None, 422, id='page-too-large'),
    pytest.param({'page': 9999, 'page_size': 1}, None, 200, id='page-just-small-enough'),
])
@pytest.mark.parametrize('http_method', ['post', 'get'])
@pytest.mark.parametrize('test_method', [
    pytest.param(perform_entries_metadata_test, id='metadata'),
    pytest.param(perform_entries_raw_test, id='raw'),
    pytest.param(perform_entries_archive_test, id='archive')])
def test_entries_pagination(client, data, pagination, response_pagination, status_code, http_method, test_method):
    response_json = test_method(
        client, pagination=pagination, status_code=status_code, http_method=http_method)

    if response_json is None:
        return

    assert_pagination(pagination, response_json['pagination'], response_json['data'], is_get=(http_method == 'get'))
