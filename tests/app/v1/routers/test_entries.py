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

from nomad.datamodel import results
from nomad.metainfo.elasticsearch_extension import entry_type
from nomad.utils.exampledata import ExampleData

from tests.test_files import example_mainfile_contents, append_raw_files  # pylint: disable=unused-import

from .common import (
    aggregation_exclude_from_search_test_parameters, assert_response, assert_base_metadata_response,
    assert_metadata_response, assert_required, assert_aggregations, assert_pagination,
    perform_metadata_test, post_query_test_parameters, get_query_test_parameters,
    perform_owner_test, owner_test_parameters, pagination_test_parameters,
    aggregation_test_parameters)
from tests.conftest import example_data as data  # pylint: disable=unused-import

'''
These are the tests for all API operations below ``entries``. The tests are organized
using the following type of methods: fixtures, ``perfrom_*_test``, ``assert_*``, and
``test_*``. While some ``test_*`` methods test individual API operations, some
test methods will test multiple API operations that use common aspects like
supporting queries, pagination, or the owner parameter. The test methods will use
``perform_*_test`` methods as an parameter. Similarely, the ``assert_*`` methods allow
to assert for certain aspects in the responses.
'''


def perform_entries_metadata_test(*args, **kwargs):
    kwargs.update(endpoint='entries')
    return perform_metadata_test(*args, **kwargs)


def perform_entries_raw_test(
        client, headers={}, query={}, owner=None, files={}, total=-1, files_per_entry=5,
        status_code=200, http_method='get'):

    if type(total) == int:
        total_entries = total_mainfiles = total
    else:
        total_entries, total_mainfiles = total

    if owner == 'all':
        # This operation is not allow for owner 'all'
        status_code = 401

    if http_method == 'post':
        body = {'query': query, 'files': files}
        if owner is not None:
            body['owner'] = owner
        response = client.post('entries/raw/query', headers=headers, json=body)

    elif http_method == 'get':
        params = dict(**query)
        params.update(**files)
        if owner is not None:
            params['owner'] = owner
        response = client.get('entries/raw?%s' % urlencode(params, doseq=True), headers=headers)

    else:
        assert False

    assert_response(response, status_code)
    if status_code == 200:
        assert_raw_zip_file(
            response, files=total_mainfiles * files_per_entry + 1, manifest_entries=total_entries,
            compressed=files.get('compress', False))


def perform_entries_rawdir_test(
        client, owner=None, headers={}, status_code=200,
        total=None, http_method='get', files_per_entry=-1, **kwargs):

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
            'entries/rawdir?%s' % urlencode(params, doseq=True), headers=headers)

    elif http_method == 'post':
        body = dict(**kwargs)
        if owner is not None:
            body['owner'] = owner
        response = client.post('entries/rawdir/query', headers=headers, json=body)

    else:
        assert False

    response_json = assert_base_metadata_response(response, status_code=status_code)

    if response_json is None:
        return None

    assert 'pagination' in response_json
    if total is not None:
        assert response_json['pagination']['total'] == total

    assert_entries_rawdir_response(response_json, files_per_entry=files_per_entry)

    return response_json


def perform_entries_archive_download_test(
        client, headers={}, query={}, owner=None, files={},
        total=-1, status_code=200, http_method='get'):

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
            response, total=total,
            compressed=files.get('compress', False))


def perform_entries_archive_test(
        client, headers={}, total=-1, status_code=200, http_method='get', **kwargs):

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
    if total >= 0:
        assert json_response['pagination']['total'] == total
    for archive_data in json_response['data']:
        required = kwargs.get('required', '*')
        archive = archive_data['archive']
        if required == '*':
            for key in ['metadata', 'run']:
                assert key in archive
        else:
            for key in required: assert key in archive
            for key in archive: assert key in required

    return json_response


def assert_raw_zip_file(
        response, files: int = -1, manifest_entries: int = -1, compressed: bool = False):

    manifest_keys = ['entry_id', 'upload_id', 'mainfile']

    assert len(response.content) > 0
    with zipfile.ZipFile(io.BytesIO(response.content)) as zip_file:
        with zip_file.open('manifest.json', 'r') as f:
            manifest = json.load(f)

        with_missing_files = any(entry['entry_id'] == 'id_02' for entry in manifest)
        with_overlapping_files = any(entry['entry_id'] == 'id_11' for entry in manifest)

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


def assert_entries_rawdir_response(response_json, files_per_entry: int = -1):
    assert 'data' in response_json
    for entry in response_json['data']:
        assert_entry_rawdir(entry, files_per_entry)


def assert_entry_rawdir_response(response_json, files_per_entry: int = -1):
    for key in ['entry_id', 'data']:
        assert key in response_json
    assert_entry_rawdir(response_json['data'], files_per_entry=files_per_entry)


def assert_entry_rawdir(data, files_per_entry: int = -1):
    for key in ['upload_id', 'entry_id', 'files']:
        assert key in data
    files = data['files']
    if files_per_entry >= 0:
        if data['entry_id'] == 'id_02':
            # missing files
            assert len(files) == 0
        elif data['entry_id'] in ['id_10', 'id_11']:
            # overlapping files
            assert len(files) == files_per_entry + 1
        else:
            assert len(files) == files_per_entry
    for file_ in files:
        assert 'size' in file_
        assert 'path' in file_


def assert_archive_zip_file(response, total: int = -1, compressed: bool = False):
    manifest_keys = ['entry_id', 'upload_id', 'path', 'parser_name']

    assert len(response.content) > 0
    with zipfile.ZipFile(io.BytesIO(response.content)) as zip_file:
        assert zip_file.testzip() is None
        with zip_file.open('manifest.json', 'r') as f:
            manifest = json.load(f)

        with_missing_files = any(entry['entry_id'] == 'id_02' for entry in manifest)

        zip_files = set(zip_file.namelist())
        if total >= 0:
            assert len(zip_files) == total + 1 - (1 if with_missing_files else 0)
            assert (zip_file.getinfo(zip_file.namelist()[0]).compress_type > 0) == compressed

        for path in zip_files:
            assert path.endswith('.json')
            with zip_file.open(path, 'r') as f:
                data = json.load(f)
                if path != 'manifest.json':
                    for key in ['entry_id', 'archive']:
                        assert key in data
                    assert_archive(data['archive'])

        if total >= 0:
            assert len(manifest) == total

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
    for key in ['entry_id', 'upload_id', 'parser_name', 'archive']:
        assert key in response_json['data']
    assert_archive(response_json['data']['archive'], required=required)


def assert_archive(archive, required=None):
    for key in ['metadata']:
        assert key in archive


n_code_names = results.Simulation.program_name.a_elasticsearch[0].default_aggregation_size
program_name = 'results.method.simulation.program_name'


def test_entries_all_metrics(client, data):
    aggregations = {
        quantity: {
            'terms': {
                'quantity': quantity, 'metrics': [metric for metric in entry_type.metrics]
            }
        }
        for quantity in entry_type.quantities if entry_type.quantities[quantity].aggregateable}
    response_json = perform_entries_metadata_test(
        client, aggregations=aggregations, status_code=200, http_method='post')
    for name, agg in aggregations.items():
        assert_aggregations(response_json, name, agg['terms'])


@pytest.mark.parametrize(
    'aggregation, total, size, status_code, user',
    aggregation_test_parameters(entity_id='entry_id', material_prefix='results.material.', entry_prefix='', total=23) + [
        pytest.param(
            {
                'terms': {
                    'quantity': 'upload_id',
                    'entries': {
                        'size': 10,
                        'required': {'exclude': ['files', 'mainfile']}
                    }
                }
            },
            8, 8, 200, 'test_user', id='entries-exclude'),
        pytest.param(
            {'terms': {'quantity': 'entry_id', 'include': '_0'}},
            9, 9, 200, None, id='filter'),
        pytest.param(
            {'terms': {'quantity': 'entry_id', 'include': '.*_0.*'}},
            -1, -1, 422, None, id='bad-filter')
    ])
def test_entries_aggregations(client, data, test_user_auth, aggregation, total, size, status_code, user):
    headers = {}
    if user == 'test_user':
        headers = test_user_auth

    aggregations = {'test_agg_name': aggregation}

    response_json = perform_entries_metadata_test(
        client, headers=headers, owner='visible', aggregations=aggregations,
        pagination=dict(page_size=0),
        status_code=status_code, http_method='post')

    if response_json is None:
        return

    for aggregation_obj in aggregation.values():
        assert_aggregations(
            response_json, 'test_agg_name', aggregation_obj, total=total, size=size,
            default_key='entry_id')


@pytest.mark.parametrize(
    'query,agg_data,total,status_code',
    aggregation_exclude_from_search_test_parameters(resource='entries', total_per_entity=1, total=23))
def test_entries_aggregations_exclude_from_search(client, data, query, agg_data, total, status_code):
    aggs, types, lengths = agg_data
    response_json = perform_entries_metadata_test(
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
    pytest.param({'include': ['entry_id', 'upload_id']}, 200, id='include'),
    pytest.param({'include': ['results.*', 'upload_id']}, 200, id='include-section'),
    pytest.param({'exclude': ['upload_id']}, 200, id='exclude'),
    pytest.param({'exclude': ['missspelled', 'upload_id']}, 422, id='bad-quantitiy'),
    pytest.param({'exclude': ['entry_id']}, 200, id='exclude-id'),
    pytest.param({'exclude': ['results.material.*']}, 200, id='exclude-sub-section'),
    pytest.param({'exclude': ['files', 'results.material.*', 'results.method.*']}, 200, id='exclude-multiple'),
    pytest.param({'include': ['upload_id']}, 200, id='include-id')
])
@pytest.mark.parametrize('http_method', ['post', 'get'])
def test_entries_required(client, data, required, status_code, http_method):
    response_json = perform_entries_metadata_test(
        client, required=required, pagination={'page_size': 1}, status_code=status_code, http_method=http_method)

    if response_json is None:
        return

    assert_required(response_json['data'][0], required, default_key='entry_id')


@pytest.mark.parametrize('entry_id, required, status_code', [
    pytest.param('id_01', {}, 200, id='id'),
    pytest.param('doesnotexist', {}, 404, id='404'),
    pytest.param('id_01', {'include': ['entry_id', 'upload_id']}, 200, id='include'),
    pytest.param('id_01', {'exclude': ['upload_id']}, 200, id='exclude'),
    pytest.param('id_01', {'exclude': ['entry_id', 'upload_id']}, 200, id='exclude-entry-id')
])
def test_entry_metadata(client, data, entry_id, required, status_code):
    response = client.get('entries/%s?%s' % (entry_id, urlencode(required, doseq=True)))
    response_json = assert_metadata_response(response, status_code=status_code)

    if response_json is None:
        return

    assert_required(response_json['data'], required, default_key='entry_id')


@pytest.mark.parametrize('query, files, total, files_per_entry, status_code', [
    pytest.param({}, {}, 23, 5, 200, id='all'),
    pytest.param({'entry_id': 'id_01'}, {}, 1, 5, 200, id='all'),
    pytest.param({program_name: 'DOESNOTEXIST'}, {}, 0, 5, 200, id='empty')
])
@pytest.mark.parametrize('http_method', ['post', 'get'])
def test_entries_rawdir(client, data, query, files, total, files_per_entry, status_code, http_method):
    perform_entries_rawdir_test(
        client, status_code=status_code, query=query, files=files, total=total,
        files_per_entry=files_per_entry, http_method=http_method)


@pytest.mark.parametrize('query, files, total, files_per_entry, status_code', [
    pytest.param({}, {}, 23, 5, 200, id='all'),
    pytest.param({program_name: 'DOESNOTEXIST'}, {}, 0, 5, 200, id='empty'),
    pytest.param({}, {'glob_pattern': '*.json'}, 23, 1, 200, id='glob'),
    pytest.param({}, {'re_pattern': '[a-z]*\\.aux'}, 23, 4, 200, id='re'),
    pytest.param({}, {'re_pattern': 'test_entry_02'}, 1, 5, 200, id='re-filter-entries'),
    pytest.param({}, {'re_pattern': 'test_entry_02/.*\\.json'}, 1, 1, 200, id='re-filter-entries-and-files'),
    pytest.param({}, {'glob_pattern': '*.json', 're_pattern': '.*\\.aux'}, 23, 4, 200, id='re-overwrites-glob'),
    pytest.param({}, {'re_pattern': '**'}, -1, -1, 422, id='bad-re-pattern'),
    pytest.param({}, {'compress': True}, 23, 5, 200, id='compress'),
    pytest.param({}, {'include_files': ['1.aux']}, 23, 1, 200, id='file'),
    pytest.param({}, {'include_files': ['1.aux', '2.aux']}, 23, 2, 200, id='files')
])
@pytest.mark.parametrize('http_method', ['post', 'get'])
def test_entries_raw(client, data, query, files, total, files_per_entry, status_code, http_method):
    perform_entries_raw_test(
        client, status_code=status_code, query=query, files=files, total=total,
        files_per_entry=files_per_entry, http_method=http_method)


@pytest.mark.parametrize('http_method', ['post', 'get'])
@pytest.mark.parametrize('test_method', [
    pytest.param(perform_entries_raw_test, id='raw'),
    pytest.param(perform_entries_archive_download_test, id='archive-download')])
def test_entries_download_max(monkeypatch, client, data, test_method, http_method):
    monkeypatch.setattr('nomad.config.max_entry_download', 20)

    test_method(client, status_code=400, http_method=http_method)


@pytest.mark.parametrize('entry_id, files_per_entry, status_code', [
    pytest.param('id_01', 5, 200, id='id'),
    pytest.param('id_embargo', -1, 404, id='404'),
    pytest.param('doesnotexist', -1, 404, id='404')])
def test_entry_rawdir(client, data, entry_id, files_per_entry, status_code):
    response = client.get('entries/%s/rawdir' % entry_id)
    assert_response(response, status_code)
    if status_code == 200:
        assert_entry_rawdir_response(response.json(), files_per_entry=files_per_entry)


@pytest.mark.parametrize('entry_id, files, files_per_entry, status_code', [
    pytest.param('id_01', {}, 5, 200, id='id'),
    pytest.param('doesnotexist', {}, -1, 404, id='404'),
    pytest.param('id_01', {'glob_pattern': '*.json'}, 1, 200, id='glob'),
    pytest.param('id_01', {'re_pattern': '[a-z]*\\.aux'}, 4, 200, id='re'),
    pytest.param('id_01', {'re_pattern': '**'}, -1, 422, id='bad-re-pattern'),
    pytest.param('id_01', {'compress': True}, 5, 200, id='compress'),
    pytest.param('id_01', {'include_files': ['1.aux']}, 1, 200, id='file'),
    pytest.param('id_01', {'include_files': ['1.aux', '2.aux']}, 2, 200, id='files')
])
def test_entry_raw(client, data, entry_id, files, files_per_entry, status_code):
    response = client.get('entries/%s/raw?%s' % (entry_id, urlencode(files, doseq=True)))
    assert_response(response, status_code)
    if status_code == 200:
        assert_raw_zip_file(
            response, files=files_per_entry + 1, manifest_entries=1,
            compressed=files.get('compress', False))


@pytest.fixture(scope='function')
def example_data_with_compressed_files(elastic_module, raw_files_module, mongo_module, test_user, other_test_user, normalized):
    data = ExampleData(main_author=test_user)

    data.create_upload(
        upload_id='with_compr_published',
        published=True)
    data.create_entry(
        upload_id='with_compr_published',
        entry_id='with_compr_published',
        mainfile='test_content/test_entry/mainfile.json')
    data.create_upload(
        upload_id='with_compr_unpublished',
        published=False)
    data.create_entry(
        upload_id='with_compr_unpublished',
        entry_id='with_compr_unpublished',
        mainfile='test_content/test_entry/mainfile.json')

    data.save()

    append_raw_files(
        'with_compr_published', 'tests/data/api/mainfile.xz',
        'test_content/test_entry/mainfile.xz')
    append_raw_files(
        'with_compr_published', 'tests/data/api/mainfile.gz',
        'test_content/test_entry/mainfile.gz')
    append_raw_files(
        'with_compr_unpublished', 'tests/data/api/mainfile.xz',
        'test_content/test_entry/mainfile.xz')
    append_raw_files(
        'with_compr_unpublished', 'tests/data/api/mainfile.gz',
        'test_content/test_entry/mainfile.gz')

    yield

    data.delete()
    from nomad.search import search
    assert search(query=dict(upload_id='with_compr_published')).pagination.total == 0


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
    pytest.param('with_compr_published', 'mainfile.xz', {'decompress': True}, 200, id='decompress-xz-published'),
    pytest.param('with_compr_published', 'mainfile.gz', {'decompress': True}, 200, id='decompress-gz-published'),
    pytest.param('with_compr_unpublished', 'mainfile.xz', {'decompress': True, 'user': 'test-user'}, 200, id='decompress-xz-unpublished'),
    pytest.param('with_compr_unpublished', 'mainfile.gz', {'decompress': True, 'user': 'test-user'}, 200, id='decompress-gz-unpublished'),
    pytest.param('id_unpublished', 'mainfile.json', {}, 404, id='404-unpublished'),
    pytest.param('id_embargo_1', 'mainfile.json', {}, 404, id='404-embargo-no-user'),
    pytest.param('id_embargo_1', 'mainfile.json', {'user': 'other-test-user'}, 404, id='404-embargo-no-access'),
    pytest.param('id_embargo_1', 'mainfile.json', {'user': 'test-user'}, 200, id='embargo-main_author'),
    pytest.param('id_embargo_w_coauthor_1', 'mainfile.json', {'user': 'other-test-user'}, 200, id='embargo-coauthor'),
    pytest.param('id_embargo_w_reviewer_1', 'mainfile.json', {'user': 'other-test-user'}, 200, id='embargo-reviewer')
])
def test_entry_raw_file(
        client, data, example_data_with_compressed_files, example_mainfile_contents, test_user_auth, other_test_user_auth,
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
        f'entries/{entry_id}/raw/{path}?{urlencode(params, doseq=True)}',
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


@pytest.mark.parametrize('query, files, total, status_code', [
    pytest.param({}, {}, 23, 200, id='all'),
    pytest.param({program_name: 'DOESNOTEXIST'}, {}, -1, 200, id='empty'),
    pytest.param({}, {'compress': True}, 23, 200, id='compress')
])
@pytest.mark.parametrize('http_method', ['post', 'get'])
def test_entries_archive_download(client, data, query, files, total, status_code, http_method):
    perform_entries_archive_download_test(
        client, status_code=status_code, query=query, files=files, total=total,
        http_method=http_method)


@pytest.mark.parametrize('required, status_code', [
    pytest.param('*', 200, id='full'),
    pytest.param({'metadata': '*'}, 200, id='partial'),
    pytest.param({'run': {'system[NOTANINT]': '*'}}, 422, id='bad-required-1'),
    pytest.param({'metadata': {'viewers[NOTANINT]': '*'}}, 422, id='bad-required-2'),
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


@pytest.mark.parametrize('entry_id, status_code', [
    pytest.param('id_01', 200, id='id'),
    pytest.param('id_02', 404, id='404-not-visible'),
    pytest.param('doesnotexist', 404, id='404-does-not-exist')])
def test_entry_archive_download(client, data, entry_id, status_code):
    response = client.get('entries/%s/archive/download' % entry_id)
    assert_response(response, status_code)
    if status_code == 200:
        archive = response.json()
        assert 'metadata' in archive
        assert 'run' in archive


@pytest.mark.parametrize('entry_id, required, status_code', [
    pytest.param('id_01', '*', 200, id='full'),
    pytest.param('id_02', '*', 404, id='404'),
    pytest.param('id_01', {'metadata': '*'}, 200, id='partial'),
    pytest.param('id_01', {'run': {'system[NOTANINT]': '*'}}, 422, id='bad-required-1'),
    pytest.param('id_01', {'metadata': {'viewers[NOTANINT]': '*'}}, 422, id='bad-required-2'),
    pytest.param('id_01', {'DOESNOTEXIST': '*'}, 422, id='bad-required-3'),
    pytest.param('id_01', {'resolve-inplace': 'NotBool', 'workflow': '*'}, 422, id='bad-required-4'),
    pytest.param('id_01', {'resolve-inplace': True, 'metadata': 'include-resolved'}, 200, id='resolve-inplace')
])
def test_entry_archive_query(client, data, entry_id, required, status_code):
    response = client.post('entries/%s/archive/query' % entry_id, json={
        'required': required
    })
    assert_response(response, status_code)
    if status_code == 200:
        assert_archive_response(response.json(), required=required)


elements = 'results.material.elements'
n_elements = 'results.material.n_elements'


@pytest.mark.parametrize('query, status_code, total', post_query_test_parameters(
    'entry_id', total=23, material_prefix='results.material.', entry_prefix='') + [
    pytest.param({'pid': '123'}, 200, 1, id='number-valued-string'),
    pytest.param({'optimade_filter': 'nelements = 2'}, 200, 23, id='optimade-filter-positive'),
    pytest.param({'optimade_filter': 'nelements < 2'}, 200, 0, id='optimade-filter-negative'),
    pytest.param({'optimade_filter': '#broken syntax'}, 422, 0, id='optimade-filter-broken-syntax'),
    pytest.param({'optimade_filter': 'doesnotexist = 1'}, 422, 0, id='optimade-filter-broken-semantics')
])
@pytest.mark.parametrize('test_method', [
    pytest.param(perform_entries_metadata_test, id='metadata'),
    pytest.param(perform_entries_raw_test, id='raw'),
    pytest.param(perform_entries_rawdir_test, id='rawdir'),
    pytest.param(perform_entries_archive_test, id='archive'),
    pytest.param(perform_entries_archive_download_test, id='archive-download')])
def test_entries_post_query(client, data, query, status_code, total, test_method):
    response_json = test_method(client, query=query, status_code=status_code, total=total, http_method='post')

    response = client.post('entries/query', json={'query': query})
    response_json = assert_metadata_response(response, status_code=status_code)

    if response_json is None:
        return

    if 'pagination' not in response_json:
        return

    pagination = response_json['pagination']
    assert pagination['total'] == total
    assert pagination['page_size'] == 10
    assert pagination['order_by'] == 'entry_id'
    assert pagination['order'] == 'asc'
    assert ('next_page_after_value' in pagination) == (total > 10)


@pytest.mark.parametrize('query, status_code, total', get_query_test_parameters(
    'entry_id', total=23, material_prefix='results.material.', entry_prefix='') + [
        pytest.param({'q': 'domain__dft'}, 200, 23, id='enum')])
@pytest.mark.parametrize('test_method', [
    pytest.param(perform_entries_metadata_test, id='metadata'),
    pytest.param(perform_entries_raw_test, id='raw'),
    pytest.param(perform_entries_rawdir_test, id='rawdir'),
    pytest.param(perform_entries_archive_test, id='archive'),
    pytest.param(perform_entries_archive_download_test, id='archive-download')])
def test_entries_get_query(client, data, query, status_code, total, test_method):
    response_json = test_method(
        client, query=query, status_code=status_code, total=total, http_method='get')

    if response_json is None:
        return

    if 'pagination' not in response_json:
        return

    response = client.get('entries?%s' % urlencode(query, doseq=True))

    response_json = assert_metadata_response(response, status_code=status_code)

    if response_json is None:
        return

    pagination = response_json['pagination']
    assert pagination['total'] == total
    assert pagination['page_size'] == 10
    assert pagination['order_by'] == 'entry_id'
    assert pagination['order'] == 'asc'
    assert ('next_page_after_value' in pagination) == (total > 10)


@pytest.mark.parametrize(
    'owner, user, status_code, total_entries, total_mainfiles, total_materials',
    owner_test_parameters())
@pytest.mark.parametrize('http_method', ['post', 'get'])
@pytest.mark.parametrize('test_method', [
    pytest.param(perform_entries_metadata_test, id='metadata'),
    pytest.param(perform_entries_raw_test, id='raw'),
    pytest.param(perform_entries_rawdir_test, id='rawdir'),
    pytest.param(perform_entries_archive_test, id='archive'),
    pytest.param(perform_entries_archive_download_test, id='archive-download')])
def test_entries_owner(
        client, data, test_user_auth, other_test_user_auth, admin_user_auth,
        owner, user, status_code, total_entries, total_mainfiles, total_materials,
        http_method, test_method):

    total = (total_entries, total_mainfiles) if test_method == perform_entries_raw_test else total_entries
    perform_owner_test(
        client, test_user_auth, other_test_user_auth, admin_user_auth,
        owner, user, status_code, total, http_method, test_method)


@pytest.mark.parametrize('pagination, response_pagination, status_code', pagination_test_parameters(
    elements='results.material.elements',
    n_elements='results.material.n_elements',
    crystal_system='results.material.symmetry.crystal_system',
    total=23))
@pytest.mark.parametrize('http_method', ['post', 'get'])
@pytest.mark.parametrize('test_method', [
    pytest.param(perform_entries_metadata_test, id='metadata'),
    pytest.param(perform_entries_rawdir_test, id='rawdir'),
    pytest.param(perform_entries_archive_test, id='archive')])
def test_entries_pagination(client, data, pagination, response_pagination, status_code, http_method, test_method):
    response_json = test_method(
        client, pagination=pagination, status_code=status_code, http_method=http_method)

    if response_json is None:
        return

    assert_pagination(pagination, response_json['pagination'], response_json['data'], is_get=(http_method == 'get'))

    if response_pagination is None:
        return
    for key in response_pagination:
        if response_pagination[key] is None:
            assert key not in response_json['pagination']
        else:
            assert response_json['pagination'][key] == response_pagination[key]
    if len(response_json['data']) > 0 and 'order_by' not in pagination:
        if response_pagination['next_page_after_value'] is not None:
            assert response_json['data'][-1]['entry_id'] == response_pagination['next_page_after_value']
