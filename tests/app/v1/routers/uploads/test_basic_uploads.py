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

import asyncio
import io
import os
import tempfile
import time
import zipfile
from collections.abc import Iterable
from datetime import datetime, timezone
from typing import Any, Literal

import pytest
import requests
from fastapi.testclient import TestClient

from nomad import files, infrastructure, processing
from nomad.bundles import BundleExporter
from nomad.common import now
from nomad.config import config
from nomad.config.models.plugins import ExampleUploadEntryPoint
from nomad.datamodel import EntryMetadata
from nomad.files import PublicUploadFiles, StagingUploadFiles, UploadFiles
from nomad.processing import Entry, ProcessStatus, Upload
from nomad.utils.exampledata import ExampleData
from tests.app.v1.routers.common import assert_response, perform_get
from tests.config.models.test_plugins import (
    mock_example_upload_entry_point,
    mock_plugin_package,
)
from tests.fixtures.infrastructure import DataciteMock, TemporalWorkerContext
from tests.processing import test_data as test_processing
from tests.processing.test_edit_metadata import (
    all_admin_metadata,
    all_coauthor_metadata,
    assert_metadata_edited,
)
from tests.test_files import (
    assert_upload_files,
    empty_file,
    example_file_aux,
    example_file_corrupt_zip,
    example_file_mainfile_different_atoms,
    example_file_unparsable,
    example_file_vasp_with_binary,
)
from tests.test_search import assert_search_upload
from tests.utils import assert_doi_name, build_url, set_upload_entry_metadata

from .common import assert_upload

"""
These are the tests for all API operations below ``uploads``. The tests are organized
using the following type of methods: fixtures, ``perform_*``, ``assert_*``, and
``test_*``. While some ``test_*`` methods test individual API operations, some
test methods will test multiple API operations that use common aspects like
supporting queries, pagination, or the owner parameter. The test methods will use
``perform_*`` methods as a parameter. Similarly, the ``assert_*`` methods allow
to assert for certain aspects in the responses.
"""


def perform_post_put_file(
    client,
    action,
    url,
    mode,
    file_paths,
    user_auth=None,
    token=None,
    accept='application/json',
    **query_args,
):
    """Posts or puts a file."""
    if isinstance(file_paths, str):
        file_paths = [file_paths]
    headers = {'Accept': accept}
    if user_auth:
        headers.update(user_auth)
    if mode == 'local_path':
        assert len(file_paths) == 1
        query_args.update(local_path=file_paths[0])
    if token:
        headers['Upload-Token'] = token
    url = build_url(url, query_args)

    if action == 'POST':
        func = client.post
    elif action == 'PUT':
        func = client.put
    else:
        assert False, f'Invalid action provided: {action}'

    if not file_paths:
        response = func(url, data='', headers=headers)
    else:
        if mode == 'multipart':
            if len(file_paths) == 1:
                with open(file_paths[0], 'rb') as f:
                    response = func(url, files={'file': f}, headers=headers)
            else:
                files_list = []
                open_files = []
                try:
                    for file_path in file_paths:
                        filename = os.path.basename(file_path)
                        f = open(file_path, 'rb')
                        open_files.append(f)
                        files_list.append(('file', (filename, f)))
                    response = func(url, files=files_list, headers=headers)
                finally:
                    for f in open_files:
                        f.close()
        elif mode == 'stream':
            assert len(file_paths) == 1
            with open(file_paths[0], 'rb') as f:
                response = func(url, data=f.read(), headers=headers)
        elif mode == 'local_path':
            response = func(url, headers=headers)
        else:
            assert False, f'Invalid value for mode provided {mode}'

    return response


def perform_post_upload_action(
    client, user_auth, upload_id, action, json=None, **query_args
):
    return client.post(
        build_url(f'uploads/{upload_id}/action/{action}', query_args),
        headers=user_auth,
        json=json,
    )


def assert_file_upload_and_processing(
    auth_headers,
    upload_tokens,
    client,
    action,
    url,
    mode,
    user,
    upload_id,
    source_paths,
    target_path,
    query_args,
    accept_json,
    use_upload_token,
    expected_status_code,
    expected_process_status,
    expected_mainfiles,
    published,
    all_entries_should_succeed,
):
    """
    Uploads a file, using the given action (POST or PUT), url, query arguments, and checks
    the results.
    """
    source_paths = source_paths or []
    if isinstance(source_paths, str):
        source_paths = [source_paths]
    user_auth = auth_headers[user]
    # Use either token or bearer token for the post operation (never both)
    user_auth_action = user_auth
    if use_upload_token:
        token = upload_tokens[user]
        user_auth_action = None
    else:
        token = None
    accept = 'application/json' if accept_json else '*'
    processed_response_data = None
    response = perform_post_put_file(
        client,
        action,
        url,
        mode,
        source_paths,
        user_auth_action,
        token,
        accept,
        **query_args,
    )

    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        if accept_json:
            response_json = response.json()
            upload_id = response_json['upload_id']
            if expected_process_status:
                assert (
                    response_json['data']['process_status'] == expected_process_status
                )
            assert_upload(response_json)
        else:
            assert 'Thanks for uploading' in response.text
            if not upload_id:
                return None, None

        if example_file_corrupt_zip in source_paths:
            processed_response_data = assert_processing_fails(
                client, upload_id, user_auth
            )
        else:
            processed_response_data = assert_processing(
                client,
                upload_id,
                user_auth,
                published=published,
                all_entries_should_succeed=all_entries_should_succeed,
            )

            # Check that files got copied as expected
            if query_args.get('auto_decompress', True):
                for source_path in source_paths:
                    upload_files = files.UploadFiles.get(upload_id)
                    file_name = os.path.basename(source_path)
                    if zipfile.is_zipfile(source_path):
                        with open(source_path, 'rb') as f:
                            zf = zipfile.ZipFile(f)
                            for path in zf.namelist():
                                if not path.endswith('/'):
                                    target_path_full = os.path.join(target_path, path)
                                    assert upload_files.raw_exists(target_path_full)
                                    assert upload_files.raw_isfile(target_path_full)
                    elif os.path.isdir(source_path):
                        for root, _, filepaths in os.walk(source_path):
                            for filepath in filepaths:
                                rel_dir = os.path.relpath(root, source_path)
                                path = (
                                    filepath
                                    if rel_dir == '.'
                                    else os.path.join(rel_dir, filepath)
                                )
                                target_path_full = os.path.join(target_path, path)
                                assert upload_files.raw_exists(target_path_full)
                                assert upload_files.raw_isfile(target_path_full)
                    else:
                        if mode == 'stream':
                            # Must specify file_name
                            file_name = query_args['file_name']
                        target_path_full = os.path.join(target_path, file_name)
                        assert upload_files.raw_exists(target_path_full)
                        assert upload_files.raw_isfile(target_path_full)
                        assert (
                            upload_files.raw_file_size(target_path_full)
                            == os.stat(source_path).st_size
                        )
            else:
                upload_files = files.UploadFiles.get(upload_id)
                file_name = os.path.basename(source_paths[0])
                target_path_full = os.path.join(target_path, file_name)
                assert upload_files.raw_exists(target_path_full)
                assert upload_files.raw_isfile(target_path_full)

        assert_expected_mainfiles(upload_id, expected_mainfiles)
    return response, processed_response_data


def assert_expected_mainfiles(upload_id, expected_mainfiles):
    if expected_mainfiles is not None:
        entries = [e.mainfile for e in Entry.objects(upload_id=upload_id)]
        assert set(entries) == set(expected_mainfiles), 'Wrong entries found'
        for entry in Entry.objects(upload_id=upload_id):
            if (
                not isinstance(expected_mainfiles, dict)
                or expected_mainfiles[entry.mainfile]
            ):
                assert entry.process_status == ProcessStatus.SUCCESS
            else:
                assert entry.process_status == ProcessStatus.FAILURE


def assert_upload_does_not_exist(client, upload_id: str, user_auth):
    block_until_completed(client, upload_id, user_auth)

    response = perform_get(client, 'uploads/{upload_id}', user_auth)
    assert_response(response, 404)

    assert Upload.objects(upload_id=upload_id).first() is None
    assert Entry.objects(upload_id=upload_id).count() is 0

    mongo_db = infrastructure.mongo_client[config.mongo.db_name]
    mongo_collection = mongo_db['archive']
    assert mongo_collection.count_documents({}) == 0

    upload_files = UploadFiles.get(upload_id)
    assert upload_files is None or isinstance(upload_files, PublicUploadFiles)


def assert_processing(
    client,
    upload_id,
    user_auth,
    check_search=True,
    check_files=True,
    published=False,
    all_entries_should_succeed=True,
):
    response_data = block_until_completed(client, upload_id, user_auth)

    assert response_data['process_status'] in (
        ProcessStatus.SUCCESS,
        ProcessStatus.READY,
    )
    assert not response_data['process_running']

    response_entries = perform_get(client, f'uploads/{upload_id}/entries', user_auth)
    assert_response(response_entries, 200)
    response_entries_json = response_entries.json()
    response_entries_data = response_entries_json['data']
    all_entries_succesful = True

    for entry in response_entries_data:
        entry_succeeded = entry['process_status'] == ProcessStatus.SUCCESS
        if not entry_succeeded:
            all_entries_succesful = False
            if all_entries_should_succeed:
                assert False, 'One or more entries failed to process'
        pagination = response_entries_json['pagination']
        assert pagination['total'] < pagination['page_size']

    entries = get_upload_entries_metadata(response_entries_data)
    if check_files:
        expected_file_class = (
            files.PublicUploadFiles if published else files.StagingUploadFiles
        )
        assert_upload_files(upload_id, entries, expected_file_class)
    if check_search and all_entries_succesful:
        assert_search_upload(
            entries,
            additional_keys=[
                'results.material.elements',
                'results.method.simulation.program_name',
            ],
            upload_id=upload_id,
        )
    return response_data


def assert_processing_fails(client, upload_id, user_auth):
    response_data = block_until_completed(client, upload_id, user_auth)

    assert response_data['process_status'] == ProcessStatus.FAILURE
    return response_data


def assert_gets_published(
    client,
    upload_id,
    user_auth,
    from_oasis=False,
    current_embargo_length=0,
    **query_args,
):
    embargo_length = query_args.get('embargo_length', current_embargo_length)

    block_until_completed(client, upload_id, user_auth)

    upload_proc = Upload.objects(upload_id=upload_id).first()
    assert upload_proc is not None
    assert upload_proc.published is True
    assert upload_proc.from_oasis == from_oasis
    assert upload_proc.embargo_length == embargo_length

    with upload_proc.entries_metadata() as entries:
        for entry in entries:
            assert entry.with_embargo == (embargo_length > 0)

    assert_upload_files(upload_id, entries, files.PublicUploadFiles, published=True)


def assert_pagination(pagination, expected_pagination):
    """Checks that the contents of `paginaion` matches what is expected."""
    for key, value in expected_pagination.items():
        if value is None:
            assert key not in pagination, (
                f'No value expected for {key}, got {pagination[key]}'
            )
        elif value is Any:
            assert pagination.get(key) is not None, (
                f'Value expected for {key}, got None'
            )
        else:
            assert pagination.get(key) == value, (
                f'For {key} we expecte {value}, but got {pagination.get(key)}'
            )


def block_until_completed(client, upload_id: str, user_auth):
    """Blocks until the processing of the given upload is finished."""
    start_time = time.time()
    while time.time() - start_time < config.tests.default_timeout:
        time.sleep(0.1)
        response = client.get(f'uploads/{upload_id}', headers=user_auth)
        if response.status_code == 200:
            response_json = response.json()
            assert_upload(response_json)
            response_data = response_json['data']
            if not response_data['process_running']:
                return response_data
        elif response.status_code == 404:
            return None
        else:
            raise Exception(
                f'unexpected status code while blocking for upload processing: {str(response.status_code)}'
            )
    raise Exception('Timed out while waiting for upload processing to finish')


def get_upload_entries_metadata(
    entries: list[dict[str, Any]],
) -> Iterable[EntryMetadata]:
    """
    Create a iterable of :class:`EntryMetadata` from a API upload json record, plus a
    with_embargo flag fetched from mongodb.
    """
    return [
        EntryMetadata(
            domain='dft',
            entry_id=entry['entry_id'],
            mainfile=entry['mainfile'],
            with_embargo=Upload.get(entry['upload_id']).with_embargo,
        )
        for entry in entries
    ]


@pytest.mark.parametrize(
    'user, upload_id_key, expected_status_code',
    [
        pytest.param('user1', 'id_unpublished_w', 200, id='valid-upload_id'),
        pytest.param('user1', 'id_published_doi', 200, id='published-with-doi'),
        pytest.param('user1', 'silly_value', 404, id='invalid-upload_id'),
        pytest.param(None, 'id_unpublished_w', 401, id='no-credentials'),
        pytest.param('invalid', 'id_unpublished_w', 401, id='invalid-credentials'),
        pytest.param('user2', 'id_unpublished_w', 403, id='no-access'),
        pytest.param('user0', 'id_unpublished_w', 200, id='admin-access'),
    ],
)
def test_get_upload(
    auth_headers,
    client,
    example_data_writeable,
    example_data_published_doi,
    user,
    upload_id_key,
    expected_status_code,
):
    """Tests the endpoint for getting an upload by upload_id."""
    if upload_id_key in example_data_writeable:
        upload_id = example_data_writeable[upload_id_key]
    else:
        upload_id = upload_id_key
    response = perform_get(client, f'uploads/{upload_id}', auth_headers[user])
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        assert_upload(response.json())


@pytest.mark.parametrize(
    'mode, user, upload_id, source_paths, target_path, query_args, accept_json, use_upload_token, expected_status_code, expected_mainfiles',
    [
        pytest.param(
            'stream',
            None,
            None,
            example_file_aux,
            '',
            {'file_name': 'blah.aux'},
            True,
            False,
            401,
            None,
            id='no-credentials',
        ),
        pytest.param(
            'stream',
            None,
            None,
            example_file_aux,
            '',
            {'file_name': 1},
            True,
            False,
            401,
            None,
            id='filename-not-str',
        ),
        pytest.param(
            'stream',
            'invalid',
            None,
            example_file_aux,
            '',
            {'file_name': 'blah.aux'},
            True,
            False,
            401,
            None,
            id='invalid-credentials',
        ),
        pytest.param(
            'stream',
            'invalid',
            None,
            example_file_aux,
            '',
            {'file_name': 'blah.aux'},
            True,
            True,
            401,
            None,
            id='invalid-credentials-token',
        ),
        pytest.param(
            'multipart',
            'user0',
            'id_published_w',
            example_file_aux,
            '',
            {},
            True,
            False,
            400,
            None,
            id='published',
        ),
        pytest.param(
            'multipart',
            'user0',
            'id_processing_w',
            example_file_aux,
            '',
            {},
            True,
            False,
            400,
            None,
            id='processing',
        ),
        pytest.param(
            'multipart',
            'user2',
            'silly_value',
            example_file_aux,
            '',
            {},
            True,
            False,
            404,
            None,
            id='bad-upload_id',
        ),
        pytest.param(
            'multipart',
            'user2',
            None,
            example_file_aux,
            '',
            {},
            True,
            False,
            403,
            None,
            id='no-access-to-upload',
        ),
        pytest.param(
            'multipart',
            'user1',
            None,
            None,
            '',
            {},
            True,
            False,
            400,
            None,
            id='no-file',
        ),
        pytest.param(
            'local_path',
            'user1',
            None,
            example_file_aux,
            '',
            {},
            True,
            False,
            403,
            None,
            id='local_path-not-admin',
        ),
        pytest.param(
            'stream',
            'user1',
            None,
            example_file_aux,
            '',
            {},
            True,
            False,
            400,
            None,
            id='stream-no-file_name',
        ),
        pytest.param(
            'stream',
            'user1',
            'id_unpublished_w',
            example_file_aux,
            'test_content/test_embargo_entry',
            {'file_name': 'mainfile.json', 'overwrite_if_exists': False},
            True,
            False,
            409,
            None,
            id='cannot-overwrite-existing',
        ),
        pytest.param(
            'stream',
            'user1',
            'id_unpublished_w',
            None,
            'test_content/test_embargo_entry',
            {
                'file_name': '2.aux',
                'copy_or_move_source_path': 'test_content/test_embargo_entry/1.aux',
                'copy_or_move': 'copy',
            },
            True,
            False,
            409,
            None,
            id='copy-file-to-rawdir-already-exists',
        ),
        pytest.param(
            'multipart',
            'user1',
            None,
            example_file_aux,
            '',
            {},
            True,
            False,
            200,
            ['examples_template/template.json'],
            id='multipart',
        ),
        pytest.param(
            'stream',
            'user1',
            None,
            None,
            '',
            {
                'file_name': 'template.json',
                'copy_or_move_source_path': 'examples_template/template.json',
                'copy_or_move': 'copy',
            },
            True,
            False,
            200,
            {'template.json': True, 'examples_template/template.json': True},
            id='copy-file-to-rawdir',
        ),
        pytest.param(
            'stream',
            'user1',
            None,
            None,
            '',
            {
                'file_name': 'template_2.json',
                'copy_or_move_source_path': 'examples_template/template.json',
                'copy_or_move': 'copy',
            },
            True,
            False,
            200,
            {'examples_template/template.json': True},
            id='copy-with-rename-file-to-rawdir',
        ),
        pytest.param(
            'stream',
            'user1',
            None,
            None,
            '',
            {
                'file_name': 'template.json',
                'copy_or_move_source_path': 'examples_template/template.json',
                'copy_or_move': 'move',
            },
            True,
            False,
            200,
            {'template.json': True},
            id='move-file-to-rawdir',
        ),
        pytest.param(
            'stream',
            'user1',
            None,
            None,
            '',
            {
                'file_name': 'template_2.json',
                'copy_or_move_source_path': 'examples_template/template.json',
                'copy_or_move': 'move',
            },
            True,
            False,
            200,
            None,
            id='move-with-rename-file-to-rawdir',
        ),
        pytest.param(
            'stream',
            'user1',
            None,
            example_file_aux,
            '',
            {'file_name': 'blah.aux'},
            True,
            False,
            200,
            ['examples_template/template.json'],
            id='stream',
        ),
        pytest.param(
            'local_path',
            'user0',
            None,
            example_file_aux,
            '',
            {},
            True,
            False,
            200,
            ['examples_template/template.json'],
            id='local_path',
        ),
        pytest.param(
            'stream',
            'user1',
            None,
            example_file_aux,
            '',
            {'file_name': 'blah.aux'},
            True,
            True,
            200,
            ['examples_template/template.json'],
            id='token-auth',
        ),
        pytest.param(
            'stream',
            'user1',
            None,
            example_file_aux,
            'dir1/dir2/dir3',
            {'file_name': 'blah.aux'},
            True,
            False,
            200,
            ['examples_template/template.json'],
            id='file-to-subfolder',
        ),
        pytest.param(
            'stream',
            'user1',
            None,
            example_file_vasp_with_binary,
            'dir1/dir2',
            {'file_name': 'tmp.zip'},
            True,
            False,
            200,
            [
                'examples_template/template.json',
                'dir1/dir2/examples_vasp/xml/Si.xml',
                'dir1/dir2/examples_vasp/xml/perovskite.xml.gz',
            ],
            id='zip-to-subfolder',
        ),
        pytest.param(
            'stream',
            'user1',
            None,
            example_file_aux,
            'examples_template',
            {'file_name': 'template.json'},
            True,
            False,
            200,
            {'examples_template/template.json': False},
            id='overwrite-and-destroy-old-mainfile',
        ),
        pytest.param(
            'stream',
            'user1',
            None,
            example_file_vasp_with_binary,
            '',
            {'file_name': 'tmp.zip'},
            True,
            False,
            200,
            [
                'examples_template/template.json',
                'examples_vasp/xml/Si.xml',
                'examples_vasp/xml/perovskite.xml.gz',
            ],
            id='unzip-and-add-new-mainfiles',
        ),
        pytest.param(
            'stream',
            'user1',
            None,
            example_file_corrupt_zip,
            '',
            {'file_name': 'tmp.zip'},
            True,
            False,
            400,
            None,
            id='bad-zip',
        ),
        pytest.param(
            'multipart',
            'user1',
            None,
            example_file_aux,
            'examples_template',
            {'wait_for_processing': True},
            True,
            False,
            200,
            ['examples_template/template.json'],
            id='wait_for_processing-auxfile-add',
        ),
        pytest.param(
            'multipart',
            'user1',
            None,
            example_file_mainfile_different_atoms,
            'dir1/dir2',
            {'wait_for_processing': True},
            True,
            False,
            200,
            ['examples_template/template.json', 'dir1/dir2/template.json'],
            id='wait_for_processing-mainfile-add',
        ),
        pytest.param(
            'multipart',
            'user1',
            None,
            example_file_mainfile_different_atoms,
            'dir1/dir2',
            {'wait_for_processing': True, 'include_archive': True},
            True,
            False,
            200,
            ['examples_template/template.json', 'dir1/dir2/template.json'],
            id='wait_for_processing-mainfile-add-include_archive',
        ),
        pytest.param(
            'multipart',
            'user1',
            None,
            example_file_mainfile_different_atoms,
            'examples_template',
            {'wait_for_processing': True},
            True,
            False,
            200,
            ['examples_template/template.json'],
            id='wait_for_processing-mainfile-overwrite',
        ),
        pytest.param(
            'multipart',
            'user1',
            None,
            example_file_unparsable,
            'examples_template',
            {'wait_for_processing': True},
            True,
            False,
            200,
            {'examples_template/template.json': False},
            id='wait_for_processing-mainfile-overwrite-destroy',
        ),
        pytest.param(
            'multipart',
            'user1',
            None,
            example_file_vasp_with_binary,
            'examples_template',
            {'wait_for_processing': True},
            True,
            False,
            400,
            None,
            id='wait_for_processing-zipfile',
        ),
        pytest.param(
            'multipart',
            'user1',
            None,
            [example_file_vasp_with_binary, example_file_aux],
            'dir1',
            {'file_name': 'tmp.zip'},
            True,
            False,
            200,
            [
                'examples_template/template.json',
                'dir1/examples_vasp/xml/Si.xml',
                'dir1/examples_vasp/xml/perovskite.xml.gz',
            ],
            id='upload-multiple-vasp-and-aux',
        ),
        pytest.param(
            'multipart',
            'user1',
            None,
            [example_file_aux, example_file_corrupt_zip],
            'dir1',
            {'file_name': 'tmp.zip'},
            True,
            False,
            400,
            ['examples_template/template.json'],
            id='upload-multiple-one-corrupted-zip',
        ),
        pytest.param(
            'stream',
            'user1',
            None,
            empty_file,
            '',
            {'file_name': 'empty.zip', 'auto_decompress': False},
            True,
            False,
            200,
            None,
            id='disable-default-decompression',
        ),
    ],
)
@pytest.mark.asyncio
async def test_put_upload_raw_path(
    auth_headers,
    upload_tokens,
    client,
    elastic_function,
    temporal_worker,
    non_empty_processed_with_temporal,
    example_data_writeable,
    mode,
    user,
    upload_id,
    source_paths,
    target_path,
    query_args,
    accept_json,
    use_upload_token,
    expected_status_code,
    expected_mainfiles,
):
    if upload_id is None:
        upload_id = non_empty_processed_with_temporal.upload_id
    elif example_data_upload_id := example_data_writeable.get(upload_id):
        upload_id = example_data_upload_id

    action = 'PUT'
    url = f'uploads/{upload_id}/raw/{target_path}'
    published = False
    all_entries_should_succeed = not (
        isinstance(expected_mainfiles, dict) and False in expected_mainfiles.values()
    )
    expected_process_status = (
        ProcessStatus.SUCCESS if 'wait_for_processing' in query_args else None
    )

    async with temporal_worker():
        response, _ = await asyncio.to_thread(
            lambda: assert_file_upload_and_processing(
                auth_headers,
                upload_tokens,
                client,
                action,
                url,
                mode,
                user,
                upload_id,
                source_paths,
                target_path,
                query_args,
                accept_json,
                use_upload_token,
                expected_status_code,
                expected_process_status,
                expected_mainfiles,
                published,
                all_entries_should_succeed,
            )
        )

    if response.status_code == 200 and accept_json:
        response_json = response.json()
        processing = response_json['processing']
        if 'wait_for_processing' in query_args:
            assert processing
            assert processing['upload_id'] == upload_id
            assert processing['path'] == os.path.join(
                target_path, os.path.basename(source_paths)
            )
            if source_paths == example_file_aux:
                # Not a mainfile
                for k in ('entry_id', 'parser_name', 'entry', 'archive'):
                    assert processing[k] is None
            else:
                # Mainfile was added
                if source_paths == example_file_unparsable:
                    expected_entry_process_status = ProcessStatus.FAILURE
                else:
                    expected_entry_process_status = ProcessStatus.SUCCESS
                assert processing['entry_id'] is not None
                assert processing['parser_name'] is not None
                assert (
                    processing['entry']['process_status']
                    == expected_entry_process_status
                )
                assert (processing['archive'] is None) == (
                    not query_args.get('include_archive')
                )
        else:
            assert not processing


@pytest.mark.parametrize(
    'user, upload_id, path, expected_status_code',
    [
        pytest.param(
            'user1', 'id_published_w', 'test_content/newdir', 400, id='published'
        ),
        pytest.param(
            None, 'id_unpublished_w', 'test_content/newdir', 401, id='no-credentials'
        ),
        pytest.param(
            'user2',
            'id_unpublished_w',
            'test_content/newdir',
            403,
            id='no-access',
        ),
        pytest.param(
            'user0',
            'id_unpublished_w',
            'test_content/newdir',
            200,
            id='admin-access',
        ),
        pytest.param(
            'user1',
            'id_unpublished_w',
            'test_content/test_embargo_entry/newdir',
            200,
            id='ok',
        ),
        pytest.param(
            'user1',
            'id_unpublished_w',
            'test_content/chars?! "\'@#$%&\\()[]{}=+`´^~*,.;:|<>',
            200,
            id='special-chars',
        ),
        pytest.param(
            'user1',
            'id_unpublished_w',
            'test_content/test_embargo_entry/mainfile.json/newdir',
            400,
            id='bad-path',
        ),
    ],
)
@pytest.mark.asyncio
async def test_post_upload_raw_create_dir_path(
    auth_headers,
    client,
    temporal_worker,
    example_data_writeable,
    user,
    upload_id,
    path,
    expected_status_code,
):
    if upload_id in example_data_writeable:
        upload_id = example_data_writeable[upload_id]
    url = f'uploads/{upload_id}/raw-create-dir/{requests.utils.quote(path)}'
    async with temporal_worker():
        response = await asyncio.to_thread(
            lambda: client.post(url, headers=auth_headers[user])
        )
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        upload = Upload.get(upload_id)
        assert upload.upload_files.raw_exists(path)
        assert not upload.upload_files.raw_isfile(path)


@pytest.mark.parametrize(
    'user, upload_id, path, use_upload_token, expected_status_code, expected_mainfiles',
    [
        pytest.param(
            'user1',
            None,
            'examples_template/1.aux',
            False,
            200,
            ['examples_template/template.json'],
            id='delete-aux-file',
        ),
        pytest.param(
            'user1',
            None,
            'examples_template/template.json',
            False,
            200,
            [],
            id='delete-main-file',
        ),
        pytest.param('user1', None, '', False, 200, [], id='delete-root'),
        pytest.param(
            'user1',
            None,
            'examples_template',
            False,
            200,
            [],
            id='delete-subfolder',
        ),
        pytest.param(
            'user1',
            None,
            'examples_template/1.aux',
            True,
            200,
            ['examples_template/template.json'],
            id='delete-token-access',
        ),
        pytest.param(
            'user0',
            None,
            'examples_template/1.aux',
            False,
            200,
            ['examples_template/template.json'],
            id='delete-admin-access',
        ),
        pytest.param(
            'user2',
            None,
            'examples_template/1.aux',
            False,
            403,
            None,
            id='no-access',
        ),
        pytest.param(
            None,
            None,
            'examples_template/1.aux',
            False,
            401,
            None,
            id='no-credentials',
        ),
        pytest.param(
            'invalid',
            None,
            'examples_template/1.aux',
            False,
            401,
            None,
            id='invalid-credentials',
        ),
        pytest.param(
            'invalid',
            None,
            'examples_template/1.aux',
            True,
            401,
            None,
            id='invalid-credentials-token',
        ),
        pytest.param(
            'user1',
            'id_published_w',
            'examples_template/1.aux',
            False,
            400,
            None,
            id='published',
        ),
        pytest.param(
            'user1',
            'id_processing_w',
            'examples_template/1.aux',
            False,
            400,
            None,
            id='processing',
        ),
    ],
)
@pytest.mark.asyncio
async def test_delete_upload_raw_path(
    auth_headers,
    client,
    temporal_worker,
    non_empty_processed_with_temporal,
    example_data_writeable,
    upload_tokens,
    user,
    upload_id,
    path,
    use_upload_token,
    expected_status_code,
    expected_mainfiles,
):
    static_upload_id = upload_id
    if upload_id is None:
        upload_id = non_empty_processed_with_temporal.upload_id
    elif example_data_upload_id := example_data_writeable.get(upload_id):
        upload_id = example_data_upload_id
    user_auth = auth_headers[user]
    # Use either token or bearer token for the post operation (never both)
    if use_upload_token:
        headers = {'Upload-Token': upload_tokens[user]}
    else:
        headers = dict(user_auth or {})

    if static_upload_id == 'id_processing_w':
        # Ensure file exists (otherwise we get 404, which is not what we want to test)
        upload_files = StagingUploadFiles(upload_id)
        upload_files.add_rawfiles(
            'tests/data/proc/examples_template/1.aux', 'examples_template'
        )

    async with temporal_worker():
        response = await asyncio.to_thread(
            lambda: client.delete(
                build_url(f'uploads/{upload_id}/raw/{path}', query_args={}),
                headers=headers,
            )
        )
        assert_response(response, expected_status_code)
        if expected_status_code == 200:
            await asyncio.to_thread(
                lambda: assert_processing(client, upload_id, user_auth)
            )
            # Check that path to remove has disappeared
            upload_files = StagingUploadFiles(upload_id)
            if path == '':
                # Deleting the root folder = the folder itself should be emptied, but not deleted.
                assert not list(upload_files.raw_listdir(''))
            else:
                # Deleting a file or folder within the raw folder - it should disappear.
                assert not upload_files.raw_exists(path)

            assert_expected_mainfiles(upload_id, expected_mainfiles)


@pytest.mark.parametrize(
    'user, upload_id, kwargs',
    [
        pytest.param(
            'user1',
            'id_unpublished_w',
            dict(metadata=all_coauthor_metadata),
            id='edit-all',
        ),
        pytest.param(
            'user1',
            'id_unpublished_w',
            dict(
                metadata=dict(coauthors='unknown'),
                expected_error_loc=('metadata', 'coauthors'),
            ),
            id='edit-coauthor-unknown-fails',
        ),
        pytest.param(
            'user1',
            'id_published_w',
            dict(metadata=dict(embargo_length=0)),
            id='lift-embargo',
        ),
        pytest.param(
            'user0',
            'id_published_w',
            dict(metadata=all_admin_metadata),
            id='protected-admin',
        ),
        pytest.param(
            'user1',
            'id_unpublished_w',
            dict(
                metadata=dict(main_author='lhofstadter'),
                expected_error_loc=('metadata', 'main_author'),
            ),
            id='protected-not-admin',
        ),
        pytest.param(
            'user1',
            'silly_value',
            dict(
                metadata=dict(upload_name='test_name'),
                expected_error_loc=('upload_id',),
            ),
            id='bad-upload_id',
        ),
        pytest.param(
            'user0',
            'id_published_w',
            dict(metadata=dict(upload_name='test_name')),
            id='published-admin',
        ),
        pytest.param(
            'user1',
            'id_published_w',
            dict(metadata=dict(upload_name='test_name')),
            id='published-not-admin',
        ),
        pytest.param(
            None,
            'id_unpublished_w',
            dict(metadata=dict(upload_name='test_name'), expected_status_code=401),
            id='no-credentials',
        ),
        pytest.param(
            'invalid',
            'id_unpublished_w',
            dict(metadata=dict(upload_name='test_name'), expected_status_code=401),
            id='invalid-credentials',
        ),
        pytest.param(
            'user2',
            'id_unpublished_w',
            dict(
                metadata=dict(upload_name='test_name'),
                expected_error_loc=('metadata', 'upload_name'),
            ),
            id='no-access',
        ),
        pytest.param(
            'user2',
            'id_unpublished_w',
            dict(metadata=dict(upload_name='test_name'), add_coauthor=True),
            id='coauthor-access',
        ),
        pytest.param(
            'user1',
            'id_empty_w',
            dict(metadata=dict(upload_name='test_name')),
            id='empty-upload-ok',
        ),
        pytest.param(
            'user1',
            'id_unpublished_w',
            dict(
                query={
                    'and': [
                        {'upload_create_time:gt': '2021-01-01'},
                        {'published': False},
                    ]
                },
                owner='user',
                metadata=dict(comment='a test comment'),
            ),
            id='query-ok',
        ),
        pytest.param(
            'user1',
            'id_unpublished_w',
            dict(
                query={
                    'and': [
                        {'upload_create_time:gt': '2021-01-01'},
                        {'published': False},
                    ]
                },
                owner='user',
                metadata=dict(upload_name='a test name'),
                expected_error_loc=('metadata', 'upload_name'),
            ),
            id='query-cannot-edit-upload-data',
        ),
        pytest.param(
            'user1',
            'id_unpublished_w',
            dict(
                query={'upload_create_time:lt': '2021-01-01'},
                owner='user',
                metadata=dict(comment='a test comment'),
                expected_error_loc=('query',),
            ),
            id='query-no-results',
        ),
    ],
)
@pytest.mark.asyncio
async def test_post_upload_edit(
    auth_headers,
    client,
    temporal_worker,
    example_data_writeable,
    example_datasets,
    users_dict,
    user,
    upload_id,
    kwargs,
):
    """
    Note, since the endpoint basically just forwards the request to
    `MetadataEditRequestHandler.edit_metadata`, we only do very simple verification here,
    the more extensive testnig is done in `tests.processing.test_edit_metadata`.
    """
    if upload_id in example_data_writeable:
        upload_id = example_data_writeable[upload_id]
    user_auth = auth_headers[user]
    user = users_dict.get(user)
    query = kwargs.get('query')
    owner = kwargs.get('owner')
    metadata = kwargs.get('metadata')
    entries = kwargs.get('entries')
    entries_key = kwargs.get('entries_key')
    verify_only = kwargs.get('verify_only', False)
    expected_error_loc = kwargs.get('expected_error_loc')
    expected_status_code = kwargs.get('expected_status_code')
    affected_upload_ids = kwargs.get('affected_upload_ids', [upload_id])
    affected_upload_ids = [
        example_data_writeable.get(uid, uid) for uid in affected_upload_ids
    ]

    expected_metadata = kwargs.get('expected_metadata', metadata)

    add_coauthor = kwargs.get('add_coauthor', False)
    async with temporal_worker() as env:
        if add_coauthor:
            upload = Upload.get(upload_id)
            await asyncio.to_thread(
                lambda: upload.edit_upload_metadata(
                    edit_request_json={'metadata': {'coauthors': user.user_id}},
                    user_id=upload.main_author,
                )
            )

        edit_request_json = dict(
            query=query,
            owner=owner,
            metadata=metadata,
            entries=entries,
            entries_key=entries_key,
            verify_only=verify_only,
        )
        url = f'uploads/{upload_id}/edit'
        edit_start = now().isoformat()[0:22]
        response = await asyncio.to_thread(
            lambda: client.post(url, headers=user_auth, json=edit_request_json)
        )
    if expected_error_loc:
        assert_response(response, 422)
        error_locs = [tuple(d['loc']) for d in response.json()['detail']]
        assert expected_error_loc in error_locs
    elif expected_status_code not in (None, 200):
        assert_response(response, expected_status_code)
    else:
        assert_response(response, 200)
        assert_metadata_edited(user, expected_metadata, affected_upload_ids, edit_start)


@pytest.mark.parametrize(
    'mode, source_paths, query_args, user, use_upload_token, test_limit, accept_json, expected_status_code',
    [
        pytest.param(
            'multipart',
            example_file_vasp_with_binary,
            dict(upload_name='test_name'),
            'user1',
            False,
            False,
            True,
            200,
            id='multipart',
        ),
        pytest.param(
            'multipart',
            example_file_vasp_with_binary,
            dict(),
            'user1',
            False,
            False,
            True,
            200,
            id='multipart-no-name',
        ),
        pytest.param(
            'multipart',
            example_file_vasp_with_binary,
            dict(upload_name='test_name'),
            'user1',
            True,
            False,
            True,
            200,
            id='multipart-token',
        ),
        pytest.param(
            'stream',
            example_file_vasp_with_binary,
            dict(embargo_length=0, upload_name='test_name'),
            'user1',
            False,
            False,
            True,
            200,
            id='stream-no-embargo',
        ),
        pytest.param(
            'stream',
            example_file_vasp_with_binary,
            dict(embargo_length=7),
            'user1',
            False,
            False,
            True,
            200,
            id='stream-no-name-embargoed',
        ),
        pytest.param(
            'stream',
            example_file_vasp_with_binary,
            dict(embargo_length=37),
            'user1',
            False,
            False,
            True,
            400,
            id='stream-invalid-embargo',
        ),
        pytest.param(
            'stream',
            example_file_vasp_with_binary,
            dict(upload_name='test_name'),
            'user1',
            True,
            False,
            True,
            200,
            id='stream-token',
        ),
        pytest.param(
            'local_path',
            example_file_vasp_with_binary,
            dict(),
            'user0',
            False,
            False,
            True,
            200,
            id='local_path_file',
        ),
        pytest.param(
            'local_path',
            'tests/data/proc/example_upload',
            dict(upload_name='test_name'),
            'user0',
            False,
            False,
            True,
            200,
            id='local_path_folder',
        ),
        pytest.param(
            'local_path',
            example_file_vasp_with_binary,
            dict(),
            'user1',
            False,
            False,
            True,
            403,
            id='local_path-not-admin',
        ),
        pytest.param(
            'stream',
            example_file_vasp_with_binary,
            dict(),
            'user1',
            False,
            False,
            False,
            200,
            id='no-accept-json',
        ),
        pytest.param(
            'multipart',
            example_file_vasp_with_binary,
            dict(),
            None,
            False,
            False,
            True,
            401,
            id='no-credentials',
        ),
        pytest.param(
            'multipart',
            example_file_vasp_with_binary,
            dict(),
            'invalid',
            False,
            False,
            True,
            401,
            id='invalid-credentials',
        ),
        pytest.param(
            'multipart',
            example_file_vasp_with_binary,
            dict(),
            'invalid',
            True,
            False,
            True,
            401,
            id='invalid-credentials-token',
        ),
        pytest.param(
            'stream',
            [],
            dict(upload_name='test_name'),
            'user1',
            False,
            False,
            True,
            200,
            id='no-file',
        ),
        pytest.param(
            'stream',
            example_file_aux,
            dict(file_name='1.aux'),
            'user1',
            False,
            False,
            True,
            200,
            id='stream-non-zip-file',
        ),
        pytest.param(
            'stream',
            example_file_aux,
            dict(),
            'user1',
            False,
            False,
            True,
            400,
            id='stream-non-zip-file-no-file_name',
        ),
        pytest.param(
            'stream',
            example_file_vasp_with_binary,
            dict(upload_name='test_name', publish_directly=True),
            'user1',
            False,
            False,
            True,
            200,
            id='publish_directly',
        ),
        pytest.param(
            'stream',
            empty_file,
            dict(upload_name='test_name', publish_directly=True),
            'user1',
            False,
            False,
            True,
            200,
            id='publish_directly-empty',
        ),
        pytest.param(
            'stream',
            example_file_vasp_with_binary,
            dict(upload_name='test_name'),
            'user1',
            False,
            True,
            True,
            400,
            id='upload-limit-exceeded',
        ),
        pytest.param(
            'multipart',
            example_file_corrupt_zip,
            dict(),
            'user1',
            False,
            False,
            True,
            200,
            id='bad-zip',
        ),
        pytest.param(
            'multipart',
            [example_file_aux, example_file_mainfile_different_atoms],
            dict(),
            'user1',
            False,
            False,
            True,
            200,
            id='upload-multiple-files',
        ),
        pytest.param(
            'multipart',
            [example_file_aux, example_file_corrupt_zip],
            dict(),
            'user1',
            False,
            False,
            True,
            200,
            id='upload-multiple-files-one-corrupt',
        ),
        pytest.param(
            None,
            [],
            dict(example_upload_id='test'),
            'user1',
            False,
            False,
            True,
            200,
            id='example-upload',
        ),
        pytest.param(
            'stream',
            empty_file,
            dict(upload_name='test_name', auto_decompress=False, file_name='empty.zip'),
            'user1',
            False,
            False,
            True,
            200,
            id='disable-default-decompression',
        ),
    ],
)
@pytest.mark.asyncio
async def test_post_upload(
    auth_headers,
    upload_tokens,
    client,
    temporal_worker,
    monkeypatch,
    empty_upload,
    non_empty_example_upload,
    mode,
    source_paths,
    query_args,
    user,
    use_upload_token,
    test_limit,
    accept_json,
    expected_status_code,
):
    """
    Posts an upload, with different arguments.
    """
    if isinstance(source_paths, str):
        source_paths = [source_paths]
    if test_limit:
        monkeypatch.setattr('nomad.config.services.upload_limit', 0)

    # Create a mocked example upload + files if testing example uploads
    is_example_upload = query_args.get('example_upload_id')
    if is_example_upload:
        temp_dir = tempfile.TemporaryDirectory()
        package_directory = temp_dir.name
        filepath = os.path.join(package_directory, 'data.txt')
        with open(filepath, 'w'):
            pass
        assert os.path.exists(filepath)
        mock_plugin_package(monkeypatch, package_directory)
        mock_example_upload_entry_point(
            monkeypatch,
            ExampleUploadEntryPoint(
                id='test',
                title='test',
                description='test',
                category='test',
                resources='data.txt',
            ),
        )
    action = 'POST'
    url = 'uploads'
    published = query_args.get('publish_directly') and not source_paths == [empty_file]
    all_entries_should_succeed = True
    target_path = ''
    expected_mainfiles = None
    upload_id = None  # Not determined yet
    expected_process_status = None

    async with temporal_worker():
        _, processed_response_data = await asyncio.to_thread(
            lambda: assert_file_upload_and_processing(
                auth_headers,
                upload_tokens,
                client,
                action,
                url,
                mode,
                user,
                upload_id,
                source_paths,
                target_path,
                query_args,
                accept_json,
                use_upload_token,
                expected_status_code,
                expected_process_status,
                expected_mainfiles,
                published,
                all_entries_should_succeed,
            )
        )

    if is_example_upload:
        temp_dir.cleanup()

    if expected_status_code == 200 and processed_response_data:
        expected_upload_name = query_args.get('upload_name')
        if not expected_upload_name:
            if is_example_upload:
                expected_upload_name = 'test'
            elif mode in ('multipart', 'local_path') and len(source_paths) == 1:
                expected_upload_name = os.path.basename(source_paths[0])
            elif mode == 'stream':
                expected_upload_name = query_args.get('file_name')

        assert processed_response_data.get('upload_name') == expected_upload_name

    if query_args.get('publish_directly'):
        upload_id = processed_response_data['upload_id']
        upload_proc = Upload.objects(upload_id=upload_id).first()
        if source_paths == [empty_file]:
            assert not upload_proc.published
        else:
            assert_gets_published(
                client, upload_id, auth_headers['user1'], **query_args
            )


@pytest.mark.parametrize(
    'upload_id, publish, user, expected_status_code',
    [
        pytest.param(None, True, 'user0', 200, id='published-admin'),
        pytest.param(None, True, 'user1', 403, id='published-not-admin'),
        pytest.param(None, False, 'user1', 200, id='not-published'),
        pytest.param(None, False, None, 401, id='no-credentials'),
        pytest.param(None, False, 'invalid', 401, id='invalid-credentials'),
        pytest.param(None, False, 'user2', 403, id='no-access'),
        pytest.param('id_processing_w', False, 'user1', 400, id='already-processing'),
        pytest.param('silly_value', False, 'user1', 404, id='invalid-upload_id'),
    ],
)
@pytest.mark.asyncio
async def test_post_upload_action_process(
    auth_headers,
    client,
    temporal_worker,
    monkeypatch,
    example_data_writeable,
    non_empty_processed_with_temporal,
    internal_example_user_metadata,
    upload_id,
    publish,
    user,
    expected_status_code,
):
    async with temporal_worker():
        if publish:
            set_upload_entry_metadata(
                non_empty_processed_with_temporal, internal_example_user_metadata
            )
            await asyncio.to_thread(non_empty_processed_with_temporal.publish_upload)

        monkeypatch.setattr('nomad.config.meta.version', 're_process_test_version')
        monkeypatch.setattr('nomad.config.meta.commit', 're_process_test_commit')
        user_auth = auth_headers[user]

        if upload_id is None:
            upload_id = non_empty_processed_with_temporal.upload_id
        elif example_data_upload_id := example_data_writeable.get(upload_id):
            upload_id = example_data_upload_id
        response = await asyncio.to_thread(
            lambda: perform_post_upload_action(client, user_auth, upload_id, 'process')
        )
        assert_response(response, expected_status_code)
        if expected_status_code == 200:
            await asyncio.to_thread(
                lambda: assert_processing(
                    client,
                    upload_id,
                    auth_headers['user1'],
                    check_files=False,
                    published=True,
                )
            )


@pytest.mark.parametrize(
    'upload_id, user, owner, query, include_parent_folders, expected_status_code, expect_exists, expect_not_exists',
    [
        pytest.param(
            'id_published_w',
            'user0',
            None,
            {'entry_id': 'id_published_w_entry'},
            False,
            400,
            [],
            [],
            id='published-admin',
        ),
        pytest.param(
            'id_unpublished_w',
            'user2',
            None,
            {'entry_id': 'id_unpublished_w_entry'},
            False,
            403,
            [],
            [],
            id='unpublished-no-access',
        ),
        pytest.param(
            'id_unpublished_w',
            'user1',
            None,
            None,
            False,
            400,
            [],
            [],
            id='no-query',
        ),
        pytest.param(
            'id_unpublished_w',
            'user1',
            None,
            {'entry_id': ['id_unpublished_w_entry', 'silly']},
            False,
            200,
            ['test_content/test_embargo_entry/1.aux'],
            ['test_content/test_embargo_entry/mainfile.json'],
            id='ok',
        ),
        pytest.param(
            'id_unpublished_w',
            'user1',
            None,
            {'entry_id': 'id_unpublished_w_entry'},
            True,
            200,
            ['test_content'],
            ['test_content/test_embargo_entry'],
            id='ok-delete-folder',
        ),
        pytest.param(
            'id_unpublished_w',
            'user0',
            'admin',
            {'entry_id': 'id_unpublished_w_entry'},
            False,
            200,
            ['test_content/test_embargo_entry/1.aux'],
            ['test_content/test_embargo_entry/mainfile.json'],
            id='ok-admin-access',
        ),
    ],
)
@pytest.mark.asyncio
async def test_post_upload_action_delete_entry_files(
    auth_headers,
    client,
    temporal_worker: TemporalWorkerContext,
    example_data_writeable,
    upload_id,
    user,
    owner,
    query,
    include_parent_folders,
    expected_status_code,
    expect_exists,
    expect_not_exists,
):
    if upload_id in example_data_writeable:
        upload_id = example_data_writeable[upload_id]
    json: dict = {}
    if include_parent_folders is not None:
        json.update(include_parent_folders=include_parent_folders)
    if owner is not None:
        json.update(owner=owner)
    if query is not None:
        json.update(query=query)

    async with temporal_worker() as env:
        response = await asyncio.to_thread(
            lambda: perform_post_upload_action(
                client, auth_headers[user], upload_id, 'delete-entry-files', json=json
            )
        )
        assert_response(response, expected_status_code)
        if expected_status_code == 200:
            upload = Upload.get(upload_id)
            upload.reload()
            # allow enough time to start processing the workflow
            while True:
                if upload.process_status == ProcessStatus.PENDING:
                    await asyncio.to_thread(lambda: time.sleep(1))
                    upload.reload()
                else:
                    break
            if workflow_ids := upload.workflow_ids:
                handle = env.client.get_workflow_handle(workflow_ids[0])
                await handle.result()
            for path in expect_exists or []:
                assert upload.upload_files.raw_exists(path), (
                    f'Missing expected path: {path}'
                )
            for path in expect_not_exists or []:
                assert not upload.upload_files.raw_exists(path), (
                    f'Expected path not to exist: {path}'
                )


@pytest.mark.parametrize(
    'upload_id, user, expected_status_code',
    [
        pytest.param('id_unpublished_w', 'user1', 200, id='delete-own'),
        pytest.param('id_unpublished_w', 'user2', 403, id='delete-others-not-admin'),
        pytest.param('id_unpublished_w', 'user0', 200, id='delete-others-admin'),
        pytest.param('id_published_w', 'user1', 403, id='delete-own-published'),
        pytest.param(
            'id_published_w', 'user0', 200, id='delete-others-published-admin'
        ),
        pytest.param('silly_value', 'user1', 404, id='invalid-upload_id'),
        pytest.param('id_unpublished_w', None, 401, id='no-credentials'),
        pytest.param('id_unpublished_w', 'invalid', 401, id='invalid-credentials'),
    ],
)
@pytest.mark.asyncio
async def test_delete_upload(
    auth_headers,
    client,
    temporal_worker,
    example_data_writeable,
    upload_id,
    user,
    expected_status_code,
):
    """Uploads a file, and then tries to delete it, with different parameters and users."""
    if upload_id in example_data_writeable:
        upload_id = example_data_writeable[upload_id]
    async with temporal_worker():
        # Run blocking call in thread pool
        response = await asyncio.to_thread(
            lambda: client.delete(f'uploads/{upload_id}', headers=auth_headers[user])
        )
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        assert_upload_does_not_exist(client, upload_id, auth_headers['user1'])


@pytest.mark.parametrize(
    'authorized, expected_status_code',
    [pytest.param(True, 200, id='ok'), pytest.param(False, 401, id='not-authorized')],
)
def test_get_command_examples(auth_headers, client, authorized, expected_status_code):
    response = perform_get(
        client,
        'uploads/command-examples',
        user_auth=auth_headers['user1'] if authorized else None,
    )
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        data = response.json()
        for k in (
            'upload_url',
            'upload_command',
            'upload_command_with_name',
            'upload_progress_command',
            'upload_command_form',
            'upload_tar_command',
        ):
            assert k in data
        assert '/api/v1/uploads' in data['upload_command']


@pytest.mark.parametrize(
    'has_write_access,is_published,upload_state,expected_status',
    [
        pytest.param(True, False, ProcessStatus.PENDING, 200, id='success-case'),
        pytest.param(False, False, ProcessStatus.PENDING, 403, id='permission-denied'),
        pytest.param(True, True, ProcessStatus.PENDING, 400, id='published-upload'),
        pytest.param(
            True, False, ProcessStatus.SUCCESS, 400, id='success-state-invalid'
        ),
        pytest.param(
            True, False, ProcessStatus.FAILURE, 400, id='failure-state-invalid'
        ),
    ],
)
def test_stop_processing_action(
    has_write_access,
    is_published,
    upload_state,
    expected_status,
    non_empty_uploaded,
    user1,
    user2,
    auth_headers,
    client,
    temporal_worker,
    monkeypatch,
):
    """Tests the endpoint for stopping the processing of an upload."""
    upload_id, _ = non_empty_uploaded

    # Create upload with appropriate owner based on access test
    upload_owner = user1 if has_write_access else user2
    upload = Upload.create(
        upload_id=upload_id,
        main_author=upload_owner,
        publish_time=datetime.now(timezone.utc) if is_published else None,
        workflow_ids=['example-workflow-id'],
    )
    upload.save()
    upload.process_status = upload_state
    upload.save()

    # Always use user1's auth headers for the request
    user_auth = auth_headers['user1']

    # Mock the stop processing workflow method
    async def mock_stop_processing_workflows(self):
        pass

    monkeypatch.setattr(
        Upload, '_stop_processing_workflows', mock_stop_processing_workflows
    )

    # Perform the request
    response = perform_post_upload_action(
        client, user_auth, upload_id, 'stop-processing'
    )

    assert_response(response, expected_status)

    if expected_status == 200:
        upload.reload()
        assert len(upload.workflow_ids) == 0
        assert upload.process_status == ProcessStatus.READY
        assert upload.last_status_message == 'Processing stopped'


@pytest.mark.parametrize(
    'kwargs',
    [
        pytest.param(dict(expected_status_code=200), id='no-args'),
        pytest.param(
            dict(query_args={'embargo_length': 12}, expected_status_code=200),
            id='non-standard-embargo',
        ),
        pytest.param(
            dict(query_args={'embargo_length': 24}, expected_status_code=200),
            id='non-standard-embargo-length-only',
        ),
        pytest.param(
            dict(query_args={'embargo_length': 100}, expected_status_code=400),
            id='illegal-embargo-length',
        ),
        pytest.param(
            dict(query_args={'embargo_length': 0}, expected_status_code=200),
            id='no-embargo',
        ),
        pytest.param(
            dict(upload_id='id_empty_w', expected_status_code=400), id='empty'
        ),
        pytest.param(
            dict(upload_id='id_processing_w', expected_status_code=400), id='processing'
        ),
        pytest.param(
            dict(upload_id='id_published_w', expected_status_code=400),
            id='already-published',
        ),
        pytest.param(dict(user=None, expected_status_code=401), id='no-credentials'),
        pytest.param(
            dict(user='invalid', expected_status_code=401), id='invalid-credentials'
        ),
        pytest.param(dict(user='user2', expected_status_code=403), id='no-access'),
    ],
)
@pytest.mark.asyncio
async def test_post_upload_action_publish(
    auth_headers, client, temporal_worker, example_data_writeable, kwargs
):
    """Tests the publish action with various arguments."""
    upload_id = kwargs.get('upload_id', 'id_unpublished_w')
    if upload_id in example_data_writeable:
        upload_id = example_data_writeable[upload_id]
    query_args = kwargs.get('query_args', {})
    expected_status_code = kwargs.get('expected_status_code', 200)
    user = kwargs.get('user', 'user1')
    user_auth = auth_headers[user]
    async with temporal_worker():
        response = await asyncio.to_thread(
            lambda: perform_post_upload_action(
                client, user_auth, upload_id, 'publish', **query_args
            )
        )

    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        upload = assert_upload(response.json())
        assert upload['process_running']

        assert_gets_published(
            client, upload_id, user_auth, current_embargo_length=12, **query_args
        )


@pytest.mark.parametrize(
    'upload_id, user, preprocess, expected_status_code',
    [
        pytest.param('id_published_w', 'user1', None, 200, id='ok'),
        pytest.param('id_published_w', 'user2', None, 403, id='no-access'),
        pytest.param('id_published_w', 'user2', 'make-coauthor', 200, id='ok-coauthor'),
        pytest.param('id_published_w', None, None, 401, id='no-credentials'),
        pytest.param('id_published_w', 'invalid', None, 401, id='invalid-credentials'),
        pytest.param('id_unpublished_w', 'user1', None, 400, id='not-published'),
        pytest.param('id_published_w', 'user1', 'lift', 400, id='already-lifted'),
    ],
)
@pytest.mark.asyncio
async def test_post_upload_action_lift_embargo(
    auth_headers,
    client,
    example_data_writeable,
    users_dict,
    upload_id,
    user,
    preprocess,
    expected_status_code,
    temporal_worker,
):
    if upload_id in example_data_writeable:
        upload_id = example_data_writeable[upload_id]
    user_auth = auth_headers[user]
    user = users_dict.get(user)

    async with temporal_worker():
        if preprocess:
            if preprocess == 'lift':
                metadata = {'embargo_length': 0}
            elif preprocess == 'make-coauthor':
                metadata = {'coauthors': user.user_id}
            upload = Upload.get(upload_id)
            await upload._start_edit_upload_metadata_workflow(
                dict(metadata=metadata), config.services.admin_user_id
            )

        response = await asyncio.to_thread(
            lambda: perform_post_upload_action(
                client, user_auth, upload_id, 'lift-embargo'
            )
        )
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        assert_metadata_edited(user, {'embargo_length': 0}, [upload_id])


@pytest.mark.parametrize(
    'upload_id, user, query_args, expected_status_code',
    [
        pytest.param('id_published_w', 'user1', dict(), 200, id='published-owner'),
        pytest.param('id_published_w', 'user0', dict(), 200, id='published-admin'),
        pytest.param('id_published_w', 'user2', dict(), 403, id='published-not-owner'),
        pytest.param(
            'id_published_w',
            'user1',
            dict(include_raw_files=False),
            200,
            id='published-owner-exclude-raw',
        ),
        pytest.param(
            'id_published_w',
            'user1',
            dict(include_archive_files=False),
            200,
            id='published-owner-exclude-archive',
        ),
        pytest.param('id_unpublished_w', 'user1', dict(), 200, id='unpublished-owner'),
        pytest.param('id_unpublished_w', 'user0', dict(), 200, id='unpublished-admin'),
        pytest.param(
            'id_unpublished_w',
            'user2',
            dict(),
            403,
            id='unpublished-not-owner',
        ),
    ],
)
@pytest.mark.asyncio
async def test_get_upload_bundle(
    auth_headers,
    client,
    temporal_worker,
    example_data_writeable,
    upload_id,
    user,
    query_args,
    expected_status_code,
):
    if upload_id in example_data_writeable:
        upload_id = example_data_writeable[upload_id]
    include_raw_files = query_args.get('include_raw_files', True)
    include_archive_files = query_args.get('include_archive_files', True)

    url = build_url(f'uploads/{upload_id}/bundle', query_args)
    response = perform_get(client, url, user_auth=auth_headers[user])
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        with zipfile.ZipFile(io.BytesIO(response.content)) as zip_file:
            upload = Upload.get(upload_id)
            upload_files = upload.upload_files
            expected_files = set(['bundle_info.json'])
            for dirpath, __, filenames in os.walk(upload_files.os_path):
                for filename in filenames:
                    os_path = os.path.join(dirpath, filename)
                    rel_path = os.path.relpath(os_path, upload_files.os_path)
                    include = False
                    include |= (
                        rel_path.startswith('raw')
                        and not rel_path.endswith('.h5')
                        and include_raw_files
                    )
                    include |= rel_path.startswith('archive') and include_archive_files
                    if include:
                        expected_files.add(rel_path)
            assert expected_files == set(zip_file.namelist())


@pytest.mark.parametrize(
    'publish, test_duplicate, user, export_args, query_args, expected_status_code',
    [
        pytest.param(True, False, 'user0', dict(), dict(), 200, id='published-admin'),
        pytest.param(
            False, False, 'user0', dict(), dict(), 200, id='unpublished-admin'
        ),
        pytest.param(True, True, 'user0', dict(), dict(), 400, id='duplicate'),
        pytest.param(True, False, 'user2', dict(), dict(), 200, id='not-oasis-admin'),
        pytest.param(True, False, None, dict(), dict(), 401, id='no-credentials'),
    ],
)
@pytest.mark.asyncio
async def test_post_upload_bundle(
    auth_headers,
    client,
    temporal_worker,
    non_empty_uploaded,
    internal_example_user_metadata,
    publish,
    test_duplicate,
    user,
    users_dict,
    export_args,
    query_args,
    expected_status_code,
):
    async with temporal_worker():
        non_empty_processed = await asyncio.to_thread(
            lambda: test_processing.run_processing(
                non_empty_uploaded, users_dict[user or 'user0']
            )
        )
        # Create the bundle
        set_upload_entry_metadata(non_empty_processed, internal_example_user_metadata)
        if publish:
            await asyncio.to_thread(non_empty_processed.publish_upload)
            await non_empty_processed.await_workflows()
        upload = non_empty_processed
        upload_id = upload.upload_id
        export_path = os.path.join(config.fs.tmp, 'bundle_' + upload_id)
        export_args_with_defaults = dict(
            export_as_stream=False,
            export_path=export_path,
            zipped=True,
            overwrite=True,
            export_settings=config.bundle_export.default_settings,
        )
        export_args_with_defaults.update(export_args)
        BundleExporter(upload, **export_args_with_defaults).export_bundle()

        if not test_duplicate:
            # Delete the upload so we can import the bundle without id collisions
            upload.delete_upload_local()
        # Finally, import the bundle
        user_auth = auth_headers[user]
        response = await asyncio.to_thread(
            lambda: perform_post_put_file(
                client,
                'POST',
                'uploads/bundle',
                'stream',
                export_path,
                user_auth,
                **query_args,
            )
        )
    assert_response(response, expected_status_code)
    if expected_status_code == 200:
        assert_processing(client, upload_id, user_auth, published=publish)
        upload = Upload.get(upload_id)
        assert upload.from_oasis and upload.oasis_deployment_url


def _raw_path_exists(upload_id: str, path: str):
    return Upload.get(upload_id).upload_files.raw_exists(path)


async def _perform_move_or_copy(
    client: TestClient,
    user,
    upload_id: str,
    source_path: str,
    new_file_name: str,
    copy_or_move: Literal['copy', 'move'],
    # This is the path of the parent folder where is supposed to end the file
    # If empty string it will be stored in the raw directory
    final_destination_folder_path: str,
    trigger_processing: bool = True,
):
    return await asyncio.to_thread(
        lambda: client.put(
            build_url(
                f'uploads/{upload_id}/raw/{final_destination_folder_path}',
                query_args={
                    'copy_or_move': copy_or_move,
                    'file_name': new_file_name,
                    'copy_or_move_source_path': source_path,
                    'trigger_processing': trigger_processing,
                },
            ),
            headers=user,
        )
    )


@pytest.mark.parametrize(
    'source_path, new_file_name, expected_status_code, expected_error_message, orignal_file_should_exist, trigger_processing',
    [
        pytest.param(
            'examples_template/0.aux',
            'random_file_name.aux',
            200,
            None,
            False,
            True,
            id='success-rename-file',
        ),
        pytest.param(
            'examples_template/0.aux',
            'random_file_name.aux',
            200,
            None,
            False,
            False,
            id='success-rename-file-without-reprocessing',
        ),
        pytest.param(
            'examples_template/0.aux',
            '1.aux',
            409,
            'The provided path already exists',
            True,
            True,
            id='conflicting-file-rename',
        ),
        pytest.param(
            'examples_template/non-existing-file.aux',
            'random-name.aux',
            409,
            'No file or folder with that source path',
            False,
            True,
            id='renaming-a-non-existing-file',
        ),
        # TODO: Add folder rename tests when folder rename is supported
    ],
)
@pytest.mark.asyncio
async def test_rename_file_or_folder(
    temporal_worker,
    non_empty_processed_with_temporal: processing.Upload,
    client: TestClient,
    elastic_function,
    auth_headers,
    source_path: str,
    new_file_name: str,
    expected_status_code: int,
    expected_error_message: None | str,
    orignal_file_should_exist: bool,
    trigger_processing: bool,
):
    upload_id: str = non_empty_processed_with_temporal.upload_id
    user = auth_headers['user1']
    async with temporal_worker() as env:
        parent_folder = '/'.join(source_path.split('/')[:-1])
        rename_result = await _perform_move_or_copy(
            client,
            user,
            upload_id,
            source_path=source_path,
            new_file_name=new_file_name,
            copy_or_move='move',
            final_destination_folder_path=parent_folder,
            trigger_processing=trigger_processing,
        )
        assert rename_result.status_code == expected_status_code
        if expected_status_code == 200:
            await _assert_trigger_reprocessing_behavior(
                env,
                client,
                upload_id,
                user,
                trigger_processing,
            )
            assert not _raw_path_exists(upload_id, source_path)
            assert _raw_path_exists(upload_id, f'{parent_folder}/{new_file_name}')
        else:
            body = rename_result.json()
            message = body.get('detail')
            assert message is not None
            if expected_error_message is not None:
                assert expected_error_message in message

            if orignal_file_should_exist:
                assert _raw_path_exists(upload_id, source_path)
            else:
                assert not _raw_path_exists(upload_id, source_path)


async def _get_list_of_started_workflows(temporal_env):
    matching_workflows = []
    async for execution in temporal_env.client.list_workflows():
        matching_workflows.append(execution.raw_info.type.name)
    return matching_workflows


async def _assert_trigger_reprocessing_behavior(
    env,
    client: TestClient,
    upload_id: str,
    user,
    trigger_processing: None | bool,
):
    """
    Waits for possible processing and asserts that the correct workflows are started based on the trigger_processing flag.
    """
    if trigger_processing is True or (trigger_processing is None):
        await asyncio.to_thread(lambda: assert_processing(client, upload_id, user))
        matching_workflows = await _get_list_of_started_workflows(env)
        assert 'UpdateUploadWorkflow' in matching_workflows
        assert 'ProcessUploadWorkflow' in matching_workflows
    else:
        await asyncio.to_thread(lambda: block_until_completed(client, upload_id, user))
        matching_workflows = await _get_list_of_started_workflows(env)
        assert 'UpdateUploadWorkflow' in matching_workflows
        assert 'ProcessUploadWorkflow' not in matching_workflows


@pytest.mark.parametrize(
    'trigger_processing',
    [
        pytest.param(True, id='trigger-processing'),
        pytest.param(False, id='no-processing'),
        pytest.param(None, id='default-processing'),
    ],
)
@pytest.mark.asyncio
async def test_delete_raw_path_trigger_processing_option(
    temporal_worker,
    non_empty_processed_with_temporal,
    client: TestClient,
    auth_headers,
    trigger_processing: None | bool,
):
    upload_id: str = non_empty_processed_with_temporal.upload_id
    user = auth_headers['user1']
    path_to_delete = 'examples_template/1.aux'
    async with temporal_worker() as env:
        url = build_url(
            f'uploads/{upload_id}/raw/{path_to_delete}',
            query_args={'trigger_processing': trigger_processing},
        )
        delete_result = await asyncio.to_thread(
            lambda: client.delete(
                url,
                headers=user,
            )
        )
        assert_response(delete_result, 200)
        await _assert_trigger_reprocessing_behavior(
            env,
            client,
            upload_id,
            user,
            trigger_processing,
        )

        upload_files = StagingUploadFiles(upload_id)
        assert not upload_files.raw_path_exists(path_to_delete)


@pytest.mark.parametrize(
    'trigger_processing',
    [
        pytest.param(True, id='trigger-reprocessing'),
        pytest.param(False, id='no-reprocessing'),
        pytest.param(None, id='default-reprocessing'),
    ],
)
@pytest.mark.asyncio
async def test_put_upload_raw_path_trigger_processing_option(
    temporal_worker,
    non_empty_processed_with_temporal,
    client: TestClient,
    auth_headers,
    trigger_processing: None | bool,
):
    upload_id: str = non_empty_processed_with_temporal.upload_id
    user = auth_headers['user1']
    path_to_upload = example_file_aux  # contains '1.aux'
    async with temporal_worker() as env:
        url = build_url(
            f'uploads/{upload_id}/raw/',
            query_args={
                'trigger_processing': trigger_processing,
            },
        )
        upload_response = await asyncio.to_thread(
            lambda: perform_post_put_file(
                client,
                'PUT',
                url,
                'multipart',
                path_to_upload,
                user,
            )
        )
        assert_response(upload_response, 200)
        await _assert_trigger_reprocessing_behavior(
            env,
            client,
            upload_id,
            user,
            trigger_processing,
        )

        upload_files = StagingUploadFiles(upload_id)
        assert upload_files.raw_path_exists('1.aux')


@pytest.fixture
def create_upload(elastic_function, raw_files_function, mongo_function, user1):
    default_upload = dict(
        upload_id='upload_id',
    )

    default_entry = dict(
        upload_id=default_upload['upload_id'],
        entry_id='entry_id',
    )

    def _create(
        *,
        upload: dict | None = None,
        entry: dict | None = None,
        skip_entry: bool = False,
    ):
        data = ExampleData(main_author=user1)
        upload = default_upload | (upload or {})
        data.create_upload(**upload)

        if not skip_entry:
            entry = default_entry | (entry or {})
            data.create_entry(**entry)

        data.save()
        return data

    return _create


@pytest.mark.parametrize(
    'upload_label, user, datacite_enabled, status_code',
    [
        pytest.param('published', 'user1', True, 200, id='plain'),
        pytest.param('published', None, True, 401, id='no-user'),
        pytest.param('published', 'user2', True, 403, id='wrong-user'),
        pytest.param('published', 'user1', False, 403, id='datacite-disabled'),
        pytest.param('with_doi', 'user1', True, 400, id='with-doi'),
        pytest.param('unpublished', 'user1', True, 400, id='unpublished'),
        pytest.param('empty', 'user1', True, 400, id='empty'),
        pytest.param(None, 'user1', True, 404, id='non-existing'),
    ],
)
def test_assign_doi_upload(
    datacite_mock: DataciteMock,
    auth_headers,
    client,
    create_upload,
    upload_label,
    user,
    datacite_enabled,
    status_code,
):
    datacite_mock.set_enabled(datacite_enabled)

    if upload_label == 'published':
        data = create_upload(upload={'publish_time': now()})
    elif upload_label == 'with_doi':
        upload = {'publish_time': now(), 'doi': {'id': '10.83696/test-doi'}}
        data = create_upload(upload=upload)
    elif upload_label == 'unpublished':
        data = create_upload()
    elif upload_label == 'empty':
        data = create_upload(skip_entry=True)

    headers = auth_headers[user]
    response = client.post(f'uploads/upload_id/action/assign-doi', headers=headers)

    assert_response(response, status_code)
    if not datacite_enabled:
        assert 'not enabled' in response.json()['detail']
    if status_code != 200:
        return

    response = response.json()
    assert_upload(response)
    doi_name = response['data']['doi']['id']
    assert_doi_name(doi_name)


def test_assign_doi_upload_datacite_error(
    datacite_mock: DataciteMock,
    auth_headers,
    client,
    create_upload,
):
    datacite_mock.set_requests(401, False, 'Bad credentials.')
    data = create_upload(upload={'publish_time': now()})

    headers = auth_headers['user1']
    response = client.post(f'uploads/upload_id/action/assign-doi', headers=headers)

    assert_response(response, 500)
    msg = response.json()['detail']
    assert 'An error occurred while creating the DOI draft at DataCite.' in msg
    upload = Upload.get('upload_id')
    assert upload.doi is None


def test_assign_doi_upload_mongo_error(
    mock_mongo_fail_save,
    datacite_mock: DataciteMock,
    auth_headers,
    client,
    create_upload,
):
    data = create_upload(upload={'publish_time': now()})
    mock_mongo_fail_save(Upload)

    headers = auth_headers['user1']
    response = client.post(f'uploads/upload_id/action/assign-doi', headers=headers)

    assert_response(response, 500)
    msg = response.json()['detail']
    assert 'An error occurred while saving the upload doi to the database.' in msg
    upload = Upload.get('upload_id')
    assert upload.doi is None
