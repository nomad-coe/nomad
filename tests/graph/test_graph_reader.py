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

import json
from datetime import datetime

import pytest
import yaml

from nomad.datamodel import EntryArchive, ServerContext
from nomad.graph.graph_reader import (
    EntryReader,
    FileSystemReader,
    GeneralReader,
    MongoReader,
    Token,
    UploadReader,
    UserReader,
)
from nomad.graph.lazy_wrapper import LazyWrapper
from nomad.utils.exampledata import ExampleData
from tests.normalizing.conftest import run_processing, simulationworkflowschema
from tests.utils import ListWithSortKey


def rprint(msg):
    print(msg)
    # try:
    #     import rich
    #     rich.print(msg)
    # except ImportError:
    #     print(msg)


def assert_time(i, j):
    try:
        datetime.fromisoformat(i)
        datetime.fromisoformat(j)
    except Exception:
        assert i == j


def assert_list(observed, expected):
    assert len(observed) == len(expected)
    if isinstance(expected, ListWithSortKey):
        observed = sorted(observed, key=expected.sort_key)
        expected = sorted(expected, key=expected.sort_key)
    for i, j in zip(observed, expected):
        if isinstance(i, LazyWrapper):
            i = i.to_json()
        if isinstance(i, dict):
            assert_dict(i, j)
        elif isinstance(i, list):
            assert_list(i, j)
        else:
            assert_time(i, j)


def assert_dict(observed, expected):
    observed.pop(GeneralReader.__CACHE__, None)
    observed.pop('m_response', None)
    observed.pop('m_def', None)
    observed.pop('m_def_id', None)
    # we do not check the definition ID in this file
    # it has been systematically tested in the other test file
    observed.pop('definition_id', None)
    expected.pop('m_def', None)
    assert set(observed.keys()) == set(expected.keys())
    for k, v in observed.items():
        if isinstance(v, LazyWrapper):
            v = v.to_json()
        if isinstance(v, dict):
            assert_dict(v, expected[k])
        elif isinstance(v, list):
            assert_list(v, expected[k])
        elif k == 'upload_files_server_path':
            continue
        else:
            assert_time(v, expected[k])


user_dict = {
    'name': 'Sheldon Cooper',
    'first_name': 'Sheldon',
    'last_name': 'Cooper',
    'email': 'sheldon.cooper@nomad-coe.eu',
    'user_id': '00000000-0000-0000-0000-000000000001',
    'username': 'scooper',
    'is_admin': False,
    'is_oasis_admin': True,
}


def increment():
    n = 0
    while True:
        n += 1
        yield n


counter = increment()


# noinspection SpellCheckingInspection,DuplicatedCode
def test_remote_reference(json_dict, example_data_with_reference, user1):
    def __user_print(msg, required, *, result: dict = None):
        with UserReader(required, user=user1) as reader:
            if result:
                assert_dict(reader.sync_read(user1.user_id), result)
            else:
                rprint(f'\n\nExample: {next(counter)} -> {msg}:')
                rprint(required)
                rprint('output:')
                rprint(reader.sync_read(user1.user_id))

    __user_print(
        'plain user',
        {'m_request': {'directive': 'plain'}},
        result={
            'name': 'Sheldon Cooper',
            'first_name': 'Sheldon',
            'last_name': 'Cooper',
            'email': 'sheldon.cooper@nomad-coe.eu',
            'user_id': '00000000-0000-0000-0000-000000000001',
            'username': 'scooper',
            'is_admin': False,
            'is_oasis_admin': True,
        },
    )
    __user_print(
        'plain user',
        '*',
        result={
            'name': 'Sheldon Cooper',
            'first_name': 'Sheldon',
            'last_name': 'Cooper',
            'email': 'sheldon.cooper@nomad-coe.eu',
            'user_id': '00000000-0000-0000-0000-000000000001',
            'username': 'scooper',
            'is_admin': False,
            'is_oasis_admin': True,
        },
    )

    __user_print(
        'link to uploads',
        {
            'm_request': {'directive': 'plain'},
            Token.UPLOADS: {
                'm_request': {'directive': 'plain'},
            },
        },
        result={
            'name': 'Sheldon Cooper',
            'first_name': 'Sheldon',
            'last_name': 'Cooper',
            'email': 'sheldon.cooper@nomad-coe.eu',
            'user_id': '00000000-0000-0000-0000-000000000001',
            'username': 'scooper',
            'is_admin': False,
            'is_oasis_admin': True,
            Token.UPLOADS: {'id_published_with_ref': 'id_published_with_ref'},
        },
    )

    __user_print(
        'link to uploads, resolve with metadata',
        {
            'm_request': {'directive': 'plain'},
            Token.UPLOADS: {
                'm_request': {'directive': 'resolved', 'resolve_type': 'upload'},
            },
        },
        result={
            'name': 'Sheldon Cooper',
            'first_name': 'Sheldon',
            'last_name': 'Cooper',
            'email': 'sheldon.cooper@nomad-coe.eu',
            'user_id': '00000000-0000-0000-0000-000000000001',
            'username': 'scooper',
            'is_admin': False,
            'is_oasis_admin': True,
            Token.UPLOADS: {
                'id_published_with_ref': {
                    'process_running': False,
                    'current_process': 'process_upload',
                    'process_status': 'SUCCESS',
                    'last_status_message': None,
                    'errors': [],
                    'warnings': [],
                    'complete_time': '2023-03-05T21:24:22.172000',
                    'upload_id': 'id_published_with_ref',
                    'upload_name': 'name_published',
                    'upload_create_time': '2023-03-05T21:24:22.171000',
                    'main_author': user_dict,
                    'coauthors': [],
                    'reviewers': [],
                    'coauthor_groups': [],
                    'reviewer_groups': [],
                    'writers': [user_dict],
                    'viewers': [user_dict],
                    'writer_groups': [],
                    'viewer_groups': [],
                    'published': False,
                    'processing_failed': 0,
                    'processing_successful': 6,
                    'published_to': [],
                    'publish_time': None,
                    'with_embargo': False,
                    'embargo_length': 0,
                    'license': 'CC BY 4.0',
                    'n_entries': 6,
                    'upload_files_server_path': 'id_published_with_ref',
                }
            },
        },
    )

    __user_print(
        'link to uploads, resolve with metadata, form 2 using dict style with explicit upload id',
        {
            'm_request': {'directive': 'plain'},
            Token.UPLOADS: {
                'id_published_with_ref': {
                    'm_request': {
                        'directive': 'plain',
                    },
                }
            },
        },
        result={
            'name': 'Sheldon Cooper',
            'first_name': 'Sheldon',
            'last_name': 'Cooper',
            'email': 'sheldon.cooper@nomad-coe.eu',
            'user_id': '00000000-0000-0000-0000-000000000001',
            'username': 'scooper',
            'is_admin': False,
            'is_oasis_admin': True,
            Token.UPLOADS: {
                'id_published_with_ref': {
                    'process_running': False,
                    'current_process': 'process_upload',
                    'process_status': 'SUCCESS',
                    'last_status_message': None,
                    'errors': [],
                    'warnings': [],
                    'complete_time': '2023-03-05T21:30:22.807000',
                    'upload_id': 'id_published_with_ref',
                    'upload_name': 'name_published',
                    'upload_create_time': '2023-03-05T21:30:22.806000',
                    'main_author': user_dict,
                    'coauthors': [],
                    'reviewers': [],
                    'coauthor_groups': [],
                    'reviewer_groups': [],
                    'writers': [user_dict],
                    'viewers': [user_dict],
                    'writer_groups': [],
                    'viewer_groups': [],
                    'published': False,
                    'processing_failed': 0,
                    'processing_successful': 6,
                    'published_to': [],
                    'publish_time': None,
                    'with_embargo': False,
                    'embargo_length': 0,
                    'license': 'CC BY 4.0',
                    'n_entries': 6,
                    'upload_files_server_path': 'id_published_with_ref',
                }
            },
        },
    )

    __user_print(
        'link to entries directly from user, resolve with metadata',
        {
            'm_request': {'directive': 'plain'},
            Token.ENTRIES: {
                'm_request': {'directive': 'resolved', 'resolve_type': 'entry'},
            },
        },
        result={
            'name': 'Sheldon Cooper',
            'first_name': 'Sheldon',
            'last_name': 'Cooper',
            'email': 'sheldon.cooper@nomad-coe.eu',
            'user_id': '00000000-0000-0000-0000-000000000001',
            'username': 'scooper',
            'is_admin': False,
            'is_oasis_admin': True,
            Token.ENTRIES: {
                'id_01': {
                    'process_running': False,
                    'current_process': None,
                    'process_status': 'SUCCESS',
                    'last_status_message': None,
                    'errors': [],
                    'warnings': [],
                    'complete_time': None,
                    'entry_id': 'id_01',
                    'entry_create_time': '2023-03-05T21:27:24.488000',
                    'mainfile_path': 'mainfile_for_id_01',
                    'mainfile_key': None,
                    'upload_id': 'id_published_with_ref',
                    'parser_name': 'parsers/vasp',
                },
                'id_02': {
                    'process_running': False,
                    'current_process': None,
                    'process_status': 'SUCCESS',
                    'last_status_message': None,
                    'errors': [],
                    'warnings': [],
                    'complete_time': None,
                    'entry_id': 'id_02',
                    'entry_create_time': '2023-03-05T21:27:24.489000',
                    'mainfile_path': 'mainfile_for_id_02',
                    'mainfile_key': None,
                    'upload_id': 'id_published_with_ref',
                    'parser_name': 'parsers/vasp',
                },
                'id_03': {
                    'process_running': False,
                    'current_process': None,
                    'process_status': 'SUCCESS',
                    'last_status_message': None,
                    'errors': [],
                    'warnings': [],
                    'complete_time': None,
                    'entry_id': 'id_03',
                    'entry_create_time': '2023-03-05T21:27:24.490000',
                    'mainfile_path': 'mainfile_for_id_03',
                    'mainfile_key': None,
                    'upload_id': 'id_published_with_ref',
                    'parser_name': 'parsers/vasp',
                },
                'id_04': {
                    'process_running': False,
                    'current_process': None,
                    'process_status': 'SUCCESS',
                    'last_status_message': None,
                    'errors': [],
                    'warnings': [],
                    'complete_time': None,
                    'entry_id': 'id_04',
                    'entry_create_time': '2023-03-05T21:27:24.491000',
                    'mainfile_path': 'mainfile_for_id_04',
                    'mainfile_key': None,
                    'upload_id': 'id_published_with_ref',
                    'parser_name': 'parsers/vasp',
                },
                'id_05': {
                    'process_running': False,
                    'current_process': None,
                    'process_status': 'SUCCESS',
                    'last_status_message': None,
                    'errors': [],
                    'warnings': [],
                    'complete_time': None,
                    'entry_id': 'id_05',
                    'entry_create_time': '2023-03-05T21:27:24.492000',
                    'mainfile_path': 'mainfile_for_id_05',
                    'mainfile_key': None,
                    'upload_id': 'id_published_with_ref',
                    'parser_name': 'parsers/vasp',
                },
                'id_06': {
                    'process_running': False,
                    'current_process': None,
                    'process_status': 'SUCCESS',
                    'last_status_message': None,
                    'errors': [],
                    'warnings': [],
                    'complete_time': None,
                    'entry_id': 'id_06',
                    'entry_create_time': '2023-03-05T21:27:24.493000',
                    'mainfile_path': 'mainfile_for_id_06',
                    'mainfile_key': None,
                    'upload_id': 'id_published_with_ref',
                    'parser_name': 'parsers/vasp',
                },
            },
        },
    )

    __user_print(
        'link to entries directly from user, resolve with metadata, form 2 using dict style with explicit entry id',
        {'m_request': {'directive': 'plain'}, Token.ENTRIES: {'id_01': '*'}},
        result={
            'name': 'Sheldon Cooper',
            'first_name': 'Sheldon',
            'last_name': 'Cooper',
            'email': 'sheldon.cooper@nomad-coe.eu',
            'user_id': '00000000-0000-0000-0000-000000000001',
            'username': 'scooper',
            'is_admin': False,
            'is_oasis_admin': True,
            Token.ENTRIES: {
                'id_01': {
                    'process_running': False,
                    'current_process': None,
                    'process_status': 'SUCCESS',
                    'last_status_message': None,
                    'errors': [],
                    'warnings': [],
                    'complete_time': None,
                    'entry_id': 'id_01',
                    'entry_create_time': '2023-03-05T21:30:22.809000',
                    'mainfile_path': 'mainfile_for_id_01',
                    'mainfile_key': None,
                    'upload_id': 'id_published_with_ref',
                    'parser_name': 'parsers/vasp',
                }
            },
        },
    )

    __user_print(
        'link to entries from uploads, resolve with metadata',
        {
            'm_request': {'directive': 'plain'},
            Token.UPLOADS: {
                'id_published_with_ref': {
                    'm_request': {
                        'directive': 'plain',
                    },
                    Token.ENTRIES: {
                        'm_request': {'directive': 'resolved', 'resolve_type': 'entry'},
                    },
                }
            },
        },
        result={
            'name': 'Sheldon Cooper',
            'first_name': 'Sheldon',
            'last_name': 'Cooper',
            'email': 'sheldon.cooper@nomad-coe.eu',
            'user_id': '00000000-0000-0000-0000-000000000001',
            'username': 'scooper',
            'is_admin': False,
            'is_oasis_admin': True,
            Token.UPLOADS: {
                'id_published_with_ref': {
                    'process_running': False,
                    'current_process': 'process_upload',
                    'process_status': 'SUCCESS',
                    'last_status_message': None,
                    'errors': [],
                    'warnings': [],
                    'complete_time': '2023-03-05T21:31:56.873000',
                    'upload_id': 'id_published_with_ref',
                    'upload_name': 'name_published',
                    'upload_create_time': '2023-03-05T21:31:56.872000',
                    'main_author': user_dict,
                    'coauthors': [],
                    'reviewers': [],
                    'coauthor_groups': [],
                    'reviewer_groups': [],
                    'writers': [user_dict],
                    'viewers': [user_dict],
                    'writer_groups': [],
                    'viewer_groups': [],
                    'published': False,
                    'processing_failed': 0,
                    'processing_successful': 6,
                    'published_to': [],
                    'publish_time': None,
                    'with_embargo': False,
                    'embargo_length': 0,
                    'license': 'CC BY 4.0',
                    'n_entries': 6,
                    'upload_files_server_path': 'id_published_with_ref',
                    Token.ENTRIES: {
                        'id_01': {
                            'process_running': False,
                            'current_process': None,
                            'process_status': 'SUCCESS',
                            'last_status_message': None,
                            'errors': [],
                            'warnings': [],
                            'complete_time': None,
                            'entry_id': 'id_01',
                            'entry_create_time': '2023-03-05T21:31:56.875000',
                            'mainfile_path': 'mainfile_for_id_01',
                            'mainfile_key': None,
                            'upload_id': 'id_published_with_ref',
                            'parser_name': 'parsers/vasp',
                        },
                        'id_02': {
                            'process_running': False,
                            'current_process': None,
                            'process_status': 'SUCCESS',
                            'last_status_message': None,
                            'errors': [],
                            'warnings': [],
                            'complete_time': None,
                            'entry_id': 'id_02',
                            'entry_create_time': '2023-03-05T21:31:56.876000',
                            'mainfile_path': 'mainfile_for_id_02',
                            'mainfile_key': None,
                            'upload_id': 'id_published_with_ref',
                            'parser_name': 'parsers/vasp',
                        },
                        'id_03': {
                            'process_running': False,
                            'current_process': None,
                            'process_status': 'SUCCESS',
                            'last_status_message': None,
                            'errors': [],
                            'warnings': [],
                            'complete_time': None,
                            'entry_id': 'id_03',
                            'entry_create_time': '2023-03-05T21:31:56.877000',
                            'mainfile_path': 'mainfile_for_id_03',
                            'mainfile_key': None,
                            'upload_id': 'id_published_with_ref',
                            'parser_name': 'parsers/vasp',
                        },
                        'id_04': {
                            'process_running': False,
                            'current_process': None,
                            'process_status': 'SUCCESS',
                            'last_status_message': None,
                            'errors': [],
                            'warnings': [],
                            'complete_time': None,
                            'entry_id': 'id_04',
                            'entry_create_time': '2023-03-05T21:31:56.878000',
                            'mainfile_path': 'mainfile_for_id_04',
                            'mainfile_key': None,
                            'upload_id': 'id_published_with_ref',
                            'parser_name': 'parsers/vasp',
                        },
                        'id_05': {
                            'process_running': False,
                            'current_process': None,
                            'process_status': 'SUCCESS',
                            'last_status_message': None,
                            'errors': [],
                            'warnings': [],
                            'complete_time': None,
                            'entry_id': 'id_05',
                            'entry_create_time': '2023-03-05T21:31:56.879000',
                            'mainfile_path': 'mainfile_for_id_05',
                            'mainfile_key': None,
                            'upload_id': 'id_published_with_ref',
                            'parser_name': 'parsers/vasp',
                        },
                        'id_06': {
                            'process_running': False,
                            'current_process': None,
                            'process_status': 'SUCCESS',
                            'last_status_message': None,
                            'errors': [],
                            'warnings': [],
                            'complete_time': None,
                            'entry_id': 'id_06',
                            'entry_create_time': '2023-03-05T21:31:56.880000',
                            'mainfile_path': 'mainfile_for_id_06',
                            'mainfile_key': None,
                            'upload_id': 'id_published_with_ref',
                            'parser_name': 'parsers/vasp',
                        },
                    },
                }
            },
        },
    )

    __user_print(
        'link to entries from uploads, resolve with metadata, dict style',
        {
            'm_request': {'directive': 'plain'},
            Token.UPLOADS: {
                'id_published_with_ref': {
                    'm_request': {
                        'directive': 'plain',
                    },
                    Token.ENTRIES: {
                        'id_01': {
                            'm_request': {
                                'directive': 'plain',
                            }
                        },
                        '*': {
                            'm_request': {
                                'directive': 'plain',
                                'include': ['entry_id', 'mainfile_path'],
                            },
                        },
                    },
                }
            },
        },
        result={
            'name': 'Sheldon Cooper',
            'first_name': 'Sheldon',
            'last_name': 'Cooper',
            'email': 'sheldon.cooper@nomad-coe.eu',
            'user_id': '00000000-0000-0000-0000-000000000001',
            'username': 'scooper',
            'is_admin': False,
            'is_oasis_admin': True,
            Token.UPLOADS: {
                'id_published_with_ref': {
                    'process_running': False,
                    'current_process': 'process_upload',
                    'process_status': 'SUCCESS',
                    'last_status_message': None,
                    'errors': [],
                    'warnings': [],
                    'complete_time': '2023-03-05T21:31:56.873000',
                    'upload_id': 'id_published_with_ref',
                    'upload_name': 'name_published',
                    'upload_create_time': '2023-03-05T21:31:56.872000',
                    'main_author': user_dict,
                    'coauthors': [],
                    'reviewers': [],
                    'coauthor_groups': [],
                    'reviewer_groups': [],
                    'writers': [user_dict],
                    'viewers': [user_dict],
                    'writer_groups': [],
                    'viewer_groups': [],
                    'published': False,
                    'processing_failed': 0,
                    'processing_successful': 6,
                    'published_to': [],
                    'publish_time': None,
                    'with_embargo': False,
                    'embargo_length': 0,
                    'license': 'CC BY 4.0',
                    'n_entries': 6,
                    'upload_files_server_path': 'id_published_with_ref',
                    Token.ENTRIES: {
                        'id_02': {
                            'entry_id': 'id_02',
                            'mainfile_path': 'mainfile_for_id_02',
                        },
                        'id_03': {
                            'entry_id': 'id_03',
                            'mainfile_path': 'mainfile_for_id_03',
                        },
                        'id_04': {
                            'entry_id': 'id_04',
                            'mainfile_path': 'mainfile_for_id_04',
                        },
                        'id_05': {
                            'entry_id': 'id_05',
                            'mainfile_path': 'mainfile_for_id_05',
                        },
                        'id_06': {
                            'entry_id': 'id_06',
                            'mainfile_path': 'mainfile_for_id_06',
                        },
                        'id_01': {
                            'process_running': False,
                            'current_process': None,
                            'process_status': 'SUCCESS',
                            'last_status_message': None,
                            'errors': [],
                            'warnings': [],
                            'complete_time': None,
                            'entry_id': 'id_01',
                            'entry_create_time': '2023-03-05T21:31:56.875000',
                            'mainfile_path': 'mainfile_for_id_01',
                            'mainfile_key': None,
                            'upload_id': 'id_published_with_ref',
                            'parser_name': 'parsers/vasp',
                        },
                    },
                }
            },
        },
    )

    __user_print(
        'uploads to entries back to uploads',
        {
            Token.UPLOADS: {
                'id_published_with_ref': {
                    Token.ENTRIES: {
                        'id_01': {
                            'm_request': {
                                'directive': 'plain',
                            },
                            'upload_id': {
                                'm_request': {
                                    'directive': 'resolved',
                                    'resolve_type': 'upload',
                                }
                            },
                        }
                    }
                }
            }
        },
        result={
            Token.UPLOADS: {
                'id_published_with_ref': {
                    Token.ENTRIES: {
                        'id_01': {
                            'process_running': False,
                            'current_process': None,
                            'process_status': 'SUCCESS',
                            'last_status_message': None,
                            'errors': [],
                            'warnings': [],
                            'complete_time': None,
                            'entry_id': 'id_01',
                            'entry_create_time': '2023-03-05T21:31:56.875000',
                            'mainfile_path': 'mainfile_for_id_01',
                            'mainfile_key': None,
                            'parser_name': 'parsers/vasp',
                            'upload_id': {
                                'process_running': False,
                                'current_process': 'process_upload',
                                'process_status': 'SUCCESS',
                                'last_status_message': None,
                                'errors': [],
                                'warnings': [],
                                'complete_time': '2023-03-05T21:31:56.873000',
                                'upload_id': 'id_published_with_ref',
                                'upload_name': 'name_published',
                                'upload_create_time': '2023-03-05T21:31:56.872000',
                                'main_author': user_dict,
                                'coauthors': [],
                                'reviewers': [],
                                'coauthor_groups': [],
                                'reviewer_groups': [],
                                'writers': [user_dict],
                                'viewers': [user_dict],
                                'writer_groups': [],
                                'viewer_groups': [],
                                'published': False,
                                'processing_failed': 0,
                                'processing_successful': 6,
                                'published_to': [],
                                'publish_time': None,
                                'with_embargo': False,
                                'embargo_length': 0,
                                'license': 'CC BY 4.0',
                                'n_entries': 6,
                                'upload_files_server_path': 'id_published_with_ref',
                            },
                        }
                    }
                }
            }
        },
    )

    __user_print(
        'uploads to entries to archive',
        {
            Token.UPLOADS: {
                'id_published_with_ref': {
                    Token.ENTRIES: {
                        'id_01': {
                            'm_request': {
                                'directive': 'plain',
                            },
                            Token.ARCHIVE: {
                                'm_request': {
                                    'directive': 'plain',
                                    'include': ['results'],
                                },
                            },
                        }
                    }
                }
            }
        },
        result={
            Token.UPLOADS: {
                'id_published_with_ref': {
                    Token.ENTRIES: {
                        'id_01': {
                            'process_running': False,
                            'current_process': None,
                            'process_status': 'SUCCESS',
                            'last_status_message': None,
                            'errors': [],
                            'warnings': [],
                            'complete_time': None,
                            'entry_id': 'id_01',
                            'entry_create_time': '2023-03-05T21:31:56.875000',
                            'mainfile_path': 'mainfile_for_id_01',
                            'mainfile_key': None,
                            'upload_id': 'id_published_with_ref',
                            'parser_name': 'parsers/vasp',
                            Token.ARCHIVE: {
                                'results': {
                                    'material': {
                                        'dimensionality': '3D',
                                        'material_id': 'test_material_id',
                                        'elements': ['H', 'O'],
                                        'symmetry': {'crystal_system': 'cubic'},
                                    },
                                    'method': {
                                        'simulation': {
                                            'program_name': 'VASP',
                                            'dft': {'xc_functional_type': 'GGA'},
                                        }
                                    },
                                    'properties': {
                                        'n_calculations': 1,
                                        'electronic': {
                                            'dos_electronic': [
                                                {
                                                    'spin_polarized': False,
                                                    'band_gap': [{'type': 'indirect'}],
                                                }
                                            ]
                                        },
                                    },
                                }
                            },
                        }
                    }
                }
            }
        },
    )

    def __upload_print(msg, required, *, result: dict = None):
        with UploadReader(required, user=user1) as reader:
            if result:
                assert_dict(reader.sync_read('id_published_with_ref'), result)
            else:
                rprint(f'\n\nExample: {next(counter)} -> {msg}:')
                rprint(required)
                rprint('output:')
                rprint(reader.sync_read('id_published_with_ref'))

    __upload_print(
        'plain upload reader',
        {
            'm_request': {
                'directive': 'plain',
            },
        },
        result={
            'process_running': False,
            'current_process': 'process_upload',
            'process_status': 'SUCCESS',
            'last_status_message': None,
            'errors': [],
            'warnings': [],
            'complete_time': '2023-03-05T22:12:08.421000',
            'upload_id': 'id_published_with_ref',
            'upload_name': 'name_published',
            'upload_create_time': '2023-03-05T22:12:08.420000',
            'main_author': user_dict,
            'coauthors': [],
            'reviewers': [],
            'coauthor_groups': [],
            'reviewer_groups': [],
            'writers': [user_dict],
            'viewers': [user_dict],
            'writer_groups': [],
            'viewer_groups': [],
            'published': False,
            'published_to': [],
            'publish_time': None,
            'with_embargo': False,
            'embargo_length': 0,
            'processing_failed': 0,
            'processing_successful': 6,
            'license': 'CC BY 4.0',
            'n_entries': 6,
            'upload_files_server_path': 'id_published_with_ref',
        },
    )

    __upload_print(
        'upload, resolve raw files',
        {
            'm_request': '*',
            Token.RAW: {
                'm_request': '*',
                'mainfile_for_id_01': {
                    'm_request': {
                        'directive': 'resolved',
                    },
                },
            },
        },
        result={
            'process_running': False,
            'current_process': 'process_upload',
            'process_status': 'SUCCESS',
            'last_status_message': None,
            'errors': [],
            'warnings': [],
            'complete_time': '2023-03-05T22:16:52.436000',
            'upload_id': 'id_published_with_ref',
            'upload_name': 'name_published',
            'upload_create_time': '2023-03-05T22:16:52.435000',
            'main_author': user_dict,
            'coauthors': [],
            'reviewers': [],
            'coauthor_groups': [],
            'reviewer_groups': [],
            'writers': [user_dict],
            'viewers': [user_dict],
            'writer_groups': [],
            'viewer_groups': [],
            'published': False,
            'published_to': [],
            'publish_time': None,
            'with_embargo': False,
            'embargo_length': 0,
            'processing_failed': 0,
            'processing_successful': 6,
            'license': 'CC BY 4.0',
            'n_entries': 6,
            'upload_files_server_path': 'id_published_with_ref',
            Token.RAW: {
                'm_is': 'Directory',
                '1.aux': {'path': '1.aux', 'm_is': 'File', 'size': 8},
                '2.aux': {'path': '2.aux', 'm_is': 'File', 'size': 8},
                '3.aux': {'path': '3.aux', 'm_is': 'File', 'size': 8},
                '4.aux': {'path': '4.aux', 'm_is': 'File', 'size': 8},
                'mainfile_for_id_02': {
                    'path': 'mainfile_for_id_02',
                    'm_is': 'File',
                    'size': 3227,
                },
                'mainfile_for_id_03': {
                    'path': 'mainfile_for_id_03',
                    'm_is': 'File',
                    'size': 3227,
                },
                'mainfile_for_id_04': {
                    'path': 'mainfile_for_id_04',
                    'm_is': 'File',
                    'size': 3227,
                },
                'mainfile_for_id_05': {
                    'path': 'mainfile_for_id_05',
                    'm_is': 'File',
                    'size': 3227,
                },
                'mainfile_for_id_06': {
                    'path': 'mainfile_for_id_06',
                    'm_is': 'File',
                    'size': 3227,
                },
                'mainfile_for_id_01': {
                    'path': 'mainfile_for_id_01',
                    'm_is': 'File',
                    'size': 3227,
                    Token.ENTRY: {
                        'process_running': False,
                        'current_process': None,
                        'process_status': 'SUCCESS',
                        'last_status_message': None,
                        'errors': [],
                        'warnings': [],
                        'complete_time': None,
                        'entry_id': 'id_01',
                        'entry_create_time': '2023-03-05T22:16:52.438000',
                        'mainfile_path': 'mainfile_for_id_01',
                        'mainfile_key': None,
                        'upload_id': 'id_published_with_ref',
                        'parser_name': 'parsers/vasp',
                    },
                },
            },
        },
    )

    __upload_print(
        'upload, resolve user',
        {
            'm_request': {
                'directive': 'plain',
            },
            'viewers': {'m_request': {'directive': 'resolved', 'resolve_type': 'user'}},
        },
        result={
            'process_running': False,
            'current_process': 'process_upload',
            'process_status': 'SUCCESS',
            'last_status_message': None,
            'errors': [],
            'warnings': [],
            'complete_time': '2023-03-05T22:16:52.436000',
            'upload_id': 'id_published_with_ref',
            'upload_name': 'name_published',
            'upload_create_time': '2023-03-05T22:16:52.435000',
            'main_author': user_dict,
            'coauthors': [],
            'reviewers': [],
            'coauthor_groups': [],
            'reviewer_groups': [],
            'writers': [user_dict],
            'writer_groups': [],
            'viewer_groups': [],
            'published': False,
            'published_to': [],
            'publish_time': None,
            'with_embargo': False,
            'embargo_length': 0,
            'license': 'CC BY 4.0',
            'n_entries': 6,
            'processing_failed': 0,
            'processing_successful': 6,
            'upload_files_server_path': 'id_published_with_ref',
            'viewers': [user_dict],
        },
    )

    __upload_print(
        'resolve itself using upload id',
        {
            'm_request': {
                'directive': 'plain',
            },
            'upload_id': {
                'm_request': {'directive': 'resolved', 'resolve_type': 'upload'}
            },
        },
        result={
            'process_running': False,
            'current_process': 'process_upload',
            'process_status': 'SUCCESS',
            'last_status_message': None,
            'errors': [],
            'warnings': [],
            'complete_time': '2023-03-05T22:16:52.436000',
            'upload_name': 'name_published',
            'upload_create_time': '2023-03-05T22:16:52.435000',
            'main_author': user_dict,
            'coauthors': [],
            'reviewers': [],
            'coauthor_groups': [],
            'reviewer_groups': [],
            'writers': [user_dict],
            'viewers': [user_dict],
            'writer_groups': [],
            'viewer_groups': [],
            'published': False,
            'published_to': [],
            'publish_time': None,
            'with_embargo': False,
            'embargo_length': 0,
            'license': 'CC BY 4.0',
            'n_entries': 6,
            'processing_failed': 0,
            'processing_successful': 6,
            'upload_files_server_path': 'id_published_with_ref',
            'upload_id': {
                'process_running': False,
                'current_process': 'process_upload',
                'process_status': 'SUCCESS',
                'last_status_message': None,
                'errors': [],
                'warnings': [],
                'complete_time': '2023-03-05T22:16:52.436000',
                'upload_id': 'id_published_with_ref',
                'upload_name': 'name_published',
                'upload_create_time': '2023-03-05T22:16:52.435000',
                'main_author': user_dict,
                'coauthors': [],
                'reviewers': [],
                'coauthor_groups': [],
                'reviewer_groups': [],
                'writers': [user_dict],
                'viewers': [user_dict],
                'writer_groups': [],
                'viewer_groups': [],
                'published': False,
                'published_to': [],
                'publish_time': None,
                'with_embargo': False,
                'embargo_length': 0,
                'license': 'CC BY 4.0',
                'n_entries': 6,
                'processing_failed': 0,
                'processing_successful': 6,
                'upload_files_server_path': 'id_published_with_ref',
            },
        },
    )

    __upload_print(
        'resolve itself twice then go to entry',
        {
            'upload_id': {
                'm_request': {
                    'directive': 'plain',
                },
                'upload_id': {
                    Token.ENTRIES: {
                        'id_01': {
                            'm_request': {
                                'directive': 'plain',
                            },
                        }
                    }
                },
            }
        },
        result={
            'upload_id': {
                'process_running': False,
                'current_process': 'process_upload',
                'process_status': 'SUCCESS',
                'last_status_message': None,
                'errors': [],
                'warnings': [],
                'complete_time': '2023-03-05T22:16:52.436000',
                'upload_name': 'name_published',
                'upload_create_time': '2023-03-05T22:16:52.435000',
                'main_author': user_dict,
                'coauthors': [],
                'reviewers': [],
                'coauthor_groups': [],
                'reviewer_groups': [],
                'writers': [user_dict],
                'viewers': [user_dict],
                'writer_groups': [],
                'viewer_groups': [],
                'published': False,
                'published_to': [],
                'publish_time': None,
                'with_embargo': False,
                'embargo_length': 0,
                'license': 'CC BY 4.0',
                'n_entries': 6,
                'processing_failed': 0,
                'processing_successful': 6,
                'upload_files_server_path': 'id_published_with_ref',
                'upload_id': {
                    Token.ENTRIES: {
                        'id_01': {
                            'process_running': False,
                            'current_process': None,
                            'process_status': 'SUCCESS',
                            'last_status_message': None,
                            'errors': [],
                            'warnings': [],
                            'complete_time': None,
                            'entry_id': 'id_01',
                            'entry_create_time': '2023-03-05T22:16:52.438000',
                            'mainfile_path': 'mainfile_for_id_01',
                            'mainfile_key': None,
                            'upload_id': 'id_published_with_ref',
                            'parser_name': 'parsers/vasp',
                        }
                    }
                },
            }
        },
    )

    __upload_print(
        'resolve itself twice then go to entry, collect other info on different levels',
        {
            'm_request': {
                'directive': 'plain',
            },
            'upload_id': {
                'upload_id': {
                    Token.ENTRIES: {
                        'id_01': {
                            'm_request': {
                                'directive': 'plain',
                            },
                        }
                    }
                }
            },
        },
        result={
            'process_running': False,
            'current_process': 'process_upload',
            'process_status': 'SUCCESS',
            'last_status_message': None,
            'errors': [],
            'warnings': [],
            'complete_time': '2023-03-05T22:16:52.436000',
            'upload_name': 'name_published',
            'upload_create_time': '2023-03-05T22:16:52.435000',
            'main_author': user_dict,
            'coauthors': [],
            'reviewers': [],
            'coauthor_groups': [],
            'reviewer_groups': [],
            'writers': [user_dict],
            'viewers': [user_dict],
            'writer_groups': [],
            'viewer_groups': [],
            'published': False,
            'published_to': [],
            'publish_time': None,
            'with_embargo': False,
            'embargo_length': 0,
            'license': 'CC BY 4.0',
            'n_entries': 6,
            'processing_failed': 0,
            'processing_successful': 6,
            'upload_files_server_path': 'id_published_with_ref',
            'upload_id': {
                'upload_id': {
                    Token.ENTRIES: {
                        'id_01': {
                            'process_running': False,
                            'current_process': None,
                            'process_status': 'SUCCESS',
                            'last_status_message': None,
                            'errors': [],
                            'warnings': [],
                            'complete_time': None,
                            'entry_id': 'id_01',
                            'entry_create_time': '2023-03-05T22:16:52.438000',
                            'mainfile_path': 'mainfile_for_id_01',
                            'mainfile_key': None,
                            'upload_id': 'id_published_with_ref',
                            'parser_name': 'parsers/vasp',
                        }
                    }
                }
            },
        },
    )

    def __entry_print(msg, required, *, to_file: bool = False, result: dict = None):
        with EntryReader(required, user=user1) as reader:
            if result:
                assert_dict(reader.sync_read('id_03'), result)
            else:
                rprint(f'\n\nExample: {next(counter)} -> {msg}:')
                rprint(required)
                if not to_file:
                    rprint('output:')
                    rprint(reader.sync_read('id_03'))
                else:
                    with open('entry_reader_test.json', 'w') as f:
                        f.write(json.dumps(reader.sync_read('id_03')))

    __entry_print(
        'plain entry reader',
        {
            'm_request': {
                'directive': 'plain',
            },
        },
        result={
            'process_running': False,
            'current_process': None,
            'process_status': 'SUCCESS',
            'last_status_message': None,
            'errors': [],
            'warnings': [],
            'complete_time': None,
            'entry_id': 'id_03',
            'entry_create_time': '2023-03-05T22:16:52.440000',
            'mainfile_path': 'mainfile_for_id_03',
            'mainfile_key': None,
            'upload_id': 'id_published_with_ref',
            'parser_name': 'parsers/vasp',
        },
    )
    __entry_print(
        'plain entry reader, resolve inplace',
        {
            Token.ARCHIVE: {
                'm_request': {
                    'directive': 'resolved',
                    'resolve_inplace': True,
                    'include': ['workflow2'],
                },
            }
        },
        result={
            Token.UPLOADS: {
                'id_published_with_ref': {
                    Token.ENTRIES: {
                        'id_01': {
                            Token.ARCHIVE: {
                                'workflow2': {
                                    'results': {
                                        'calculation_result_ref': 'uploads/id_published_with_ref/entries/id_01/archive/run/0/calculation/1'
                                    }
                                },
                                'run': [
                                    {
                                        'calculation': [
                                            None,
                                            {
                                                'system_ref': 'uploads/id_published_with_ref/entries/id_01/archive/run/0/system/1',
                                                'energy': {'total': {'value': 0.2}},
                                                'dos_electronic': [
                                                    {
                                                        'energies': [
                                                            0.0,
                                                            0.1,
                                                            1.0,
                                                            1.0,
                                                            1.0,
                                                            1.0,
                                                            1.0,
                                                            1.0,
                                                            1.0,
                                                            1.0,
                                                            1.0,
                                                            1.0,
                                                            1.0,
                                                            1.0,
                                                        ]
                                                    }
                                                ],
                                            },
                                        ],
                                        'system': [
                                            None,
                                            {
                                                'atoms': {'labels': ['H']},
                                                'symmetry': [
                                                    {'space_group_number': 221}
                                                ],
                                            },
                                        ],
                                    }
                                ],
                            }
                        }
                    }
                }
            },
            Token.ARCHIVE: {
                'workflow2': {
                    'tasks': [
                        {
                            'task': 'uploads/id_published_with_ref/entries/id_01/archive/workflow2'
                        }
                    ]
                }
            },
        },
    )
    __entry_print(
        'plain entry reader, resolve to root',
        {
            Token.ARCHIVE: {
                'm_request': {
                    'directive': 'resolved',
                    'resolve_inplace': False,
                    'include': ['workflow2'],
                },
            }
        },
        result={
            Token.UPLOADS: {
                'id_published_with_ref': {
                    Token.ENTRIES: {
                        'id_01': {
                            Token.ARCHIVE: {
                                'workflow2': {
                                    'results': {
                                        'calculation_result_ref': 'uploads/id_published_with_ref/entries/id_01/archive/run/0/calculation/1'
                                    }
                                },
                                'run': [
                                    {
                                        'calculation': [
                                            None,
                                            {
                                                'system_ref': 'uploads/id_published_with_ref/entries/id_01/archive/run/0/system/1',
                                                'energy': {'total': {'value': 0.2}},
                                                'dos_electronic': [
                                                    {
                                                        'energies': [
                                                            0.0,
                                                            0.1,
                                                            1.0,
                                                            1.0,
                                                            1.0,
                                                            1.0,
                                                            1.0,
                                                            1.0,
                                                            1.0,
                                                            1.0,
                                                            1.0,
                                                            1.0,
                                                            1.0,
                                                            1.0,
                                                        ]
                                                    }
                                                ],
                                            },
                                        ],
                                        'system': [
                                            None,
                                            {
                                                'atoms': {'labels': ['H']},
                                                'symmetry': [
                                                    {'space_group_number': 221}
                                                ],
                                            },
                                        ],
                                    }
                                ],
                            }
                        }
                    }
                }
            },
            Token.ARCHIVE: {
                'workflow2': {
                    'tasks': [
                        {
                            'task': 'uploads/id_published_with_ref/entries/id_01/archive/workflow2'
                        }
                    ]
                }
            },
        },
    )
    __entry_print(
        'plain entry reader, resolve to root',
        {
            Token.ARCHIVE: {
                'metadata': {
                    'm_request': {'directive': 'plain', 'depth': 1, 'max_list_size': 1},
                }
            }
        },
        result={
            'archive': {
                'metadata': {
                    'domain': 'dft',
                    'embargo_length': 0,
                    'entry_create_time': '2024-05-28T19:14:10.754059+00:00',
                    'entry_hash': 'dummy_hash_id_03',
                    'entry_id': 'id_03',
                    'entry_references': '__INTERNAL__:../uploads/id_published_with_ref/archive/id_03#/metadata/entry_references',
                    'license': 'CC BY 4.0',
                    'main_author': '00000000-0000-0000-0000-000000000001',
                    'mainfile': 'mainfile_for_id_03',
                    'n_quantities': 66,
                    'parser_name': 'parsers/vasp',
                    'processed': True,
                    'published': False,
                    'quantities': '__INTERNAL__:../uploads/id_published_with_ref/archive/id_03#/metadata/quantities',
                    'section_defs': '__INTERNAL__:../uploads/id_published_with_ref/archive/id_03#/metadata/section_defs',
                    'sections': '__INTERNAL__:../uploads/id_published_with_ref/archive/id_03#/metadata/sections',
                    'text_search_contents': '__INTERNAL__:../uploads/id_published_with_ref/archive/id_03#/metadata/text_search_contents',
                    'upload_create_time': '2024-05-28T19:14:10.749059+00:00',
                    'upload_id': 'id_published_with_ref',
                    'upload_name': 'name_published',
                    'with_embargo': False,
                }
            }
        },
    )
    if simulationworkflowschema is not None:
        __entry_print(
            'entry reader to definition reader',
            {
                Token.ARCHIVE: {
                    'workflow2': {
                        'm_def': {
                            'm_request': {
                                'directive': 'plain',
                            }
                        }
                    }
                }
            },
            result={
                'metainfo': {
                    'simulationworkflowschema.general': {
                        'section_definitions': [
                            None,
                            None,
                            {
                                'name': 'SimulationWorkflow',
                                'base_sections': [
                                    'metainfo/nomad.datamodel.metainfo.workflow/section_definitions/3'
                                ],
                                'sub_sections': [
                                    {
                                        'name': 'method',
                                        'sub_section': 'metainfo/simulationworkflowschema.general/section_definitions/0',
                                    },
                                    {
                                        'name': 'results',
                                        'categories': ['/category_definitions/0'],
                                        'sub_section': 'metainfo/simulationworkflowschema.general/section_definitions/1',
                                    },
                                ],
                            },
                        ]
                    }
                },
                'archive': {
                    'workflow2': {
                        'm_def': {
                            'm_def': 'metainfo/simulationworkflowschema.general/section_definitions/2'
                        }
                    }
                },
            },
        )
        __entry_print(
            'plain without rewriting references',
            {
                Token.ARCHIVE: {
                    'workflow2': {
                        'm_request': {
                            'directive': 'plain',
                        }
                    }
                }
            },
            result={
                'archive': {
                    'workflow2': {
                        'tasks': [
                            {
                                'm_def': 'nomad.datamodel.metainfo.workflow.TaskReference',
                                'task': '../entries/id_01/archive#/workflow2',
                            }
                        ]
                    }
                }
            },
        )
        __entry_print(
            'plain with rewriting references',
            {
                Token.ARCHIVE: {
                    'workflow2': {
                        'm_request': {
                            'directive': 'plain',
                            'always_rewrite_references': True,
                        }
                    }
                }
            },
            result={
                'archive': {
                    'workflow2': {
                        'tasks': [
                            {
                                'task': 'uploads/id_published_with_ref/entries/id_01/archive/workflow2'
                            }
                        ]
                    }
                }
            },
        )

    __entry_print(
        'go to upload, resolve explicitly',
        {
            'm_request': {
                'directive': 'plain',
            },
            'upload_id': {
                'm_request': {'directive': 'resolved', 'resolve_type': 'upload'}
            },
        },
        result={
            'process_running': False,
            'current_process': None,
            'process_status': 'SUCCESS',
            'last_status_message': None,
            'errors': [],
            'warnings': [],
            'complete_time': None,
            'entry_id': 'id_03',
            'entry_create_time': '2023-03-05T22:16:52.440000',
            'mainfile_path': 'mainfile_for_id_03',
            'mainfile_key': None,
            'parser_name': 'parsers/vasp',
            'upload_id': {
                'process_running': False,
                'current_process': 'process_upload',
                'process_status': 'SUCCESS',
                'last_status_message': None,
                'errors': [],
                'warnings': [],
                'complete_time': '2023-03-05T22:16:52.436000',
                'upload_id': 'id_published_with_ref',
                'upload_name': 'name_published',
                'upload_create_time': '2023-03-05T22:16:52.435000',
                'main_author': user_dict,
                'coauthors': [],
                'reviewers': [],
                'coauthor_groups': [],
                'reviewer_groups': [],
                'writers': [user_dict],
                'viewers': [user_dict],
                'writer_groups': [],
                'viewer_groups': [],
                'published': False,
                'published_to': [],
                'publish_time': None,
                'with_embargo': False,
                'embargo_length': 0,
                'processing_failed': 0,
                'processing_successful': 6,
                'license': 'CC BY 4.0',
                'n_entries': 6,
                'upload_files_server_path': 'id_published_with_ref',
            },
        },
    )

    __entry_print(
        'go to upload, resolve implicitly, resolve main author explicitly',
        {
            'm_request': {
                'directive': 'plain',
            },
            'upload_id': {
                'm_request': {
                    'directive': 'plain',
                },
                'main_author': {
                    'm_request': {'directive': 'resolved', 'resolve_type': 'user'}
                },
            },
        },
        result={
            'process_running': False,
            'current_process': None,
            'process_status': 'SUCCESS',
            'last_status_message': None,
            'errors': [],
            'warnings': [],
            'complete_time': None,
            'entry_id': 'id_03',
            'entry_create_time': '2023-03-05T22:16:52.440000',
            'mainfile_path': 'mainfile_for_id_03',
            'mainfile_key': None,
            'parser_name': 'parsers/vasp',
            'upload_id': {
                'process_running': False,
                'current_process': 'process_upload',
                'process_status': 'SUCCESS',
                'last_status_message': None,
                'errors': [],
                'warnings': [],
                'complete_time': '2023-03-05T22:16:52.436000',
                'upload_id': 'id_published_with_ref',
                'upload_name': 'name_published',
                'upload_create_time': '2023-03-05T22:16:52.435000',
                'coauthors': [],
                'reviewers': [],
                'coauthor_groups': [],
                'reviewer_groups': [],
                'writers': [user_dict],
                'viewers': [user_dict],
                'writer_groups': [],
                'viewer_groups': [],
                'published': False,
                'published_to': [],
                'publish_time': None,
                'with_embargo': False,
                'embargo_length': 0,
                'license': 'CC BY 4.0',
                'n_entries': 6,
                'processing_failed': 0,
                'processing_successful': 6,
                'upload_files_server_path': 'id_published_with_ref',
                'main_author': user_dict,
            },
        },
    )

    def __fs_print(msg, required, *, result: dict = None):
        with FileSystemReader(required, user=user1) as reader:
            if result:
                assert_dict(reader.sync_read('id_published_with_ref'), result)
            else:
                rprint(f'\n\nExample: {next(counter)} -> {msg}:')
                rprint(required)
                rprint('output:')
                rprint(reader.sync_read('id_published_with_ref'))

    __fs_print(
        'plain file system reader',
        {
            'm_request': {
                'directive': 'plain',
            },
        },
        result={
            'm_is': 'Directory',
            '1.aux': {'path': '1.aux', 'm_is': 'File', 'size': 8},
            '2.aux': {'path': '2.aux', 'm_is': 'File', 'size': 8},
            '3.aux': {'path': '3.aux', 'm_is': 'File', 'size': 8},
            '4.aux': {'path': '4.aux', 'm_is': 'File', 'size': 8},
            'mainfile_for_id_01': {
                'path': 'mainfile_for_id_01',
                'm_is': 'File',
                'size': 3227,
            },
            'mainfile_for_id_02': {
                'path': 'mainfile_for_id_02',
                'm_is': 'File',
                'size': 3227,
            },
            'mainfile_for_id_03': {
                'path': 'mainfile_for_id_03',
                'm_is': 'File',
                'size': 3227,
            },
            'mainfile_for_id_04': {
                'path': 'mainfile_for_id_04',
                'm_is': 'File',
                'size': 3227,
            },
            'mainfile_for_id_05': {
                'path': 'mainfile_for_id_05',
                'm_is': 'File',
                'size': 3227,
            },
            'mainfile_for_id_06': {
                'path': 'mainfile_for_id_06',
                'm_is': 'File',
                'size': 3227,
            },
        },
    )

    __fs_print(
        'go to entry',
        {
            'm_request': {
                'directive': 'resolved',
            },
        },
        result={
            'm_is': 'Directory',
            '1.aux': {'path': '1.aux', 'm_is': 'File', 'size': 8},
            '2.aux': {'path': '2.aux', 'm_is': 'File', 'size': 8},
            '3.aux': {'path': '3.aux', 'm_is': 'File', 'size': 8},
            '4.aux': {'path': '4.aux', 'm_is': 'File', 'size': 8},
            'mainfile_for_id_01': {
                'path': 'mainfile_for_id_01',
                'm_is': 'File',
                'size': 3227,
                Token.ENTRY: {
                    'process_running': False,
                    'current_process': None,
                    'process_status': 'SUCCESS',
                    'last_status_message': None,
                    'errors': [],
                    'warnings': [],
                    'complete_time': None,
                    'entry_id': 'id_01',
                    'entry_create_time': '2023-03-05T22:16:52.438000',
                    'mainfile_path': 'mainfile_for_id_01',
                    'mainfile_key': None,
                    'upload_id': 'id_published_with_ref',
                    'parser_name': 'parsers/vasp',
                },
            },
            'mainfile_for_id_02': {
                'path': 'mainfile_for_id_02',
                'm_is': 'File',
                'size': 3227,
                Token.ENTRY: {
                    'process_running': False,
                    'current_process': None,
                    'process_status': 'SUCCESS',
                    'last_status_message': None,
                    'errors': [],
                    'warnings': [],
                    'complete_time': None,
                    'entry_id': 'id_02',
                    'entry_create_time': '2023-03-05T22:16:52.439000',
                    'mainfile_path': 'mainfile_for_id_02',
                    'mainfile_key': None,
                    'upload_id': 'id_published_with_ref',
                    'parser_name': 'parsers/vasp',
                },
            },
            'mainfile_for_id_03': {
                'path': 'mainfile_for_id_03',
                'm_is': 'File',
                'size': 3227,
                Token.ENTRY: {
                    'process_running': False,
                    'current_process': None,
                    'process_status': 'SUCCESS',
                    'last_status_message': None,
                    'errors': [],
                    'warnings': [],
                    'complete_time': None,
                    'entry_id': 'id_03',
                    'entry_create_time': '2023-03-05T22:16:52.440000',
                    'mainfile_path': 'mainfile_for_id_03',
                    'mainfile_key': None,
                    'upload_id': 'id_published_with_ref',
                    'parser_name': 'parsers/vasp',
                },
            },
            'mainfile_for_id_04': {
                'path': 'mainfile_for_id_04',
                'm_is': 'File',
                'size': 3227,
                Token.ENTRY: {
                    'process_running': False,
                    'current_process': None,
                    'process_status': 'SUCCESS',
                    'last_status_message': None,
                    'errors': [],
                    'warnings': [],
                    'complete_time': None,
                    'entry_id': 'id_04',
                    'entry_create_time': '2023-03-05T22:16:52.441000',
                    'mainfile_path': 'mainfile_for_id_04',
                    'mainfile_key': None,
                    'upload_id': 'id_published_with_ref',
                    'parser_name': 'parsers/vasp',
                },
            },
            'mainfile_for_id_05': {
                'path': 'mainfile_for_id_05',
                'm_is': 'File',
                'size': 3227,
                Token.ENTRY: {
                    'process_running': False,
                    'current_process': None,
                    'process_status': 'SUCCESS',
                    'last_status_message': None,
                    'errors': [],
                    'warnings': [],
                    'complete_time': None,
                    'entry_id': 'id_05',
                    'entry_create_time': '2023-03-05T22:16:52.442000',
                    'mainfile_path': 'mainfile_for_id_05',
                    'mainfile_key': None,
                    'upload_id': 'id_published_with_ref',
                    'parser_name': 'parsers/vasp',
                },
            },
            'mainfile_for_id_06': {
                'path': 'mainfile_for_id_06',
                'm_is': 'File',
                'size': 3227,
                Token.ENTRY: {
                    'process_running': False,
                    'current_process': None,
                    'process_status': 'SUCCESS',
                    'last_status_message': None,
                    'errors': [],
                    'warnings': [],
                    'complete_time': None,
                    'entry_id': 'id_06',
                    'entry_create_time': '2023-03-05T22:16:52.443000',
                    'mainfile_path': 'mainfile_for_id_06',
                    'mainfile_key': None,
                    'upload_id': 'id_published_with_ref',
                    'parser_name': 'parsers/vasp',
                },
            },
        },
    )

    __fs_print(
        'go to selected entry',
        {
            'm_request': {
                'directive': 'plain',
            },
            'mainfile_for_id_01': {
                'm_request': {
                    'directive': 'resolved',
                },
            },
        },
        result={
            'm_is': 'Directory',
            '1.aux': {'path': '1.aux', 'm_is': 'File', 'size': 8},
            '2.aux': {'path': '2.aux', 'm_is': 'File', 'size': 8},
            '3.aux': {'path': '3.aux', 'm_is': 'File', 'size': 8},
            '4.aux': {'path': '4.aux', 'm_is': 'File', 'size': 8},
            'mainfile_for_id_02': {
                'path': 'mainfile_for_id_02',
                'm_is': 'File',
                'size': 3227,
            },
            'mainfile_for_id_03': {
                'path': 'mainfile_for_id_03',
                'm_is': 'File',
                'size': 3227,
            },
            'mainfile_for_id_04': {
                'path': 'mainfile_for_id_04',
                'm_is': 'File',
                'size': 3227,
            },
            'mainfile_for_id_05': {
                'path': 'mainfile_for_id_05',
                'm_is': 'File',
                'size': 3227,
            },
            'mainfile_for_id_06': {
                'path': 'mainfile_for_id_06',
                'm_is': 'File',
                'size': 3227,
            },
            'mainfile_for_id_01': {
                'path': 'mainfile_for_id_01',
                'm_is': 'File',
                'size': 3227,
                Token.ENTRY: {
                    'process_running': False,
                    'current_process': None,
                    'process_status': 'SUCCESS',
                    'last_status_message': None,
                    'errors': [],
                    'warnings': [],
                    'complete_time': None,
                    'entry_id': 'id_01',
                    'entry_create_time': '2023-03-05T22:16:52.438000',
                    'mainfile_path': 'mainfile_for_id_01',
                    'mainfile_key': None,
                    'upload_id': 'id_published_with_ref',
                    'parser_name': 'parsers/vasp',
                },
            },
        },
    )

    __fs_print(
        'go to selected entry then to upload',
        {
            'm_request': {
                'directive': 'plain',
            },
            'mainfile_for_id_01': {
                'm_request': {
                    'directive': 'plain',
                },
                Token.ENTRY: {
                    'upload_id': {
                        'm_request': {
                            'directive': 'resolved',
                            'resolve_type': 'upload',
                        },
                    }
                },
            },
        },
        result={
            'm_is': 'Directory',
            '1.aux': {'path': '1.aux', 'm_is': 'File', 'size': 8},
            '2.aux': {'path': '2.aux', 'm_is': 'File', 'size': 8},
            '3.aux': {'path': '3.aux', 'm_is': 'File', 'size': 8},
            '4.aux': {'path': '4.aux', 'm_is': 'File', 'size': 8},
            'mainfile_for_id_02': {
                'path': 'mainfile_for_id_02',
                'm_is': 'File',
                'size': 3227,
            },
            'mainfile_for_id_03': {
                'path': 'mainfile_for_id_03',
                'm_is': 'File',
                'size': 3227,
            },
            'mainfile_for_id_04': {
                'path': 'mainfile_for_id_04',
                'm_is': 'File',
                'size': 3227,
            },
            'mainfile_for_id_05': {
                'path': 'mainfile_for_id_05',
                'm_is': 'File',
                'size': 3227,
            },
            'mainfile_for_id_06': {
                'path': 'mainfile_for_id_06',
                'm_is': 'File',
                'size': 3227,
            },
            'mainfile_for_id_01': {
                'path': 'mainfile_for_id_01',
                'm_is': 'File',
                'size': 3227,
                Token.ENTRY: {
                    'upload_id': {
                        'process_running': False,
                        'current_process': 'process_upload',
                        'process_status': 'SUCCESS',
                        'last_status_message': None,
                        'errors': [],
                        'warnings': [],
                        'complete_time': '2023-03-05T22:16:52.436000',
                        'upload_id': 'id_published_with_ref',
                        'upload_name': 'name_published',
                        'upload_create_time': '2023-03-05T22:16:52.435000',
                        'main_author': user_dict,
                        'coauthors': [],
                        'reviewers': [],
                        'coauthor_groups': [],
                        'reviewer_groups': [],
                        'writers': [user_dict],
                        'viewers': [user_dict],
                        'writer_groups': [],
                        'viewer_groups': [],
                        'published': False,
                        'published_to': [],
                        'publish_time': None,
                        'with_embargo': False,
                        'embargo_length': 0,
                        'license': 'CC BY 4.0',
                        'processing_failed': 0,
                        'processing_successful': 6,
                        'n_entries': 6,
                        'upload_files_server_path': 'id_published_with_ref',
                    }
                },
            },
        },
    )

    __fs_print(
        'go to selected entry then to upload, skipping file info',
        {
            'm_request': {
                'directive': 'plain',
            },
            'mainfile_for_id_01': {
                Token.ENTRY: {
                    'upload_id': {
                        'm_request': {
                            'directive': 'resolved',
                            'resolve_type': 'upload',
                        },
                    }
                }
            },
        },
        result={
            'm_is': 'Directory',
            '1.aux': {'path': '1.aux', 'm_is': 'File', 'size': 8},
            '2.aux': {'path': '2.aux', 'm_is': 'File', 'size': 8},
            '3.aux': {'path': '3.aux', 'm_is': 'File', 'size': 8},
            '4.aux': {'path': '4.aux', 'm_is': 'File', 'size': 8},
            'mainfile_for_id_02': {
                'path': 'mainfile_for_id_02',
                'm_is': 'File',
                'size': 3227,
            },
            'mainfile_for_id_03': {
                'path': 'mainfile_for_id_03',
                'm_is': 'File',
                'size': 3227,
            },
            'mainfile_for_id_04': {
                'path': 'mainfile_for_id_04',
                'm_is': 'File',
                'size': 3227,
            },
            'mainfile_for_id_05': {
                'path': 'mainfile_for_id_05',
                'm_is': 'File',
                'size': 3227,
            },
            'mainfile_for_id_06': {
                'path': 'mainfile_for_id_06',
                'm_is': 'File',
                'size': 3227,
            },
            'mainfile_for_id_01': {
                'path': 'mainfile_for_id_01',
                'm_is': 'File',
                'size': 3227,
                Token.ENTRY: {
                    'upload_id': {
                        'process_running': False,
                        'current_process': 'process_upload',
                        'process_status': 'SUCCESS',
                        'last_status_message': None,
                        'errors': [],
                        'warnings': [],
                        'complete_time': '2023-03-05T22:16:52.436000',
                        'upload_id': 'id_published_with_ref',
                        'upload_name': 'name_published',
                        'upload_create_time': '2023-03-05T22:16:52.435000',
                        'main_author': user_dict,
                        'coauthors': [],
                        'reviewers': [],
                        'coauthor_groups': [],
                        'reviewer_groups': [],
                        'writers': [user_dict],
                        'viewers': [user_dict],
                        'writer_groups': [],
                        'viewer_groups': [],
                        'published': False,
                        'published_to': [],
                        'publish_time': None,
                        'with_embargo': False,
                        'embargo_length': 0,
                        'license': 'CC BY 4.0',
                        'n_entries': 6,
                        'processing_failed': 0,
                        'processing_successful': 6,
                        'upload_files_server_path': 'id_published_with_ref',
                    }
                },
            },
        },
    )

    __fs_print(
        'go to selected entry then to upload, resolve user',
        {'mainfile_for_id_01': {Token.ENTRY: {'upload_id': {'main_author': '*'}}}},
        result={
            'm_is': 'Directory',
            'mainfile_for_id_01': {
                Token.ENTRY: {
                    'upload_id': {
                        'main_author': {
                            'name': 'Sheldon Cooper',
                            'first_name': 'Sheldon',
                            'last_name': 'Cooper',
                            'email': 'sheldon.cooper@nomad-coe.eu',
                            'user_id': '00000000-0000-0000-0000-000000000001',
                            'username': 'scooper',
                            'is_admin': False,
                            'is_oasis_admin': True,
                        }
                    }
                }
            },
        },
    )


# noinspection DuplicatedCode,SpellCheckingInspection
def test_group_reader(groups_function, user1):
    def __ge_print(msg, required, *, to_file: bool = False, result: dict = None):
        with MongoReader(required, user=user1) as reader:
            if result:
                assert_dict(reader.sync_read(), result)
            else:
                rprint(f'\n\nExample: {next(counter)} -> {msg}:')
                rprint(required)
                if not to_file:
                    rprint('output:')
                    rprint(reader.sync_read())
                else:
                    with open('archive_reader_test.json', 'w') as f:
                        f.write(json.dumps(reader.sync_read()))

    __ge_print(
        'general start from group; id: *',
        {
            Token.GROUP: {
                'GGGGGGGGGGGGGGGGGGGG14': '*',
            }
        },
        result={
            'group': {
                'GGGGGGGGGGGGGGGGGGGG14': {
                    'group_id': 'GGGGGGGGGGGGGGGGGGGG14',
                    'group_name': 'Group 14',
                    'owner': {
                        'name': 'Sheldon Cooper',
                        'first_name': 'Sheldon',
                        'last_name': 'Cooper',
                        'email': 'sheldon.cooper@nomad-coe.eu',
                        'user_id': '00000000-0000-0000-0000-000000000001',
                        'username': 'scooper',
                        'is_admin': False,
                        'is_oasis_admin': True,
                    },
                    'members': ListWithSortKey(
                        (
                            {
                                'name': 'Rajesh Koothrappali',
                                'first_name': 'Rajesh',
                                'last_name': 'Koothrappali',
                                'email': 'rajesh.koothrappali@nomad-fairdi.tests.de',
                                'user_id': '00000000-0000-0000-0000-000000000004',
                                'username': 'rkoothrappali',
                                'is_admin': False,
                            },
                            {
                                'name': 'Sheldon Cooper',
                                'first_name': 'Sheldon',
                                'last_name': 'Cooper',
                                'email': 'sheldon.cooper@nomad-coe.eu',
                                'user_id': '00000000-0000-0000-0000-000000000001',
                                'username': 'scooper',
                                'is_admin': False,
                                'is_oasis_admin': True,
                            },
                        ),
                        sort_key=lambda x: x['user_id'],
                    ),
                }
            }
        },
    )
    __ge_print(
        'general start from group; id: owner.email',
        {
            Token.GROUP: {
                'GGGGGGGGGGGGGGGGGGGG14': {
                    'owner': {
                        'email': '*',
                    },
                },
            }
        },
        result={
            'group': {
                'GGGGGGGGGGGGGGGGGGGG14': {
                    'owner': {'email': 'sheldon.cooper@nomad-coe.eu'}
                }
            }
        },
    )
    __ge_print(
        'general start from group; *: group_name',
        {
            Token.GROUP: {
                'm_request': {
                    'pagination': {'page_size': 20},
                },
                '*': {'group_name': '*'},
            }
        },
        result={
            'group': {
                'GGGGGGGGGGGGGGGGGGGGG0': {'group_name': 'Group 0'},
                'GGGGGGGGGGGGGGGGGGGGG1': {'group_name': 'Group 1'},
                'GGGGGGGGGGGGGGGGGGGGG2': {'group_name': 'Group 2'},
                'GGGGGGGGGGGGGGGGGGGGG3': {'group_name': 'Group 3'},
                'GGGGGGGGGGGGGGGGGGGGG6': {'group_name': 'Group 6'},
                'GGGGGGGGGGGGGGGGGGGGG8': {'group_name': 'Group 8'},
                'GGGGGGGGGGGGGGGGGGGGG9': {'group_name': 'Group 9'},
                'GGGGGGGGGGGGGGGGGGGG14': {'group_name': 'Group 14'},
                'GGGGGGGGGGGGGGGGGGGG15': {'group_name': 'Group 15'},
                'GGGGGGGGGGGGGGGGGGGG18': {'group_name': 'Group 18'},
                'GGGGGGGGGGGGGGGGGGGG19': {'group_name': 'Group 19'},
                'GGGGGGGGGGGGGGGGGGG123': {'group_name': 'Group 123'},
                'GGGGGGGGGGGGGGGGGGUniq': {'group_name': 'Group Uniq'},
                'GGGGGGGGGGGGGGTwin One': {'group_name': 'Group Twin One'},
                'GGGGGGGGGGGGGGTwin Two': {'group_name': 'Group Twin Two'},
                'GGGGGGGGGOne Two Three': {'group_name': 'Group One Two Three'},
            }
        },
    )
    __ge_print(
        'general start from group; query: [group_id]; *: group_name',
        {
            Token.GROUP: {
                'm_request': {'query': {'group_id': ['GGGGGGGGGGGGGGGGGGG123']}},
                '*': {'group_name': '*'},
            }
        },
        result={
            'group': {
                'GGGGGGGGGGGGGGGGGGG123': {'group_name': 'Group 123'},
            }
        },
    )
    __ge_print(
        'general start from group; query: user_id; *: group_name',
        {
            Token.GROUP: {
                'm_request': {
                    'query': {'user_id': '00000000-0000-0000-0000-000000000008'}
                },
                '*': {'group_name': '*'},
            }
        },
        result={
            'group': {
                'GGGGGGGGGGGGGGGGGGGGG8': {'group_name': 'Group 8'},
                'GGGGGGGGGGGGGGGGGGGG18': {'group_name': 'Group 18'},
            }
        },
    )
    __ge_print(
        'general start from group; query: search_terms; *: group_name',
        {
            Token.GROUP: {
                'm_request': {'query': {'search_terms': 'win'}},
                '*': {'group_name': '*'},
            }
        },
        result={
            'group': {
                'GGGGGGGGGGGGGGTwin One': {'group_name': 'Group Twin One'},
                'GGGGGGGGGGGGGGTwin Two': {'group_name': 'Group Twin Two'},
            }
        },
    )
    __ge_print(
        'general start from group; query: user_id, search_terms; *: group_name',
        {
            Token.GROUP: {
                'm_request': {
                    'query': {
                        'user_id': '00000000-0000-0000-0000-000000000008',
                        'search_terms': '1',
                    }
                },
                '*': {'group_name': '*'},
            }
        },
        result={
            'group': {
                'GGGGGGGGGGGGGGGGGGGG18': {'group_name': 'Group 18'},
            }
        },
    )


# noinspection DuplicatedCode,SpellCheckingInspection
def test_general_reader(json_dict, example_data_with_reference, user1):
    def __ge_print(msg, required, *, to_file: bool = False, result: dict = None):
        with MongoReader(required, user=user1) as reader:
            if result:
                assert_dict(reader.sync_read(), result)
            else:
                rprint(f'\n\nExample: {next(counter)} -> {msg}:')
                rprint(required)
                if not to_file:
                    rprint('output:')
                    rprint(reader.sync_read())
                else:
                    with open('archive_reader_test.json', 'w') as f:
                        f.write(json.dumps(reader.sync_read()))

    __ge_print(
        'general start from entry',
        {
            Token.ENTRIES: {
                'm_request': {
                    'directive': 'resolved',
                    'resolve_type': 'entry',
                    'pagination': {'page_size': 2, 'page_after_value': 'id_03'},
                },
            }
        },
        result={
            Token.ENTRIES: {
                'id_04': {
                    'process_running': False,
                    'current_process': None,
                    'process_status': 'SUCCESS',
                    'last_status_message': None,
                    'errors': [],
                    'warnings': [],
                    'complete_time': None,
                    'entry_id': 'id_04',
                    'entry_create_time': '2023-03-05T22:29:55.842000',
                    'mainfile_path': 'mainfile_for_id_04',
                    'mainfile_key': None,
                    'upload_id': 'id_published_with_ref',
                    'parser_name': 'parsers/vasp',
                },
                'id_05': {
                    'process_running': False,
                    'current_process': None,
                    'process_status': 'SUCCESS',
                    'last_status_message': None,
                    'errors': [],
                    'warnings': [],
                    'complete_time': None,
                    'entry_id': 'id_05',
                    'entry_create_time': '2023-03-05T22:29:55.843000',
                    'mainfile_path': 'mainfile_for_id_05',
                    'mainfile_key': None,
                    'upload_id': 'id_published_with_ref',
                    'parser_name': 'parsers/vasp',
                },
            }
        },
    )
    __ge_print(
        'general start from entry with wildcard',
        {
            Token.ENTRIES: {
                'id_01': {
                    'm_request': {
                        'directive': 'plain',
                    },
                },
                '*': {
                    'm_request': {
                        'directive': 'plain',
                        'include': ['entry_id', 'mainfile_path'],
                    },
                    'upload_id': '*',
                },
            }
        },
        result={
            Token.ENTRIES: {
                'id_01': {
                    'process_running': False,
                    'current_process': None,
                    'process_status': 'SUCCESS',
                    'last_status_message': None,
                    'errors': [],
                    'warnings': [],
                    'complete_time': None,
                    'entry_id': 'id_01',
                    'entry_create_time': '2023-03-05T22:20:46.586000',
                    'mainfile_path': 'mainfile_for_id_01',
                    'mainfile_key': None,
                    'upload_id': 'id_published_with_ref',
                    'parser_name': 'parsers/vasp',
                },
                'id_02': {
                    'entry_id': 'id_02',
                    'mainfile_path': 'mainfile_for_id_02',
                    'upload_id': 'id_published_with_ref',
                },
                'id_03': {
                    'entry_id': 'id_03',
                    'mainfile_path': 'mainfile_for_id_03',
                    'upload_id': 'id_published_with_ref',
                },
                'id_04': {
                    'entry_id': 'id_04',
                    'mainfile_path': 'mainfile_for_id_04',
                    'upload_id': 'id_published_with_ref',
                },
                'id_05': {
                    'entry_id': 'id_05',
                    'mainfile_path': 'mainfile_for_id_05',
                    'upload_id': 'id_published_with_ref',
                },
                'id_06': {
                    'entry_id': 'id_06',
                    'mainfile_path': 'mainfile_for_id_06',
                    'upload_id': 'id_published_with_ref',
                },
            }
        },
    )
    __ge_print(
        'general start from entry to metadata',
        {
            Token.ENTRIES: {
                'id_01': {
                    'metadata': {
                        'results': {
                            'm_request': {
                                'directive': 'plain',
                            }
                        }
                    }
                }
            }
        },
        result={
            'entries': {
                'id_01': {
                    'metadata': {
                        'results': {
                            'material': {
                                'dimensionality': '3D',
                                'symmetry': {'crystal_system': 'cubic'},
                                'elements': ['H', 'O'],
                                'elements_exclusive': 'H O',
                                'material_id': 'test_material_id',
                                'structural_type': 'not processed',
                                'n_elements': 2,
                            },
                            'method': {
                                'simulation': {
                                    'program_version': 'not processed',
                                    'program_version_internal': 'not processed',
                                    'dft': {
                                        'basis_set_type': 'unavailable',
                                        'core_electron_treatment': 'unavailable',
                                        'xc_functional_type': 'GGA',
                                        'xc_functional_names': [],
                                        'jacobs_ladder': 'not processed',
                                    },
                                    'program_name': 'VASP',
                                }
                            },
                            'properties': {
                                'available_properties': ['dos_electronic'],
                                'n_calculations': 1,
                                'electronic': {
                                    'dos_electronic': [
                                        {
                                            'spin_polarized': False,
                                            'band_gap': [{'type': 'indirect'}],
                                        }
                                    ]
                                },
                            },
                        }
                    }
                }
            }
        },
    )
    # only check if those keys work
    # the result order is not checked
    for order_by in ['entry_create_time', 'mainfile_path']:
        __ge_print(
            'general start from entry WITHOUT retrieval of metadata (just listing)',
            {
                Token.SEARCH: {
                    'm_request': {
                        'directive': 'plain',
                        'pagination': {'page_size': 2, 'page': 2, 'order_by': order_by},
                        'query': {'owner': 'user'},
                    },
                }
            },
            result={
                Token.SEARCH: {
                    'id_03': 'id_03',
                    'id_04': 'id_04',
                }
            },
        )
    __ge_print(
        'general start from upload',
        {
            Token.UPLOADS: {
                'm_request': {
                    'directive': 'resolved',
                    'resolve_type': 'upload',
                    'pagination': {
                        'page_size': 10,
                    },
                },
            }
        },
        result={
            Token.UPLOADS: {
                'id_published_with_ref': {
                    'process_running': False,
                    'current_process': 'process_upload',
                    'process_status': 'SUCCESS',
                    'last_status_message': None,
                    'errors': [],
                    'warnings': [],
                    'complete_time': '2023-03-05T22:20:46.584000',
                    'upload_id': 'id_published_with_ref',
                    'upload_name': 'name_published',
                    'upload_create_time': '2023-03-05T22:20:46.583000',
                    'main_author': user_dict,
                    'coauthors': [],
                    'reviewers': [],
                    'coauthor_groups': [],
                    'reviewer_groups': [],
                    'writers': [user_dict],
                    'viewers': [user_dict],
                    'writer_groups': [],
                    'viewer_groups': [],
                    'published': False,
                    'published_to': [],
                    'publish_time': None,
                    'with_embargo': False,
                    'embargo_length': 0,
                    'license': 'CC BY 4.0',
                    'n_entries': 6,
                    'processing_failed': 0,
                    'processing_successful': 6,
                    'upload_files_server_path': 'id_published_with_ref',
                },
            }
        },
    )
    # only check if those keys work
    # the result order is not checked
    for order_by in ['entry_create_time', 'mainfile_path']:
        __ge_print(
            'general start from entry with query and pagination',
            {
                Token.ENTRIES: {
                    'm_request': {
                        'directive': 'plain',
                        'pagination': {'page_size': 10, 'order_by': order_by},
                    },
                }
            },
            result={
                'entries': {
                    'id_01': 'id_01',
                    'id_02': 'id_02',
                    'id_03': 'id_03',
                    'id_04': 'id_04',
                    'id_05': 'id_05',
                    'id_06': 'id_06',
                }
            },
        )
    __ge_print(
        'general start from upload with query and pagination',
        {
            Token.UPLOADS: {
                'm_request': {
                    'directive': 'resolved',
                    'resolve_type': 'upload',
                    'pagination': {'page_size': 10, 'order_by': 'upload_create_time'},
                    'query': {'is_processing': False},
                },
            }
        },
        result={
            Token.UPLOADS: {
                'id_published_with_ref': {
                    'process_running': False,
                    'current_process': 'process_upload',
                    'process_status': 'SUCCESS',
                    'last_status_message': None,
                    'errors': [],
                    'warnings': [],
                    'complete_time': '2023-03-05T22:20:46.584000',
                    'upload_id': 'id_published_with_ref',
                    'upload_name': 'name_published',
                    'upload_create_time': '2023-03-05T22:20:46.583000',
                    'main_author': user_dict,
                    'coauthors': [],
                    'reviewers': [],
                    'coauthor_groups': [],
                    'reviewer_groups': [],
                    'writers': [user_dict],
                    'viewers': [user_dict],
                    'writer_groups': [],
                    'viewer_groups': [],
                    'published': False,
                    'published_to': [],
                    'publish_time': None,
                    'with_embargo': False,
                    'embargo_length': 0,
                    'license': 'CC BY 4.0',
                    'n_entries': 6,
                    'processing_failed': 0,
                    'processing_successful': 6,
                    'upload_files_server_path': 'id_published_with_ref',
                },
            }
        },
    )
    __ge_print(
        'general start from user, does NOT perform search from security',
        {
            Token.USER: {
                'me': {
                    'm_request': {'directive': 'plain'},
                }
            }
        },
        result={
            Token.USER: {
                'me': {
                    'name': 'Sheldon Cooper',
                    'first_name': 'Sheldon',
                    'last_name': 'Cooper',
                    'email': 'sheldon.cooper@nomad-coe.eu',
                    'user_id': '00000000-0000-0000-0000-000000000001',
                    'username': 'scooper',
                    'is_admin': False,
                    'is_oasis_admin': True,
                }
            }
        },
    )
    __ge_print(
        'general start from me, with its user id like fields resolved',
        {
            Token.USER: {
                'me': {
                    'm_request': {'directive': 'resolved', 'resolve_type': 'user'},
                }
            }
        },
        result={
            Token.USER: {
                'me': {
                    'name': 'Sheldon Cooper',
                    'first_name': 'Sheldon',
                    'last_name': 'Cooper',
                    'email': 'sheldon.cooper@nomad-coe.eu',
                    'user_id': {
                        'name': 'Sheldon Cooper',
                        'first_name': 'Sheldon',
                        'last_name': 'Cooper',
                        'email': 'sheldon.cooper@nomad-coe.eu',
                        'user_id': '00000000-0000-0000-0000-000000000001',
                        'username': 'scooper',
                        'is_admin': False,
                        'is_oasis_admin': True,
                    },
                    'username': 'scooper',
                    'is_admin': False,
                    'is_oasis_admin': True,
                }
            }
        },
    )


# noinspection DuplicatedCode,SpellCheckingInspection
def test_metainfo_reader(mongo_infra, user1):
    def __ge_print(msg, required, *, to_file: bool = False, result: dict = None):
        with MongoReader(required, user=user1) as reader:
            if result:
                assert_dict(reader.sync_read(), result)
            else:
                rprint(f'\n\nExample: {next(counter)} -> {msg}:')
                rprint(required)
                if not to_file:
                    rprint('output:')
                    rprint(reader.sync_read())
                else:
                    with open('archive_reader_test.json', 'w') as f:
                        f.write(json.dumps(reader.sync_read()))

    __ge_print(
        'general start from metainfo',
        {
            Token.METAINFO: {
                'nomad.datamodel.metainfo.simulation.run': {
                    'section_definitions[2]': {
                        'm_request': {'directive': 'plain'},
                    }
                }
            }
        },
        result={
            'metainfo': {
                'nomad.datamodel.metainfo.simulation.run': {
                    'section_definitions': [
                        None,
                        None,
                        {
                            'name': 'MessageRun',
                            'description': 'Contains warning, error, and info messages of the run.',
                            'quantities': [
                                {
                                    'name': 'type',
                                    'description': 'Type of the message. Can be one of warning, error, info, debug.',
                                    'type': {'type_kind': 'python', 'type_data': 'str'},
                                    'shape': [],
                                },
                                {
                                    'name': 'value',
                                    'description': 'Value of the message of the computational program, given by type.',
                                    'type': {'type_kind': 'python', 'type_data': 'str'},
                                    'shape': [],
                                },
                            ],
                        },
                    ]
                }
            },
        },
    )

    __ge_print(
        'general start from metainfo',
        {
            Token.METAINFO: {
                'm_request': {
                    'include': ['*nomad.datamodel.metainfo.simulation.run'],
                    'pagination': {'page_size': 50},
                },
                '*': {'m_request': {'index': [2]}},
            }
        },
        result={
            'metainfo': {
                'nomad.datamodel.metainfo.simulation.run': {
                    'name': 'nomad.datamodel.metainfo.simulation.run',
                    'category_definitions': [
                        {
                            'description': 'Information that *in theory* should not affect the results of the calculations (e.g.,\ntiming).',
                            'name': 'AccessoryInfo',
                        },
                        {
                            'categories': ['/category_definitions/0'],
                            'description': 'Contains information on the program that generated the data, i.e. the program_name,\nprogram_version, program_compilation_host and program_compilation_datetime as direct\nchildren of this field.',
                            'name': 'ProgramInfo',
                        },
                    ],
                    'section_definitions': [
                        {
                            'description': 'Contains the specifications of the program.',
                            'name': 'Program',
                            'quantities': [
                                {
                                    'categories': [
                                        '/category_definitions/0',
                                        '/category_definitions/1',
                                    ],
                                    'description': 'Contains the program compilation date and time from *Unix epoch* (00:00:00 UTC on\n1 January 1970) in seconds. For date and times without a timezone, the default\ntimezone GMT is used.',
                                    'name': 'compilation_datetime',
                                    'shape': [],
                                    'type': {
                                        'type_kind': 'numpy',
                                        'type_data': 'float64',
                                    },
                                    'unit': 'second',
                                },
                                {
                                    'categories': [
                                        '/category_definitions/0',
                                        '/category_definitions/1',
                                    ],
                                    'description': 'Specifies the host on which the program was compiled.',
                                    'name': 'compilation_host',
                                    'shape': [],
                                    'type': {'type_kind': 'python', 'type_data': 'str'},
                                },
                                {
                                    'categories': [
                                        '/category_definitions/0',
                                        '/category_definitions/1',
                                    ],
                                    'description': 'Specifies the name of the program that generated the data.',
                                    'name': 'name',
                                    'shape': [],
                                    'type': {'type_kind': 'python', 'type_data': 'str'},
                                },
                                {
                                    'categories': [
                                        '/category_definitions/0',
                                        '/category_definitions/1',
                                    ],
                                    'description': 'Specifies the official release version of the program that was used.',
                                    'name': 'version',
                                    'shape': [],
                                    'type': {'type_kind': 'python', 'type_data': 'str'},
                                },
                                {
                                    'categories': ['/category_definitions/1'],
                                    'description': 'Specifies a program version tag used internally for development purposes.\nAny kind of tagging system is supported, including git commit hashes.',
                                    'name': 'version_internal',
                                    'type': {'type_kind': 'python', 'type_data': 'str'},
                                },
                            ],
                        },
                        {
                            'description': 'Contains information on timing information of the run.',
                            'name': 'TimeRun',
                            'quantities': [
                                {
                                    'description': 'Stores the end time of the run on CPU 1.',
                                    'name': 'cpu1_end',
                                    'shape': [],
                                    'type': {
                                        'type_kind': 'numpy',
                                        'type_data': 'float64',
                                    },
                                    'unit': 'second',
                                },
                                {
                                    'description': 'Stores the start time of the run on CPU 1.',
                                    'name': 'cpu1_start',
                                    'shape': [],
                                    'type': {
                                        'type_kind': 'numpy',
                                        'type_data': 'float64',
                                    },
                                    'unit': 'second',
                                },
                                {
                                    'description': 'Stores the end date of the run as time since the *Unix epoch* (00:00:00 UTC on 1\nJanuary 1970) in seconds. For date and times without a timezone, the default\ntimezone GMT is used.',
                                    'name': 'date_end',
                                    'shape': [],
                                    'type': {
                                        'type_kind': 'numpy',
                                        'type_data': 'float64',
                                    },
                                    'unit': 'second',
                                },
                                {
                                    'description': 'Stores the start date of the run as time since the *Unix epoch* (00:00:00 UTC on 1\nJanuary 1970) in seconds. For date and times without a timezone, the default\ntimezone GMT is used.',
                                    'name': 'date_start',
                                    'shape': [],
                                    'type': {
                                        'type_kind': 'numpy',
                                        'type_data': 'float64',
                                    },
                                    'unit': 'second',
                                },
                                {
                                    'description': 'Stores the internal wall-clock time at the end of the run.',
                                    'name': 'wall_end',
                                    'shape': [],
                                    'type': {
                                        'type_kind': 'numpy',
                                        'type_data': 'float64',
                                    },
                                    'unit': 'second',
                                },
                                {
                                    'description': 'Stores the internal wall-clock time from the start of the run.',
                                    'name': 'wall_start',
                                    'shape': [],
                                    'type': {
                                        'type_kind': 'numpy',
                                        'type_data': 'float64',
                                    },
                                    'unit': 'second',
                                },
                            ],
                        },
                        {
                            'description': 'Contains warning, error, and info messages of the run.',
                            'name': 'MessageRun',
                            'quantities': [
                                {
                                    'description': 'Type of the message. Can be one of warning, error, info, debug.',
                                    'name': 'type',
                                    'shape': [],
                                    'type': {'type_kind': 'python', 'type_data': 'str'},
                                },
                                {
                                    'description': 'Value of the message of the computational program, given by type.',
                                    'name': 'value',
                                    'shape': [],
                                    'type': {'type_kind': 'python', 'type_data': 'str'},
                                },
                            ],
                        },
                        {
                            'description': 'Every section run represents a single call of a program.',
                            'name': 'Run',
                            'base_sections': [
                                'metainfo/nomad.datamodel.data/section_definitions/0'
                            ],
                            'quantities': [
                                {
                                    'description': 'Contains the nomad uri of a raw the data file connected to the current run. There\nshould be an value for the main_file_uri and all ancillary files.',
                                    'name': 'calculation_file_uri',
                                    'shape': [],
                                    'type': {'type_kind': 'python', 'type_data': 'str'},
                                },
                                {
                                    'description': 'Indicates whether this run terminated properly (true), or if it was killed or\nexited with an error code unequal to zero (false).',
                                    'name': 'clean_end',
                                    'shape': [],
                                    'type': {
                                        'type_kind': 'python',
                                        'type_data': 'bool',
                                    },
                                },
                                {
                                    'description': 'Number of references to the current section calculation.',
                                    'name': 'n_references',
                                    'shape': [],
                                    'type': {
                                        'type_kind': 'numpy',
                                        'type_data': 'int32',
                                    },
                                },
                                {
                                    'description': 'An optional calculation id, if one is found in the code input/output files.',
                                    'name': 'raw_id',
                                    'shape': [],
                                    'type': {'type_kind': 'python', 'type_data': 'str'},
                                },
                                {
                                    'categories': ['/category_definitions/0'],
                                    'description': 'Links the the current section to other run sections. Such a link is necessary for\nexample for workflows that may contain a series of runs.',
                                    'name': 'runs_ref',
                                    'shape': ['n_references'],
                                    'type': {
                                        'type_kind': 'reference',
                                        'type_data': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/3',
                                    },
                                },
                                {
                                    'categories': ['/category_definitions/0'],
                                    'description': 'Links the current section run to a section run containing the calculations from\nwhich the current section starts.',
                                    'name': 'starting_run_ref',
                                    'shape': [],
                                    'type': {
                                        'type_kind': 'reference',
                                        'type_data': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/3',
                                    },
                                },
                            ],
                            'sub_sections': [
                                {
                                    'name': 'calculation',
                                    'repeats': True,
                                    'sub_section': 'metainfo/nomad.datamodel.metainfo.simulation.calculation/section_definitions/36',
                                },
                                {
                                    'name': 'message',
                                    'sub_section': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/2',
                                },
                                {
                                    'name': 'method',
                                    'repeats': True,
                                    'sub_section': 'metainfo/nomad.datamodel.metainfo.simulation.method/section_definitions/44',
                                },
                                {
                                    'name': 'program',
                                    'sub_section': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/0',
                                },
                                {
                                    'name': 'system',
                                    'repeats': True,
                                    'sub_section': 'metainfo/nomad.datamodel.metainfo.simulation.system/section_definitions/8',
                                },
                                {
                                    'name': 'time_run',
                                    'sub_section': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/1',
                                },
                            ],
                        },
                    ],
                    'all_quantities': {
                        'Program.compilation_datetime': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/0/quantities/0',
                        'Program.compilation_host': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/0/quantities/1',
                        'Program.name': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/0/quantities/2',
                        'Program.version': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/0/quantities/3',
                        'Program.version_internal': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/0/quantities/4',
                        'TimeRun.cpu1_end': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/1/quantities/0',
                        'TimeRun.cpu1_start': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/1/quantities/1',
                        'TimeRun.date_end': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/1/quantities/2',
                        'TimeRun.date_start': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/1/quantities/3',
                        'TimeRun.wall_end': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/1/quantities/4',
                        'TimeRun.wall_start': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/1/quantities/5',
                        'MessageRun.type': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/2/quantities/0',
                        'MessageRun.value': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/2/quantities/1',
                        'Run.calculation_file_uri': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/3/quantities/0',
                        'Run.clean_end': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/3/quantities/1',
                        'Run.n_references': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/3/quantities/2',
                        'Run.raw_id': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/3/quantities/3',
                        'Run.runs_ref': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/3/quantities/4',
                        'Run.starting_run_ref': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/3/quantities/5',
                    },
                    'all_sub_sections': {
                        'Calculation': 'metainfo/nomad.datamodel.metainfo.simulation.calculation/section_definitions/36',
                        'Method': 'metainfo/nomad.datamodel.metainfo.simulation.method/section_definitions/44',
                        'Program': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/0',
                        'TimeRun': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/1',
                        'MessageRun': 'metainfo/nomad.datamodel.metainfo.simulation.run/section_definitions/2',
                        'System': 'metainfo/nomad.datamodel.metainfo.simulation.system/section_definitions/8',
                    },
                    'all_base_sections': {
                        'ArchiveSection': 'metainfo/nomad.datamodel.data/section_definitions/0'
                    },
                },
            }
        },
    )

    __ge_print(
        'general start from metainfo',
        {
            Token.METAINFO: {
                'm_request': {
                    'include': ['*test_data'],
                    'pagination': {'page_size': 50},
                }
            }
        },
        result={
            'metainfo': {
                'tests.processing.test_data': {
                    'name': 'tests.processing.test_data',
                    'section_definitions': [
                        {
                            'name': 'TestBatchSample',
                            'base_sections': [
                                'metainfo/nomad.datamodel.data/section_definitions/1'
                            ],
                            'quantities': [
                                {
                                    'description': 'Id for the batch',
                                    'name': 'batch_id',
                                    'type': {'type_kind': 'python', 'type_data': 'str'},
                                },
                                {
                                    'm_annotations': {
                                        'eln': [{'component': 'RichTextEditQuantity'}]
                                    },
                                    'description': 'Comments',
                                    'name': 'comments',
                                    'type': {'type_kind': 'python', 'type_data': 'str'},
                                },
                                {
                                    'description': 'Sample index',
                                    'name': 'sample_number',
                                    'type': {'type_kind': 'python', 'type_data': 'int'},
                                },
                            ],
                        },
                        {
                            'name': 'TestBatch',
                            'base_sections': [
                                'metainfo/nomad.datamodel.data/section_definitions/1'
                            ],
                            'quantities': [
                                {
                                    'm_annotations': {
                                        'eln': [{'component': 'StringEditQuantity'}]
                                    },
                                    'description': 'Id for the batch',
                                    'name': 'batch_id',
                                    'type': {'type_kind': 'python', 'type_data': 'str'},
                                },
                                {
                                    'm_annotations': {
                                        'eln': [{'component': 'NumberEditQuantity'}]
                                    },
                                    'description': 'Number of samples in batch',
                                    'name': 'n_samples',
                                    'type': {'type_kind': 'python', 'type_data': 'int'},
                                },
                                {
                                    'more': {
                                        'descriptions': 'The samples in the batch.',
                                        'type_data': 'metainfo/tests.processing.test_data/section_definitions/0',
                                    },
                                    'name': 'sample_refs',
                                    'shape': ['*'],
                                    'type': {
                                        'type_kind': 'reference',
                                        'type_data': 'metainfo/tests.processing.test_data/section_definitions/0',
                                    },
                                },
                            ],
                        },
                        {
                            'name': 'TestSection',
                            'base_sections': [
                                'metainfo/nomad.datamodel.data/section_definitions/0'
                            ],
                        },
                        {
                            'name': 'TestReferenceSection',
                            'base_sections': [
                                'metainfo/nomad.datamodel.data/section_definitions/0'
                            ],
                            'quantities': [
                                {
                                    'name': 'reference',
                                    'type': {
                                        'type_kind': 'reference',
                                        'type_data': 'metainfo/tests.processing.test_data/section_definitions/2',
                                    },
                                }
                            ],
                        },
                        {
                            'name': 'TestData',
                            'base_sections': [
                                'metainfo/nomad.datamodel.data/section_definitions/1'
                            ],
                            'sub_sections': [
                                {
                                    'name': 'reference_section',
                                    'sub_section': 'metainfo/tests.processing.test_data/section_definitions/3',
                                },
                                {
                                    'name': 'test_section',
                                    'repeats': True,
                                    'sub_section': 'metainfo/tests.processing.test_data/section_definitions/2',
                                },
                            ],
                        },
                    ],
                    'all_quantities': {
                        'TestBatchSample.batch_id': 'metainfo/tests.processing.test_data/section_definitions/0/quantities/0',
                        'TestBatchSample.comments': 'metainfo/tests.processing.test_data/section_definitions/0/quantities/1',
                        'TestBatchSample.sample_number': 'metainfo/tests.processing.test_data/section_definitions/0/quantities/2',
                        'TestBatch.batch_id': 'metainfo/tests.processing.test_data/section_definitions/1/quantities/0',
                        'TestBatch.n_samples': 'metainfo/tests.processing.test_data/section_definitions/1/quantities/1',
                        'TestBatch.sample_refs': 'metainfo/tests.processing.test_data/section_definitions/1/quantities/2',
                        'TestReferenceSection.reference': 'metainfo/tests.processing.test_data/section_definitions/3/quantities/0',
                    },
                    'all_sub_sections': {
                        'TestSection': 'metainfo/tests.processing.test_data/section_definitions/2',
                        'TestReferenceSection': 'metainfo/tests.processing.test_data/section_definitions/3',
                    },
                    'all_base_sections': {
                        'ArchiveSection': 'metainfo/nomad.datamodel.data/section_definitions/0',
                        'EntryData': 'metainfo/nomad.datamodel.data/section_definitions/1',
                    },
                },
            }
        },
    )


# noinspection DuplicatedCode,SpellCheckingInspection
def test_general_reader_search(json_dict, example_data_with_reference, user1):
    def __ge_print(msg, required, *, to_file: bool = False, result: dict = None):
        with MongoReader(required, user=user1) as reader:
            if result:
                assert_dict(reader.sync_read(), result)
            else:
                rprint(f'\n\nExample: {next(counter)} -> {msg}:')
                rprint(required)
                if not to_file:
                    rprint('output:')
                    rprint(reader.sync_read())
                else:
                    with open('archive_reader_test.json', 'w') as f:
                        f.write(json.dumps(reader.sync_read()))

    __ge_print(
        'general start from elastic search',
        {
            Token.SEARCH: {
                'm_request': {'query': {}, 'exclude': ['*']},
                'id_01': {Token.ENTRIES: {'mainfile': {'..': '*'}}},
            }
        },
        result={
            'search': {
                'id_01': {
                    'entries': {
                        'mainfile': {
                            'mainfile_for_id_01': {
                                '..': {
                                    'm_is': 'Directory',
                                    '1.aux': {
                                        'm_is': 'File',
                                        'path': '1.aux',
                                        'size': 8,
                                    },
                                    '2.aux': {
                                        'm_is': 'File',
                                        'path': '2.aux',
                                        'size': 8,
                                    },
                                    '3.aux': {
                                        'm_is': 'File',
                                        'path': '3.aux',
                                        'size': 8,
                                    },
                                    '4.aux': {
                                        'm_is': 'File',
                                        'path': '4.aux',
                                        'size': 8,
                                    },
                                    'mainfile_for_id_01': {
                                        'm_is': 'File',
                                        'path': 'mainfile_for_id_01',
                                        'size': 3227,
                                    },
                                    'mainfile_for_id_02': {
                                        'm_is': 'File',
                                        'path': 'mainfile_for_id_02',
                                        'size': 3227,
                                    },
                                    'mainfile_for_id_03': {
                                        'm_is': 'File',
                                        'path': 'mainfile_for_id_03',
                                        'size': 3227,
                                    },
                                    'mainfile_for_id_04': {
                                        'm_is': 'File',
                                        'path': 'mainfile_for_id_04',
                                        'size': 3227,
                                    },
                                    'mainfile_for_id_05': {
                                        'm_is': 'File',
                                        'path': 'mainfile_for_id_05',
                                        'size': 3227,
                                    },
                                    'mainfile_for_id_06': {
                                        'm_is': 'File',
                                        'path': 'mainfile_for_id_06',
                                        'size': 3227,
                                    },
                                }
                            }
                        }
                    }
                },
            }
        },
    )


@pytest.fixture(scope='function')
def custom_data(user1, proc_infra):
    yaml_archive = yaml.safe_load(
        """
---
definitions:
  name: test_package_name
  section_definitions:
  - name: MySection
    base_sections:
    - nomad.datamodel.data.EntryData
    quantities:
    - name: my_quantity
      type:
        type_kind: python
        type_data: str
    - name: datetime_list
      type:
        type_kind: custom
        type_data: nomad.metainfo.data_type.Datetime
      shape:
      - "*"
data:
  m_def: "/definitions/section_definitions/0"
  my_quantity: test_value
  datetime_list:
  - '2022-04-01'
  - '2022-04-02'
"""
    )
    archive = EntryArchive.m_from_dict(yaml_archive, m_context=ServerContext())
    data = ExampleData(main_author=user1)

    data.create_upload(
        upload_id='id_custom', upload_name='name_published', published=True
    )
    data.create_entry(
        upload_id='id_custom', entry_id='id_example', entry_archive=archive
    )
    data.save(
        with_files=True,
        with_es=True,
        with_mongo=True,
        additional_files_path='tests/data/proc/nested.zip',
    )

    yield data

    data.delete()


def test_custom_schema_archive_and_definition(user1, custom_data):
    def __entry_print(msg, required, *, to_file: bool = False, result: dict = None):
        with EntryReader(required, user=user1) as reader:
            response = reader.sync_read('id_example')
            if result:
                assert_dict(response, result)
            else:
                rprint(f'\n\nExample: {next(counter)} -> {msg}:')
                rprint(required)
                if not to_file:
                    rprint('output:')
                    rprint(response)
                else:
                    with open('entry_reader_test.json', 'w') as f:
                        f.write(json.dumps(response))

    __entry_print(
        'custom',
        {
            'm_request': {
                'directive': 'plain',
            },
            Token.ARCHIVE: {
                'data': {
                    'm_request': {
                        'directive': 'plain',
                    },
                    'm_def': {
                        'm_request': {
                            'directive': 'plain',
                            'export_whole_package': True,
                        },
                    },
                }
            },
        },
        result={
            'process_running': False,
            'current_process': None,
            'process_status': 'SUCCESS',
            'last_status_message': None,
            'errors': [],
            'warnings': [],
            'complete_time': None,
            'entry_id': 'id_example',
            'entry_create_time': '2025-08-20T10:46:57.535000',
            'mainfile_key': None,
            'upload_id': 'id_custom',
            'parser_name': 'parsers/vasp',
            'mainfile_path': 'mainfile_for_id_example',
            'uploads': {
                'id_custom': {
                    'entries': {
                        'id_example': {
                            'archive': {
                                'definitions': {
                                    'name': 'test_package_name',
                                    'section_definitions': [
                                        {
                                            'name': 'MySection',
                                            'base_sections': [
                                                'metainfo/nomad.datamodel.data/section_definitions/1'
                                            ],
                                            'quantities': [
                                                {
                                                    'name': 'my_quantity',
                                                    'type': {
                                                        'type_kind': 'python',
                                                        'type_data': 'str',
                                                    },
                                                },
                                                {
                                                    'name': 'datetime_list',
                                                    'shape': ['*'],
                                                    'type': {
                                                        'type_kind': 'custom',
                                                        'type_data': 'nomad.metainfo.data_type.Datetime',
                                                    },
                                                },
                                            ],
                                        }
                                    ],
                                    'all_quantities': {
                                        'MySection.my_quantity': 'uploads/id_custom/entries/id_example/archive/definitions/section_definitions/0/quantities/0',
                                        'MySection.datetime_list': 'uploads/id_custom/entries/id_example/archive/definitions/section_definitions/0/quantities/1',
                                    },
                                    'all_sub_sections': {},
                                    'all_base_sections': {
                                        'ArchiveSection': 'metainfo/nomad.datamodel.data/section_definitions/0',
                                        'EntryData': 'metainfo/nomad.datamodel.data/section_definitions/1',
                                    },
                                }
                            }
                        }
                    }
                }
            },
            'archive': {
                'data': {
                    'my_quantity': 'test_value',
                    'datetime_list': [
                        '2022-04-01T00:00:00+00:00',
                        '2022-04-02T00:00:00+00:00',
                    ],
                    'm_def': {
                        'm_def': 'uploads/id_custom/entries/id_example/archive/definitions/section_definitions/0',
                    },
                }
            },
        },
    )

    __entry_print(
        'custom',
        {
            Token.ARCHIVE: {
                'data': {
                    'm_def': {
                        'm_request': {
                            'directive': 'resolved',
                            'export_whole_package': True,
                            'depth': 1,
                        },
                    },
                }
            },
        },
        result={
            'uploads': {
                'id_custom': {
                    'entries': {
                        'id_example': {
                            'archive': {
                                'definitions': {
                                    'name': 'test_package_name',
                                    'section_definitions': [
                                        {
                                            'name': 'MySection',
                                            'base_sections': [
                                                'metainfo/nomad.datamodel.data/section_definitions/1'
                                            ],
                                            'quantities': [
                                                {
                                                    'name': 'my_quantity',
                                                    'type': {
                                                        'type_kind': 'python',
                                                        'type_data': 'str',
                                                    },
                                                },
                                                {
                                                    'name': 'datetime_list',
                                                    'shape': ['*'],
                                                    'type': {
                                                        'type_kind': 'custom',
                                                        'type_data': 'nomad.metainfo.data_type.Datetime',
                                                    },
                                                },
                                            ],
                                        }
                                    ],
                                    'all_quantities': {
                                        'MySection.my_quantity': 'uploads/id_custom/entries/id_example/archive/definitions/section_definitions/0/quantities/0',
                                        'MySection.datetime_list': 'uploads/id_custom/entries/id_example/archive/definitions/section_definitions/0/quantities/1',
                                    },
                                    'all_sub_sections': {},
                                    'all_base_sections': {
                                        'ArchiveSection': 'metainfo/nomad.datamodel.data/section_definitions/0',
                                        'EntryData': 'metainfo/nomad.datamodel.data/section_definitions/1',
                                    },
                                    'base_sections': [
                                        'metainfo/nomad.datamodel.data/section_definitions/1'
                                    ],
                                }
                            }
                        }
                    }
                }
            },
            'metainfo': {
                'nomad.datamodel.data': {
                    'name': 'nomad.datamodel.data',
                    'category_definitions': [
                        {
                            'name': 'EntryDataCategory',
                        },
                        {
                            'categories': ['/category_definitions/0'],
                            'label': 'Third-party ELN Integration',
                            'name': 'ElnIntegrationCategory',
                        },
                        {
                            'categories': ['/category_definitions/0'],
                            'label': 'Basic ELN',
                            'name': 'BasicElnCategory',
                        },
                        {
                            'categories': ['/category_definitions/0'],
                            'label': 'Example ELNs',
                            'name': 'ElnExampleCategory',
                        },
                        {
                            'categories': ['/category_definitions/0'],
                            'label': 'Use-cases',
                            'name': 'UseCaseElnCategory',
                        },
                        {
                            'categories': ['/category_definitions/0'],
                            'label': 'Workflows',
                            'name': 'WorkflowsElnCategory',
                        },
                    ],
                    'section_definitions': [
                        {
                            'description': 'Base class for sections in a NOMAD archive. Provides a framework for custom section normalization via the `normalize` function.',
                            'name': 'ArchiveSection',
                        },
                        {
                            'description': 'An empty base section definition. This can be used to add new top-level sections to an entry.',
                            'name': 'EntryData',
                            'base_sections': [
                                'metainfo/nomad.datamodel.data/section_definitions/0'
                            ],
                        },
                        {
                            'description': 'A person that is author of data in NOMAD or references by NOMAD.',
                            'name': 'Author',
                            'quantities': [
                                {
                                    'description': 'The name of the company and institutes the user identifies with',
                                    'name': 'affiliation',
                                    'type': {'type_kind': 'python', 'type_data': 'str'},
                                },
                                {
                                    'description': 'The address of the given affiliation',
                                    'name': 'affiliation_address',
                                    'type': {'type_kind': 'python', 'type_data': 'str'},
                                },
                                {
                                    'name': 'email',
                                    'type': {'type_kind': 'python', 'type_data': 'str'},
                                },
                                {
                                    'description': 'The users first name (including all other given names)',
                                    'name': 'first_name',
                                    'type': {
                                        'type_kind': 'custom',
                                        'type_data': 'nomad.metainfo.data_type.Capitalized',
                                    },
                                },
                                {
                                    'description': 'The users last name',
                                    'name': 'last_name',
                                    'type': {
                                        'type_kind': 'custom',
                                        'type_data': 'nomad.metainfo.data_type.Capitalized',
                                    },
                                },
                                {
                                    'm_annotations': {
                                        'elasticsearch': [
                                            'writers.name',
                                            'writers.name.text',
                                            'writers.name__suggestion',
                                        ]
                                    },
                                    'name': 'name',
                                    'type': {'type_kind': 'python', 'type_data': 'str'},
                                    'virtual': True,
                                },
                            ],
                        },
                        {
                            'm_annotations': {'pydantic': ['PydanticModel']},
                            'description': 'A NOMAD user. Typically a NOMAD user has a NOMAD account. The user related data is managed by\nNOMAD keycloak user-management system. Users are used to denote authors,\nreviewers, and owners of datasets.',
                            'name': 'User',
                            'base_sections': [
                                'metainfo/nomad.datamodel.data/section_definitions/2'
                            ],
                            'quantities': [
                                {
                                    'description': 'The time the account was created',
                                    'name': 'created',
                                    'type': {
                                        'type_kind': 'custom',
                                        'type_data': 'nomad.metainfo.data_type.Datetime',
                                    },
                                },
                                {
                                    'description': 'Bool that indicated, if the user is the admin',
                                    'name': 'is_admin',
                                    'type': {
                                        'type_kind': 'python',
                                        'type_data': 'bool',
                                    },
                                    'virtual': True,
                                },
                                {
                                    'name': 'is_oasis_admin',
                                    'default': False,
                                    'type': {
                                        'type_kind': 'python',
                                        'type_data': 'bool',
                                    },
                                },
                                {
                                    'description': 'Optional, legacy user id from the old NOMAD CoE repository.',
                                    'name': 'repo_user_id',
                                    'type': {'type_kind': 'python', 'type_data': 'str'},
                                },
                                {
                                    'm_annotations': {
                                        'elasticsearch': ['writers.user_id']
                                    },
                                    'description': 'The unique, persistent keycloak UUID',
                                    'name': 'user_id',
                                    'type': {'type_kind': 'python', 'type_data': 'str'},
                                },
                                {
                                    'description': 'The unique, persistent, user chosen username',
                                    'name': 'username',
                                    'type': {'type_kind': 'python', 'type_data': 'str'},
                                },
                            ],
                        },
                    ],
                    'all_quantities': {
                        'Author.affiliation': 'metainfo/nomad.datamodel.data/section_definitions/2/quantities/0',
                        'Author.affiliation_address': 'metainfo/nomad.datamodel.data/section_definitions/2/quantities/1',
                        'Author.email': 'metainfo/nomad.datamodel.data/section_definitions/2/quantities/2',
                        'Author.first_name': 'metainfo/nomad.datamodel.data/section_definitions/2/quantities/3',
                        'Author.last_name': 'metainfo/nomad.datamodel.data/section_definitions/2/quantities/4',
                        'Author.name': 'metainfo/nomad.datamodel.data/section_definitions/2/quantities/5',
                        'User.created': 'metainfo/nomad.datamodel.data/section_definitions/3/quantities/0',
                        'User.is_admin': 'metainfo/nomad.datamodel.data/section_definitions/3/quantities/1',
                        'User.is_oasis_admin': 'metainfo/nomad.datamodel.data/section_definitions/3/quantities/2',
                        'User.repo_user_id': 'metainfo/nomad.datamodel.data/section_definitions/3/quantities/3',
                        'User.user_id': 'metainfo/nomad.datamodel.data/section_definitions/3/quantities/4',
                        'User.username': 'metainfo/nomad.datamodel.data/section_definitions/3/quantities/5',
                    },
                    'all_sub_sections': {},
                    'all_base_sections': {
                        'ArchiveSection': 'metainfo/nomad.datamodel.data/section_definitions/0',
                        'Author': 'metainfo/nomad.datamodel.data/section_definitions/2',
                    },
                }
            },
            'archive': {
                'data': {
                    'm_def': {
                        'm_def': 'uploads/id_custom/entries/id_example/archive/definitions/section_definitions/0',
                    }
                }
            },
        },
    )

    __entry_print(
        'custom',
        {
            Token.ARCHIVE: {
                'data': {
                    'm_request': {
                        'directive': 'plain',
                    },
                    'm_def': {
                        'm_request': {
                            'directive': 'plain',
                        },
                    },
                }
            },
        },
        result={
            'uploads': {
                'id_custom': {
                    'entries': {
                        'id_example': {
                            'archive': {
                                'definitions': {
                                    'section_definitions': [
                                        {
                                            'name': 'MySection',
                                            'base_sections': [
                                                'metainfo/nomad.datamodel.data/section_definitions/1'
                                            ],
                                            'quantities': [
                                                {
                                                    'name': 'my_quantity',
                                                    'type': {
                                                        'type_kind': 'python',
                                                        'type_data': 'str',
                                                    },
                                                },
                                                {
                                                    'name': 'datetime_list',
                                                    'shape': ['*'],
                                                    'type': {
                                                        'type_kind': 'custom',
                                                        'type_data': 'nomad.metainfo.data_type.Datetime',
                                                    },
                                                },
                                            ],
                                        }
                                    ]
                                }
                            }
                        }
                    }
                }
            },
            'archive': {
                'data': {
                    'my_quantity': 'test_value',
                    'datetime_list': [
                        '2022-04-01T00:00:00+00:00',
                        '2022-04-02T00:00:00+00:00',
                    ],
                    'm_def': {
                        'm_def': 'uploads/id_custom/entries/id_example/archive/definitions/section_definitions/0',
                    },
                }
            },
        },
    )

    def __fs_print(msg, required, *, result: dict = None):
        with FileSystemReader(required, user=user1) as reader:
            if result:
                assert_dict(reader.sync_read('id_custom'), result)
            else:
                rprint(f'\n\nExample: {next(counter)} -> {msg}:')
                rprint(required)
                rprint('output:')
                rprint(reader.sync_read('id_custom'))

    __fs_print(
        'one level deep second page',
        {
            'm_request': {
                'directive': 'plain',
                'depth': 1,
                'pagination': {'page_size': 10, 'page': 2},
            },
        },
        result={
            'm_is': 'Directory',
            'mainfile_for_id_example': {
                'm_is': 'File',
                'path': 'mainfile_for_id_example',
                'size': 3227,
            },
        },
    )

    __fs_print(
        'two levels',
        {
            'm_request': {
                'directive': 'plain',
                'depth': 2,
                'pagination': {'page_size': 10, 'page': 2},
            },
        },
        result={
            'm_is': 'Directory',
            '3.aux': {'m_is': 'File', 'path': '3.aux', 'size': 8},
            '4.aux': {'m_is': 'File', 'path': '4.aux', 'size': 8},
            'edge_names': {
                '!"┬º$%&()=?.txt': {
                    'm_is': 'File',
                    'path': 'edge_names/!"┬º$%&()=?.txt',
                    'size': 0,
                },
                'suuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuper-long.txt': {
                    'm_is': 'File',
                    'path': 'edge_names/suuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuper-long.txt',
                    'size': 0,
                },
            },
            'entry.archive.json': {
                'm_is': 'File',
                'path': 'entry.archive.json',
                'size': 185,
            },
            'file.txt': {'m_is': 'File', 'path': 'file.txt', 'size': 0},
            'mainfile_for_id_example': {
                'm_is': 'File',
                'path': 'mainfile_for_id_example',
                'size': 3227,
            },
            'many_files': {
                'file1.txt': {
                    'm_is': 'File',
                    'path': 'many_files/file1.txt',
                    'size': 0,
                },
                'file10.txt': {
                    'm_is': 'File',
                    'path': 'many_files/file10.txt',
                    'size': 0,
                },
                'file100.txt': {
                    'm_is': 'File',
                    'path': 'many_files/file100.txt',
                    'size': 0,
                },
            },
        },
    )

    __fs_print(
        'different configs',
        {
            'm_request': {
                'directive': 'plain',
                'depth': 1,
                'pagination': {'page_size': 2, 'page': 2},
            },
            'many_files': {
                'm_request': {
                    'directive': 'plain',
                    'depth': 1,
                    'pagination': {'page_size': 3, 'page': 2},
                },
            },
        },
        result={
            'm_is': 'Directory',
            'preview': {'m_is': 'Directory'},
            'subdirs': {'m_is': 'Directory'},
            'many_files': {
                'm_is': 'Directory',
                'file11.txt': {
                    'm_is': 'File',
                    'path': 'many_files/file11.txt',
                    'size': 0,
                },
                'file12.txt': {
                    'm_is': 'File',
                    'path': 'many_files/file12.txt',
                    'size': 0,
                },
                'file13.txt': {
                    'm_is': 'File',
                    'path': 'many_files/file13.txt',
                    'size': 0,
                },
            },
        },
    )


@pytest.fixture(scope='function')
def example_data_with_reference(
    elastic_function, raw_files_module, mongo_function, user1, json_dict
):
    """
    Provides a couple of entries with references.

    Only used in test_required_reader_with_remote_reference.
    """
    data = ExampleData(main_author=user1)

    data.create_upload(
        upload_id='id_published_with_ref', upload_name='name_published', published=False
    )

    ref_list = [
        {
            'results': {'calculation_result_ref': '/run/0/calculation/1'}
        },  # plain direct reference
        {
            'results': {'calculation_result_ref': '#/run/0/calculation/1'}
        },  # new-style reference
        {
            'tasks': [
                {
                    'm_def': 'nomad.datamodel.metainfo.workflow.TaskReference',
                    'task': '../entries/id_01/archive#/workflow2',
                }
            ]
        },  # reference to another archive
        {
            'tasks': [
                {
                    'm_def': 'nomad.datamodel.metainfo.workflow.TaskReference',
                    'task': '../entries/id_05/archive#/workflow2',
                }
            ]
        },  # circular reference
        {
            'tasks': [
                {
                    'm_def': 'nomad.datamodel.metainfo.workflow.TaskReference',
                    'task': '../entries/id_04/archive#/workflow2',
                }
            ]
        },  # circular reference
        {
            'tasks': [
                {
                    'm_def': 'nomad.datamodel.metainfo.workflow.TaskReference',
                    'task': 'https://another.domain/entries/id_03/archive#/workflow2',
                }
            ]
        },  # remote reference
    ]

    del json_dict['results']

    for index, ref in enumerate(ref_list):
        ref['m_def'] = 'simulationworkflowschema.SimulationWorkflow'
        json_dict['workflow2'] = ref
        data.create_entry(
            upload_id='id_published_with_ref',
            entry_id=f'id_{index + 1:02d}',
            entry_archive=EntryArchive.m_from_dict(json_dict),
        )

    for archive in data.archives.values():
        archive.metadata.apply_archive_metadata(archive)

    data.save(with_files=True, with_es=True, with_mongo=True)

    yield data
    data.delete()


@pytest.fixture(scope='function')
def json_dict():
    return {
        'metadata': {'entry_id': 'test_id', 'upload_id': 'id_published_with_ref'},
        'results': {
            'properties': {
                'electronic': {
                    'dos_electronic': [
                        {'energies': '/run/0/calculation/1/dos_electronic/0/energies'}
                    ]
                }
            }
        },
        'run': [
            {
                'm_def': 'runschema.run.Run',
                'system': [
                    {
                        'atoms': {'labels': ['He']},
                        'symmetry': [{'space_group_number': 221}],
                    },
                    {
                        'atoms': {'labels': ['H']},
                        'symmetry': [{'space_group_number': 221}],
                    },
                ],
                'calculation': [
                    {
                        'system_ref': '/run/0/system/1',
                        'energy': {'total': {'value': 0.1}},
                    },
                    {
                        'system_ref': '/run/0/system/1',
                        'energy': {'total': {'value': 0.2}},
                        'dos_electronic': [
                            {'energies': [0.0, 0.1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]}
                        ],
                        'eigenvalues': [],
                    },
                    {
                        'system_ref': '/run/0/system/1',
                        'energy': {'total': {'value': 0.1}},
                    },
                ],
            }
        ],
        'workflow2': {
            'm_def': 'simulationworkflowschema.SimulationWorkflow',
            'results': {'calculation_result_ref': '/run/0/calculation/1'},
        },
    }


@pytest.fixture(scope='function')
def example_data_with_figure(proc_infra, user1):
    data = ExampleData(main_author=user1)

    data.create_upload(
        upload_id='id_published_with_ref', upload_name='name_published', published=False
    )

    directory = 'tests/data/datamodel/metainfo/plotly'
    mainfile = 'plotly.schema.archive.yaml'
    data.create_entry(
        upload_id='id_published_with_ref',
        entry_id='id_plotly',
        entry_archive=run_processing(directory, mainfile),
    )

    for archive in data.archives.values():
        archive.metadata.apply_archive_metadata(archive)

    data.save(with_files=True, with_es=True, with_mongo=True)

    yield data

    data.delete()


@pytest.mark.parametrize(
    'query,result',
    [
        # plain get default quantities
        # the references are not resolved
        pytest.param(
            {
                Token.ARCHIVE: {
                    'data': {
                        'figures[0]': {
                            'm_request': {
                                'directive': 'plain',
                            },
                        }
                    }
                },
            },
            {
                'archive': {
                    'data': {
                        'figures': [
                            {
                                'label': 'graph object 1',
                                'figure': {
                                    'data': {'x': '#xArr', 'y': '#xArr'},
                                    'layout': {
                                        'title': {'text': 'Plot in section level'},
                                        'xaxis': {'title': {'text': 'x data'}},
                                        'yaxis': {'title': {'text': 'y data'}},
                                    },
                                },
                            }
                        ]
                    }
                }
            },
            id='plain-read',
        ),
        pytest.param(
            {
                Token.ARCHIVE: {
                    'data': {
                        'figures[0]': {
                            'm_request': {
                                'directive': 'resolved',
                            },
                        }
                    }
                },
            },
            {
                'archive': {
                    'data': {
                        'xArr': [1.1, 2.0, 3.0, 4.0, 5.0],
                        'figures': [
                            {
                                'label': 'graph object 1',
                                'figure': {
                                    'data': {'x': '#xArr', 'y': '#xArr'},
                                    'layout': {
                                        'title': {'text': 'Plot in section level'},
                                        'xaxis': {'title': {'text': 'x data'}},
                                        'yaxis': {'title': {'text': 'y data'}},
                                    },
                                },
                            }
                        ],
                    }
                }
            },
            id='read-resolved',
        ),
    ],
)
def test_figure_resolution(user1, example_data_with_figure, query, result):
    def __entry_print(required, *, result=None):
        with EntryReader(required, user=user1) as reader:
            response = reader.sync_read('id_plotly')
            if result:
                assert_dict(response, result)

    __entry_print(query, result=result)
