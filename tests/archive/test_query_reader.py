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

from nomad.archive.query_reader import (
    EntryReader, UploadReader, UserReader, FileSystemReader, MongoReader, GeneralReader
)
from nomad.datamodel import EntryArchive
from nomad.utils.exampledata import ExampleData


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


def assert_list(l1, l2):
    assert len(l1) == len(l2)
    for i, j in zip(l1, l2):
        if isinstance(i, dict):
            assert_dict(i, j)
        elif isinstance(i, list):
            assert_list(i, j)
        else:
            assert_time(i, j)


def assert_dict(d1, d2):
    if GeneralReader.__CACHE__ in d1:
        del d1[GeneralReader.__CACHE__]
    if 'pagination' in d1:
        del d1['pagination']
    if 'm_def' in d1:
        del d1['m_def']
    assert set(d1.keys()) == set(d2.keys())
    for k, v in d1.items():
        if isinstance(v, dict):
            assert_dict(v, d2[k])
        elif isinstance(v, list):
            assert_list(v, d2[k])
        elif k == 'upload_files_server_path':
            continue
        else:
            assert_time(v, d2[k])


# noinspection SpellCheckingInspection,DuplicatedCode
def test_remote_reference(json_dict, example_data_with_reference, test_user):
    def increment():
        n = 0
        while True:
            n += 1
            yield n

    counter = increment()

    def __user_print(msg, required, *, result: dict = None):
        with UserReader(required, user=test_user) as reader:
            if result:
                assert_dict(reader.read(test_user.user_id), result)
            else:
                rprint(f'\n\nExample: {next(counter)} -> {msg}:')
                rprint(required)
                rprint('output:')
                rprint(reader.read(test_user.user_id))

    __user_print('plain user', {
        'm_request': {
            'directive': 'plain'
        }
    }, result={
        'name': 'Sheldon Cooper', 'first_name': 'Sheldon', 'last_name': 'Cooper',
        'email': 'sheldon.cooper@nomad-coe.eu', 'user_id': '00000000-0000-0000-0000-000000000001',
        'username': 'scooper', 'is_admin': False, 'is_oasis_admin': True})
    __user_print('plain user', '*', result={
        'name': 'Sheldon Cooper', 'first_name': 'Sheldon', 'last_name': 'Cooper',
        'email': 'sheldon.cooper@nomad-coe.eu', 'user_id': '00000000-0000-0000-0000-000000000001',
        'username': 'scooper', 'is_admin': False, 'is_oasis_admin': True})

    __user_print('link to uploads', {
        'm_request': {
            'directive': 'plain'
        },
        'm_uploads': {
            'm_request': {
                'directive': 'plain'
            },
        }
    }, result={
        'name': 'Sheldon Cooper', 'first_name': 'Sheldon', 'last_name': 'Cooper',
        'email': 'sheldon.cooper@nomad-coe.eu', 'user_id': '00000000-0000-0000-0000-000000000001',
        'username': 'scooper', 'is_admin': False, 'is_oasis_admin': True,
        'm_uploads': {'id_published_with_ref': 'id_published_with_ref'}})

    __user_print('link to uploads, resolve with metadata', {
        'm_request': {
            'directive': 'plain'
        },
        'm_uploads': {
            'm_request': {
                'directive': 'resolved',
                'resolve_type': 'upload'
            },
        }
    }, result={
        'name': 'Sheldon Cooper', 'first_name': 'Sheldon', 'last_name': 'Cooper',
        'email': 'sheldon.cooper@nomad-coe.eu', 'user_id': '00000000-0000-0000-0000-000000000001',
        'username': 'scooper', 'is_admin': False, 'is_oasis_admin': True, 'm_uploads': {
            'id_published_with_ref': {
                'process_running': False, 'current_process': 'process_upload',
                'process_status': 'SUCCESS', 'last_status_message': None, 'errors': [],
                'warnings': [], 'complete_time': '2023-03-05T21:24:22.172000',
                'upload_id': 'id_published_with_ref', 'upload_name': 'name_published',
                'upload_create_time': '2023-03-05T21:24:22.171000',
                'main_author': '00000000-0000-0000-0000-000000000001', 'coauthors': [],
                'reviewers': [], 'viewers': ['00000000-0000-0000-0000-000000000001'],
                'writers': ['00000000-0000-0000-0000-000000000001'], 'published': False,
                'published_to': [], 'publish_time': None, 'with_embargo': False,
                'embargo_length': 0, 'license': 'CC BY 4.0', 'entries': 6,
                'upload_files_server_path': 'id_published_with_ref'}}})

    __user_print('link to uploads, resolve with metadata, form 2 using dict style with explicit upload id', {
        'm_request': {
            'directive': 'plain'
        },
        'm_uploads': {
            'id_published_with_ref': {
                'm_request': {
                    'directive': 'plain',
                },
            }
        }
    }, result={
        'name': 'Sheldon Cooper', 'first_name': 'Sheldon', 'last_name': 'Cooper',
        'email': 'sheldon.cooper@nomad-coe.eu', 'user_id': '00000000-0000-0000-0000-000000000001',
        'username': 'scooper', 'is_admin': False, 'is_oasis_admin': True, 'm_uploads': {
            'id_published_with_ref': {
                'process_running': False, 'current_process': 'process_upload',
                'process_status': 'SUCCESS', 'last_status_message': None, 'errors': [],
                'warnings': [], 'complete_time': '2023-03-05T21:30:22.807000',
                'upload_id': 'id_published_with_ref', 'upload_name': 'name_published',
                'upload_create_time': '2023-03-05T21:30:22.806000',
                'main_author': '00000000-0000-0000-0000-000000000001', 'coauthors': [],
                'reviewers': [], 'viewers': ['00000000-0000-0000-0000-000000000001'],
                'writers': ['00000000-0000-0000-0000-000000000001'], 'published': False,
                'published_to': [], 'publish_time': None, 'with_embargo': False,
                'embargo_length': 0, 'license': 'CC BY 4.0', 'entries': 6,
                'upload_files_server_path': 'id_published_with_ref'}}})

    __user_print('link to entries directly from user, resolve with metadata', {
        'm_request': {
            'directive': 'plain'
        },
        'm_entries': {
            'm_request': {
                'directive': 'resolved',
                'resolve_type': 'entry'
            },
        }
    }, result={
        'name': 'Sheldon Cooper', 'first_name': 'Sheldon', 'last_name': 'Cooper',
        'email': 'sheldon.cooper@nomad-coe.eu', 'user_id': '00000000-0000-0000-0000-000000000001',
        'username': 'scooper', 'is_admin': False, 'is_oasis_admin': True, 'm_entries': {
            'id_01': {
                'process_running': False, 'current_process': None, 'process_status': 'SUCCESS',
                'last_status_message': None, 'errors': [], 'warnings': [], 'complete_time': None,
                'entry_id': 'id_01', 'entry_create_time': '2023-03-05T21:27:24.488000',
                'mainfile': 'mainfile_for_id_01', 'mainfile_key': None, 'upload_id': 'id_published_with_ref',
                'parser_name': 'parsers/vasp', 'entry_metadata': None},
            'id_02': {
                'process_running': False, 'current_process': None, 'process_status': 'SUCCESS',
                'last_status_message': None, 'errors': [], 'warnings': [], 'complete_time': None,
                'entry_id': 'id_02', 'entry_create_time': '2023-03-05T21:27:24.489000',
                'mainfile': 'mainfile_for_id_02', 'mainfile_key': None, 'upload_id': 'id_published_with_ref',
                'parser_name': 'parsers/vasp', 'entry_metadata': None},
            'id_03': {
                'process_running': False, 'current_process': None, 'process_status': 'SUCCESS',
                'last_status_message': None, 'errors': [], 'warnings': [], 'complete_time': None,
                'entry_id': 'id_03', 'entry_create_time': '2023-03-05T21:27:24.490000',
                'mainfile': 'mainfile_for_id_03', 'mainfile_key': None, 'upload_id': 'id_published_with_ref',
                'parser_name': 'parsers/vasp', 'entry_metadata': None},
            'id_04': {
                'process_running': False, 'current_process': None, 'process_status': 'SUCCESS',
                'last_status_message': None, 'errors': [], 'warnings': [], 'complete_time': None,
                'entry_id': 'id_04', 'entry_create_time': '2023-03-05T21:27:24.491000',
                'mainfile': 'mainfile_for_id_04', 'mainfile_key': None, 'upload_id': 'id_published_with_ref',
                'parser_name': 'parsers/vasp', 'entry_metadata': None},
            'id_05': {
                'process_running': False, 'current_process': None, 'process_status': 'SUCCESS',
                'last_status_message': None, 'errors': [], 'warnings': [], 'complete_time': None,
                'entry_id': 'id_05', 'entry_create_time': '2023-03-05T21:27:24.492000',
                'mainfile': 'mainfile_for_id_05', 'mainfile_key': None, 'upload_id': 'id_published_with_ref',
                'parser_name': 'parsers/vasp', 'entry_metadata': None},
            'id_06': {
                'process_running': False, 'current_process': None, 'process_status': 'SUCCESS',
                'last_status_message': None, 'errors': [], 'warnings': [], 'complete_time': None,
                'entry_id': 'id_06', 'entry_create_time': '2023-03-05T21:27:24.493000',
                'mainfile': 'mainfile_for_id_06', 'mainfile_key': None, 'upload_id': 'id_published_with_ref',
                'parser_name': 'parsers/vasp', 'entry_metadata': None}}})

    __user_print(
        'link to entries directly from user, resolve with metadata, form 2 using dict style with explicit entry id', {
            'm_request': {
                'directive': 'plain'
            },
            'm_entries': {
                'id_01': '*'
            }
        }, result={
            'name': 'Sheldon Cooper', 'first_name': 'Sheldon', 'last_name': 'Cooper',
            'email': 'sheldon.cooper@nomad-coe.eu', 'user_id': '00000000-0000-0000-0000-000000000001',
            'username': 'scooper', 'is_admin': False, 'is_oasis_admin': True, 'm_entries': {
                'id_01': {
                    'process_running': False, 'current_process': None, 'process_status': 'SUCCESS',
                    'last_status_message': None, 'errors': [], 'warnings': [], 'complete_time': None,
                    'entry_id': 'id_01', 'entry_create_time': '2023-03-05T21:30:22.809000',
                    'mainfile': 'mainfile_for_id_01', 'mainfile_key': None, 'upload_id': 'id_published_with_ref',
                    'parser_name': 'parsers/vasp', 'entry_metadata': None}}})

    __user_print('link to entries from uploads, resolve with metadata', {
        'm_request': {
            'directive': 'plain'
        },
        'm_uploads': {
            'id_published_with_ref': {
                'm_request': {
                    'directive': 'plain',
                },
                'm_entries': {
                    'm_request': {
                        'directive': 'resolved',
                        'resolve_type': 'entry'
                    },
                }
            }

        }
    }, result={
        'name': 'Sheldon Cooper', 'first_name': 'Sheldon', 'last_name': 'Cooper',
        'email': 'sheldon.cooper@nomad-coe.eu', 'user_id': '00000000-0000-0000-0000-000000000001',
        'username': 'scooper', 'is_admin': False, 'is_oasis_admin': True, 'm_uploads': {
            'id_published_with_ref': {
                'process_running': False, 'current_process': 'process_upload',
                'process_status': 'SUCCESS', 'last_status_message': None, 'errors': [],
                'warnings': [], 'complete_time': '2023-03-05T21:31:56.873000',
                'upload_id': 'id_published_with_ref', 'upload_name': 'name_published',
                'upload_create_time': '2023-03-05T21:31:56.872000',
                'main_author': '00000000-0000-0000-0000-000000000001', 'coauthors': [],
                'reviewers': [], 'viewers': ['00000000-0000-0000-0000-000000000001'],
                'writers': ['00000000-0000-0000-0000-000000000001'], 'published': False,
                'published_to': [], 'publish_time': None, 'with_embargo': False,
                'embargo_length': 0, 'license': 'CC BY 4.0', 'entries': 6,
                'upload_files_server_path': 'id_published_with_ref',
                'm_entries': {
                    'id_01': {
                        'process_running': False, 'current_process': None,
                        'process_status': 'SUCCESS', 'last_status_message': None,
                        'errors': [], 'warnings': [], 'complete_time': None,
                        'entry_id': 'id_01',
                        'entry_create_time': '2023-03-05T21:31:56.875000',
                        'mainfile': 'mainfile_for_id_01', 'mainfile_key': None,
                        'upload_id': 'id_published_with_ref',
                        'parser_name': 'parsers/vasp', 'entry_metadata': None},
                    'id_02': {
                        'process_running': False, 'current_process': None,
                        'process_status': 'SUCCESS', 'last_status_message': None,
                        'errors': [], 'warnings': [], 'complete_time': None,
                        'entry_id': 'id_02',
                        'entry_create_time': '2023-03-05T21:31:56.876000',
                        'mainfile': 'mainfile_for_id_02', 'mainfile_key': None,
                        'upload_id': 'id_published_with_ref',
                        'parser_name': 'parsers/vasp', 'entry_metadata': None},
                    'id_03': {
                        'process_running': False, 'current_process': None,
                        'process_status': 'SUCCESS', 'last_status_message': None,
                        'errors': [], 'warnings': [], 'complete_time': None,
                        'entry_id': 'id_03',
                        'entry_create_time': '2023-03-05T21:31:56.877000',
                        'mainfile': 'mainfile_for_id_03', 'mainfile_key': None,
                        'upload_id': 'id_published_with_ref',
                        'parser_name': 'parsers/vasp', 'entry_metadata': None},
                    'id_04': {
                        'process_running': False, 'current_process': None,
                        'process_status': 'SUCCESS', 'last_status_message': None,
                        'errors': [], 'warnings': [], 'complete_time': None,
                        'entry_id': 'id_04',
                        'entry_create_time': '2023-03-05T21:31:56.878000',
                        'mainfile': 'mainfile_for_id_04', 'mainfile_key': None,
                        'upload_id': 'id_published_with_ref',
                        'parser_name': 'parsers/vasp', 'entry_metadata': None},
                    'id_05': {
                        'process_running': False, 'current_process': None,
                        'process_status': 'SUCCESS', 'last_status_message': None,
                        'errors': [], 'warnings': [], 'complete_time': None,
                        'entry_id': 'id_05',
                        'entry_create_time': '2023-03-05T21:31:56.879000',
                        'mainfile': 'mainfile_for_id_05', 'mainfile_key': None,
                        'upload_id': 'id_published_with_ref',
                        'parser_name': 'parsers/vasp', 'entry_metadata': None},
                    'id_06': {
                        'process_running': False, 'current_process': None,
                        'process_status': 'SUCCESS', 'last_status_message': None,
                        'errors': [], 'warnings': [], 'complete_time': None,
                        'entry_id': 'id_06',
                        'entry_create_time': '2023-03-05T21:31:56.880000',
                        'mainfile': 'mainfile_for_id_06', 'mainfile_key': None,
                        'upload_id': 'id_published_with_ref',
                        'parser_name': 'parsers/vasp',
                        'entry_metadata': None}}}}})

    __user_print('link to entries from uploads, resolve with metadata, dict style', {
        'm_request': {
            'directive': 'plain'
        },
        'm_uploads': {
            'id_published_with_ref': {
                'm_request': {
                    'directive': 'plain',
                },
                'm_entries': {
                    'id_01': {
                        'm_request': {
                            'directive': 'plain',
                        }
                    },
                    '*': {
                        'm_request': {
                            'directive': 'plain',
                            'include': ['entry_id', 'mainfile']
                        },
                    }
                }
            }

        }
    }, result={
        'name': 'Sheldon Cooper', 'first_name': 'Sheldon', 'last_name': 'Cooper',
        'email': 'sheldon.cooper@nomad-coe.eu', 'user_id': '00000000-0000-0000-0000-000000000001',
        'username': 'scooper', 'is_admin': False, 'is_oasis_admin': True, 'm_uploads': {
            'id_published_with_ref': {
                'process_running': False, 'current_process': 'process_upload',
                'process_status': 'SUCCESS', 'last_status_message': None, 'errors': [],
                'warnings': [], 'complete_time': '2023-03-05T21:31:56.873000',
                'upload_id': 'id_published_with_ref', 'upload_name': 'name_published',
                'upload_create_time': '2023-03-05T21:31:56.872000',
                'main_author': '00000000-0000-0000-0000-000000000001', 'coauthors': [],
                'reviewers': [], 'viewers': ['00000000-0000-0000-0000-000000000001'],
                'writers': ['00000000-0000-0000-0000-000000000001'], 'published': False,
                'published_to': [], 'publish_time': None, 'with_embargo': False,
                'embargo_length': 0, 'license': 'CC BY 4.0', 'entries': 6,
                'upload_files_server_path': 'id_published_with_ref',
                'm_entries': {
                    'id_02': {'entry_id': 'id_02', 'mainfile': 'mainfile_for_id_02'},
                    'id_03': {'entry_id': 'id_03', 'mainfile': 'mainfile_for_id_03'},
                    'id_04': {'entry_id': 'id_04', 'mainfile': 'mainfile_for_id_04'},
                    'id_05': {'entry_id': 'id_05', 'mainfile': 'mainfile_for_id_05'},
                    'id_06': {'entry_id': 'id_06', 'mainfile': 'mainfile_for_id_06'},
                    'id_01': {
                        'process_running': False, 'current_process': None,
                        'process_status': 'SUCCESS', 'last_status_message': None,
                        'errors': [], 'warnings': [], 'complete_time': None,
                        'entry_id': 'id_01',
                        'entry_create_time': '2023-03-05T21:31:56.875000',
                        'mainfile': 'mainfile_for_id_01', 'mainfile_key': None,
                        'upload_id': 'id_published_with_ref',
                        'parser_name': 'parsers/vasp',
                        'entry_metadata': None}}}}})

    __user_print('uploads to entries back to uploads', {
        'm_uploads': {
            'id_published_with_ref': {
                'm_entries': {
                    'id_01': {
                        'm_request': {
                            'directive': 'plain',
                        },
                        'upload_id': {
                            'm_request': {
                                'directive': 'resolved',
                                'resolve_type': 'upload'
                            }
                        }
                    }
                }
            }

        }
    }, result={'m_uploads': {'id_published_with_ref': {'m_entries': {
        'id_01': {
            'process_running': False, 'current_process': None, 'process_status': 'SUCCESS',
            'last_status_message': None, 'errors': [], 'warnings': [], 'complete_time': None, 'entry_id': 'id_01',
            'entry_create_time': '2023-03-05T21:31:56.875000', 'mainfile': 'mainfile_for_id_01',
            'mainfile_key': None, 'parser_name': 'parsers/vasp', 'entry_metadata': None,
            'upload_id': {
                'process_running': False, 'current_process': 'process_upload',
                'process_status': 'SUCCESS', 'last_status_message': None, 'errors': [], 'warnings': [],
                'complete_time': '2023-03-05T21:31:56.873000', 'upload_id': 'id_published_with_ref',
                'upload_name': 'name_published', 'upload_create_time': '2023-03-05T21:31:56.872000',
                'main_author': '00000000-0000-0000-0000-000000000001', 'coauthors': [], 'reviewers': [],
                'viewers': ['00000000-0000-0000-0000-000000000001'],
                'writers': ['00000000-0000-0000-0000-000000000001'], 'published': False,
                'published_to': [], 'publish_time': None, 'with_embargo': False, 'embargo_length': 0,
                'license': 'CC BY 4.0', 'entries': 6,
                'upload_files_server_path': 'id_published_with_ref'}}}}}})

    __user_print('uploads to entries back to uploads, resolve viewers', {
        'm_uploads': {
            'id_published_with_ref': {
                'm_entries': {
                    'id_01': {
                        'm_request': {
                            'directive': 'plain',
                        },
                        'upload_id': {
                            'm_request': {
                                'directive': 'plain',
                            },
                            'viewers': {
                                'username': {
                                    'm_request': '*'
                                }
                            }
                        }
                    }
                }
            }

        }
    }, result={'m_uploads': {'id_published_with_ref': {'m_entries': {
        'id_01': {
            'process_running': False, 'current_process': None, 'process_status': 'SUCCESS',
            'last_status_message': None, 'errors': [], 'warnings': [], 'complete_time': None, 'entry_id': 'id_01',
            'entry_create_time': '2023-03-05T21:31:56.875000', 'mainfile': 'mainfile_for_id_01',
            'mainfile_key': None, 'parser_name': 'parsers/vasp', 'entry_metadata': None,
            'upload_id': {
                'process_running': False, 'current_process': 'process_upload',
                'process_status': 'SUCCESS', 'last_status_message': None, 'errors': [], 'warnings': [],
                'complete_time': '2023-03-05T21:31:56.873000', 'upload_id': 'id_published_with_ref',
                'upload_name': 'name_published', 'upload_create_time': '2023-03-05T21:31:56.872000',
                'main_author': '00000000-0000-0000-0000-000000000001', 'coauthors': [], 'reviewers': [],
                'writers': ['00000000-0000-0000-0000-000000000001'], 'published': False,
                'published_to': [], 'publish_time': None, 'with_embargo': False, 'embargo_length': 0,
                'license': 'CC BY 4.0', 'entries': 6,
                'upload_files_server_path': 'id_published_with_ref',
                'viewers': [{'username': 'scooper'}]}}}}}})

    __user_print('uploads to entries to archive', {
        'm_uploads': {
            'id_published_with_ref': {
                'm_entries': {
                    'id_01': {
                        'm_request': {
                            'directive': 'plain',
                        },
                        'm_archive': {
                            'm_request': {
                                'directive': 'plain',
                                'include': ['results']
                            },
                        }
                    }
                }
            }

        }
    }, result={'m_uploads': {'id_published_with_ref': {'m_entries': {
        'id_01': {
            'process_running': False, 'current_process': None, 'process_status': 'SUCCESS',
            'last_status_message': None, 'errors': [], 'warnings': [], 'complete_time': None, 'entry_id': 'id_01',
            'entry_create_time': '2023-03-05T21:31:56.875000', 'mainfile': 'mainfile_for_id_01',
            'mainfile_key': None, 'upload_id': 'id_published_with_ref', 'parser_name': 'parsers/vasp',
            'entry_metadata': None, 'm_archive': {'results': {
                'material': {
                    'material_id': 'test_material_id', 'elements': ['H', 'O'],
                    'symmetry': {'crystal_system': 'cubic'}},
                'method': {'simulation': {'program_name': 'VASP', 'dft': {'xc_functional_type': 'GGA'}}},
                'properties': {'n_calculations': 1, 'electronic': {
                    'dos_electronic': [{'spin_polarized': False, 'band_gap': [{'type': 'indirect'}]}]}}}}}}}}})

    def __upload_print(msg, required, *, result: dict = None):
        with UploadReader(required, user=test_user) as reader:
            if result:
                assert_dict(reader.read('id_published_with_ref'), result)
            else:
                rprint(f'\n\nExample: {next(counter)} -> {msg}:')
                rprint(required)
                rprint('output:')
                rprint(reader.read('id_published_with_ref'))

    __upload_print('plain upload reader', {
        'm_request': {
            'directive': 'plain',
        },
    }, result={
        'process_running': False, 'current_process': 'process_upload', 'process_status': 'SUCCESS',
        'last_status_message': None, 'errors': [], 'warnings': [], 'complete_time': '2023-03-05T22:12:08.421000',
        'upload_id': 'id_published_with_ref', 'upload_name': 'name_published',
        'upload_create_time': '2023-03-05T22:12:08.420000',
        'main_author': '00000000-0000-0000-0000-000000000001', 'coauthors': [], 'reviewers': [],
        'viewers': ['00000000-0000-0000-0000-000000000001'], 'writers': ['00000000-0000-0000-0000-000000000001'],
        'published': False, 'published_to': [], 'publish_time': None, 'with_embargo': False, 'embargo_length': 0,
        'license': 'CC BY 4.0', 'entries': 6,
        'upload_files_server_path': 'id_published_with_ref'})

    __upload_print('upload, resolve raw files', {
        'm_request': '*',
        'm_raw': {
            'm_request': '*',
            'mainfile_for_id_01': {
                'm_request': {
                    'directive': 'resolved',
                },
            }
        }
    }, result={
        'process_running': False, 'current_process': 'process_upload', 'process_status': 'SUCCESS',
        'last_status_message': None, 'errors': [], 'warnings': [], 'complete_time': '2023-03-05T22:16:52.436000',
        'upload_id': 'id_published_with_ref', 'upload_name': 'name_published',
        'upload_create_time': '2023-03-05T22:16:52.435000',
        'main_author': '00000000-0000-0000-0000-000000000001', 'coauthors': [], 'reviewers': [],
        'viewers': ['00000000-0000-0000-0000-000000000001'], 'writers': ['00000000-0000-0000-0000-000000000001'],
        'published': False, 'published_to': [], 'publish_time': None, 'with_embargo': False, 'embargo_length': 0,
        'license': 'CC BY 4.0', 'entries': 6,
        'upload_files_server_path': 'id_published_with_ref',
        'm_raw': {
            'mainfile_for_id_02': {
                'path': 'mainfile_for_id_02', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
            'mainfile_for_id_03': {
                'path': 'mainfile_for_id_03', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
            'mainfile_for_id_04': {
                'path': 'mainfile_for_id_04', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
            'mainfile_for_id_05': {
                'path': 'mainfile_for_id_05', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
            'mainfile_for_id_06': {
                'path': 'mainfile_for_id_06', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
            'mainfile_for_id_01': {
                'path': 'mainfile_for_id_01', 'is_file': True, 'size': 3231, 'access': 'unpublished',
                'm_entries': {
                    'process_running': False, 'current_process': None,
                    'process_status': 'SUCCESS', 'last_status_message': None,
                    'errors': [], 'warnings': [], 'complete_time': None,
                    'entry_id': 'id_01',
                    'entry_create_time': '2023-03-05T22:16:52.438000',
                    'mainfile': 'mainfile_for_id_01', 'mainfile_key': None,
                    'upload_id': 'id_published_with_ref',
                    'parser_name': 'parsers/vasp', 'entry_metadata': None}}}})

    __upload_print('upload, resolve user', {
        'm_request': {
            'directive': 'plain',
        },
        'viewers': {
            'm_request': {
                'directive': 'resolved',
                'resolve_type': 'user'
            }
        }
    }, result={
        'process_running': False, 'current_process': 'process_upload', 'process_status': 'SUCCESS',
        'last_status_message': None, 'errors': [], 'warnings': [], 'complete_time': '2023-03-05T22:16:52.436000',
        'upload_id': 'id_published_with_ref', 'upload_name': 'name_published',
        'upload_create_time': '2023-03-05T22:16:52.435000',
        'main_author': '00000000-0000-0000-0000-000000000001', 'coauthors': [], 'reviewers': [],
        'writers': ['00000000-0000-0000-0000-000000000001'], 'published': False, 'published_to': [],
        'publish_time': None, 'with_embargo': False, 'embargo_length': 0, 'license': 'CC BY 4.0', 'entries': 6,
        'upload_files_server_path': 'id_published_with_ref',
        'viewers': [{
            'name': 'Sheldon Cooper', 'first_name': 'Sheldon', 'last_name': 'Cooper',
            'email': 'sheldon.cooper@nomad-coe.eu', 'user_id': '00000000-0000-0000-0000-000000000001',
            'username': 'scooper', 'is_admin': False, 'is_oasis_admin': True}]})

    __upload_print('resolve itself using upload id', {
        'm_request': {
            'directive': 'plain',
        },
        'upload_id': {
            'm_request': {
                'directive': 'resolved',
                'resolve_type': 'upload'
            }
        }
    }, result={
        'process_running': False, 'current_process': 'process_upload', 'process_status': 'SUCCESS',
        'last_status_message': None, 'errors': [], 'warnings': [], 'complete_time': '2023-03-05T22:16:52.436000',
        'upload_name': 'name_published', 'upload_create_time': '2023-03-05T22:16:52.435000',
        'main_author': '00000000-0000-0000-0000-000000000001', 'coauthors': [], 'reviewers': [],
        'viewers': ['00000000-0000-0000-0000-000000000001'], 'writers': ['00000000-0000-0000-0000-000000000001'],
        'published': False, 'published_to': [], 'publish_time': None, 'with_embargo': False, 'embargo_length': 0,
        'license': 'CC BY 4.0', 'entries': 6,
        'upload_files_server_path': 'id_published_with_ref',
        'upload_id': {
            'process_running': False, 'current_process': 'process_upload', 'process_status': 'SUCCESS',
            'last_status_message': None, 'errors': [], 'warnings': [],
            'complete_time': '2023-03-05T22:16:52.436000', 'upload_id': 'id_published_with_ref',
            'upload_name': 'name_published', 'upload_create_time': '2023-03-05T22:16:52.435000',
            'main_author': '00000000-0000-0000-0000-000000000001', 'coauthors': [], 'reviewers': [],
            'viewers': ['00000000-0000-0000-0000-000000000001'],
            'writers': ['00000000-0000-0000-0000-000000000001'], 'published': False,
            'published_to': [], 'publish_time': None, 'with_embargo': False, 'embargo_length': 0,
            'license': 'CC BY 4.0', 'entries': 6,
            'upload_files_server_path': 'id_published_with_ref'}})

    __upload_print('resolve itself twice then go to entry', {
        'upload_id': {
            'm_request': {
                'directive': 'plain',
            },
            'upload_id': {
                'm_entries': {
                    'id_01': {
                        'm_request': {
                            'directive': 'plain',
                        },
                    }
                }
            }
        }
    }, result={'upload_id': {
        'process_running': False, 'current_process': 'process_upload', 'process_status': 'SUCCESS',
        'last_status_message': None, 'errors': [], 'warnings': [],
        'complete_time': '2023-03-05T22:16:52.436000', 'upload_name': 'name_published',
        'upload_create_time': '2023-03-05T22:16:52.435000',
        'main_author': '00000000-0000-0000-0000-000000000001', 'coauthors': [], 'reviewers': [],
        'viewers': ['00000000-0000-0000-0000-000000000001'],
        'writers': ['00000000-0000-0000-0000-000000000001'], 'published': False,
        'published_to': [], 'publish_time': None, 'with_embargo': False, 'embargo_length': 0,
        'license': 'CC BY 4.0', 'entries': 6,
        'upload_files_server_path': 'id_published_with_ref',
        'upload_id': {
            'm_entries': {
                'id_01': {
                    'process_running': False, 'current_process': None,
                    'process_status': 'SUCCESS',
                    'last_status_message': None, 'errors': [],
                    'warnings': [], 'complete_time': None,
                    'entry_id': 'id_01',
                    'entry_create_time': '2023-03-05T22:16:52.438000',
                    'mainfile': 'mainfile_for_id_01',
                    'mainfile_key': None,
                    'upload_id': 'id_published_with_ref',
                    'parser_name': 'parsers/vasp',
                    'entry_metadata': None}}}}})

    __upload_print('resolve itself twice then go to entry, collect other info on different levels', {
        'm_request': {
            'directive': 'plain',
        },
        'upload_id': {
            'upload_id': {
                'm_entries': {
                    'id_01': {
                        'm_request': {
                            'directive': 'plain',
                        },
                    }
                }
            }
        }
    }, result={
        'process_running': False, 'current_process': 'process_upload', 'process_status': 'SUCCESS',
        'last_status_message': None, 'errors': [], 'warnings': [], 'complete_time': '2023-03-05T22:16:52.436000',
        'upload_name': 'name_published', 'upload_create_time': '2023-03-05T22:16:52.435000',
        'main_author': '00000000-0000-0000-0000-000000000001', 'coauthors': [], 'reviewers': [],
        'viewers': ['00000000-0000-0000-0000-000000000001'], 'writers': ['00000000-0000-0000-0000-000000000001'],
        'published': False, 'published_to': [], 'publish_time': None, 'with_embargo': False, 'embargo_length': 0,
        'license': 'CC BY 4.0', 'entries': 6,
        'upload_files_server_path': 'id_published_with_ref',
        'upload_id': {
            'upload_id': {
                'm_entries': {
                    'id_01': {
                        'process_running': False, 'current_process': None, 'process_status': 'SUCCESS',
                        'last_status_message': None, 'errors': [], 'warnings': [], 'complete_time': None,
                        'entry_id': 'id_01', 'entry_create_time': '2023-03-05T22:16:52.438000',
                        'mainfile': 'mainfile_for_id_01', 'mainfile_key': None,
                        'upload_id': 'id_published_with_ref', 'parser_name': 'parsers/vasp',
                        'entry_metadata': None}}}}})

    def __entry_print(msg, required, *, to_file: bool = False, result: dict = None):
        with EntryReader(required, user=test_user) as reader:
            if result:
                assert_dict(reader.read('id_03'), result)
            else:
                rprint(f'\n\nExample: {next(counter)} -> {msg}:')
                rprint(required)
                if not to_file:
                    rprint('output:')
                    rprint(reader.read('id_03'))
                else:
                    with open('entry_reader_test.json', 'w') as f:
                        f.write(json.dumps(reader.read('id_03')))

    __entry_print('plain entry reader', {
        'm_request': {
            'directive': 'plain',
        },
    }, result={
        'process_running': False, 'current_process': None, 'process_status': 'SUCCESS',
        'last_status_message': None, 'errors': [], 'warnings': [], 'complete_time': None, 'entry_id': 'id_03',
        'entry_create_time': '2023-03-05T22:16:52.440000', 'mainfile': 'mainfile_for_id_03',
        'mainfile_key': None, 'upload_id': 'id_published_with_ref', 'parser_name': 'parsers/vasp',
        'entry_metadata': None})
    __entry_print('plain entry reader, resolve inplace', {
        'm_archive': {
            'm_request': {
                'directive': 'resolved',
                'resolve_inplace': True,
                'include': ['workflow2']
            },
        }
    }, result={'m_uploads': {'id_published_with_ref': {'m_entries': {'id_01': {'m_archive': {
        'workflow2': {
            'results': {'calculation_result_ref': 'm_uploads/id_published_with_ref/m_entries/id_01/m_archive/run/0/calculation/1'}},
        'run': [{
            'calculation': [
                None, {
                    'system_ref': 'm_uploads/id_published_with_ref/m_entries/id_01/m_archive/run/0/system/1',
                    'energy': {
                        'total': {
                            'value': 0.2}},
                    'dos_electronic': [
                        {'energies': [
                            0.0, 0.1,
                            1.0, 1.0,
                            1.0, 1.0,
                            1.0, 1.0,
                            1.0, 1.0,
                            1.0, 1.0,
                            1.0,
                            1.0]}]}],
            'system': [None, {
                'atoms': {
                    'labels': [
                        'H']},
                'symmetry': [{
                    'space_group_number': 221}]}]}]}}}}},
        'm_archive': {'workflow2': {'tasks': [
            {'task': 'm_uploads/id_published_with_ref/m_entries/id_01/m_archive/workflow2'}]}}})
    __entry_print('plain entry reader, resolve to root', {
        'm_archive': {
            'm_request': {
                'directive': 'resolved',
                'resolve_inplace': False,
                'include': ['workflow2']
            },
        }
    }, result={
        'm_uploads': {'id_published_with_ref': {'m_entries': {'id_01': {'m_archive': {
            'workflow2': {
                'results': {'calculation_result_ref': 'm_uploads/id_published_with_ref/m_entries/id_01/m_archive/run/0/calculation/1'}},
            'run': [{
                'calculation': [
                    None, {
                        'system_ref': 'm_uploads/id_published_with_ref/m_entries/id_01/m_archive/run/0/system/1',
                        'energy': {
                            'total': {
                                'value': 0.2}},
                        'dos_electronic': [
                            {'energies': [
                                0.0, 0.1,
                                1.0, 1.0,
                                1.0, 1.0,
                                1.0, 1.0,
                                1.0, 1.0,
                                1.0, 1.0,
                                1.0,
                                1.0]}]}],
                'system': [None, {
                    'atoms': {
                        'labels': [
                            'H']},
                    'symmetry': [{
                        'space_group_number': 221}]}]}]}}}}},
        'm_archive': {'workflow2': {
            'tasks': [{'task': 'm_uploads/id_published_with_ref/m_entries/id_01/m_archive/workflow2'}]}}})
    # __entry_print('plain entry reader, resolve inplace', {
    #     'm_archive': {
    #         'workflow[0]': {
    #             'm_def': {
    #                 'm_request': {
    #                     'directive': 'resolved',
    #                     'resolve_inplace': False,
    #                 },
    #                 'inherited_sections': '*'
    #             },
    #             'm_request': {
    #                 'directive': 'resolved',
    #                 'resolve_inplace': True,
    #                 'resolve_depth': 2,
    #                 'max_list_size': 5
    #             },
    #         }
    #     }
    # }, to_file=True)

    __entry_print('go to upload, resolve explicitly', {
        'm_request': {
            'directive': 'plain',
        },
        'upload_id': {
            'm_request': {
                'directive': 'resolved',
                'resolve_type': 'upload'
            }
        }
    }, result={
        'process_running': False, 'current_process': None, 'process_status': 'SUCCESS',
        'last_status_message': None, 'errors': [], 'warnings': [], 'complete_time': None, 'entry_id': 'id_03',
        'entry_create_time': '2023-03-05T22:16:52.440000', 'mainfile': 'mainfile_for_id_03',
        'mainfile_key': None, 'parser_name': 'parsers/vasp', 'entry_metadata': None,
        'upload_id': {
            'process_running': False, 'current_process': 'process_upload', 'process_status': 'SUCCESS',
            'last_status_message': None, 'errors': [], 'warnings': [],
            'complete_time': '2023-03-05T22:16:52.436000', 'upload_id': 'id_published_with_ref',
            'upload_name': 'name_published', 'upload_create_time': '2023-03-05T22:16:52.435000',
            'main_author': '00000000-0000-0000-0000-000000000001', 'coauthors': [], 'reviewers': [],
            'viewers': ['00000000-0000-0000-0000-000000000001'],
            'writers': ['00000000-0000-0000-0000-000000000001'], 'published': False,
            'published_to': [], 'publish_time': None, 'with_embargo': False, 'embargo_length': 0,
            'license': 'CC BY 4.0', 'entries': 6,
            'upload_files_server_path': 'id_published_with_ref'}})

    __entry_print('go to upload, resolve implicitly, resolve main author explicitly', {
        'm_request': {
            'directive': 'plain',
        },
        'upload_id': {
            'm_request': {
                'directive': 'plain',
            },
            'main_author': {
                'm_request': {
                    'directive': 'resolved',
                    'resolve_type': 'user'
                }
            }
        }
    }, result={
        'process_running': False, 'current_process': None, 'process_status': 'SUCCESS',
        'last_status_message': None, 'errors': [], 'warnings': [], 'complete_time': None, 'entry_id': 'id_03',
        'entry_create_time': '2023-03-05T22:16:52.440000', 'mainfile': 'mainfile_for_id_03',
        'mainfile_key': None, 'parser_name': 'parsers/vasp', 'entry_metadata': None,
        'upload_id': {
            'process_running': False, 'current_process': 'process_upload', 'process_status': 'SUCCESS',
            'last_status_message': None, 'errors': [], 'warnings': [],
            'complete_time': '2023-03-05T22:16:52.436000', 'upload_id': 'id_published_with_ref',
            'upload_name': 'name_published', 'upload_create_time': '2023-03-05T22:16:52.435000',
            'coauthors': [], 'reviewers': [], 'viewers': ['00000000-0000-0000-0000-000000000001'],
            'writers': ['00000000-0000-0000-0000-000000000001'], 'published': False,
            'published_to': [], 'publish_time': None, 'with_embargo': False, 'embargo_length': 0,
            'license': 'CC BY 4.0', 'entries': 6,
            'upload_files_server_path': 'id_published_with_ref',
            'main_author': {
                'name': 'Sheldon Cooper', 'first_name': 'Sheldon', 'last_name': 'Cooper',
                'email': 'sheldon.cooper@nomad-coe.eu',
                'user_id': '00000000-0000-0000-0000-000000000001', 'username': 'scooper',
                'is_admin': False, 'is_oasis_admin': True}}})

    def __fs_print(msg, required, *, result: dict = None):
        with FileSystemReader(required, user=test_user) as reader:
            if result:
                assert_dict(reader.read('id_published_with_ref'), result)
            else:
                rprint(f'\n\nExample: {next(counter)} -> {msg}:')
                rprint(required)
                rprint('output:')
                rprint(reader.read('id_published_with_ref'))

    __fs_print('plain file system reader', {
        'm_request': {
            'directive': 'plain',
        },
    }, result={
        'mainfile_for_id_01': {'path': 'mainfile_for_id_01', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
        'mainfile_for_id_02': {'path': 'mainfile_for_id_02', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
        'mainfile_for_id_03': {'path': 'mainfile_for_id_03', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
        'mainfile_for_id_04': {'path': 'mainfile_for_id_04', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
        'mainfile_for_id_05': {'path': 'mainfile_for_id_05', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
        'mainfile_for_id_06': {'path': 'mainfile_for_id_06', 'is_file': True, 'size': 3231, 'access': 'unpublished'}})

    __fs_print('go to entry', {
        'm_request': {
            'directive': 'resolved',
        },
    }, result={
        'mainfile_for_id_01': {
            'path': 'mainfile_for_id_01', 'is_file': True, 'size': 3231, 'access': 'unpublished',
            'm_entries': {
                'process_running': False, 'current_process': None,
                'process_status': 'SUCCESS', 'last_status_message': None, 'errors': [],
                'warnings': [], 'complete_time': None, 'entry_id': 'id_01',
                'entry_create_time': '2023-03-05T22:16:52.438000',
                'mainfile': 'mainfile_for_id_01', 'mainfile_key': None,
                'upload_id': 'id_published_with_ref', 'parser_name': 'parsers/vasp',
                'entry_metadata': None}},
        'mainfile_for_id_02': {
            'path': 'mainfile_for_id_02', 'is_file': True, 'size': 3231, 'access': 'unpublished',
            'm_entries': {
                'process_running': False, 'current_process': None,
                'process_status': 'SUCCESS', 'last_status_message': None, 'errors': [],
                'warnings': [], 'complete_time': None, 'entry_id': 'id_02',
                'entry_create_time': '2023-03-05T22:16:52.439000',
                'mainfile': 'mainfile_for_id_02', 'mainfile_key': None,
                'upload_id': 'id_published_with_ref', 'parser_name': 'parsers/vasp',
                'entry_metadata': None}},
        'mainfile_for_id_03': {
            'path': 'mainfile_for_id_03', 'is_file': True, 'size': 3231, 'access': 'unpublished',
            'm_entries': {
                'process_running': False, 'current_process': None,
                'process_status': 'SUCCESS', 'last_status_message': None, 'errors': [],
                'warnings': [], 'complete_time': None, 'entry_id': 'id_03',
                'entry_create_time': '2023-03-05T22:16:52.440000',
                'mainfile': 'mainfile_for_id_03', 'mainfile_key': None,
                'upload_id': 'id_published_with_ref', 'parser_name': 'parsers/vasp',
                'entry_metadata': None}},
        'mainfile_for_id_04': {
            'path': 'mainfile_for_id_04', 'is_file': True, 'size': 3231, 'access': 'unpublished',
            'm_entries': {
                'process_running': False, 'current_process': None,
                'process_status': 'SUCCESS', 'last_status_message': None, 'errors': [],
                'warnings': [], 'complete_time': None, 'entry_id': 'id_04',
                'entry_create_time': '2023-03-05T22:16:52.441000',
                'mainfile': 'mainfile_for_id_04', 'mainfile_key': None,
                'upload_id': 'id_published_with_ref', 'parser_name': 'parsers/vasp',
                'entry_metadata': None}},
        'mainfile_for_id_05': {
            'path': 'mainfile_for_id_05', 'is_file': True, 'size': 3231, 'access': 'unpublished',
            'm_entries': {
                'process_running': False, 'current_process': None,
                'process_status': 'SUCCESS', 'last_status_message': None, 'errors': [],
                'warnings': [], 'complete_time': None, 'entry_id': 'id_05',
                'entry_create_time': '2023-03-05T22:16:52.442000',
                'mainfile': 'mainfile_for_id_05', 'mainfile_key': None,
                'upload_id': 'id_published_with_ref', 'parser_name': 'parsers/vasp',
                'entry_metadata': None}},
        'mainfile_for_id_06': {
            'path': 'mainfile_for_id_06', 'is_file': True, 'size': 3231, 'access': 'unpublished',
            'm_entries': {
                'process_running': False, 'current_process': None,
                'process_status': 'SUCCESS', 'last_status_message': None, 'errors': [],
                'warnings': [], 'complete_time': None, 'entry_id': 'id_06',
                'entry_create_time': '2023-03-05T22:16:52.443000',
                'mainfile': 'mainfile_for_id_06', 'mainfile_key': None,
                'upload_id': 'id_published_with_ref', 'parser_name': 'parsers/vasp',
                'entry_metadata': None}}})

    __fs_print('go to selected entry', {
        'm_request': {
            'directive': 'plain',
        },
        'mainfile_for_id_01': {
            'm_request': {
                'directive': 'resolved',
            },
        }
    }, result={
        'mainfile_for_id_02': {'path': 'mainfile_for_id_02', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
        'mainfile_for_id_03': {'path': 'mainfile_for_id_03', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
        'mainfile_for_id_04': {'path': 'mainfile_for_id_04', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
        'mainfile_for_id_05': {'path': 'mainfile_for_id_05', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
        'mainfile_for_id_06': {'path': 'mainfile_for_id_06', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
        'mainfile_for_id_01': {
            'path': 'mainfile_for_id_01', 'is_file': True, 'size': 3231, 'access': 'unpublished',
            'm_entries': {
                'process_running': False, 'current_process': None,
                'process_status': 'SUCCESS', 'last_status_message': None, 'errors': [],
                'warnings': [], 'complete_time': None, 'entry_id': 'id_01',
                'entry_create_time': '2023-03-05T22:16:52.438000',
                'mainfile': 'mainfile_for_id_01', 'mainfile_key': None,
                'upload_id': 'id_published_with_ref', 'parser_name': 'parsers/vasp',
                'entry_metadata': None}}})

    __fs_print('go to selected entry then to upload', {
        'm_request': {
            'directive': 'plain',
        },
        'mainfile_for_id_01': {
            'm_request': {
                'directive': 'plain',
            },
            'm_entries': {
                'upload_id': {
                    'm_request': {
                        'directive': 'resolved',
                        'resolve_type': 'upload'
                    },
                }
            }
        }
    }, result={
        'mainfile_for_id_02': {'path': 'mainfile_for_id_02', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
        'mainfile_for_id_03': {'path': 'mainfile_for_id_03', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
        'mainfile_for_id_04': {'path': 'mainfile_for_id_04', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
        'mainfile_for_id_05': {'path': 'mainfile_for_id_05', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
        'mainfile_for_id_06': {'path': 'mainfile_for_id_06', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
        'mainfile_for_id_01': {
            'path': 'mainfile_for_id_01', 'is_file': True, 'size': 3231, 'access': 'unpublished',
            'm_entries': {
                'upload_id': {
                    'process_running': False, 'current_process': 'process_upload',
                    'process_status': 'SUCCESS', 'last_status_message': None, 'errors': [],
                    'warnings': [], 'complete_time': '2023-03-05T22:16:52.436000',
                    'upload_id': 'id_published_with_ref', 'upload_name': 'name_published',
                    'upload_create_time': '2023-03-05T22:16:52.435000',
                    'main_author': '00000000-0000-0000-0000-000000000001', 'coauthors': [],
                    'reviewers': [], 'viewers': ['00000000-0000-0000-0000-000000000001'],
                    'writers': ['00000000-0000-0000-0000-000000000001'],
                    'published': False, 'published_to': [], 'publish_time': None,
                    'with_embargo': False, 'embargo_length': 0, 'license': 'CC BY 4.0',
                    'entries': 6,
                    'upload_files_server_path': 'id_published_with_ref'}}}})

    __fs_print('go to selected entry then to upload, skipping file info', {
        'm_request': {
            'directive': 'plain',
        },
        'mainfile_for_id_01': {
            'm_entries': {
                'upload_id': {
                    'm_request': {
                        'directive': 'resolved',
                        'resolve_type': 'upload'
                    },
                }
            }
        }
    }, result={
        'mainfile_for_id_02': {'path': 'mainfile_for_id_02', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
        'mainfile_for_id_03': {'path': 'mainfile_for_id_03', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
        'mainfile_for_id_04': {'path': 'mainfile_for_id_04', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
        'mainfile_for_id_05': {'path': 'mainfile_for_id_05', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
        'mainfile_for_id_06': {'path': 'mainfile_for_id_06', 'is_file': True, 'size': 3231, 'access': 'unpublished'},
        'mainfile_for_id_01': {'m_entries': {
            'upload_id': {
                'process_running': False, 'current_process': 'process_upload', 'process_status': 'SUCCESS',
                'last_status_message': None, 'errors': [], 'warnings': [],
                'complete_time': '2023-03-05T22:16:52.436000', 'upload_id': 'id_published_with_ref',
                'upload_name': 'name_published', 'upload_create_time': '2023-03-05T22:16:52.435000',
                'main_author': '00000000-0000-0000-0000-000000000001', 'coauthors': [], 'reviewers': [],
                'viewers': ['00000000-0000-0000-0000-000000000001'],
                'writers': ['00000000-0000-0000-0000-000000000001'], 'published': False, 'published_to': [],
                'publish_time': None, 'with_embargo': False, 'embargo_length': 0, 'license': 'CC BY 4.0',
                'entries': 6,
                'upload_files_server_path': 'id_published_with_ref'}}}})

    __fs_print('go to selected entry then to upload, resolve user', {
        'mainfile_for_id_01': {
            'm_entries': {
                'upload_id': {
                    'main_author': {
                        'm_request': {
                            'directive': 'resolved',
                            'resolve_type': 'user'
                        },
                    }
                }
            }
        }
    }, result={'mainfile_for_id_01': {'m_entries': {'upload_id': {
        'main_author': {
            'name': 'Sheldon Cooper', 'first_name': 'Sheldon', 'last_name': 'Cooper',
            'email': 'sheldon.cooper@nomad-coe.eu', 'user_id': '00000000-0000-0000-0000-000000000001',
            'username': 'scooper', 'is_admin': False, 'is_oasis_admin': True}}}}})


# noinspection DuplicatedCode,SpellCheckingInspection
def test_general_reader(json_dict, example_data_with_reference, test_user):
    def increment():
        n = 0
        while True:
            n += 1
            yield n

    counter = increment()

    def __ge_print(msg, required, *, to_file: bool = False, result: dict = None):
        with MongoReader(required, user=test_user) as reader:
            if result:
                assert_dict(reader.read(), result)
            else:
                rprint(f'\n\nExample: {next(counter)} -> {msg}:')
                rprint(required)
                if not to_file:
                    rprint('output:')
                    rprint(reader.read())
                else:
                    with open('archive_reader_test.json', 'w') as f:
                        f.write(json.dumps(reader.read()))

    __ge_print('general start from entry', {
        'm_entries': {
            'm_request': {
                'directive': 'resolved',
                'resolve_type': 'entry',
                'pagination': {
                    'page_size': 2,
                    'page_after_value': 'id_03'
                }
            },
        }
    }, result={'m_entries': {
        'id_04': {
            'process_running': False, 'current_process': None, 'process_status': 'SUCCESS',
            'last_status_message': None, 'errors': [], 'warnings': [], 'complete_time': None,
            'entry_id': 'id_04', 'entry_create_time': '2023-03-05T22:29:55.842000',
            'mainfile': 'mainfile_for_id_04', 'mainfile_key': None,
            'upload_id': 'id_published_with_ref', 'parser_name': 'parsers/vasp',
            'entry_metadata': None},
        'id_05': {
            'process_running': False, 'current_process': None, 'process_status': 'SUCCESS',
            'last_status_message': None, 'errors': [], 'warnings': [], 'complete_time': None,
            'entry_id': 'id_05', 'entry_create_time': '2023-03-05T22:29:55.843000',
            'mainfile': 'mainfile_for_id_05', 'mainfile_key': None,
            'upload_id': 'id_published_with_ref', 'parser_name': 'parsers/vasp',
            'entry_metadata': None}}})
    __ge_print('general start from entry with wildcard', {
        'm_entries': {
            'id_01': {
                'm_request': {
                    'directive': 'plain',
                },
            },
            '*': {
                'm_request': {
                    'directive': 'plain',
                    'include': ['entry_id', 'mainfile']
                },
                "upload_id": "*"
            }
        }
    }, result={
        'm_entries': {
            'id_01': {
                'process_running': False, 'current_process': None, 'process_status': 'SUCCESS',
                'last_status_message': None, 'errors': [], 'warnings': [], 'complete_time': None,
                'entry_id': 'id_01', 'entry_create_time': '2023-03-05T22:20:46.586000',
                'mainfile': 'mainfile_for_id_01', 'mainfile_key': None,
                'upload_id': 'id_published_with_ref', 'parser_name': 'parsers/vasp',
                'entry_metadata': None},
            'id_02': {'entry_id': 'id_02', 'mainfile': 'mainfile_for_id_02', 'upload_id': 'id_published_with_ref'},
            'id_03': {'entry_id': 'id_03', 'mainfile': 'mainfile_for_id_03', 'upload_id': 'id_published_with_ref'},
            'id_04': {'entry_id': 'id_04', 'mainfile': 'mainfile_for_id_04', 'upload_id': 'id_published_with_ref'},
            'id_05': {'entry_id': 'id_05', 'mainfile': 'mainfile_for_id_05', 'upload_id': 'id_published_with_ref'},
            'id_06': {'entry_id': 'id_06', 'mainfile': 'mainfile_for_id_06', 'upload_id': 'id_published_with_ref'}}})
    __ge_print('general start from entry WITHOUT retrieval of metadata (just listing)', {
        'm_entries': {
            'm_request': {
                'directive': 'plain',
                'pagination': {
                    'page_size': 2,
                    'page': 2
                },
                'query': {
                    'owner': 'user'
                }
            },
        }
    }, result={'m_entries': {'id_03': 'id_03', 'id_04': 'id_04'}})
    __ge_print('general start from upload', {
        'm_uploads': {
            'm_request': {
                'directive': 'resolved',
                'resolve_type': 'upload',
                'pagination': {
                    'page_size': 10,
                }
            },
        }
    }, result={
        'm_uploads': {
            'id_published_with_ref': {
                'process_running': False, 'current_process': 'process_upload',
                'process_status': 'SUCCESS', 'last_status_message': None,
                'errors': [], 'warnings': [],
                'complete_time': '2023-03-05T22:20:46.584000',
                'upload_id': 'id_published_with_ref',
                'upload_name': 'name_published',
                'upload_create_time': '2023-03-05T22:20:46.583000',
                'main_author': '00000000-0000-0000-0000-000000000001',
                'coauthors': [], 'reviewers': [],
                'viewers': ['00000000-0000-0000-0000-000000000001'],
                'writers': ['00000000-0000-0000-0000-000000000001'],
                'published': False, 'published_to': [], 'publish_time': None,
                'with_embargo': False, 'embargo_length': 0,
                'license': 'CC BY 4.0', 'entries': 6,
                'upload_files_server_path': 'id_published_with_ref'}}})
    __ge_print('general start from upload with query and pagination', {
        'm_uploads': {
            'm_request': {
                'directive': 'resolved',
                'resolve_type': 'upload',
                'pagination': {
                    'page_size': 10,
                    'order_by': 'upload_create_time'
                },
                'query': {
                    'is_processing': False
                }
            },
        }
    }, result={
        'm_uploads': {
            'id_published_with_ref': {
                'process_running': False, 'current_process': 'process_upload',
                'process_status': 'SUCCESS', 'last_status_message': None,
                'errors': [], 'warnings': [],
                'complete_time': '2023-03-05T22:20:46.584000',
                'upload_id': 'id_published_with_ref',
                'upload_name': 'name_published',
                'upload_create_time': '2023-03-05T22:20:46.583000',
                'main_author': '00000000-0000-0000-0000-000000000001',
                'coauthors': [], 'reviewers': [],
                'viewers': ['00000000-0000-0000-0000-000000000001'],
                'writers': ['00000000-0000-0000-0000-000000000001'],
                'published': False, 'published_to': [], 'publish_time': None,
                'with_embargo': False, 'embargo_length': 0,
                'license': 'CC BY 4.0', 'entries': 6,
                'upload_files_server_path': 'id_published_with_ref'}}})
    __ge_print('general start from user, does NOT perform search from security', {
        'm_users': {
            'me': {
                'm_request': {
                    'directive': 'plain'
                },
            }
        }
    }, result={'m_users': {'me': {
        'name': 'Sheldon Cooper', 'first_name': 'Sheldon', 'last_name': 'Cooper',
        'email': 'sheldon.cooper@nomad-coe.eu',
        'user_id': '00000000-0000-0000-0000-000000000001', 'username': 'scooper',
        'is_admin': False, 'is_oasis_admin': True}}})
    __ge_print('general start from me, with its user id like fields resolved', {
        'm_users': {
            'me': {
                'm_request': {
                    'directive': 'resolved',
                    'resolve_type': 'user'
                },
            }
        }
    }, result={'m_users': {'me': {
        'name': 'Sheldon Cooper', 'first_name': 'Sheldon', 'last_name': 'Cooper',
        'email': 'sheldon.cooper@nomad-coe.eu',
        'user_id': {
            'name': 'Sheldon Cooper', 'first_name': 'Sheldon', 'last_name': 'Cooper',
            'email': 'sheldon.cooper@nomad-coe.eu',
            'user_id': '00000000-0000-0000-0000-000000000001', 'username': 'scooper',
            'is_admin': False, 'is_oasis_admin': True}, 'username': 'scooper',
        'is_admin': False, 'is_oasis_admin': True}}})
    # __ge_print('general start from entry to archive', {
    #     'm_entries': {
    #         'id_03': {
    #             'm_archive': {
    #                 'workflow[0]': {
    #                     'm_request': {
    #                         'directive': 'resolved',
    #                         'resolve_inplace': False,
    #                         'include_definition': 'both',
    #                     }
    #                 }
    #             }
    #         }
    #     }
    # })
    # __ge_print('general start from entry to archive', {
    #     'm_entries': {
    #         'id_01': {
    #             'm_archive': '*'
    #         }
    #     }
    # }, to_file=True)


@pytest.fixture(scope='function')
def example_data_with_reference(elastic_function, raw_files_module, mongo_function, test_user, json_dict):
    '''
    Provides a couple of entries with references.

    Only used in test_required_reader_with_remote_reference.
    '''
    data = ExampleData(main_author=test_user)

    data.create_upload(upload_id='id_published_with_ref', upload_name='name_published', published=False)

    ref_list = [
        {'results': {'calculation_result_ref': '/run/0/calculation/1'}},  # plain direct reference
        {'results': {'calculation_result_ref': '#/run/0/calculation/1'}},  # new-style reference
        {'tasks': [{
            'm_def': 'nomad.datamodel.metainfo.workflow.TaskReference',
            'task': '../entries/id_01/archive#/workflow2'}]},  # reference to another archive
        {'tasks': [{
            'm_def': 'nomad.datamodel.metainfo.workflow.TaskReference',
            'task': '../entries/id_05/archive#/workflow2'}]},  # circular reference
        {'tasks': [{
            'm_def': 'nomad.datamodel.metainfo.workflow.TaskReference',
            'task': '../entries/id_04/archive#/workflow2'}]},  # circular reference
        {'tasks': [{
            'm_def': 'nomad.datamodel.metainfo.workflow.TaskReference',
            'task': 'https://another.domain/entries/id_03/archive#/workflow2'}]}  # remote reference
    ]

    del json_dict['results']

    for index, ref in enumerate(ref_list):
        ref['m_def'] = 'nomad.datamodel.metainfo.simulation.workflow.SimulationWorkflow'
        json_dict['workflow2'] = ref
        data.create_entry(
            upload_id='id_published_with_ref',
            entry_id=f'id_{index + 1:02d}',
            entry_archive=EntryArchive.m_from_dict(json_dict))

    for archive in data.archives.values():
        archive.metadata.apply_archive_metadata(archive)

    data.save(with_files=True, with_es=True, with_mongo=True)

    yield data
    data.delete()


@pytest.fixture(scope='function')
def json_dict():
    return {
        "metadata": {
            "entry_id": "test_id",
            'upload_id': 'id_published_with_ref'
        },
        "results": {
            "properties": {
                "electronic": {
                    "dos_electronic": [{
                        "energies": "/run/0/calculation/1/dos_electronic/0/energies"
                    }]
                }
            }
        },
        "run": [
            {
                "system": [
                    {
                        "atoms": {
                            "labels": [
                                "He"
                            ]
                        },
                        "symmetry": [
                            {
                                "space_group_number": 221
                            }
                        ]
                    },
                    {
                        "atoms": {
                            "labels": [
                                "H"
                            ]
                        },
                        "symmetry": [
                            {
                                "space_group_number": 221
                            }
                        ]
                    }
                ],
                "calculation": [
                    {
                        "system_ref": "/run/0/system/1",
                        "energy": {
                            "total": {
                                "value": 0.1
                            }
                        }
                    },
                    {
                        "system_ref": "/run/0/system/1",
                        "energy": {
                            "total": {
                                "value": 0.2
                            }
                        },
                        "dos_electronic": [
                            {
                                "energies": [0.0, 0.1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
                            }
                        ],
                        "eigenvalues": [
                        ]
                    },
                    {
                        "system_ref": "/run/0/system/1",
                        "energy": {
                            "total": {
                                "value": 0.1
                            }
                        }
                    }
                ]
            }
        ],
        "workflow2": {
            "m_def": "nomad.datamodel.metainfo.simulation.workflow.SimulationWorkflow",
            "results": {
                "calculation_result_ref": "/run/0/calculation/1"
            }
        }
    }
