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

from nomad.graph.graph_reader import EntryReader, Token
from nomad.datamodel import EntryArchive
from nomad.metainfo import Section, Quantity
from nomad.utils.exampledata import ExampleData
from tests.archive.test_archive import assert_dict
from nomad.datamodel.metainfo import simulationworkflowschema, SCHEMA_IMPORT_ERROR

# try:
#     from rich.pretty import pprint
# except ImportError:
#     def pprint(x):
#         print(x)


def assert_path_exists(path, response):
    for segment in path:
        response = response[segment]
    # empty container implies wrong id
    if isinstance(response, dict):
        response.pop('m_errors', None)
    if not response:
        raise KeyError


def test_graph_query_random(auth_headers, client, example_data):
    user_auth = auth_headers['user1']
    response = client.post(
        'graph/raw_query',
        json={
            Token.UPLOADS: {
                'id_embargo': {
                    Token.RAW: {
                        'test_content': {
                            'id_embargo_1': {'m_request': {'directive': 'plain'}}
                        }
                    },
                    Token.ENTRIES: {
                        'id_embargo_1': {
                            'mainfile': {'m_request': {'directive': 'plain'}}
                        }
                    },
                },
            }
        },
        headers={'Accept': 'application/json'} | (user_auth if user_auth else {}),
    )

    print(response.json())


@pytest.mark.parametrize(
    'upload_id,entry_id,user,status_code',
    [
        pytest.param('id_embargo', 'id_embargo_1', 'user1', 200, id='ok'),
        pytest.param(
            'id_child_entries',
            'id_child_entries_child1',
            'user1',
            200,
            id='child-entry',
        ),
        pytest.param('id_embargo', 'id_embargo_1', 'user0', 200, id='admin-access'),
        pytest.param('id_embargo', 'id_embargo_1', None, 401, id='no-credentials'),
        pytest.param(
            'id_embargo', 'id_embargo_1', 'invalid', 401, id='invalid-credentials'
        ),
        pytest.param('id_embargo', 'id_embargo_1', 'user2', 404, id='no-access'),
        pytest.param(
            'silly_value', 'id_embargo_1', 'user1', 404, id='invalid-upload_id'
        ),
        pytest.param('id_embargo', 'silly_value', 'user1', 404, id='invalid-entry_id'),
    ],
)
def test_graph_query(
    auth_headers, client, example_data, upload_id, entry_id, user, status_code
):
    user_auth = auth_headers[user]
    response = client.post(
        'graph/query',
        json={Token.UPLOADS: {upload_id: {Token.ENTRIES: {entry_id: '*'}}}},
        headers={'Accept': 'application/json'} | (user_auth if user_auth else {}),
    )
    target_path = (Token.UPLOADS, upload_id, Token.ENTRIES, entry_id)
    if 200 == status_code:
        assert_path_exists(target_path, response.json())
    elif 401 == status_code:
        assert response.status_code == 401
    else:
        with pytest.raises(KeyError):
            assert_path_exists(target_path, response.json())


@pytest.mark.parametrize(
    'required,status_code',
    [
        pytest.param('*', 200, id='full'),
        pytest.param({'metadata': '*'}, 200, id='partial'),
        pytest.param({'run': {'system[NOTANINT]': '*'}}, 422, id='bad-required-1'),
        pytest.param(
            {'metadata': {'viewers[NOTANINT]': '*'}}, 422, id='bad-required-2'
        ),
        pytest.param({'DOESNOTEXIST': '*'}, 422, id='bad-required-3'),
    ],
)
def test_graph_archive_query(auth_headers, client, example_data, required, status_code):
    response = client.post(
        'graph/archive/query',
        headers=auth_headers['user1'],
        json={'owner': 'user', 'required': required},
    )

    json_response = response.json()

    m_response = json_response['search'].pop('m_response', None)

    assert m_response['pagination']['total'] == 32

    for key, value in json_response[Token.SEARCH].items():
        value = value[Token.ENTRIES]
        if key == 'id_02':
            # since id_02 is missing (see example_data)
            assert 'm_errors' in value
        elif status_code != 200:
            # nonexistent required field
            assert 'm_errors' in value[Token.ARCHIVE]
        else:
            archive = value[Token.ARCHIVE]
            if required == '*':
                for item in ['metadata', 'run']:
                    assert item in archive
            else:
                for item in required:
                    assert item in archive


@pytest.fixture(scope='function')
def example_archive():
    json_dict = {
        'metadata': {'entry_id': 'test_id', 'upload_id': 'test_id'},
        'results': {
            'material': {'material_id': 'random_id'},
            'properties': {
                'electronic': {
                    'dos_electronic': [
                        {'energies': '/run/0/calculation/1/dos_electronic/0/energies'}
                    ]
                }
            },
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
                        'dos_electronic': [{'energies': [0.0, 0.1]}],
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

    archive = EntryArchive.m_from_dict(json_dict)
    assert archive.run is not None
    assert len(archive.run) == 1
    return archive


@pytest.fixture(scope='function')
def example_upload(example_archive, user1, mongo_function, elastic_function):
    data = ExampleData(main_author=user1)
    data.create_upload(
        upload_id='test_id', upload_name='name_published', published=True
    )
    data.create_entry(
        upload_id='test_id', entry_id='test_id', entry_archive=example_archive
    )
    data.save(with_files=True, with_es=True, with_mongo=True)

    yield data

    data.delete()


@pytest.mark.parametrize(
    'kwargs',
    [
        pytest.param(
            dict(
                expected_upload_ids=[
                    'id_embargo',
                    'id_embargo_w_coauthor',
                    'id_embargo_w_reviewer',
                    'id_unpublished',
                    'id_unpublished_w_coauthor',
                    'id_unpublished_w_reviewer',
                    'id_published',
                    'id_child_entries',
                    'id_processing',
                    'id_empty',
                ],
            ),
            id='no-args',
        ),
        pytest.param(
            dict(
                user='user2',
                expected_upload_ids=[
                    'id_embargo_w_coauthor',
                    'id_embargo_w_reviewer',
                    'id_unpublished_w_coauthor',
                    'id_unpublished_w_reviewer',
                ],
            ),
            id='user2',
        ),
        pytest.param(dict(user=None, expected_status_code=401), id='no-credentials'),
        pytest.param(
            dict(user='invalid', expected_status_code=401), id='invalid-credentials'
        ),
        pytest.param(
            dict(
                query_params={'is_processing': True},
                expected_upload_ids=['id_processing'],
            ),
            id='filter-is_processing-True',
        ),
        pytest.param(
            dict(
                query_params={'is_processing': False},
                expected_upload_ids=[
                    'id_embargo',
                    'id_embargo_w_coauthor',
                    'id_embargo_w_reviewer',
                    'id_unpublished',
                    'id_unpublished_w_coauthor',
                    'id_unpublished_w_reviewer',
                    'id_published',
                    'id_child_entries',
                    'id_empty',
                ],
            ),
            id='filter-is_processing-False',
        ),
        pytest.param(
            dict(
                query_params={'is_published': True},
                expected_upload_ids=[
                    'id_embargo',
                    'id_embargo_w_coauthor',
                    'id_embargo_w_reviewer',
                    'id_published',
                ],
            ),
            id='filter-is_published-True',
        ),
        pytest.param(
            dict(
                query_params={'is_published': False},
                expected_upload_ids=[
                    'id_unpublished',
                    'id_unpublished_w_coauthor',
                    'id_unpublished_w_reviewer',
                    'id_child_entries',
                    'id_processing',
                    'id_empty',
                ],
            ),
            id='filter-is_published-False',
        ),
        pytest.param(
            dict(
                query_params={'upload_id': ['id_published']},
                expected_upload_ids=['id_published'],
            ),
            id='filter-upload_id-single',
        ),
        pytest.param(
            dict(
                query_params={'upload_id': ['id_published', 'id_embargo']},
                expected_upload_ids=['id_embargo', 'id_published'],
            ),
            id='filter-upload_id-multiple',
        ),
        pytest.param(
            dict(
                query_params={'upload_name': ['name_published']},
                expected_upload_ids=['id_published'],
            ),
            id='filter-upload_name-single',
        ),
        pytest.param(
            dict(
                query_params={'upload_name': ['name_published', 'name_embargo']},
                expected_upload_ids=['id_embargo', 'id_published'],
            ),
            id='filter-upload_name-multiple',
        ),
        pytest.param(
            dict(
                pagination={'page_size': 2},
                expected_upload_ids=['id_embargo', 'id_embargo_w_coauthor'],
            ),
            id='pag-page-1',
        ),
        pytest.param(
            dict(
                pagination={'page_size': 2, 'page': 2},
                expected_upload_ids=['id_embargo_w_reviewer', 'id_unpublished'],
            ),
            id='pag-page-2',
        ),
        pytest.param(
            dict(
                pagination={'page_size': 4, 'page': 3},
                expected_upload_ids=['id_processing', 'id_empty'],
            ),
            id='pag-page-3',
        ),
        pytest.param(
            dict(
                pagination={'page_size': 5, 'page': 3},
                expected_upload_ids=[],
            ),
            id='pag-page-out-of-range',
        ),
        pytest.param(
            dict(
                pagination={'page_size': 2, 'order': 'desc'},
                expected_upload_ids=['id_empty', 'id_processing'],
            ),
            id='pag-page-order-desc',
        ),
        pytest.param(
            dict(pagination={'order_by': 'upload_id'}, expected_status_code=422),
            id='pag-invalid-order_by',
        ),
    ],
)
@pytest.mark.skipif(simulationworkflowschema is None, reason=SCHEMA_IMPORT_ERROR)
def test_get_uploads_graph(auth_headers, client, example_data, kwargs):
    user = kwargs.get('user', 'user1')
    query_params = kwargs.get('query_params', None)
    pagination = kwargs.get('pagination', None)
    expected_status_code = kwargs.get('expected_status_code', 200)
    expected_upload_ids = kwargs.get('expected_upload_ids', None)
    user_auth = auth_headers[user]

    query_body = {Token.UPLOADS: {'m_request': {}}}

    if query_params is None and pagination is None:
        # noinspection PyTypedDict
        query_body[Token.UPLOADS]['m_request'] = {'directive': '*'}
    else:
        if query_params is not None:
            query_body[Token.UPLOADS]['m_request'].setdefault('query', query_params)
        if pagination is not None:
            query_body[Token.UPLOADS]['m_request'].setdefault('pagination', pagination)

    response = client.post(
        'graph/query',
        json=query_body,
        headers={'Accept': 'application/json'} | (user_auth if user_auth else {}),
    )
    result = response.json()

    def assert_upload_ids(a, b):
        assert set(a.keys()) == set(b)

    if expected_status_code == 200:
        result[Token.UPLOADS].pop('m_response', None)
        if expected_upload_ids:
            assert_upload_ids(result[Token.UPLOADS], expected_upload_ids)
        else:
            assert result[Token.UPLOADS] == {}
    else:
        assert response.status_code == expected_status_code


@pytest.mark.parametrize(
    'depth,result',
    [
        pytest.param(
            1,
            {
                'uploads': {
                    'id_unpublished': {
                        'files': {
                            'm_is': 'Directory',
                            'm_response': {
                                'depth': 1,
                                'include': ['*'],
                                'pagination': {
                                    'order': 'asc',
                                    'page': 1,
                                    'page_size': 10,
                                    'total': 1,
                                },
                            },
                            'test_content': {
                                'm_is': 'Directory',
                            },
                        }
                    }
                }
            },
            id='depth-1',
        ),
        pytest.param(
            2,
            {
                'uploads': {
                    'id_unpublished': {
                        'files': {
                            'm_is': 'Directory',
                            'm_response': {
                                'depth': 2,
                                'include': ['*'],
                                'pagination': {
                                    'order': 'asc',
                                    'page': 1,
                                    'page_size': 10,
                                    'total': 2,
                                },
                            },
                            'test_content': {
                                'id_unpublished_1': {
                                    'm_is': 'Directory',
                                },
                                'm_is': 'Directory',
                            },
                        }
                    }
                }
            },
            id='depth-2',
        ),
        pytest.param(
            3,
            {
                'uploads': {
                    'id_unpublished': {
                        'files': {
                            'm_is': 'Directory',
                            'm_response': {
                                'depth': 3,
                                'include': ['*'],
                                'pagination': {
                                    'order': 'asc',
                                    'page': 1,
                                    'page_size': 10,
                                    'total': 7,
                                },
                            },
                            'test_content': {
                                'id_unpublished_1': {
                                    '1.aux': {
                                        'm_is': 'File',
                                        'path': 'test_content/id_unpublished_1/1.aux',
                                        'size': 8,
                                    },
                                    '2.aux': {
                                        'm_is': 'File',
                                        'path': 'test_content/id_unpublished_1/2.aux',
                                        'size': 8,
                                    },
                                    '3.aux': {
                                        'm_is': 'File',
                                        'path': 'test_content/id_unpublished_1/3.aux',
                                        'size': 8,
                                    },
                                    '4.aux': {
                                        'm_is': 'File',
                                        'path': 'test_content/id_unpublished_1/4.aux',
                                        'size': 8,
                                    },
                                    'm_is': 'Directory',
                                    'mainfile.json': {
                                        'm_is': 'File',
                                        'path': 'test_content/id_unpublished_1/mainfile.json',
                                        'size': 3227,
                                    },
                                },
                                'm_is': 'Directory',
                            },
                        }
                    }
                }
            },
            id='depth-3',
        ),
    ],
)
def test_fs_graph(auth_headers, client, example_data, depth, result):
    user_auth = auth_headers['user1']

    response = client.post(
        'graph/query',
        json={
            Token.UPLOADS: {
                'id_unpublished': {
                    Token.RAW: {
                        'm_request': {
                            'directive': 'plain',
                            'depth': depth,
                        }
                    }
                }
            }
        },
        headers={'Accept': 'application/json'} | (user_auth if user_auth else {}),
    )

    assert response.json() == result


@pytest.mark.parametrize(
    'required,error',
    [
        pytest.param('include', False, id='include-all'),
        pytest.param('*', False, id='include-all-alias'),
        # pytest.param({'metadata': '*'}, False, id='include-sub-section'),
        # pytest.param({'metadata': {'entry_id': '*'}}, False, id='include-quantity'),
        pytest.param(
            {
                'workflow2': {
                    'results': {'calculation_result_ref': {'energy': {'total': '*'}}}
                }
            },
            None,
            id='resolve-with-required',
        ),
        pytest.param(
            {'workflow2': {'results': {'calculation_result_ref': 'include-resolved'}}},
            False,
            id='resolve-with-directive',
        ),
        pytest.param(
            {'workflow2': 'include-resolved', 'results': 'include-resolved'},
            False,
            id='include-resolved',
        ),
        pytest.param(
            {
                'results': {
                    'properties': {
                        'electronic': {
                            'dos_electronic': {'energies': 'include-resolved'}
                        }
                    }
                }
            },
            False,
            id='resolve-quantity-ref',
        ),
        pytest.param(
            {'metadata': {'entry_id': {'doesnotexist': '*'}}}, True, id='not-a-section'
        ),
    ],
)
@pytest.mark.skipif(simulationworkflowschema is None, reason=SCHEMA_IMPORT_ERROR)
def test_entry_reader_with_reference(
    example_archive, required, error, user1, example_upload
):
    with EntryReader({Token.ARCHIVE: required}, user=user1) as reader:
        results = reader.sync_read('test_id')

    if error:
        assert 'm_errors' in results[Token.ARCHIVE]
    elif required in ('include', '*'):
        assert_dict(example_archive.m_to_dict(), results[Token.ARCHIVE])
    else:
        full_archive = example_archive.m_to_dict()
        assert isinstance(required, dict)

        def _goto(path, root):
            if isinstance(path, str) and '/' in path:
                current = root
                for item in [i for i in path.split('/') if i]:
                    current = current[int(item) if item.isdigit() else item]
                return current
            return path

        def assert_equal(a, b):
            if isinstance(a, list):
                assert len(a) == len(b)
                for i in range(len(a)):
                    assert_equal(a[i], b[i])
            elif isinstance(a, dict):
                a.pop('m_def', None)
                assert set(a.keys()) == set(b.keys())
                for key in a.keys():
                    assert_equal(a[key], b[key])
            elif not isinstance(a, str) and not isinstance(b, str):
                assert a == b

        def _walk(a, b, root):
            for key, value in root.items():
                key = int(key) if key.isdigit() else key
                # simplified, just consider workflow[0]
                a_child = _goto(
                    a[key] if isinstance(a, dict) else a[0][key], full_archive
                )
                b_child = _goto(b[key] if isinstance(b, dict) else b[0][key], results)

                if not isinstance(value, str):
                    _walk(a_child, b_child, value)
                else:
                    assert_equal(a_child, b_child)

        _walk(full_archive, results[Token.ARCHIVE], required)


@pytest.mark.parametrize(
    'required,result',
    [
        pytest.param(
            {'run[0]': {'m_request': {'directive': 'plain', 'include': ['c*']}}},
            {
                'run': [
                    {
                        'calculation': [
                            {
                                'energy': {'total': {'value': 0.1}},
                                'system_ref': '/run/0/system/1',
                            },
                            {
                                'dos_electronic': [{'energies': [0.0, 0.1]}],
                                'energy': {'total': {'value': 0.2}},
                                'system_ref': '/run/0/system/1',
                            },
                            {
                                'energy': {'total': {'value': 0.1}},
                                'system_ref': '/run/0/system/1',
                            },
                        ]
                    }
                ]
            },
            id='include-pattern',
        ),
    ],
)
def test_graph_query_archive_functionality(
    auth_headers, client, example_upload, required, result
):
    user_auth = auth_headers['user1']
    response = client.post(
        'graph/query',
        json={'uploads': {'test_id': {'entries': {'test_id': {'archive': required}}}}},
        headers={'Accept': 'application/json'} | (user_auth if user_auth else {}),
    )

    assert (
        response.json()['uploads']['test_id']['entries']['test_id']['archive'] == result
    )
