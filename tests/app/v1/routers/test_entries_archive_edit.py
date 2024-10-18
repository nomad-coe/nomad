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
import pytest

from nomad.datamodel.datamodel import EntryArchive, EntryMetadata
from nomad.datamodel.metainfo.basesections import BaseSection
from nomad.utils.exampledata import ExampleData
from tests.test_files import create_test_upload_files


@pytest.mark.parametrize(
    'edit, result, user',
    [
        pytest.param(
            {'changes': [{'path': 'data/name', 'new_value': 'NewName'}]},
            {'data': {'name': 'NewName'}},
            'user1',
            id='quantity',
        ),
        pytest.param(
            {'changes': [{'path': 'data/sub', 'new_value': {'name': 'NewName'}}]},
            {'data': {'name': 'TestName', 'sub': {'name': 'NewName'}}},
            'user1',
            id='sub-section',
        ),
        pytest.param(
            {'changes': [{'path': 'data/sub/0', 'new_value': {'name': 'NewName'}}]},
            {'data': {'name': 'TestName', 'sub': [{'name': 'NewName'}]}},
            'user1',
            id='repeated-sub-section',
        ),
        pytest.param(
            {'changes': [{'path': 'data/sub/name', 'new_value': 'NewName'}]},
            {'data': {'name': 'TestName', 'sub': {'name': 'NewName'}}},
            'user1',
            id='missing-sub-section',
        ),
        pytest.param(
            {'changes': [{'path': 'data/sub/0/name', 'new_value': 'NewName'}]},
            {'data': {'name': 'TestName', 'sub': [{'name': 'NewName'}]}},
            'user1',
            id='missing-repeated-sub-section',
        ),
        pytest.param(
            {'changes': [{'path': 'data/name', 'action': 'remove'}]},
            {'data': {}},
            'user1',
            id='remove-quantity',
        ),
        pytest.param(
            {
                'changes': [
                    {'path': 'data/sub/name', 'new_value': 'NewName'},
                    {'path': 'data/sub', 'action': 'remove'},
                ]
            },
            {
                'data': {
                    'name': 'TestName',
                }
            },
            'user1',
            id='remove-sub-section',
        ),
        pytest.param(
            {
                'changes': [
                    {'path': 'data/sub/1/name', 'new_value': 'NewName'},
                    {'path': 'data/sub/1', 'action': 'remove'},
                ]
            },
            {'data': {'name': 'TestName', 'sub': [None]}},
            'user1',
            id='remove-repeated-sub-section',
        ),
        pytest.param(
            {
                'changes': [
                    {'path': 'data/sub/0', 'action': 'upsert', 'new_value': {}},
                    {
                        'path': 'data/sub/0/name',
                        'action': 'upsert',
                        'new_value': 'NewName1',
                    },
                    {'path': 'data/sub/1', 'action': 'upsert', 'new_value': {}},
                    {
                        'path': 'data/sub/1/name',
                        'action': 'upsert',
                        'new_value': 'NewName2',
                    },
                ]
            },
            {
                'data': {
                    'name': 'TestName',
                    'sub': [{'name': 'NewName1'}, {'name': 'NewName2'}],
                }
            },
            'user1',
            id='add-multiple-repeated-sub-section',
        ),
    ],
)
def test_post_entry_edit(
    edit,
    result,
    user,
    client,
    auth_headers,
    users_dict,
    elastic_function,
    mongo_function,
    raw_files_function,
):
    mainfile = 'mainfile.archive.json'
    data = ExampleData(main_author=users_dict[user])
    data.create_upload(upload_id='upload_id', published=False)
    data.create_entry(entry_id='entry_id', upload_id='upload_id', mainfile=mainfile)
    data.save(with_files=False)

    upload_files = create_test_upload_files('upload_id', published=False, archives=[])
    with upload_files.raw_file(mainfile, 'wt') as f:
        json.dump(
            EntryArchive(
                metadata=EntryMetadata(
                    entry_id='entry_id',
                    mainfile=mainfile,
                ),
                data=BaseSection(name='TestName'),
            ).m_to_dict(),
            f,
        )

    user_auth = auth_headers[user]
    url = 'entries/entry_id/edit'
    response = client.post(url, headers=user_auth, json=edit)

    assert response.status_code == 200, response.text
    archive_data = None
    with upload_files.raw_file(mainfile, 'rt') as f:
        archive_data = json.load(f)

    assert json.dumps(
        {key: value for key, value in archive_data['data'].items() if key != 'm_def'}
    ) == json.dumps(result['data'])
