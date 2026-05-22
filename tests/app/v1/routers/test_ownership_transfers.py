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

import pytest

from nomad.mongo.users import OwnershipTransferRecord
from nomad.processing import Upload
from tests.app.v1.routers.common import assert_response


@pytest.mark.asyncio
async def test_upload_transfer_create_and_list_basic(
    auth_headers,
    client,
    example_data_writeable,
    users_dict,
):
    upload_id = example_data_writeable['id_unpublished_w']

    create_response = client.post(
        'ownership-transfers',
        headers=auth_headers['user1'],
        json={
            'resource_type': 'upload',
            'resource_id': upload_id,
            'target_user': users_dict['user2'].username,
            'target_user_type': 'username',
        },
    )
    assert_response(create_response, 200)
    transfer_id = create_response.json()['transfer_id']
    assert transfer_id

    list_incoming_response = client.get(
        'ownership-transfers?direction=incoming&resource_type=upload',
        headers=auth_headers['user2'],
    )
    assert_response(list_incoming_response, 200)
    incoming_items = list_incoming_response.json()['transfers']
    assert any(item['transfer_id'] == transfer_id for item in incoming_items)

    get_response = client.get(
        f'ownership-transfers/{transfer_id}?resource_type=upload',
        headers=auth_headers['user2'],
    )
    assert_response(get_response, 200)
    assert get_response.json()['resource_id'] == upload_id


@pytest.mark.parametrize('action', ['accept', 'refuse'])
@pytest.mark.asyncio
async def test_respond_to_upload_transfer_by_transfer_id(
    auth_headers,
    client,
    temporal_worker,
    example_data_writeable,
    users_dict,
    action,
):
    upload_id = example_data_writeable['id_unpublished_w']

    create_response = client.post(
        'ownership-transfers',
        headers=auth_headers['user1'],
        json={
            'resource_type': 'upload',
            'resource_id': upload_id,
            'target_user': users_dict['user2'].username,
            'target_user_type': 'username',
        },
    )
    assert_response(create_response, 200)
    transfer_id = create_response.json()['transfer_id']

    async with temporal_worker():
        response = await asyncio.to_thread(
            lambda: client.post(
                f'ownership-transfers/{transfer_id}/respond?resource_type=upload',
                headers=auth_headers['user2'],
                json={'action': action},
            )
        )

    assert_response(response, 200)
    assert (
        OwnershipTransferRecord.objects(
            resource_type=OwnershipTransferRecord.RESOURCE_TYPE_UPLOAD,
            resource_id=upload_id,
            state=OwnershipTransferRecord.STATE_PENDING,
        ).first()
        is None
    )


@pytest.mark.asyncio
async def test_cancel_upload_transfer_by_transfer_id(
    auth_headers,
    client,
    temporal_worker,
    example_data_writeable,
    users_dict,
):
    upload_id = example_data_writeable['id_unpublished_w']

    create_response = client.post(
        'ownership-transfers',
        headers=auth_headers['user1'],
        json={
            'resource_type': 'upload',
            'resource_id': upload_id,
            'target_user': users_dict['user2'].username,
            'target_user_type': 'username',
        },
    )
    assert_response(create_response, 200)
    transfer_id = create_response.json()['transfer_id']

    async with temporal_worker():
        response = await asyncio.to_thread(
            lambda: client.post(
                f'ownership-transfers/{transfer_id}/cancel?resource_type=upload',
                headers=auth_headers['user1'],
            )
        )

    assert_response(response, 200)
    assert (
        OwnershipTransferRecord.objects(
            resource_type=OwnershipTransferRecord.RESOURCE_TYPE_UPLOAD,
            resource_id=upload_id,
            state=OwnershipTransferRecord.STATE_PENDING,
        ).first()
        is None
    )


@pytest.mark.asyncio
async def test_upload_transfer_create_and_list_payload(
    auth_headers,
    client,
    example_data_writeable,
    users_dict,
):
    upload_id = example_data_writeable['id_unpublished_w']

    create_response = client.post(
        'ownership-transfers',
        headers=auth_headers['user1'],
        json={
            'resource_type': 'upload',
            'resource_id': upload_id,
            'target_user': users_dict['user2'].username,
            'target_user_type': 'username',
        },
    )
    assert_response(create_response, 200)
    created = create_response.json()
    transfer_id = created['transfer_id']
    assert created['resource_type'] == 'upload'
    assert created['resource_id'] == upload_id

    list_incoming_response = client.get(
        'ownership-transfers?direction=incoming&resource_type=upload',
        headers=auth_headers['user2'],
    )
    assert_response(list_incoming_response, 200)
    incoming_items = list_incoming_response.json()['transfers']
    assert any(item['transfer_id'] == transfer_id for item in incoming_items)

    get_response = client.get(
        f'ownership-transfers/{transfer_id}?resource_type=upload',
        headers=auth_headers['user2'],
    )
    assert_response(get_response, 200)
    assert get_response.json()['resource_id'] == upload_id


def test_transfer_unsupported_type(auth_headers, client, users_dict):
    resource_type = 'dataset'

    create_response = client.post(
        'ownership-transfers',
        headers=auth_headers['user1'],
        json={
            'resource_type': resource_type,
            'resource_id': 'dummy-dataset-id',
            'target_user': users_dict['user2'].username,
            'target_user_type': 'username',
        },
    )
    assert_response(create_response, 400)

    list_response = client.get(
        f'ownership-transfers?resource_type={resource_type}',
        headers=auth_headers['user1'],
    )
    assert_response(list_response, 400)

    get_response = client.get(
        f'ownership-transfers/dummy-transfer-id?resource_type={resource_type}',
        headers=auth_headers['user1'],
    )
    assert_response(get_response, 400)

    respond_response = client.post(
        f'ownership-transfers/dummy-transfer-id/respond?resource_type={resource_type}',
        headers=auth_headers['user1'],
        json={'action': 'refuse'},
    )
    assert_response(respond_response, 400)

    cancel_response = client.post(
        f'ownership-transfers/dummy-transfer-id/cancel?resource_type={resource_type}',
        headers=auth_headers['user1'],
    )
    assert_response(cancel_response, 400)


def test_respond_to_upload_transfer_by_transfer_id_refuse(
    auth_headers,
    client,
    example_data_writeable,
    users_dict,
):
    upload_id = example_data_writeable['id_unpublished_w']

    create_response = client.post(
        'ownership-transfers',
        headers=auth_headers['user1'],
        json={
            'resource_type': 'upload',
            'resource_id': upload_id,
            'target_user': users_dict['user2'].username,
            'target_user_type': 'username',
        },
    )
    assert_response(create_response, 200)
    transfer_id = create_response.json()['transfer_id']

    response = client.post(
        f'ownership-transfers/{transfer_id}/respond?resource_type=upload',
        headers=auth_headers['user2'],
        json={'action': 'refuse'},
    )
    assert_response(response, 200)
    assert response.json()['resource_type'] == 'upload'
    assert (
        OwnershipTransferRecord.objects(
            resource_type=OwnershipTransferRecord.RESOURCE_TYPE_UPLOAD,
            resource_id=upload_id,
            state=OwnershipTransferRecord.STATE_PENDING,
        ).first()
        is None
    )


def test_respond_to_upload_transfer_accept_releases_claim_on_workflow_failure(
    auth_headers,
    client,
    example_data_writeable,
    users_dict,
    monkeypatch,
):
    upload_id = example_data_writeable['id_unpublished_w']

    create_response = client.post(
        'ownership-transfers',
        headers=auth_headers['user1'],
        json={
            'resource_type': 'upload',
            'resource_id': upload_id,
            'target_user': users_dict['user2'].username,
            'target_user_type': 'username',
        },
    )
    assert_response(create_response, 200)
    transfer_id = create_response.json()['transfer_id']

    def _raise_transfer_failure(self, new_owner_user_id, previous_owner_user_id):
        raise RuntimeError('simulated transfer workflow failure')

    monkeypatch.setattr(Upload, 'transfer_ownership', _raise_transfer_failure)

    response = client.post(
        f'ownership-transfers/{transfer_id}/respond?resource_type=upload',
        headers=auth_headers['user2'],
        json={'action': 'accept'},
    )
    assert_response(response, 500)
    assert 'Failed to execute transfer workflow' in response.json()['detail']

    pending_record = OwnershipTransferRecord.objects(
        id=transfer_id,
        state=OwnershipTransferRecord.STATE_PENDING,
    ).first()
    assert pending_record is not None
