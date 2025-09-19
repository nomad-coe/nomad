from typing import Literal
from unittest.mock import MagicMock, Mock

import pytest
from fastapi.testclient import TestClient

from nomad.config import config
from nomad.processing import Upload
from nomad.processing.base import ProcessStatus
from tests.app.v1.routers.uploads.test_basic_uploads import (
    assert_processing,
    block_until_completed,
)

from .common import assert_upload


@pytest.fixture(autouse=True)
def setup_for_transfer_bundle(request, monkeypatch, mongo_function):
    disable_health_check_patch = request.node.get_closest_marker(
        'disable_health_check_patch'
    )
    if not disable_health_check_patch:
        monkeypatch.setattr(
            'nomad.app.v1.routers.uploads._check_external_deployment_status',
            MagicMock(),
        )
    enable_target_deployment_url_validation = request.node.get_closest_marker(
        'enable_target_deployment_url_validation'
    )
    if not enable_target_deployment_url_validation:
        monkeypatch.setattr(
            'nomad.app.v1.routers.uploads._validate_target_deployment_url',
            MagicMock(),
        )


def _perform_transfer_request(
    upload_id,
    client: TestClient,
    request_auth,
    embargo_length: int | None = None,
    target_deployment_url: str | None = None,
    target_deployment_token: str | None = None,
):
    transfer_config = {
        'auth_token': target_deployment_token,
        'target_deployment_url': target_deployment_url,
        'embargo_length': embargo_length,
    }
    transfer_config = {
        key: value for key, value in transfer_config.items() if value is not None
    }
    response = client.post(
        f'uploads/{upload_id}/action/transfer',
        headers=request_auth,
        json=transfer_config,
    )
    body = response.json()
    return response, body


def _get_token(auth_headers, user):
    token = auth_headers[user]['Authorization'].split(' ')[1]
    return token


def _compare_entries_meta_info(old_upload, new_upload, embargo_length):
    old_entry = old_upload.successful_entries[0]
    new_entry = new_upload.successful_entries[0]
    old_entry_metadata_dict = old_entry.full_entry_metadata(old_upload).m_to_dict()
    new_entry_metadata_dict = new_entry.full_entry_metadata(new_upload).m_to_dict()
    for k, v in old_entry_metadata_dict.items():
        if k == 'with_embargo':
            assert new_entry_metadata_dict[k] == (embargo_length > 0)
        elif k not in (
            'upload_id',
            'entry_id',
            'upload_create_time',
            'entry_create_time',
            'last_processing_time',
            'publish_time',
            'embargo_length',
            'n_quantities',
            'quantities',
        ):
            assert new_entry_metadata_dict[k] == v, f'Metadata not matching: {k}'


def _check_success_transfer_upload(
    response, client, upload_id, suffix, user_auth, embargo_length
):
    old_upload = Upload.get(upload_id)
    expected_status_code = 200
    assert response.status_code == expected_status_code
    upload = assert_upload(response.json())
    assert upload['current_process'] == '_publish_externally'
    assert upload['process_running']

    assert_processing(client, upload_id, user_auth, published=old_upload.published)
    assert_processing(
        client, upload_id + suffix, user_auth, published=old_upload.published
    )
    old_upload = Upload.get(upload_id)
    new_upload = Upload.get(upload_id + suffix)
    assert len(old_upload.successful_entries) == len(new_upload.successful_entries) == 1

    _compare_entries_meta_info(old_upload, new_upload, embargo_length)
    assert old_upload.published_to[0] == config.oasis.central_nomad_deployment_url
    assert new_upload.from_oasis and new_upload.oasis_deployment_url
    assert new_upload.embargo_length == embargo_length
    assert (
        new_upload.upload_files.access == 'restricted'  # type: ignore
        if embargo_length > 0
        else 'public'
    )


def _request_transfer_start(
    auth_headers,
    client: TestClient,
    oasis_publishable_upload: tuple[str, Literal['_v2']],
    embargo_length: int = 0,
    user: str = 'user0',
    target_deployment_user: str | None = None,
    check_success: bool = False,
    target_deployment_url: str | None = None,
):
    """
    Hit the endpoint to start a transfer and check success if required
    params:
        user: Used to authorize the request.
        target_deployment_user: Will be used to authorize the transfer in the target deployment. If not provided, it will use the same user of the request.
        check_success: Verify if the whole transfer process is successfull. This includes waiting for the workflow to finish and compare internal variables to check the integrity of the transfer
    """
    upload_id, suffix = oasis_publishable_upload
    user_auth = auth_headers[user]
    target_deployment_token = _get_token(auth_headers, target_deployment_user or user)
    response, body = _perform_transfer_request(
        upload_id,
        client,
        request_auth=user_auth,
        embargo_length=embargo_length,
        target_deployment_url=target_deployment_url,
        target_deployment_token=target_deployment_token,
    )
    if check_success:
        _check_success_transfer_upload(
            response, client, upload_id, suffix, user_auth, embargo_length
        )
    return response, body


@pytest.mark.parametrize(
    'embargo_length, expected_response_code',
    [
        pytest.param(-10, 422, id='embargo_length=-10'),
        pytest.param(0, 200, id='embargo_length=0'),
        pytest.param(5, 200, id='embargo_length=5'),
        pytest.param(36, 200, id='embargo_length=36'),
        pytest.param(40, 422, id='embargo_length=40'),
    ],
)
def test_embargo_length(
    auth_headers,
    client: TestClient,
    oasis_publishable_upload: tuple[str, Literal['_v2']],
    embargo_length: int,
    expected_response_code: int,
):
    response, body = _request_transfer_start(
        auth_headers,
        client,
        oasis_publishable_upload,
        embargo_length,
        check_success=expected_response_code < 400,
    )

    assert response.status_code == expected_response_code
    if expected_response_code >= 400:
        assert len(body['detail']) > 0  # Check error message info


def _check_workflow_failure(
    response, client, upload_id, user_auth, error_messages: list[str]
):
    """
    Waits until the workflow fails and check that the error messages
    are being stored in the upload
    params:
        error_messages: The messages to be checked if exist in the upload
    """
    # The workflow should successfully start
    assert response.status_code == 200

    # Check that the workflow endup failing
    old_upload_data = block_until_completed(client, upload_id, user_auth)
    assert old_upload_data['process_status'] == ProcessStatus.FAILURE
    assert len(old_upload_data['errors']) > 0
    upload_error = old_upload_data['errors'][0]
    for expected_error in error_messages:
        assert expected_error in upload_error


def test_bad_formatted_token(
    auth_headers,
    client: TestClient,
    oasis_publishable_upload: tuple[str, Literal['_v2']],
):
    upload_id, _ = oasis_publishable_upload
    auth = auth_headers['user0']

    response, _ = _perform_transfer_request(
        upload_id, client, request_auth=auth, target_deployment_token='abcdef'
    )
    _check_workflow_failure(
        response,
        client,
        upload_id,
        user_auth=auth,
        error_messages=[
            'Error message from external deployment',
            'user does not exist',
        ],
    )


def test_invalid_token(
    auth_headers,
    client: TestClient,
    oasis_publishable_upload: tuple[str, Literal['_v2']],
):
    upload_id, _ = oasis_publishable_upload
    user = 'user1'
    user_auth = auth_headers[user]

    response, _ = _request_transfer_start(
        auth_headers,
        client,
        oasis_publishable_upload,
        0,
        user,
        target_deployment_user='invalid',
        check_success=False,
    )
    _check_workflow_failure(
        response,
        client,
        upload_id,
        user_auth,
        error_messages=[
            'Error message from external deployment',
            'user does not exist',
        ],
    )


def test_workflow_failed(
    auth_headers,
    client: TestClient,
    oasis_publishable_upload: tuple[str, Literal['_v2']],
    monkeypatch,
):
    """
    Test a success transfer start and there is an internal
    failure in the workflow.
    The expected behavior is the following:
    1. Transfer start request answer with 200 (OK).
    2. Wait until the process is finished.
    3. The internal process fails with some message.
    4. The information about the failure is stored in the upload.
    """

    error_message = 'test error message'

    monkeypatch.setattr(
        'nomad.processing.data.Upload._publish_externally_local',
        Mock(side_effect=Exception(error_message)),
    )
    upload_id, _ = oasis_publishable_upload
    user = 'user0'
    response, _ = _request_transfer_start(
        auth_headers,
        client,
        oasis_publishable_upload,
        user=user,
        check_success=False,
    )
    _check_workflow_failure(
        response,
        client,
        upload_id,
        user_auth=auth_headers[user],
        error_messages=[error_message],
    )


@pytest.mark.parametrize(
    'user',
    [
        'user0',  # admin
        'user1',  # oasis admin
        'user2',  # normal user
    ],
)
def test_different_user_roles(
    auth_headers,
    client: TestClient,
    oasis_publishable_upload: tuple[str, Literal['_v2']],
    user,
):
    _request_transfer_start(
        auth_headers, client, oasis_publishable_upload, 0, user, check_success=True
    )


def test_token_not_provided(
    auth_headers,
    client: TestClient,
    oasis_publishable_upload: tuple[str, Literal['_v2']],
):
    upload_id, _ = oasis_publishable_upload
    response, body = _perform_transfer_request(
        upload_id,
        client,
        request_auth=auth_headers['user0'],
        target_deployment_token=None,
    )
    assert response.status_code == 422
    assert len(body['detail']) > 0


@pytest.mark.disable_health_check_patch(True)
def test_external_deployment_health_failed(
    auth_headers,
    client: TestClient,
    oasis_publishable_upload: tuple[str, Literal['_v2']],
    monkeypatch,
):
    error_message = 'generic message'

    monkeypatch.setattr(
        'nomad.app.v1.routers.uploads._perform_status_check',
        Mock(side_effect=Exception(error_message)),
    )
    response, body = _request_transfer_start(
        auth_headers, client, oasis_publishable_upload, check_success=False
    )
    assert 'detail' in body
    assert 'Failed to check external deployment health' in body['detail']
    assert error_message in body['detail']
    assert response.status_code == 400


@pytest.mark.enable_target_deployment_url_validation(True)
@pytest.mark.parametrize(
    'target_url, expected_message',
    [
        pytest.param('abcde', 'URL must start with http:// or https://', id='no-http'),
        pytest.param('http://', 'URL must contain a valid host', id='no-hostname'),
        pytest.param(
            'http://google.com', "URL path must end with '/api'", id='bad-ending'
        ),
    ],
)
def test_invalid_target_url(
    auth_headers,
    client: TestClient,
    oasis_publishable_upload: tuple[str, Literal['_v2']],
    target_url,
    expected_message,
):
    response, body = _request_transfer_start(
        auth_headers,
        client,
        oasis_publishable_upload,
        check_success=False,
        target_deployment_url=target_url,
    )
    assert expected_message in body['detail']
    assert response.status_code == 422


def test_default_target_url(
    auth_headers,
    client: TestClient,
    oasis_publishable_upload: tuple[str, Literal['_v2']],
):
    upload_id, _ = oasis_publishable_upload
    _request_transfer_start(
        auth_headers,
        client,
        oasis_publishable_upload,
        check_success=True,
    )
    old_upload = Upload.get(upload_id)
    assert len(old_upload.published_to) == 1
    assert old_upload.published_to[0] == config.oasis.central_nomad_deployment_url


def test_transfer_processing_upload(
    auth_headers,
    client: TestClient,
    oasis_publishable_upload: tuple[str, Literal['_v2']],
):
    """
    Transfer an upload that is being processed should fail
    """
    upload_id, _ = oasis_publishable_upload
    upload = Upload.get(upload_id)
    upload.process_upload()
    response, body = _request_transfer_start(
        auth_headers,
        client,
        oasis_publishable_upload,
        check_success=False,
    )

    assert 'detail' in body
    assert (
        body['detail']
        == 'The upload is currently being processed, operation not allowed.'
    )
    assert response.status_code == 400


def test_non_published_upload(
    auth_headers, client: TestClient, non_empty_processed: Upload
):
    """
    Non published uploads should not be able to transferred
    """
    upload_id = non_empty_processed.upload_id
    response, body = _perform_transfer_request(
        upload_id,
        client,
        request_auth=auth_headers['user0'],
        target_deployment_token=_get_token(auth_headers, 'user0'),
    )
    assert response.status_code == 400
    assert body['detail'] == 'The upload should be published first.'


def test_transfer_duplicated_upload(
    auth_headers,
    client: TestClient,
    oasis_publishable_upload: tuple[str, Literal['_v2']],
):
    upload_id, _ = oasis_publishable_upload
    user = 'user0'
    user_auth = auth_headers[user]
    _request_transfer_start(
        auth_headers,
        client,
        oasis_publishable_upload,
        user=user,
        check_success=True,
    )

    # Trying to transfer again should fail
    response, _ = _request_transfer_start(
        auth_headers,
        client,
        oasis_publishable_upload,
        user=user,
        check_success=False,
    )
    _check_workflow_failure(
        response,
        client,
        upload_id,
        user_auth,
        error_messages=[
            'Error message from external deployment',
            'Failed to import bundle: Upload with id',
            'already exists',
        ],
    )
