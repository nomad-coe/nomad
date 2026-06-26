from __future__ import annotations

from collections.abc import Iterable

from nomad.files import StagingUploadFiles, PublicUploadFiles
from nomad.infrastructure import mongo_client, setup_mongo
from nomad.processing.data import Upload


def _as_unique_user_ids(user_ids: str | Iterable[str]) -> list[str]:
    if isinstance(user_ids, str):
        return [user_ids]
    return list(dict.fromkeys(user_ids))


def _resolve_upload(
    *, upload: Upload | None = None, upload_id: str | None = None
) -> Upload:
    if upload is not None:
        return upload
    if upload_id is None:
        raise ValueError('Either upload or upload_id must be provided.')

    from nomad.processing import Upload as UploadModel

    resolved_upload = UploadModel.get(upload_id)
    if resolved_upload is None:
        raise KeyError(f'Upload {upload_id} was not found.')
    return resolved_upload


def add_upload_reviewers(
    reviewer_user_ids: str | Iterable[str],
    *,
    upload: Upload | None = None,
    upload_id: str | None = None,
) -> int:
    """Add reviewer user ids to an upload and persist only when changed."""
    resolved_upload = _resolve_upload(upload=upload, upload_id=upload_id)
    desired_user_ids = _as_unique_user_ids(reviewer_user_ids)
    if not desired_user_ids:
        return 0

    reviewers = list(resolved_upload.reviewers or [])
    reviewers_set = set(reviewers)
    added = 0

    for user_id in desired_user_ids:
        if user_id in reviewers_set:
            continue
        reviewers.append(user_id)
        reviewers_set.add(user_id)
        added += 1

    if added > 0:
        resolved_upload.reviewers = reviewers
        resolved_upload.save()

    return added


def remove_upload_reviewers(
    reviewer_user_ids: str | Iterable[str],
    *,
    upload: Upload | None = None,
    upload_id: str | None = None,
) -> int:
    """Remove reviewer user ids from an upload and persist only when changed."""
    resolved_upload = _resolve_upload(upload=upload, upload_id=upload_id)
    remove_user_ids = set(_as_unique_user_ids(reviewer_user_ids))
    if not remove_user_ids:
        return 0

    reviewers = list(resolved_upload.reviewers or [])
    updated_reviewers = [
        reviewer for reviewer in reviewers if reviewer not in remove_user_ids
    ]
    removed = len(reviewers) - len(updated_reviewers)

    if removed > 0:
        resolved_upload.reviewers = updated_reviewers
        resolved_upload.save()

    return removed


def _get_upload(upload_id: str) -> Upload:
    """
    Retrieve an upload document without applying any access checks.

    Args:
        upload_id: The unique identifier for the upload.

    Returns:
        The matching upload document.

    Raises:
        AttributeError: If no upload exists for the given id.
    """

    if mongo_client is None:
        setup_mongo()

    upload = Upload.get(upload_id)

    if upload is None:
        raise AttributeError(f'No upload found with id: {upload_id}')

    return upload


def _check_upload_access(upload_id, user_id: str) -> bool:
    """
    Check whether the given user can access the upload.

    Access is currently limited to the upload's main author and direct coauthors.

    Args:
        upload_id: The unique identifier for the upload.
        user_id: The unique identifier for the user requesting access.

    Returns:
        True if the user is authorized to access the upload, otherwise False.
    """

    upload = _get_upload(upload_id)

    is_coauthor = isinstance(upload.coauthors, list) and user_id in upload.coauthors
    is_authorized = upload.main_author == user_id or is_coauthor

    if not is_authorized:
        return False

    return True


def get_upload(upload_id: str, user_id: str) -> Upload:
    """
    Retrieve an upload after verifying user authorization.

    Args:
        upload_id: The unique identifier for the upload.
        user_id: The unique identifier for the user.

    Returns:
        The matching upload document.

    Raises:
        PermissionError: If the upload exists but the user is not authorized.
    """
    if not _check_upload_access(upload_id, user_id):
        raise PermissionError(
            f'User {user_id} is not authorized to access upload {upload_id}.'
        )

    upload = _get_upload(upload_id)

    return upload


def get_upload_files(
    upload_id: str, user_id: str
) -> StagingUploadFiles | PublicUploadFiles:
    """
    Retrieve upload files after verifying user authorization.

    If access is granted, staging files are preferred when present;
    otherwise public files are returned.

    Args:
        upload_id: The unique identifier for the upload.
        user_id: The unique identifier for the user.

    Returns:
        The upload files object for the upload.

    Raises:
        PermissionError: If the upload exists but the user is not authorized.
        AttributeError: If the upload exists but neither staging nor public files
            are available.
    """
    if not _check_upload_access(upload_id, user_id):
        raise PermissionError(
            f'User {user_id} is not authorized to access upload {upload_id}.'
        )

    # User is authorized, retrieve and return files
    if StagingUploadFiles.exists_for(upload_id):
        return StagingUploadFiles(upload_id)
    if PublicUploadFiles.exists_for(upload_id):
        return PublicUploadFiles(upload_id)

    raise AttributeError(f'No public or staging upload files not found for {upload_id}')
