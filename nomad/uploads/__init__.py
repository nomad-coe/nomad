from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from nomad.processing import Upload


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
