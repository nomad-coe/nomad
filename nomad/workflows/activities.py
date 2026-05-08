import json
import os
import random
import shutil
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

from temporalio import activity
from temporalio.exceptions import ApplicationError

from nomad.actions.heartbeat import activity_heartbeat
from nomad.config import config
from nomad.files import PublicUploadFiles, StagingUploadFiles
from nomad.parsing.parsers import parsers
from nomad.processing.base import ProcessFailure, ProcessStatus
from nomad.processing.data import Entry, Upload
from nomad.search import delete_upload
from nomad.workflows.shared_objects import (
    CleanupEntriesBatchActivityInput,
    CleanupEntriesResult,
    CleanupEntryBatchFromFileInput,
    DeleteUploadWorkflowInput,
    EditUploadMetadataWorkflowInput,
    EntriesToBeProcessedResult,
    EntryBatchFromFileInput,
    FinalizeUploadProcessingInput,
    ImportBundleWorkflowInput,
    ProcessEntryActivityInput,
    ProcessExampleUploadWorkflowInput,
    PublishExternallyWorkflowInput,
    PublishUploadWorkflowInput,
    UpdatedFilesResult,
    UploadProcessingWorkflowInput,
    UploadWorkflowIdInput,
)
from nomad.workflows.utils import CLEANUP_ENTRY_BATCH_SIZE, generate_batches

parser_min_level = min([parser.level for parser in parsers])
# If the heartbeat timeout is 10 mins, this would send a heartbeat every 60 seconds.
HEARTBEAT_FREQUENCY = (
    config.temporal.processing_timeouts.internal_processing_heartbeat_timeout / 10
)
MAX_IN_MEMORY_ENTRIES = 1000
CLEANUP_FAST_PATH_ENTRY_THRESHOLD = 100


@activity.defn
def delete_upload_search_activity(input: DeleteUploadWorkflowInput):
    with activity_heartbeat(HEARTBEAT_FREQUENCY):
        # Delete from search index
        delete_upload(input.upload_id, refresh=True)


@activity.defn
def delete_upload_files_activity(input: DeleteUploadWorkflowInput):
    with activity_heartbeat(HEARTBEAT_FREQUENCY):
        # Delete staging and public files
        for cls in (StagingUploadFiles, PublicUploadFiles):
            if cls.exists_for(input.upload_id):
                cls(input.upload_id).delete()


@activity.defn
def delete_upload_entries_activity(input: DeleteUploadWorkflowInput):
    with activity_heartbeat(HEARTBEAT_FREQUENCY):
        # Delete all entries for this upload
        Entry.objects(upload_id=input.upload_id).delete()  # type: ignore


@activity.defn
def delete_upload_record_activity(input: DeleteUploadWorkflowInput):
    with activity_heartbeat(HEARTBEAT_FREQUENCY):
        # Delete the upload itself
        upload = Upload.get(input.upload_id)
        upload.delete()


def _process_single_entry(input: ProcessEntryActivityInput):
    """Process one entry and map permanent processing failures to non-retryable errors."""
    entry = Entry.get(input.entry_id)
    try:
        entry.errors = []
        entry._process_entry_local()
        entry.on_success()
        entry.process_status = ProcessStatus.SUCCESS
        entry.complete_time = datetime.now(timezone.utc)
        entry.save()
    except Exception as e:
        entry.fail(*[e])
        entry.save()
        if isinstance(e, ProcessFailure):
            # ProcessFailure represents permanent failures (data validation, business logic errors)
            # that cannot be resolved through retries.
            raise ApplicationError(str(e), non_retryable=True) from e
        raise e
    activity.heartbeat()


@activity.defn
def process_entry_activity(input: ProcessEntryActivityInput):
    with activity_heartbeat(HEARTBEAT_FREQUENCY):
        _process_single_entry(input)


@activity.defn
def process_entry_batch_activity(inputs: list[ProcessEntryActivityInput]):
    """
    Process a batch of entries in one Temporal activity invocation.

    Non-retryable entry failures (`ProcessFailure`) are isolated to the affected
    entry and do not abort the batch. Retryable failures are raised only after
    all entries in the batch have been attempted.
    """
    with activity_heartbeat(HEARTBEAT_FREQUENCY):
        retryable_error: Exception | None = None
        for entry_input in inputs:
            try:
                _process_single_entry(entry_input)
            except ApplicationError as e:
                if getattr(e, 'non_retryable', False):
                    continue
                if retryable_error is None:
                    retryable_error = e
            except Exception as e:
                if retryable_error is None:
                    retryable_error = e

        if retryable_error is not None:
            raise retryable_error


@activity.defn
def update_files_activity(
    input: UploadProcessingWorkflowInput,
) -> UpdatedFilesResult:
    with activity_heartbeat(HEARTBEAT_FREQUENCY):
        upload = Upload.get(input.upload_id)
        file_operations = input.file_operations or []
        only_updated_files = (
            input.only_updated_files if input.only_updated_files is not None else False
        )
        updated_files = upload.update_files(file_operations, only_updated_files)
        if not updated_files:
            return UpdatedFilesResult()

        # Temporal has a 1.5MB limit on serialized activity results. For large file sets,
        # we store the data on disk and pass the file path instead of the full dataset.
        if len(updated_files) < 1000:
            return UpdatedFilesResult(files=updated_files)

        # If 1000+ files, save to JSON file and return the path
        updated_files_path = os.path.join(input.workflow_tmp_dir, 'updated_files.json')

        with open(updated_files_path, 'w') as f:
            json.dump(list(updated_files), f)

        return UpdatedFilesResult(file_path=str(updated_files_path))


@activity.defn
def match_all_activity(input: UploadProcessingWorkflowInput):
    with activity_heartbeat(HEARTBEAT_FREQUENCY):
        reprocess_settings = input.reprocess_settings or {}
        reprocess_obj = config.reprocess.customize(reprocess_settings)
        upload = Upload.get(input.upload_id)
        upload.match_all(
            reprocess_settings=reprocess_obj,
            path_filter=input.path_filter,
            updated_files=input.updated_files.get_files(),
        )


@activity.defn
def next_level_entries(
    input: UploadProcessingWorkflowInput,
) -> EntriesToBeProcessedResult | None:
    with activity_heartbeat(HEARTBEAT_FREQUENCY):
        upload = Upload.get(input.upload_id)
        next_entries = upload.next_level_entries(
            min_level=input.min_level,
            path_filter=input.path_filter,
            updated_files=input.updated_files.get_files(),
        )

        # If no entries exist at this parser level, return None
        # This signals the workflow that we're completely done (no more parser levels)
        if not next_entries:
            return None

        # Split all entries into manageable batches
        # Temporal imposes a limit of 1.5MB for the serialized result, this helps us stay within those limits.
        entry_batches = generate_batches(next_entries)

        # When dealing with multiple large batches, storing all entries in workflow state
        # would exceed Temporal's serialization limits. Instead, we persist batches to disk
        # and process them sequentially via file references.
        if len(next_entries) < MAX_IN_MEMORY_ENTRIES:
            entries_list = [
                ProcessEntryActivityInput(
                    upload_id=input.upload_id,
                    entry_id=str(entry.entry_id),
                    workflow_id=f'process-entry-workflow-child-id-{str(entry.entry_id)}-{str(upload.upload_id)}-{uuid.uuid4()}',
                )
                for entry in next_entries
            ]
            return EntriesToBeProcessedResult(
                entries=entries_list,
                next_parser_level=upload.parser_level,
                upload_id=input.upload_id,
            )

        # Persist entry IDs as one file per batch to avoid large workflow payloads
        # without introducing repeated scans through a shared file.
        batch_dir = os.path.join(input.workflow_tmp_dir, f'level_{input.min_level}')
        os.makedirs(batch_dir, exist_ok=True)
        entry_batches = generate_batches(
            next_entries,
            max_desired_batch_size=MAX_IN_MEMORY_ENTRIES,
        )

        for batch_idx, batch in enumerate(entry_batches):
            batch_file = os.path.join(batch_dir, f'entry_batch_{batch_idx}.json')
            with open(batch_file, 'w') as f:
                json.dump([str(entry.entry_id) for entry in batch], f)

        return EntriesToBeProcessedResult(
            directory=str(batch_dir),
            total_batches=len(entry_batches),
            next_parser_level=upload.parser_level,
            upload_id=input.upload_id,
        )


@activity.defn
def get_entry_batch_from_file(
    input: EntryBatchFromFileInput,
) -> list[ProcessEntryActivityInput]:
    """Load entries from a specific batch file."""
    batch_file = Path(input.batch_dir_path) / f'entry_batch_{input.batch_id}.json'
    if not batch_file.exists():
        return []

    with open(batch_file) as f:
        batch_entry_ids = json.load(f)

    return [
        ProcessEntryActivityInput(
            upload_id=input.upload_id,
            entry_id=entry_id,
            workflow_id=f'process-entry-workflow-child-id-{entry_id}-{input.upload_id}-{uuid.uuid4()}',
        )
        for entry_id in batch_entry_ids
    ]


@activity.defn
def setup_upload_for_workflow_process(input: UploadWorkflowIdInput):
    upload = Upload.get(input.upload_id)
    assert len(upload.workflow_ids) == 0, (  # type: ignore
        'Upload is currently being processed by another workflow'
    )
    upload.workflow_ids.append(input.workflow_id)  # type: ignore
    upload.errors = []
    upload.current_process = input.process_name
    upload.save()


@activity.defn
def prepare_cleanup_activity(
    input: UploadProcessingWorkflowInput,
) -> CleanupEntriesResult | None:
    """Prepare cleanup work and finish small uploads inline."""
    with activity_heartbeat(HEARTBEAT_FREQUENCY):
        upload = Upload.get(input.upload_id)
        entry_ids = [
            str(entry.entry_id)
            for entry in Entry.objects(upload_id=input.upload_id)  # type: ignore
        ]
        # Small uploads are cheaper to finish here than to route through batch orchestration.
        if len(entry_ids) < CLEANUP_FAST_PATH_ENTRY_THRESHOLD:
            upload.cleanup()
            return None

        # Larger uploads continue in the workflow via in-memory ids or batch files.
        upload.cleanup_prepare()
        if len(entry_ids) <= MAX_IN_MEMORY_ENTRIES:
            return CleanupEntriesResult(
                upload_id=input.upload_id,
                entry_ids=entry_ids,
            )

        batch_dir = os.path.join(input.workflow_tmp_dir, 'cleanup_batches')
        os.makedirs(batch_dir, exist_ok=True)
        entry_batches = generate_batches(
            entry_ids,
            max_desired_batch_size=CLEANUP_ENTRY_BATCH_SIZE,
        )

        for batch_idx, batch in enumerate(entry_batches):
            batch_file = os.path.join(batch_dir, f'cleanup_batch_{batch_idx}.json')
            with open(batch_file, 'w') as f:
                json.dump(batch, f)

        return CleanupEntriesResult(
            upload_id=input.upload_id,
            directory=str(batch_dir),
            total_batches=len(entry_batches),
        )


@activity.defn
def get_cleanup_entry_batch_from_file(
    input: CleanupEntryBatchFromFileInput,
) -> list[str]:
    batch_file = Path(input.batch_dir_path) / f'cleanup_batch_{input.batch_id}.json'
    if not batch_file.exists():
        return []

    with open(batch_file) as f:
        return json.load(f)


@activity.defn
def cleanup_entries_batch_activity(input: CleanupEntriesBatchActivityInput):
    with activity_heartbeat(HEARTBEAT_FREQUENCY):
        attempt = activity.info().attempt
        if attempt > 1:
            time.sleep(random.uniform(0, min(5.0, float(attempt))))
        upload = Upload.get(input.upload_id)
        upload.cleanup_entries_batch(input.entry_ids, refresh=input.refresh)


@activity.defn
def finalize_cleanup_activity(input: UploadProcessingWorkflowInput):
    with activity_heartbeat(HEARTBEAT_FREQUENCY):
        upload = Upload.get(input.upload_id)
        upload.cleanup_finalize()


@activity.defn
def setup_example_upload_activity(input: ProcessExampleUploadWorkflowInput):
    with activity_heartbeat(HEARTBEAT_FREQUENCY):
        upload = Upload.get(input.upload_id)
        upload.setup_example_upload(entry_point_id=input.example_upload_id)


@activity.defn
def edit_upload_metadata_activity(input: EditUploadMetadataWorkflowInput):
    with activity_heartbeat(HEARTBEAT_FREQUENCY):
        upload = Upload.get(input.upload_id)
        upload._edit_upload_metadata_local(input.edit_request_json, input.user_id)


@activity.defn
def import_bundle_activity(input: ImportBundleWorkflowInput):
    with activity_heartbeat(HEARTBEAT_FREQUENCY):
        upload = Upload.get(input.upload_id)
        upload._import_bundle_local(
            input.bundle_path, input.import_settings, input.embargo_length
        )


@activity.defn
def publish_upload_activity(input: PublishUploadWorkflowInput):
    with activity_heartbeat(HEARTBEAT_FREQUENCY):
        upload = Upload.get(input.upload_id)
        upload._publish_upload_local(input.embargo_length)


@activity.defn
def finalize_upload_processing_activity(input: FinalizeUploadProcessingInput):
    upload = Upload.get(input.upload_id)
    if input.result == 'success':
        # When processing is not triggered, set READY so processing can be started manually.
        upload.process_status = (
            ProcessStatus.SUCCESS if input.trigger_processing else ProcessStatus.READY
        )
        upload.set_last_status_message('Process completed successfully')
    else:
        upload.last_status_message = (
            input.failure_message if input.failure_message else 'Process upload failed'
        )
        errors = [input.error_details] if input.error_details else []
        upload.fail(*errors)

    if input.workflow_id in upload.workflow_ids:  # type: ignore
        upload.workflow_ids.remove(input.workflow_id)  # type: ignore
    upload.save()

    if input.workflow_tmp_dir and os.path.exists(input.workflow_tmp_dir):
        shutil.rmtree(input.workflow_tmp_dir, ignore_errors=True)


@activity.defn
def publish_externally_activity(input: PublishExternallyWorkflowInput):
    with activity_heartbeat(HEARTBEAT_FREQUENCY):
        upload = Upload.get(input.upload_id)
        upload._publish_externally_local(
            target_deployment_url=input.target_deployment_url,
            auth_token=input.auth_token,
            embargo_length=input.embargo_length,
        )


@activity.defn
def handle_heartbeat_failure_activity(input: ProcessEntryActivityInput):
    entry = Entry.get(input.entry_id)
    # A later retryable failure in the same batch can trigger heartbeat recovery
    # after this entry has already been persisted as successful.
    if entry.process_status == ProcessStatus.SUCCESS:
        return
    entry.fail(
        *[
            'Process entry failed due to a heartbeat timeout. '
            'If this keeps happening contact NOMAD/ your oasis admin for support.'
        ]
    )
    entry.save()
