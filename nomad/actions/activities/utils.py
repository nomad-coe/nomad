from collections.abc import Callable

from nomad.actions import TaskQueue
from nomad.actions.action import get_actions
from nomad.actions.manager import request_signal_input_activity
from nomad.workflows.activities import (
    cleanup_entries_batch_activity,
    complete_upload_ownership_transfer_activity,
    delete_upload_files_activity,
    delete_upload_record_activity,
    delete_upload_search_activity,
    edit_upload_metadata_activity,
    finalize_cleanup_activity,
    finalize_upload_processing_activity,
    get_cleanup_entry_batch_from_file,
    handle_batch_heartbeat_failure_activity,
    handle_heartbeat_failure_activity,
    import_bundle_activity,
    match_all_activity,
    prepare_cleanup_activity,
    prepare_next_level_entry_batches,
    process_entry_activity,
    process_entry_batch_activity,
    process_entry_batch_from_file_activity,
    publish_externally_activity,
    publish_upload_activity,
    setup_example_upload_activity,
    update_files_activity,
)


def get_nomad_internal_activities() -> list[Callable]:
    return [
        prepare_cleanup_activity,
        get_cleanup_entry_batch_from_file,
        cleanup_entries_batch_activity,
        finalize_cleanup_activity,
        finalize_upload_processing_activity,
        delete_upload_search_activity,
        delete_upload_files_activity,
        delete_upload_record_activity,
        process_entry_activity,
        process_entry_batch_activity,
        match_all_activity,
        update_files_activity,
        setup_example_upload_activity,
        edit_upload_metadata_activity,
        import_bundle_activity,
        publish_upload_activity,
        publish_externally_activity,
        handle_heartbeat_failure_activity,
        prepare_next_level_entry_batches,
        process_entry_batch_from_file_activity,
        handle_batch_heartbeat_failure_activity,
        complete_upload_ownership_transfer_activity,
    ]


def get_all_activities(task_queue: TaskQueue) -> list[Callable]:
    activities = []
    if task_queue != TaskQueue.NOMAD_INTERNAL_WORKFLOWS:
        activities.append(request_signal_input_activity)
    for action_entry_point in get_actions().values():
        if action_entry_point.task_queue == task_queue:
            action = action_entry_point.load()
            activities.extend(action.activities)
    if task_queue == TaskQueue.NOMAD_INTERNAL_WORKFLOWS:
        activities.extend(get_nomad_internal_activities())
    return list(set(activities))
