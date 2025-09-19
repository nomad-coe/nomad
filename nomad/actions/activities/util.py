from collections.abc import Callable

from nomad.actions import TaskQueue
from nomad.actions.action import get_actions
from nomad.workflows.activities import (
    cleanup_activity,
    cleanup_workflow_tmp_dir_activity,
    delete_upload_entries_activity,
    delete_upload_files_activity,
    delete_upload_record_activity,
    delete_upload_search_activity,
    edit_upload_metadata_activity,
    get_entry_batch_from_file,
    import_bundle_activity,
    match_all_activity,
    next_level_entries,
    process_entry_activity,
    process_entry_failure_activity,
    process_entry_success,
    process_upload_failure_activity,
    process_upload_success,
    publish_externally_activity,
    publish_upload_activity,
    remove_workflow_id_activity,
    setup_example_upload_activity,
    setup_upload_for_workflow_process,
    update_files_activity,
)


def get_nomad_internal_activities() -> list[Callable]:
    return [
        cleanup_workflow_tmp_dir_activity,
        delete_upload_search_activity,
        delete_upload_files_activity,
        delete_upload_entries_activity,
        delete_upload_record_activity,
        process_entry_activity,
        process_entry_success,
        process_entry_failure_activity,
        get_entry_batch_from_file,
        process_upload_failure_activity,
        process_upload_success,
        cleanup_activity,
        next_level_entries,
        match_all_activity,
        update_files_activity,
        setup_example_upload_activity,
        setup_upload_for_workflow_process,
        remove_workflow_id_activity,
        edit_upload_metadata_activity,
        import_bundle_activity,
        publish_upload_activity,
        publish_externally_activity,
    ]


def get_all_activities(task_queue: TaskQueue) -> list[Callable]:
    activities = []
    for action_entry_point in get_actions().values():
        if action_entry_point.task_queue == task_queue:
            action = action_entry_point.load()
            activities.extend(action.activities)
    if task_queue == TaskQueue.NOMAD_INTERNAL_WORKFLOWS:
        activities.extend(get_nomad_internal_activities())
    return activities
