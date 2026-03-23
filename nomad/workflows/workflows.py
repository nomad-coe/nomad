"""
All workflow class definitions for NOMAD workflows.
"""

import asyncio
from datetime import timedelta

from temporalio import workflow
from temporalio.common import Priority, RetryPolicy
from temporalio.exceptions import ActivityError

EDIT_UPLOAD_METADATA_PRIORITY = Priority(priority_key=1)
PUBLISH_UPLOAD_PRIORITY = Priority(priority_key=2)
PUBLISH_EXTERNALLY_PRIORITY = Priority(priority_key=2)
DELETE_UPLOAD_PRIORITY = Priority(priority_key=3)
IMPORT_BUNDLE_PRIORITY = Priority(priority_key=4)
PROCESS_EXAMPLE_UPLOAD_PRIORITY = Priority(priority_key=4)
BATCH_PROCESS_ENTRIES_PRIORITY = Priority(priority_key=4)
PROCESS_UPLOAD_PRIORITY = Priority(priority_key=4)
UPDATE_UPLOAD_PRIORITY = Priority(priority_key=4)
BATCH_PROCESS_ENTRY_PRIORITY = Priority(priority_key=4)
PROCESS_ENTRY_PRIORITY = Priority(priority_key=5)

with workflow.unsafe.imports_passed_through():
    from nomad.config import config
    from nomad.workflows.activities import (
        cleanup_activity,
        cleanup_workflow_tmp_dir_activity,
        delete_upload_entries_activity,
        delete_upload_files_activity,
        delete_upload_record_activity,
        delete_upload_search_activity,
        edit_upload_metadata_activity,
        get_entry_batch_from_file,
        handle_heartbeat_failure_activity,
        import_bundle_activity,
        match_all_activity,
        next_level_entries,
        parser_min_level,
        process_entry_activity,
        process_entry_batch_activity,
        process_upload_failure_activity,
        process_upload_success,
        publish_externally_activity,
        publish_upload_activity,
        remove_workflow_id_activity,
        setup_example_upload_activity,
        setup_upload_for_workflow_process,
        update_files_activity,
    )
    from nomad.workflows.shared_objects import (
        DeleteUploadWorkflowInput,
        EditUploadMetadataWorkflowInput,
        EntriesToBeProcessedResult,
        EntryBatchFromFileInput,
        ImportBundleWorkflowInput,
        ProcessEntryActivityInput,
        ProcessExampleUploadWorkflowInput,
        PublishExternallyWorkflowInput,
        PublishUploadWorkflowInput,
        UploadProcessingWorkflowInput,
        UploadWorkflowIdInput,
    )
    from nomad.workflows.utils import generate_batches


@workflow.defn
class DeleteUploadWorkflow:
    @workflow.run
    async def run(self, input: DeleteUploadWorkflowInput):
        retry_policy = RetryPolicy(
            maximum_attempts=2,
        )
        timeout = timedelta(
            seconds=config.temporal.processing_timeouts.delete_upload_timeout
        )
        heartbeat_timeout = timedelta(
            seconds=config.temporal.processing_timeouts.internal_processing_heartbeat_timeout
        )
        await workflow.execute_activity(
            delete_upload_search_activity,
            input,
            schedule_to_close_timeout=timeout,
            heartbeat_timeout=heartbeat_timeout,
            retry_policy=retry_policy,
            priority=DELETE_UPLOAD_PRIORITY,
        )
        await workflow.execute_activity(
            delete_upload_files_activity,
            input,
            schedule_to_close_timeout=timeout,
            heartbeat_timeout=heartbeat_timeout,
            retry_policy=retry_policy,
            priority=DELETE_UPLOAD_PRIORITY,
        )
        await workflow.execute_activity(
            delete_upload_entries_activity,
            input,
            schedule_to_close_timeout=timeout,
            heartbeat_timeout=heartbeat_timeout,
            retry_policy=retry_policy,
            priority=DELETE_UPLOAD_PRIORITY,
        )
        await workflow.execute_activity(
            delete_upload_record_activity,
            input,
            schedule_to_close_timeout=timeout,
            heartbeat_timeout=heartbeat_timeout,
            retry_policy=retry_policy,
            priority=DELETE_UPLOAD_PRIORITY,
        )


@workflow.defn
class ProcessEntryWorkflow:
    @workflow.run
    async def run(self, input: ProcessEntryActivityInput):
        retry_policy = RetryPolicy(
            maximum_attempts=2,
        )
        try:
            result = await workflow.execute_activity(
                process_entry_activity,
                input,
                schedule_to_close_timeout=timedelta(
                    seconds=config.temporal.processing_timeouts.process_entry_timeout
                ),
                heartbeat_timeout=timedelta(
                    seconds=config.temporal.processing_timeouts.internal_processing_heartbeat_timeout
                ),
                retry_policy=retry_policy,
                priority=PROCESS_ENTRY_PRIORITY,
            )
        except ActivityError as e:
            if 'heartbeat timeout' in str(e.cause):
                await workflow.execute_activity(
                    handle_heartbeat_failure_activity,
                    input,
                    schedule_to_close_timeout=timedelta(
                        seconds=config.temporal.processing_timeouts.process_entry_timeout
                    ),
                    priority=PROCESS_ENTRY_PRIORITY,
                )
            raise e

        return result


@workflow.defn
class BatchProcessEntriesWorkflow:
    """
    Handles processing of entry batches.

    Architecture:
    - Uses continue-as-new to process batches sequentially, preventing history buildup
    - Within each batch, processes entries in configurable micro-batches
    - Handles both file-based storage (large datasets) and in-memory storage (small datasets)

    Note: Workflow-level splitting keeps in-memory batches bounded to 1000 entries.
    """

    @workflow.run
    async def run(self, next_level_entries_result: EntriesToBeProcessedResult):
        retry_policy = RetryPolicy(
            maximum_attempts=2,
        )
        max_entries_per_batch_workflow = 1000
        # Handle file-based entry storage (used for very large uploads).
        # Entries are persisted as per-batch files and loaded on demand.
        if entry_batch_directory := next_level_entries_result.directory:
            file_semaphore = asyncio.Semaphore(
                config.temporal.entry_workflow_batch_concurrency
            )  # Max concurrent batches

            async def process_file_batch(batch_id):
                async with file_semaphore:
                    entries_to_be_processed = await workflow.execute_activity(
                        get_entry_batch_from_file,
                        EntryBatchFromFileInput(
                            upload_id=next_level_entries_result.upload_id,
                            batch_dir_path=entry_batch_directory,
                            batch_id=batch_id,
                        ),
                        schedule_to_close_timeout=timedelta(
                            seconds=config.temporal.processing_timeouts.next_level_entries_timeout
                        ),
                        retry_policy=retry_policy,
                        priority=BATCH_PROCESS_ENTRIES_PRIORITY,
                    )
                    # Recursively process this batch (which may further subdivide if
                    # it exceeds the entry limit for a single workflow run).
                    await workflow.execute_child_workflow(
                        BatchProcessEntriesWorkflow.run,
                        EntriesToBeProcessedResult(
                            entries=entries_to_be_processed,
                            upload_id=next_level_entries_result.upload_id,
                        ),
                        id=f'{workflow.info().workflow_id}-file-batch-{batch_id}',
                        parent_close_policy=workflow.ParentClosePolicy.TERMINATE,
                        retry_policy=retry_policy,
                        priority=BATCH_PROCESS_ENTRIES_PRIORITY,
                    )

            # Each sub-batch is bounded by the per-workflow entry limit.
            await asyncio.gather(
                *[
                    process_file_batch(batch_id)
                    for batch_id in range(next_level_entries_result.total_batches)
                ]
            )
        # Handle in-memory entry processing (from small uploads or loaded file batches)
        elif entries_to_be_processed := next_level_entries_result.entries:
            # Two-tier processing strategy based on batch size:
            # 1. Large batches: Split into smaller batches and process sequentially
            # 2. Small batches: Process entries directly as micro-batched activities
            if len(entries_to_be_processed) > max_entries_per_batch_workflow:
                entry_batches = list(
                    generate_batches(
                        entries_to_be_processed,
                        max_desired_batch_size=max_entries_per_batch_workflow,
                    )
                )
                current_sub_batch_index = (
                    next_level_entries_result.current_sub_batch_index
                )

                # Process current sub-batch
                current_batch = entry_batches[current_sub_batch_index]
                await self._process_entries_batch(current_batch, retry_policy)

                # Continue to next sub-batch using continue-as-new if more remain
                next_sub_batch_index = current_sub_batch_index + 1
                if next_sub_batch_index < len(entry_batches):
                    workflow.continue_as_new(
                        EntriesToBeProcessedResult(
                            entries=entries_to_be_processed,  # Keep original list for batching
                            upload_id=next_level_entries_result.upload_id,
                            current_sub_batch_index=next_sub_batch_index,
                        )
                    )
            else:
                # Process entries directly as activities when within one workflow run.
                await self._process_entries_batch(entries_to_be_processed, retry_policy)

    async def _process_entries_batch(
        self, entries: list[ProcessEntryActivityInput], retry_policy: RetryPolicy
    ):
        """
        Process a batch of entries concurrently as activities.

        Args:
            entries: List of entry inputs to process.
            retry_policy: Retry policy for activity execution
        """
        entry_activity_batch_size = max(1, config.temporal.entry_activity_batch_size)
        micro_batches = generate_batches(
            entries, max_desired_batch_size=entry_activity_batch_size
        )
        batch_activity_concurrency = max(
            1, config.temporal.entry_concurrency_target // entry_activity_batch_size
        )
        tasks = []
        self.entry_batch_semaphore = asyncio.Semaphore(batch_activity_concurrency)

        for entry_batch in micro_batches:
            task = self._process_single_entry_batch(entry_batch, retry_policy)
            tasks.append(task)

        # Use return_exceptions=True to allow individual activities to fail
        # without stopping the entire batch or failing the parent workflow
        await asyncio.gather(*tasks, return_exceptions=True)

    async def _process_single_entry_batch(
        self, inputs: list[ProcessEntryActivityInput], retry_policy: RetryPolicy
    ):
        """
        Process a micro-batch of entries with error handling for heartbeat timeouts.

        Args:
            inputs: Entry inputs to process in a single activity
            retry_policy: Retry policy for activity execution

        Returns:
            Result from the process_entry_batch_activity
        """
        timeout_seconds = (
            config.temporal.processing_timeouts.process_entry_timeout * len(inputs)
        )
        try:
            async with self.entry_batch_semaphore:
                result = await workflow.execute_activity(
                    process_entry_batch_activity,
                    inputs,
                    schedule_to_close_timeout=timedelta(seconds=timeout_seconds),
                    heartbeat_timeout=timedelta(
                        seconds=config.temporal.processing_timeouts.internal_processing_heartbeat_timeout
                    ),
                    retry_policy=retry_policy,
                    priority=BATCH_PROCESS_ENTRY_PRIORITY,
                )
                return result

        except ActivityError as e:
            # Handle heartbeat timeout failures with a dedicated recovery activity
            if 'heartbeat timeout' in str(e.cause):
                for input in inputs:
                    await workflow.execute_activity(
                        handle_heartbeat_failure_activity,
                        input,
                        schedule_to_close_timeout=timedelta(
                            seconds=config.temporal.processing_timeouts.process_entry_timeout
                        ),
                        priority=BATCH_PROCESS_ENTRIES_PRIORITY,
                    )
            raise e


@workflow.defn
class ProcessUploadWorkflow:
    """
    Specialized workflow to process an upload through multiple steps:
    1. Match all files to parsers
    2. Parse entries level by level
    3. Cleanup temporary data
    """

    @workflow.run
    async def run(self, input: UploadProcessingWorkflowInput):
        retry_policy = RetryPolicy(
            maximum_attempts=2,
        )
        heartbeat_timeout = timedelta(
            seconds=config.temporal.processing_timeouts.internal_processing_heartbeat_timeout
        )
        workflow_info = workflow.info()
        # Step 2: Match all, pass updated_files as set
        await workflow.execute_activity(
            match_all_activity,
            input,
            schedule_to_close_timeout=timedelta(
                seconds=config.temporal.processing_timeouts.match_all_timeout
            ),
            heartbeat_timeout=heartbeat_timeout,
            retry_policy=retry_policy,
            priority=PROCESS_UPLOAD_PRIORITY,
        )

        # Step 3: Parse next level
        while True:  # Outer loop: Continue until no more parser levels to process
            next_level_entries_result = await workflow.execute_activity(
                next_level_entries,
                input,
                schedule_to_close_timeout=timedelta(
                    seconds=config.temporal.processing_timeouts.next_level_entries_timeout
                ),
                heartbeat_timeout=heartbeat_timeout,
                retry_policy=retry_policy,
                priority=PROCESS_UPLOAD_PRIORITY,
            )

            # If None returned: no entries exist for this parser level at all
            # then we're done with all parser levels.
            if not next_level_entries_result:
                break

            # Delegate all batch processing complexity to BatchProcessEntriesWorkflow
            await workflow.execute_child_workflow(
                BatchProcessEntriesWorkflow.run,
                next_level_entries_result,
                id=f'{workflow_info.workflow_id}-{input.min_level}-batch-processor',
                parent_close_policy=workflow.ParentClosePolicy.TERMINATE,
                retry_policy=retry_policy,
                priority=PROCESS_UPLOAD_PRIORITY,
            )

            next_parser_level = (
                next_level_entries_result.next_parser_level or input.min_level
            )
            input.min_level = next_parser_level + 1

        # Step 4: Cleanup
        await workflow.execute_activity(
            cleanup_activity,
            input,
            schedule_to_close_timeout=timedelta(
                seconds=config.temporal.processing_timeouts.cleanup_timeout
            ),
            heartbeat_timeout=heartbeat_timeout,
            retry_policy=retry_policy,
            priority=PROCESS_UPLOAD_PRIORITY,
        )


@workflow.defn
class UpdateUploadWorkflow:
    """
    Workflow to update an upload's files and optionally reprocess them.
    1. Update files
    2. (Optional) Reprocess updated files through ProcessUploadWorkflow
    3. Mark upload as successful or failed
    By default, reprocessing is triggered unless specified otherwise.
    """

    @workflow.run
    async def run(self, input: UploadProcessingWorkflowInput):
        retry_policy = RetryPolicy(
            maximum_attempts=2,
        )
        timeout = timedelta(
            seconds=config.temporal.processing_timeouts.process_upload_timeout
        )
        heartbeat_timeout = timedelta(
            seconds=config.temporal.processing_timeouts.internal_processing_heartbeat_timeout
        )
        workflow_info = workflow.info()
        upload_workflow_input = UploadWorkflowIdInput(
            upload_id=input.upload_id,
            workflow_id=workflow_info.workflow_id,
            process_name='_process_upload',
            trigger_processing=input.trigger_processing,
        )
        try:
            # Step 0: Add workflow id to upload
            await workflow.execute_activity(
                setup_upload_for_workflow_process,
                upload_workflow_input,
                schedule_to_close_timeout=timedelta(
                    seconds=config.temporal.processing_timeouts.setup_upload_timeout
                ),
                retry_policy=retry_policy,
                priority=UPDATE_UPLOAD_PRIORITY,
            )

            # Step 1: Update files
            updated_files = await workflow.execute_activity(
                update_files_activity,
                input,
                schedule_to_close_timeout=timedelta(
                    seconds=config.temporal.processing_timeouts.update_files_timeout
                ),
                heartbeat_timeout=heartbeat_timeout,
                retry_policy=retry_policy,
                priority=UPDATE_UPLOAD_PRIORITY,
            )

            if input.trigger_processing:
                parse_all_input = UploadProcessingWorkflowInput(
                    upload_id=input.upload_id,
                    file_operations=input.file_operations,
                    reprocess_settings=input.reprocess_settings,
                    path_filter=input.path_filter,
                    only_updated_files=input.only_updated_files,
                    publish_directly_after_processing=input.publish_directly_after_processing,
                    updated_files=updated_files,
                    min_level=parser_min_level,
                    workflow_id=input.workflow_id,
                    workflow_tmp_dir=input.workflow_tmp_dir,
                )
                # Here we excecute steps:
                # 2: Match all, pass updated_files
                # 3: Parse next level(s)
                # 4: Cleanup
                await workflow.execute_child_workflow(
                    ProcessUploadWorkflow.run,
                    parse_all_input,
                    id=f'{workflow_info.workflow_id}-reprocess-{input.upload_id}',
                    parent_close_policy=workflow.ParentClosePolicy.TERMINATE,
                    # Disable retries for the child workflow; it handles its own activity failures.
                    retry_policy=RetryPolicy(maximum_attempts=1),
                    priority=UPDATE_UPLOAD_PRIORITY,
                )

            # Step 5: Mark as successful if the processing was triggered, otherwise will mark as READY
            await workflow.execute_activity(
                process_upload_success,
                upload_workflow_input,
                schedule_to_close_timeout=timedelta(
                    seconds=config.temporal.processing_timeouts.process_upload_success_timeout
                ),
                retry_policy=retry_policy,
                priority=UPDATE_UPLOAD_PRIORITY,
            )

        except Exception as e:
            # Set upload to failure status
            upload_workflow_input.failure_message = 'Process upload failed'
            if isinstance(e, ActivityError):
                upload_workflow_input.error_details = str(e.cause)
            else:
                upload_workflow_input.error_details = str(e)

            await workflow.execute_activity(
                process_upload_failure_activity,
                upload_workflow_input,
                schedule_to_close_timeout=timeout,
                retry_policy=retry_policy,
                priority=UPDATE_UPLOAD_PRIORITY,
            )
            raise e

        finally:
            # Always remove workflow id, even if processing failed
            await workflow.execute_activity(
                remove_workflow_id_activity,
                upload_workflow_input,
                schedule_to_close_timeout=timedelta(
                    seconds=config.temporal.processing_timeouts.remove_workflow_id_timeout
                ),
                retry_policy=retry_policy,
                priority=UPDATE_UPLOAD_PRIORITY,
            )
            await workflow.execute_activity(
                cleanup_workflow_tmp_dir_activity,
                input.workflow_tmp_dir,
                schedule_to_close_timeout=timedelta(
                    seconds=config.temporal.processing_timeouts.cleanup_workflow_tmp_dir_timeout
                ),
                retry_policy=retry_policy,
                priority=UPDATE_UPLOAD_PRIORITY,
            )


@workflow.defn
class ProcessExampleUploadWorkflow:
    @workflow.run
    async def run(self, input: ProcessExampleUploadWorkflowInput):
        # Step 1: Setup example upload
        timeout = timedelta(
            seconds=config.temporal.processing_timeouts.process_example_upload_timeout
        )
        heartbeat_timeout = timedelta(
            seconds=config.temporal.processing_timeouts.internal_processing_heartbeat_timeout
        )
        await workflow.execute_activity(
            setup_example_upload_activity,
            input,
            schedule_to_close_timeout=timeout,
            heartbeat_timeout=heartbeat_timeout,
            priority=PROCESS_EXAMPLE_UPLOAD_PRIORITY,
        )
        current_workflow_id = workflow.info().workflow_id

        # Step 2: Process upload using the standard workflow
        process_upload_input = UploadProcessingWorkflowInput(
            upload_id=input.upload_id,
            file_operations=input.file_operations,
            publish_directly_after_processing=input.publish_directly,
            workflow_id=current_workflow_id,
            workflow_tmp_dir=input.workflow_tmp_dir,
        )

        await workflow.execute_child_workflow(
            UpdateUploadWorkflow.run,
            process_upload_input,
            id=f'process-upload-workflow-{current_workflow_id}-{input.upload_id}',
            parent_close_policy=workflow.ParentClosePolicy.TERMINATE,
            priority=PROCESS_EXAMPLE_UPLOAD_PRIORITY,
        )


@workflow.defn
class EditUploadMetadataWorkflow:
    @workflow.run
    async def run(self, input: EditUploadMetadataWorkflowInput):
        retry_policy = RetryPolicy(
            maximum_attempts=2,
        )
        timeout = timedelta(
            seconds=config.temporal.processing_timeouts.edit_upload_metadata_timeout
        )
        heartbeat_timeout = timedelta(
            seconds=config.temporal.processing_timeouts.internal_processing_heartbeat_timeout
        )
        workflow_info = workflow.info()
        upload_workflow_input = UploadWorkflowIdInput(
            upload_id=input.upload_id,
            workflow_id=workflow_info.workflow_id,
            process_name='_edit_upload_metadata',
        )

        try:
            # Add workflow id to upload
            await workflow.execute_activity(
                setup_upload_for_workflow_process,
                upload_workflow_input,
                schedule_to_close_timeout=timeout,
                retry_policy=retry_policy,
                priority=EDIT_UPLOAD_METADATA_PRIORITY,
            )

            # Edit upload metadata
            await workflow.execute_activity(
                edit_upload_metadata_activity,
                input,
                schedule_to_close_timeout=timeout,
                heartbeat_timeout=heartbeat_timeout,
                retry_policy=retry_policy,
                priority=EDIT_UPLOAD_METADATA_PRIORITY,
            )

            # Mark as successful
            await workflow.execute_activity(
                process_upload_success,
                upload_workflow_input,
                schedule_to_close_timeout=timeout,
                retry_policy=retry_policy,
                priority=EDIT_UPLOAD_METADATA_PRIORITY,
            )
        except Exception as e:
            # Set upload to failure status
            upload_workflow_input.failure_message = 'Edit metadata failed'
            if isinstance(e, ActivityError):
                upload_workflow_input.error_details = str(e.cause)
            else:
                upload_workflow_input.error_details = str(e)

            await workflow.execute_activity(
                process_upload_failure_activity,
                upload_workflow_input,
                schedule_to_close_timeout=timeout,
                retry_policy=retry_policy,
                priority=EDIT_UPLOAD_METADATA_PRIORITY,
            )
            raise e

        finally:
            # Always remove workflow id, even if processing failed
            await workflow.execute_activity(
                remove_workflow_id_activity,
                upload_workflow_input,
                schedule_to_close_timeout=timeout,
                retry_policy=retry_policy,
                priority=EDIT_UPLOAD_METADATA_PRIORITY,
            )


@workflow.defn
class ImportBundleWorkflow:
    @workflow.run
    async def run(self, input: ImportBundleWorkflowInput):
        retry_policy = RetryPolicy(
            maximum_attempts=2,
        )
        timeout = timedelta(
            seconds=config.temporal.processing_timeouts.import_bundle_timeout
        )
        heartbeat_timeout = timedelta(
            seconds=config.temporal.processing_timeouts.internal_processing_heartbeat_timeout
        )
        workflow_info = workflow.info()
        upload_workflow_input = UploadWorkflowIdInput(
            upload_id=input.upload_id,
            workflow_id=workflow_info.workflow_id,
            process_name='_import_bundle',
        )

        try:
            # Add workflow id to upload
            await workflow.execute_activity(
                setup_upload_for_workflow_process,
                upload_workflow_input,
                schedule_to_close_timeout=timeout,
                retry_policy=retry_policy,
                priority=IMPORT_BUNDLE_PRIORITY,
            )

            # Import bundle
            await workflow.execute_activity(
                import_bundle_activity,
                input,
                schedule_to_close_timeout=timeout,
                heartbeat_timeout=heartbeat_timeout,
                retry_policy=retry_policy,
                priority=IMPORT_BUNDLE_PRIORITY,
            )

            # Mark as successful
            await workflow.execute_activity(
                process_upload_success,
                upload_workflow_input,
                schedule_to_close_timeout=timeout,
                retry_policy=retry_policy,
                priority=IMPORT_BUNDLE_PRIORITY,
            )
        except Exception as e:
            # Set upload to failure status
            upload_workflow_input.failure_message = 'Import bundle failed'
            if isinstance(e, ActivityError):
                upload_workflow_input.error_details = str(e.cause)
            else:
                upload_workflow_input.error_details = str(e)

            await workflow.execute_activity(
                process_upload_failure_activity,
                upload_workflow_input,
                schedule_to_close_timeout=timeout,
                retry_policy=retry_policy,
                priority=IMPORT_BUNDLE_PRIORITY,
            )
            raise e

        finally:
            # Always remove workflow id, even if processing failed
            await workflow.execute_activity(
                remove_workflow_id_activity,
                upload_workflow_input,
                schedule_to_close_timeout=timeout,
                retry_policy=retry_policy,
                priority=IMPORT_BUNDLE_PRIORITY,
            )


@workflow.defn
class PublishUploadWorkflow:
    @workflow.run
    async def run(self, input: PublishUploadWorkflowInput):
        retry_policy = RetryPolicy(
            maximum_attempts=2,
        )
        timeout = timedelta(
            seconds=config.temporal.processing_timeouts.publish_upload_timeout
        )
        heartbeat_timeout = timedelta(
            seconds=config.temporal.processing_timeouts.internal_processing_heartbeat_timeout
        )
        workflow_info = workflow.info()
        upload_workflow_input = UploadWorkflowIdInput(
            upload_id=input.upload_id,
            workflow_id=workflow_info.workflow_id,
            process_name='_publish_upload',
        )

        try:
            # Add workflow id to upload
            await workflow.execute_activity(
                setup_upload_for_workflow_process,
                upload_workflow_input,
                schedule_to_close_timeout=timeout,
                retry_policy=retry_policy,
                priority=PUBLISH_UPLOAD_PRIORITY,
            )

            # Publish upload
            await workflow.execute_activity(
                publish_upload_activity,
                input,
                schedule_to_close_timeout=timeout,
                heartbeat_timeout=heartbeat_timeout,
                retry_policy=retry_policy,
                priority=PUBLISH_UPLOAD_PRIORITY,
            )

            # Mark as successful
            await workflow.execute_activity(
                process_upload_success,
                upload_workflow_input,
                schedule_to_close_timeout=timeout,
                retry_policy=retry_policy,
                priority=PUBLISH_UPLOAD_PRIORITY,
            )

        except Exception as e:
            # Set upload to failure status
            upload_workflow_input.failure_message = 'Publish upload failed'
            if isinstance(e, ActivityError):
                upload_workflow_input.error_details = str(e.cause)
            else:
                upload_workflow_input.error_details = str(e)

            await workflow.execute_activity(
                process_upload_failure_activity,
                upload_workflow_input,
                schedule_to_close_timeout=timeout,
                retry_policy=retry_policy,
                priority=PUBLISH_UPLOAD_PRIORITY,
            )
            raise e

        finally:
            # Always remove workflow id, even if processing failed
            await workflow.execute_activity(
                remove_workflow_id_activity,
                upload_workflow_input,
                schedule_to_close_timeout=timeout,
                retry_policy=retry_policy,
                priority=PUBLISH_UPLOAD_PRIORITY,
            )


@workflow.defn
class PublishExternallyWorkflow:
    @workflow.run
    async def run(self, input: PublishExternallyWorkflowInput):
        retry_policy = RetryPolicy(
            maximum_attempts=2,
        )
        timeout = timedelta(
            seconds=config.temporal.processing_timeouts.publish_externally_timeout
        )
        heartbeat_timeout = timedelta(
            seconds=config.temporal.processing_timeouts.internal_processing_heartbeat_timeout
        )
        workflow_info = workflow.info()
        upload_workflow_input = UploadWorkflowIdInput(
            upload_id=input.upload_id,
            workflow_id=workflow_info.workflow_id,
            process_name='_publish_externally',
        )

        try:
            # Add workflow id to upload
            await workflow.execute_activity(
                setup_upload_for_workflow_process,
                upload_workflow_input,
                schedule_to_close_timeout=timeout,
                retry_policy=retry_policy,
                priority=PUBLISH_EXTERNALLY_PRIORITY,
            )

            # Publish externally
            await workflow.execute_activity(
                publish_externally_activity,
                input,
                schedule_to_close_timeout=timeout,
                heartbeat_timeout=heartbeat_timeout,
                retry_policy=retry_policy,
                priority=PUBLISH_EXTERNALLY_PRIORITY,
            )

            # Mark as successful
            await workflow.execute_activity(
                process_upload_success,
                upload_workflow_input,
                schedule_to_close_timeout=timeout,
                retry_policy=retry_policy,
                priority=PUBLISH_EXTERNALLY_PRIORITY,
            )

        except Exception as e:
            # Set upload to failure status
            upload_workflow_input.failure_message = 'Publish externally failed'
            if isinstance(e, ActivityError):
                upload_workflow_input.error_details = str(e.cause)
            else:
                upload_workflow_input.error_details = str(e)

            await workflow.execute_activity(
                process_upload_failure_activity,
                upload_workflow_input,
                schedule_to_close_timeout=timeout,
                retry_policy=retry_policy,
                priority=PUBLISH_EXTERNALLY_PRIORITY,
            )
            raise e

        finally:
            # Always remove workflow id, even if processing failed
            await workflow.execute_activity(
                remove_workflow_id_activity,
                upload_workflow_input,
                schedule_to_close_timeout=timeout,
                retry_policy=retry_policy,
                priority=PUBLISH_EXTERNALLY_PRIORITY,
            )
