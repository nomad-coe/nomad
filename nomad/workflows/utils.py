import math
from typing import TypeVar

T = TypeVar('T')
MAX_CONCURRENT_ACTIVITIES_PER_WORKFLOW = 1000


def get_max_entries_per_batch_workflow(entry_activity_batch_size: int) -> int:
    """
    Maximum entries one batch-processing workflow should hold before splitting.

    Temporal limits concurrent activities per workflow to 1000. Since each activity
    can process multiple entries, the effective entry ceiling scales with
    `entry_activity_batch_size`.
    """
    return MAX_CONCURRENT_ACTIVITIES_PER_WORKFLOW * max(1, entry_activity_batch_size)


def generate_batches(
    items: list[T], max_desired_batch_size=1000, max_batches=1000
) -> list[list[T]]:
    """
    Splits a list of items into batches, trying to keep batch size <= max_desired_batch_size,
    while minimizing the number of batches. The number of batches is capped at max_batches.
    """
    total_items = len(items)

    if total_items == 0:
        return []

    # Find the smallest num_batches such that batch_size <= max_desired_batch_size
    for num_batches in range(1, max_batches + 1):
        batch_size = math.ceil(total_items / num_batches)
        if batch_size <= max_desired_batch_size:
            break
    else:
        # If we never found a small enough batch_size, use max_batches
        num_batches = max_batches
        batch_size = math.ceil(total_items / num_batches)

    item_batches = []
    for i in range(num_batches):
        start_idx = i * batch_size
        end_idx = min(start_idx + batch_size, total_items)
        item_batches.append(items[start_idx:end_idx])

    return item_batches
