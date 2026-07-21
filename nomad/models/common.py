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

import datetime
from typing import Annotated

from pydantic import AfterValidator

from nomad.utils import normalize_datetime_utc

UTCDateTime = Annotated[datetime.datetime, AfterValidator(normalize_datetime_utc)]


class ProcessStatus:
    """
    Class holding constants related to the possible process statuses.

    Attributes:
        READY: The process is ready to start
        PENDING: The process has been called, but still waiting for a celery worker to start running.
        RUNNING: Currently running the main process function.
        WAITING_FOR_RESULT: Waiting for the result from some other process.
        SUCCESS: The last process completed successfully.
        FAILURE: The last process completed with a fatal failure.
        DELETED: Used to signal that the process results in the deletion of the object.

        STATUSES_PROCESSING: List of statuses where the process is still incomplete (no other
            process can be started).
        STATUSES_NOT_PROCESSING: The opposite of the above - statuses from which a new
            process can be started.
    """

    READY = 'READY'
    PENDING = 'PENDING'
    RUNNING = 'RUNNING'
    WAITING_FOR_RESULT = 'WAITING_FOR_RESULT'
    SUCCESS = 'SUCCESS'
    FAILURE = 'FAILURE'
    DELETED = 'DELETED'

    STATUSES_PROCESSING = (PENDING, RUNNING, WAITING_FOR_RESULT)
    STATUSES_NOT_PROCESSING = (READY, SUCCESS, FAILURE)
    STATUSES_COMPLETED = (SUCCESS, FAILURE)
    STATUSES_VALID_IN_DB = STATUSES_NOT_PROCESSING + STATUSES_PROCESSING
