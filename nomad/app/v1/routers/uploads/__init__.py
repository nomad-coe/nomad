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

from fastapi import APIRouter

from . import actions, bundles, default
from .utils import (
    _check_upload_not_processing,
    _get_upload_with_write_access,
    entry_to_pydantic,
    get_upload_with_read_access,
    upload_to_pydantic,
)
from .models import (
    EntryProcData,
    EntryProcDataPagination,
    PaginationResponse,
    RawDirPagination,
    UploadProcData,
    UploadProcDataPagination,
    UploadProcDataQuery,
    UploadProcDataResponse,
)

router = APIRouter()
router.include_router(actions.router)
router.include_router(bundles.router)

for route in default.router.routes:
    router.routes.append(route)
