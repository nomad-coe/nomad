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

from fastapi import FastAPI, status, Request, Depends
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
import traceback
import re

from h5grove import fastapi_utils as h5grove_router

from nomad import config, utils
from nomad.app.v1.models import User
from nomad.app.v1.routers.auth import create_user_dependency
from nomad.app.v1.routers.uploads import get_upload_with_read_access

logger = utils.get_logger(__name__)

h5grove_router.settings.base_dir = config.fs.staging


async def check_user_access(upload_id: str, user: User = Depends(create_user_dependency(required=True))):
    get_upload_with_read_access(upload_id, user)

app = FastAPI(dependencies=[Depends(check_user_access)])

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.middleware("http")
async def add_upload_folder_path(request: Request, call_next):
    upload_path = f"{request.query_params['upload_id'][0:config.fs.prefix_size]}/{request.query_params['upload_id']}/raw/"
    scope = request.scope
    file = "file=" + upload_path + request.query_params["file"]
    extension = re.search(".(nxs|h5|hdf5|hd5|hdf)", scope["query_string"].decode("utf-8"))
    scope["query_string"] = file.encode("utf-8") + scope["query_string"][extension.span()[1]:]
    response = await call_next(Request(scope))
    return response


@app.exception_handler(Exception)
async def unicorn_exception_handler(request: Request, e: Exception):
    logger.error('unexpected exception in API', url=request.url, exc_info=e)
    return JSONResponse(
        status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        content={
            'detail': {
                'reason': 'Unexpected exception while handling your request',
                'exception': str(e),
                'exception_class': e.__class__.__name__,
                'exception_traceback': traceback.format_exc()
            }
        }
    )

app.include_router(h5grove_router.router)
