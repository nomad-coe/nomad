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
from __future__ import annotations

from fastapi import FastAPI, status, Request, Depends
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
import traceback
import re
import urllib.parse
import h5py
from typing import Callable, Dict, Any, IO

from h5grove import fastapi_utils as h5grove_router, utils as h5grove_utils

from nomad import utils
from nomad.files import UploadFiles, PublicUploadFiles
from nomad.app.v1.models import User
from nomad.app.v1.routers.auth import create_user_dependency
from nomad.app.v1.routers.uploads import get_upload_with_read_access

logger = utils.get_logger(__name__)


def open_zipped_h5_file(
    filepath: str,
    create_error: Callable[[int, str], Exception],
    h5py_options: Dict[str, Any] = {},
) -> h5py.File:
    import re
    import io
    from nomad import files

    """
    Patched h5grove utils function open_file_with_error_fallback in order to open h5 file
    in zipped folder.
    """
    match = re.match(
        r'.*?/uploads/(?P<upload_id>.+?)/(?P<directory>.+?)/(?P<path_or_id>.+)',
        filepath,
    )
    if not match:
        raise create_error(404, 'File not found!')

    upload_files = files.UploadFiles.get(match['upload_id'])
    path_or_id = match['path_or_id']
    try:
        file_object: IO | str
        if match['directory'] == 'raw':
            file_object = upload_files.raw_file(path_or_id, 'rb')
        else:
            file_object = upload_files.archive_hdf5_location(path_or_id)
    except Exception:
        raise create_error(404, 'File not found!')

    try:
        f = h5py.File(file_object, **h5py_options)
    except OSError as e:
        if isinstance(e, FileNotFoundError) or 'No such file or directory' in str(e):
            raise create_error(404, 'File not found!')
        if isinstance(e, PermissionError) or 'Permission denied' in str(e):
            raise create_error(403, 'Cannot read file: Permission denied!')
        if isinstance(e, io.UnsupportedOperation):
            raise create_error(404, 'File not found!')
        raise e

    return f


h5grove_utils.open_file_with_error_fallback.__code__ = open_zipped_h5_file.__code__


async def check_user_access(
    upload_id: str, user: User = Depends(create_user_dependency(required=True))
):
    get_upload_with_read_access(upload_id, user)


app = FastAPI(dependencies=[Depends(check_user_access)])

app.add_middleware(
    CORSMiddleware,
    allow_origins=['*'],
    allow_credentials=True,
    allow_methods=['*'],
    allow_headers=['*'],
)


@app.middleware('http')
async def add_upload_folder_path(request: Request, call_next):
    upload_id = request.query_params['upload_id']
    file = request.query_params['file']
    path = request.query_params['path']
    source = request.query_params['source']

    upload_path = f'/uploads/{upload_id}/{source}/'
    if source == 'archive' and isinstance(
        UploadFiles.get(upload_id), PublicUploadFiles
    ):
        path = f'{file}{path}'

    scope = request.scope
    old_file = urllib.parse.quote(request.query_params['file'], safe=[])
    new_file = urllib.parse.quote(upload_path + request.query_params['file'], safe=[])
    scope['query_string'] = scope['query_string'].replace(
        old_file.encode('utf-8'), new_file.encode('utf-8')
    )
    query_string = scope['query_string'].decode('utf-8')
    query_string = re.sub(
        r'file=.+?(?:&|\Z)', f'file={upload_path}{file}&', query_string
    )
    query_string = re.sub(
        r'path=.+?(?:&|\Z)', f'path={urllib.parse.quote(path)}&', query_string
    )
    scope['query_string'] = query_string.encode('utf-8')

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
                'exception_traceback': traceback.format_exc(),
            }
        },
    )


app.include_router(h5grove_router.router)
