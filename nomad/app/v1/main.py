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

from typing import Any
from fastapi import FastAPI, status, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, RedirectResponse
import traceback
import orjson

from nomad import config, utils

from .common import root_path
from .routers import (
    users, entries, materials, auth, info, datasets, uploads, suggestions, metainfo,
    north, systems, graph
)

logger = utils.get_logger(__name__)


class ORJSONResponse(JSONResponse):
    media_type = "application/json"

    def render(self, content: Any) -> bytes:
        return orjson.dumps(content, option=orjson.OPT_INDENT_2 | orjson.OPT_NON_STR_KEYS)  # type: ignore


app = FastAPI(
    root_path=root_path,
    openapi_url='/openapi.json',
    docs_url='/extensions/docs',
    redoc_url='/extensions/redoc',
    swagger_ui_oauth2_redirect_url='/extensions/docs/oauth2-redirect',
    title='NOMAD API',
    version=f'v1, NOMAD {config.meta.version}',
    description=utils.strip(f'''
        Please visit the [API section of the NOMAD documentation]({config.api_url(True, 'docs/api.html')})
        for a introduction and examples.
    '''),
    default_response_class=ORJSONResponse)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


async def redirect_to_docs(req: Request):
    return RedirectResponse(f'{root_path}/extensions/docs')


# app.add_route(f'{root_path}', redirect_to_docs, include_in_schema=False)
app.add_route('/', redirect_to_docs, include_in_schema=False)


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

app.include_router(info.router, prefix='/info')
app.include_router(auth.router, prefix='/auth')
app.include_router(materials.router, prefix='/materials')
app.include_router(entries.router, prefix='/entries')
app.include_router(datasets.router, prefix='/datasets')
app.include_router(uploads.router, prefix='/uploads')
app.include_router(metainfo.router, prefix='/metainfo')
app.include_router(users.router, prefix='/users')
app.include_router(suggestions.router, prefix='/suggestions')
if config.north.enabled:
    app.include_router(north.router, prefix='/north')
app.include_router(systems.router, prefix='/systems')
app.include_router(graph.router, prefix='/graph')
