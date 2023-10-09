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

import re
import os
import hashlib
import json

from fastapi import FastAPI, Request, Response, status
from fastapi.exception_handlers import http_exception_handler as default_http_exception_handler
from starlette.exceptions import HTTPException as StarletteHTTPException
from fastapi.responses import HTMLResponse, JSONResponse
from starlette.staticfiles import StaticFiles as StarletteStaticFiles, NotModifiedResponse
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.responses import PlainTextResponse
from starlette.datastructures import Headers
from pydantic import BaseModel

from nomad import config, infrastructure
from .v1.main import app as v1_app
from .dcat.main import app as dcat_app
from .optimade import optimade_app
from .h5grove_app import app as h5grove_app


class OasisAuthenticationMiddleware(BaseHTTPMiddleware):
    async def dispatch(self, request, call_next):
        path = request.url.path
        if 'extensions' in path or 'info' in path or 'versions' in path:
            return await call_next(request)

        if 'Authorization' not in request.headers:
            return Response(
                status_code=status.HTTP_401_UNAUTHORIZED,
                content='You have to authenticate to use this Oasis endpoint.')
        else:
            token = request.headers['Authorization'].split(' ')[1]
            user, _ = infrastructure.keycloak.tokenauth(token)
            if user is None or user.email not in config.oasis.allowed_users:
                return Response(
                    status_code=status.HTTP_401_UNAUTHORIZED,
                    content='You are not authorized to access this Oasis endpoint.')

        return await call_next(request)


app = FastAPI()

if config.oasis.allowed_users is not None:
    optimade_app.add_middleware(OasisAuthenticationMiddleware)

app_base = config.services.api_base_path
app.mount(f'{app_base}/api/v1', v1_app)
app.mount(f'{app_base}/dcat', dcat_app)
app.mount(f'{app_base}/optimade', optimade_app)
app.mount(f'{app_base}/h5grove', h5grove_app)

if config.resources.enabled:
    from .resources.main import app as resources_app
    app.mount(f'{app_base}/resources', resources_app)

dist_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../dist'))
docs_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), 'static/docs'))
gui_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), 'static/gui'))
if not os.path.exists(gui_folder):
    os.makedirs(gui_folder)

configured_gui_folder = os.path.join(gui_folder, '../.gui_configured')
if os.path.exists(configured_gui_folder):
    gui_folder = configured_gui_folder


class StaticFiles(StarletteStaticFiles):

    etag_re = r'^(W/)?"?([^"]*)"?$'

    def is_not_modified(
        self, response_headers: Headers, request_headers: Headers
    ) -> bool:
        # The starlette etag implementation is not considering the "..." and W/"..." etag
        # RFC syntax used by browsers.
        try:
            if_none_match = request_headers["if-none-match"]
            match = re.match(StaticFiles.etag_re, if_none_match)
            if_none_match = match.group(2)
            etag = response_headers["etag"]
            if if_none_match == etag:
                return True
        except KeyError:
            pass

        return super().is_not_modified(response_headers, request_headers)


class GuiFiles(StaticFiles):

    gui_artifacts_data = None
    gui_env_data = None
    gui_data_etag = None

    async def get_response(self, path: str, scope) -> Response:
        if path not in ['env.js', 'artifacts.js']:
            response = await super().get_response(path, scope)
        else:
            assert GuiFiles.gui_data_etag is not None, 'Etag for gui data was not initialized'
            response = PlainTextResponse(
                GuiFiles.gui_env_data if path == 'env.js' else GuiFiles.gui_artifacts_data,
                media_type='application/javascript',
                headers=dict(etag=GuiFiles.gui_data_etag))

        request_headers = Headers(scope=scope)
        if self.is_not_modified(response.headers, request_headers):
            return NotModifiedResponse(response.headers)
        return response


app.mount(f'{app_base}/dist', StaticFiles(directory=dist_folder, check_dir=False), name='dist', )
app.mount(f'{app_base}/docs', StaticFiles(directory=docs_folder, check_dir=False), name='docs')
app.mount(f'{app_base}/gui', GuiFiles(directory=gui_folder, check_dir=False), name='gui')


@app.on_event('startup')
async def startup_event():
    from nomad.parsing.parsers import import_all_parsers
    import_all_parsers()

    from nomad import infrastructure
    # each subprocess is supposed disconnect and
    # connect again: https://jira.mongodb.org/browse/PYTHON-2090
    try:
        from mongoengine import disconnect
        disconnect()
    except Exception:
        pass

    from nomad.cli.dev import get_gui_artifacts_js
    GuiFiles.gui_artifacts_data = get_gui_artifacts_js()

    from nomad.cli.dev import get_gui_config
    GuiFiles.gui_env_data = get_gui_config()

    config_data = [
        item.json()
        for item in config.__dict__.values()
        if isinstance(item, BaseModel)]
    GuiFiles.gui_data_etag = hashlib.md5(
        json.dumps(config_data).encode(), usedforsecurity=False).hexdigest()

    infrastructure.setup()


@app.exception_handler(StarletteHTTPException)
async def http_exception_handler(request, exc):
    if exc.status_code != 404:
        return await default_http_exception_handler(request, exc)

    try:
        accept = request.headers['accept']
    except Exception:
        accept = None

    if accept is not None and 'html' in accept:
        return HTMLResponse(
            status_code=404, content=f'''
        <html>
            <head><title>{config.meta.name}</title></head>
            <body>
                <h1>NOMAD app</h1>
                <h2>info</h2>
                {'<br/>'.join(f'{key}: {value}' for key, value in config.meta.dict().items())}
                <h2>apis</h2>
                <a href="{app_base}/api/v1/extensions/docs">NOMAD API v1</a><br/>
                <a href="{app_base}/api/v1.2/extensions/docs">NOMAD API v1.2</a><br/>
                <a href="{app_base}/optimade/v1/extensions/docs">Optimade API</a><br/>
                <a href="{app_base}/dcat/extensions/docs">DCAT API</a><br/>
            </body>
        </html>
        ''')

    return JSONResponse(
        status_code=404, content={
            'detail': 'Not found',
            'info': {
                'app': config.meta.dict(),
                'apis': {
                    'v1': {
                        'root': f'{app_base}/api/v1',
                        'dashboard': f'{app_base}/api/v1/extensions/docs',
                        'documentation': f'{app_base}/api/v1/extensions/redoc',
                    },
                    'optimade': {
                        'root': f'{app_base}/optimade/v1',
                        'dashboard': f'{app_base}/optimade/v1/extensions/docs'
                    },
                    'dcat': {
                        'root': f'{app_base}/dcat',
                        'dashboard': f'{app_base}/dcat/extensions/docs'
                    }
                }
            }
        })


@app.get(f'{app_base}/alive')
async def alive():
    ''' Simple endpoint to utilize kubernetes liveness/readiness probing. '''
    return "I am, alive!"


@app.get('/-/health', status_code=status.HTTP_200_OK)
async def health():
    return {'healthcheck': 'ok'}


max_cache_ages = {
    r'\.[a-f0-9]+\.chunk\.(js|css)$': 3600 * 24 * 7,
    r'\.(html|js|css)$': config.services.html_resource_http_max_age,
    r'\.(png|jpg|gif|jpeg|ico)$': config.services.image_resource_http_max_age,
}


@app.middleware('http')
async def add_header(request: Request, call_next):
    response = await call_next(request)

    max_age = None
    for key, value in max_cache_ages.items():
        if re.search(key, str(request.url)):
            max_age = value
            break

    if max_age is not None:
        response.headers['Cache-Control'] = f'max-age={max_age}, must-revalidate'
    else:
        response.headers['Cache-Control'] = f'max-age=0, no-cache, no-store, must-revalidate'

    # The etags that we and starlette produce do not follow the RFC, because they do not
    # start with a " as the RFC specifies. Nginx considers them weak etags and will strip
    # these if gzip is enabled.
    if not response.headers.get('etag', '"').startswith('"'):
        response.headers['etag'] = f'"{response.headers.get("etag")}"'

    return response
