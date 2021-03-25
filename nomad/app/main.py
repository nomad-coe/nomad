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

from fastapi import FastAPI, status, Response
from fastapi.responses import JSONResponse, HTMLResponse
# We use a2wsgi. It is an alternative to the fastapi provided WSGIMiddleware that manages
# to stream requests instead of buffering them.
from a2wsgi import WSGIMiddleware
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.exceptions import HTTPException as StarletteHTTPException
from fastapi.exception_handlers import http_exception_handler as default_http_exception_handler

from nomad import config, infrastructure

from .optimade import optimade_app
from .flask import app as flask_app
from .v1.main import app as v1_app


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
            user, _ = infrastructure.keycloak.auth(request.headers)
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
app.mount(f'{app_base}/optimade', optimade_app)
app.mount(app_base, WSGIMiddleware(flask_app))


@app.on_event('startup')
async def startup_event():
    from nomad import infrastructure
    # each subprocess is supposed disconnect connect again: https://jira.mongodb.org/browse/PYTHON-2090
    try:
        from mongoengine import disconnect
        disconnect()
    except Exception:
        pass

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
        return HTMLResponse(status_code=404, content=f'''
        <html>
            <head><title>{config.meta.name}</title></head>
            <body>
                <h1>NOMAD app</h1>
                <h2>info</h2>
                {'<br/>'.join(f'{key}: {value}' for key, value in config.meta.items())}
                <h2>apis</h2>
                <a href="{app_base}/api">NOMAD API v0</a><br/>
                <a href="{app_base}/api/v1/extensions/docs">NOMAD API v1</a><br/>
                <a href="{app_base}/optimade/v1/extensions/docs">Optimade API</a><br/>
                <a href="{app_base}/dcat">DCAT API</a><br/>
            </body>
        </html>
        ''')

    return JSONResponse(status_code=404, content={
        'detail': 'Not found',
        'info': {
            'app': config.meta,
            'apis': {
                'v0': {
                    'root': f'{app_base}/api',
                    'dashboard': f'{app_base}/api',
                },
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
                    'dashboard': f'{app_base}/dcat'
                }
            }
        }
    })
