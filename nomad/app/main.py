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

import hashlib
import json
import re
from contextlib import asynccontextmanager

from fastapi import FastAPI, HTTPException, Response, status
from fastapi.exception_handlers import (
    http_exception_handler as default_http_exception_handler,
)
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi_cache import FastAPICache
from starlette.exceptions import HTTPException as StarletteHTTPException
from starlette.middleware.base import BaseHTTPMiddleware
from temporalio.client import Client

from nomad.actions.client import get_client
from nomad.config import config
from nomad.config.models.plugins import APIEntryPoint
from nomad.mongo.cache import MongoBackend
from nomad.utils.structlogging import get_logger

from .static import GuiFiles
from .static import app as static_files_app
from .v1.main import app as v1_app
from .v1.routers.auth import resolve_user


class OasisAuthenticationMiddleware(BaseHTTPMiddleware):
    def __init__(
        self,
        app,
        whitelist: set[str] | None = None,
    ) -> None:
        """
        Middleware to enforce authentication on protected routes.

        Args:
            app: The ASGI application.
            whitelist (Iterable[str], optional): A list of regex strings
                for URL path patterns that are exempt from authentication.
        """
        super().__init__(app)
        self.whitelist_patterns = [re.compile(pat) for pat in (whitelist or [])]

    async def dispatch(self, request, call_next):
        # Skip if global auth is off or route is whitelisted
        router_path = request.url.path.removeprefix(
            f'{config.services.api_base_path}/api/v1'
        )
        if not config.oasis.require_authentication or any(
            pat.search(router_path) for pat in self.whitelist_patterns
        ):
            return await call_next(request)

        # Extract tokens (dependency injection isn’t available now)
        # Here any token would be allowed
        bearer_token = None
        if 'Authorization' in request.headers:
            parts = request.headers['Authorization'].split()
            if len(parts) == 2 and parts[0].lower() == 'bearer':
                bearer_token = parts[1]

        try:
            # NOTE: cannot handle `form_data` as the stream would be consumed
            _user = resolve_user(
                request=request,
                bearer_token=bearer_token,
                upload_token=request.query_params.get('token'),
                signature_token=request.query_params.get('signature_token'),
                required=True,
            )
        except HTTPException as exc:
            return Response(status_code=exc.status_code, content=exc.detail)

        return await call_next(request)


OASIS_AUTH_WHITELIST: dict[str, set[str]] = {
    'v1_app': {'^/auth', '^/info$', '^/extensions/', '^/openapi.json$'},
    'optimade_app': {'/extensions', '/info', '^/versions$'},
    'dcat_app': {'^/extensions/docs', '^/openapi.json$'},
    'h5grove_app': {'^/docs', '^/redoc$', '^/openapi.json$'},
}


@asynccontextmanager
async def lifespan(app: FastAPI):
    from nomad import infrastructure
    from nomad.cli.dev import get_gui_artifacts_js, get_gui_config
    from nomad.metainfo.elasticsearch_extension import entry_type
    from nomad.parsing.parsers import import_all_parsers

    import_all_parsers()

    # each subprocess is supposed disconnect and
    # connect again: https://jira.mongodb.org/browse/PYTHON-2090
    try:
        from mongoengine import disconnect

        disconnect()
    except Exception:
        pass

    entry_type.reload_quantities_dynamic()
    GuiFiles.gui_artifacts_data = get_gui_artifacts_js()
    GuiFiles.gui_env_data = get_gui_config()

    data = {
        'artifacts': GuiFiles.gui_artifacts_data,
        'gui_config': GuiFiles.gui_env_data,
    }
    GuiFiles.gui_data_etag = hashlib.md5(
        json.dumps(data).encode(), usedforsecurity=False
    ).hexdigest()

    infrastructure.setup()

    FastAPICache.init(backend=MongoBackend())

    if config.temporal.enabled:
        try:
            app.state.temporal_client = await get_client()
            yield
        except Exception as e:
            logger = get_logger(__name__)

            logger.error(f'Failed to connect to temporal', exc_info=e)
            raise
    else:
        yield


app = FastAPI(lifespan=lifespan)

app_base = config.services.api_base_path


def temporal_client() -> Client:
    return app.state.temporal_client


@app.get(f'{app_base}/alive')
async def alive():
    return 'I am, alive!'


@app.get('/-/health', status_code=status.HTTP_200_OK)
async def health():
    return {'healthcheck': 'ok'}


app.mount(f'{app_base}/api/v1', v1_app)
v1_app.add_middleware(
    OasisAuthenticationMiddleware,
    whitelist=OASIS_AUTH_WHITELIST['v1_app'],
)
v1_app.add_middleware(
    CORSMiddleware,  # CORS has to be the first to act on request
    allow_origins=['*'],
    allow_credentials=True,
    allow_methods=['*'],
    allow_headers=['*'],
    expose_headers=['Content-Disposition'],
)

if config.services.optimade_enabled:
    from .optimade import optimade_app

    app.mount(f'{app_base}/optimade', optimade_app)

    optimade_app.add_middleware(
        OasisAuthenticationMiddleware,
        whitelist=OASIS_AUTH_WHITELIST['optimade_app'],
    )

if config.services.dcat_enabled:
    from .dcat.main import app as dcat_app

    app.mount(f'{app_base}/dcat', dcat_app)
    dcat_app.add_middleware(
        OasisAuthenticationMiddleware,
        whitelist=OASIS_AUTH_WHITELIST['dcat_app'],
    )

if config.services.h5grove_enabled:
    from .h5grove_app import app as h5grove_app

    app.mount(f'{app_base}/h5grove', h5grove_app)
    h5grove_app.add_middleware(
        OasisAuthenticationMiddleware,
        whitelist=OASIS_AUTH_WHITELIST['h5grove_app'],
    )

# Add API plugins
for entry_point in config.plugins.entry_points.filtered_values():
    if isinstance(entry_point, APIEntryPoint):
        api_app = entry_point.load()
        assert isinstance(api_app, FastAPI), (
            f'Error loading entry point "{entry_point.id}": The load method of an API entry point must return a FastAPI instance'
        )
        app.mount(f'{app_base}/{entry_point.prefix}', api_app)

# Make sure to mount this last, as it is a catch-all routes that are not yet mounted.
app.mount(app_base, static_files_app)


@app.exception_handler(StarletteHTTPException)
async def http_exception_handler(request, exc):
    if exc.status_code != status.HTTP_404_NOT_FOUND:
        return await default_http_exception_handler(request, exc)

    try:
        accept = request.headers['accept']
    except Exception:
        accept = None

    if accept is not None and 'html' in accept:
        return HTMLResponse(
            status_code=status.HTTP_404_NOT_FOUND,
            content=f"""
        <html>
            <head><title>{config.meta.name}</title></head>
            <body>
                <h1>NOMAD app</h1>
                <h2>info</h2>
                {'<br/>'.join(f'{key}: {value}' for key, value in config.meta.model_dump().items())}
                <h2>apis</h2>
                <a href="{app_base}/api/v1/extensions/docs">NOMAD API v1</a><br/>
                <a href="{app_base}/optimade/v1/extensions/docs">Optimade API</a><br/>
                <a href="{app_base}/dcat/extensions/docs">DCAT API</a><br/>
            </body>
        </html>
        """,
        )

    return JSONResponse(
        status_code=status.HTTP_404_NOT_FOUND,
        content={
            'detail': 'Not found',
            'info': {
                'app': config.meta.model_dump(),
                'apis': {
                    'v1': {
                        'root': f'{app_base}/api/v1',
                        'dashboard': f'{app_base}/api/v1/extensions/docs',
                    },
                    'optimade': {
                        'root': f'{app_base}/optimade/v1',
                        'dashboard': f'{app_base}/optimade/v1/extensions/docs',
                    },
                    'dcat': {
                        'root': f'{app_base}/dcat',
                        'dashboard': f'{app_base}/dcat/extensions/docs',
                    },
                },
            },
        },
    )
