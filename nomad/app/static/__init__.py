import hashlib
import re
import os

from fastapi import FastAPI, Request, Response
from fastapi.exception_handlers import (
    http_exception_handler as default_http_exception_handler,
)
from starlette.staticfiles import (
    StaticFiles as StarletteStaticFiles,
    NotModifiedResponse,
)
from starlette.responses import PlainTextResponse, FileResponse
from starlette.datastructures import Headers

from nomad.config import config


app = FastAPI()

docs_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), 'docs'))
gui_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), 'gui'))
if not os.path.exists(gui_folder):
    os.makedirs(gui_folder)

configured_gui_folder = os.path.join(
    config.fs.working_directory, 'run', 'gui_configured'
)
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
            if_none_match = request_headers['if-none-match']
            match = re.match(StaticFiles.etag_re, if_none_match)
            if match is not None:
                if_none_match = match.group(2)
                etag = response_headers['etag']
                if if_none_match == etag:
                    return True
        except KeyError:
            pass

        return super().is_not_modified(response_headers, request_headers)


app.mount(
    '/docs', StaticFiles(directory=docs_folder, check_dir=False, html=True), name='docs'
)


class GuiFiles(StaticFiles):
    gui_artifacts_path: str | None = None
    gui_env_data: str | None = None
    gui_data_etag: str | None = None

    async def get_response(self, path: str, scope) -> Response:
        if path not in ['env.js', 'artifacts.js']:
            response = await super().get_response(path, scope)
        else:
            assert GuiFiles.gui_data_etag is not None, (
                'Etag for gui data was not initialized'
            )
            if path == 'env.js':
                response = PlainTextResponse(
                    GuiFiles.gui_env_data,
                    media_type='application/javascript',
                    headers=dict(etag=GuiFiles.gui_data_etag),
                )
            else:
                response = FileResponse(
                    GuiFiles.gui_artifacts_path,
                    media_type='application/javascript',
                    headers=dict(etag=GuiFiles.gui_data_etag),
                )

        request_headers = Headers(scope=scope)
        if self.is_not_modified(response.headers, request_headers):
            return NotModifiedResponse(response.headers)
        return response

    @classmethod
    def bootstrap(cls, artifacts_path, env_data):
        cls.gui_artifacts_path = artifacts_path
        cls.gui_env_data = env_data

        hash_obj = hashlib.md5(env_data.encode(), usedforsecurity=False)
        with open(cls.gui_artifacts_path, 'rb') as f:
            while chunk := f.read(2**20):
                hash_obj.update(chunk)
        cls.gui_data_etag = hash_obj.hexdigest()


gui_files_app = GuiFiles(directory=gui_folder, check_dir=False)
app.mount('/gui', gui_files_app, name='gui')

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
        response.headers['Cache-Control'] = (
            'max-age=0, no-cache, no-store, must-revalidate'
        )

    # The etags that we and starlette produce do not follow the RFC, because they do not
    # start with a " as the RFC specifies. Nginx considers them weak etags and will strip
    # these if gzip is enabled.
    if not response.headers.get('etag', '"').startswith('"'):
        response.headers['etag'] = f'"{response.headers.get("etag")}"'

    return response
