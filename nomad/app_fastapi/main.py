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

from fastapi import FastAPI, status, Request
from fastapi.responses import JSONResponse, RedirectResponse
from fastapi.middleware.wsgi import WSGIMiddleware
import traceback

from nomad import config, utils
from nomad.app import app as flask_app

from .routers import users, entries, auth, datasets
from .optimade import optimade_app


logger = utils.get_logger(__name__)

app = FastAPI(
    openapi_url='%s/api/v1/openapi.json' % config.services.api_base_path,
    docs_url='%s/api/v1/extensions/docs' % config.services.api_base_path,
    redoc_url='%s/api/v1/extensions/redoc' % config.services.api_base_path,
    swagger_ui_oauth2_redirect_url='%s/api/v1/docs/oauth2-redirect' % config.services.api_base_path,

    title='NOMAD API',
    version='v1, NOMAD %s@%s' % (config.meta.version, config.meta.commit),
    description=utils.strip('''
        **Disclaimer!** This is the new NOMAD API. It is still under development and only includes a
        part of the NOMAD API functionality. You can still use the old flask-based API
        as `/api` and the optimade API as `/optimade/v1`.

        ## Getting started

        ... TODO put the examples and tutorial here ...

        ## Conventions

        ### Paths

        The various API operations are organized with the following path scheme. The first
        part of the path, describes the data entity that is covered by
        the operations below (e.g. `entries`, `users`, `datasets`, `uploads`). For example
        everything below `entries` will be about searching entries, getting
        an entry, editing entries, etc.

        The second (optional and variable) path segment allows to denote a specific entity instance,
        e.g. a specific entry or dataset, usually by id. With out such a variable second
        path segment, its about all instances, e.g. searching entries or listing all datasets.

        Optional (if available) further path segments will determine the variety and format
        of data. This is mostly for entries to distinguish the metadata, raw, and archive
        data or distinguish between listing (i.e. paginated json) and downloading
        (i.e. streaming a zip-file)

        Further, we try to adhere to the paradim of getting and posting resources. Therefore,
        when you post a complex query, you will not post it to `/entries` (a query is not an entry),
        but `/entries/query`. Here *query* being a kind of virtual resource.

        ### Parameters and bodies for GET and POST operations

        We offer **GET** and **POST** versions for many complex operations. The idea is that
        **GET** is easy to use, e.g. via curl or simply in the browser, while **POST**
        allows to provide more complex parameters (i.e. a JSON body). For example to
        search for entries, you can use the **GET** operation `/entries` to specify simple
        queries via URL, e.g. `/entries?code_name=VASP&atoms=Ti`, but you would use
        **POST** `/entries/query` to provide a complex nested queries, e.g. with logical
        operators.

        Typicall the **POST** version is a super-set of the functionality of the **GET**
        version. But, most top-level parameters in the **POST** body, will be available
        in the **GET** version as URL parameters with the same name and meaning. This
        is especially true for reoccuring parameters for general API concepts like pagination
        or specifying required result fields.

        ### Response layout

        Typically a response will mirror all input parameters in the normalized form that
        was used to perform the operation.

        Some of these will be augmented with result values. For example the pagination
        section of a request will be augmented with the total available number.

        The actual requested data, will be placed under the key `data`.

        ## About Authentication

        NOMAD is an open datasharing platform, and most of the API operations do not require
        any authorization and can be freely used without a user or credentials. However,
        to upload data, edit data, or view your own and potentially unpublished data,
        the API needs to authenticate you.

        The NOMAD API uses OAuth and tokens to authenticate users. We provide simple operations
        that allow you to acquire an *access token* via username and password based
        authentication (`/auth/token`). The resulting access token can then be used on all operations
        (e.g. that support or require authentication).

        To use authentication in the dashboard, simply use the Authorize button. The
        dashboard GUI will manage the access token and use it while you try out the various
        operations.
    '''))


async def redirect_to_docs(req: Request):
    return RedirectResponse('%s/api/v1/extensions/docs' % config.services.api_base_path)


app.add_route('%s/api/v1' % config.services.api_base_path, redirect_to_docs, include_in_schema=False)
app.add_route('%s/api/v1/' % config.services.api_base_path, redirect_to_docs, include_in_schema=False)


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

app.include_router(auth.router, prefix='%s/api/v1/auth' % config.services.api_base_path)
app.include_router(users.router, prefix='%s/api/v1/users' % config.services.api_base_path)
app.include_router(entries.router, prefix='%s/api/v1/entries' % config.services.api_base_path)
app.include_router(datasets.router, prefix='%s/api/v1/datasets' % config.services.api_base_path)

app.mount('%s/optimade' % config.services.api_base_path, optimade_app)
app.mount(config.services.api_base_path, WSGIMiddleware(flask_app))
