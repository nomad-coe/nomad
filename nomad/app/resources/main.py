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
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
import traceback
from celery.signals import worker_process_init

from nomad import config, utils

from .routers import resources


logger = utils.get_logger(__name__)

mongo_client_resources = None

app = FastAPI(
    openapi_url='/openapi.json',
    docs_url='/extensions/docs',
    redoc_url='/extensions/redoc',
    swagger_ui_oauth2_redirect_url='/extensions/docs/oauth2-redirect',
    title='Resources API',
    version='v1, NOMAD %s@%s' % (config.meta.version, config.meta.commit),
    description='NOMAD\'s API for serving related external resources')

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


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


@app.middleware('http')
async def setup_fastapi(request: Request, callback):

    setup_mongo()

    response = await callback(request)
    return response


@worker_process_init.connect
def setup_celery(**kwargs):
    setup_mongo()


def setup_mongo():
    global mongo_client_resources

    if mongo_client_resources is None:
        from mongoengine import connect
        mongo_client_resources = connect(
            db=config.resources.db_name, alias='resources', host=config.mongo.host, port=config.mongo.port)


def remove_mongo():
    setup_mongo()
    mongo_client_resources.drop_database(config.resources.db_name)


app.include_router(resources.router)
