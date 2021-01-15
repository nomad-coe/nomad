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

from fastapi import FastAPI
from fastapi.middleware.wsgi import WSGIMiddleware

from nomad import config

from .optimade import optimade_app
from .flask import app as flask_app
from .v1.main import app as v1_app


app = FastAPI()

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
