# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

'''
This module comprises the nomad@FAIRDI APIs. Currently there is NOMAD's official api, and
we will soon at the optimade api. The app module also servers documentation, gui, and
alive.
'''
from flask import Flask, Blueprint, jsonify, url_for, abort, request, make_response
from flask_restplus import Api, representations
from flask_cors import CORS
from werkzeug.exceptions import HTTPException
from werkzeug.wsgi import DispatcherMiddleware  # pylint: disable=E0611
import os.path
import random
from structlog import BoundLogger
import collections
from mongoengine.base.datastructures import BaseList
import orjson

from nomad import config, utils as nomad_utils

from .api import blueprint as api_blueprint, api
from .optimade import blueprint as optimade_blueprint, api as optimade
from .docs import blueprint as docs_blueprint
from .dist import blueprint as dist_blueprint
from .gui import blueprint as gui_blueprint
from .encyclopedia import blueprint as encyclopedia_blueprint
from . import common


def dump_json(data):
    def default(data):
        if isinstance(data, collections.OrderedDict):
            return dict(data)

        if data.__class__.__name__ == 'BaseList':
            return list(data)

        raise TypeError

    return orjson.dumps(
        data, default=default,
        option=orjson.OPT_INDENT_2 | orjson.OPT_NON_STR_KEYS)


# replace the json implementation of flask_restplus
def output_json(data, code, headers=None):
    dumped = dump_json(data) + b'\n'

    resp = make_response(dumped, code)
    resp.headers.extend(headers or {})
    return resp


api.representation('application/json')(output_json)
optimade.representation('application/json')(output_json)


@property  # type: ignore
def specs_url(self):
    '''
    Fixes issue where swagger-ui makes a call to swagger.json over HTTP.
    This can ONLY be used on servers that actually use HTTPS.  On servers that use HTTP,
    this code should not be used at all.
    '''
    return url_for(self.endpoint('specs'), _external=True, _scheme='https')


if config.services.https:
    Api.specs_url = specs_url


app = Flask(__name__)
''' The Flask app that serves all APIs. '''

app.config.APPLICATION_ROOT = common.base_path  # type: ignore
app.config.RESTPLUS_MASK_HEADER = False  # type: ignore
app.config.RESTPLUS_MASK_SWAGGER = False  # type: ignore
app.config.SWAGGER_UI_OPERATION_ID = True  # type: ignore
app.config.SWAGGER_UI_REQUEST_DURATION = True  # type: ignore

app.config['SECRET_KEY'] = config.services.api_secret


def api_base_path_response(env, resp):
    resp('200 OK', [('Content-Type', 'text/plain')])
    return [
        ('Development nomad api server. Api is served under %s/.' %
            config.services.api_base_path).encode('utf-8')]


app.wsgi_app = DispatcherMiddleware(  # type: ignore
    api_base_path_response, {config.services.api_base_path: app.wsgi_app})

CORS(app)

app.register_blueprint(api_blueprint, url_prefix='/api')
app.register_blueprint(optimade_blueprint, url_prefix='/optimade')
app.register_blueprint(docs_blueprint, url_prefix='/docs')
app.register_blueprint(dist_blueprint, url_prefix='/dist')
app.register_blueprint(gui_blueprint, url_prefix='/gui')
app.register_blueprint(encyclopedia_blueprint, url_prefix='/encyclopedia')


@app.errorhandler(Exception)
def handle(error: Exception):
    status_code = getattr(error, 'code', 500)
    if not isinstance(status_code, int):
        status_code = 500
    if status_code < 100:
        status_code = 500

    name = getattr(error, 'name', 'Internal Server Error')
    description = getattr(error, 'description', 'No description available')
    data = dict(
        code=status_code,
        name=name,
        description=description)
    data.update(getattr(error, 'data', []))
    response = jsonify(data)
    response.status_code = status_code
    if status_code == 500:
        local_logger = common.logger
        # the logger is created in before_request, if the error was created before that
        # logger can be None
        if local_logger is None:
            local_logger = nomad_utils.get_logger(__name__)

        # TODO the error seems not to be the actual exception, therefore
        # there might be no stacktrace. Maybe there is a way to get the actual
        # exception/stacktrace
        local_logger.error('internal server error', error=str(error), exc_info=error)

    return response


@app.route('/alive')
def alive():
    ''' Simple endpoint to utilize kubernetes liveness/readiness probing. '''
    return "I am, alive!"


@app.before_request
def before_request():
    # api logger
    args = getattr(request, 'view_args')
    if args is None:
        args = {}
    else:
        args = dict(**args)

    args.update(
        name=__name__,
        blueprint=str(request.blueprint),
        endpoint=request.endpoint,
        method=request.method,
        url=request.url,
        json=request.json,
        args=request.args)

    common.logger = nomad_utils.get_logger(**args)

    # chaos monkey
    if config.services.api_chaos > 0:
        if random.randint(0, 100) <= config.services.api_chaos:
            abort(random.choice([400, 404, 500]), 'With best wishes from the chaos monkey.')


@app.before_first_request
def setup():
    from nomad import infrastructure

    if not app.config['TESTING']:
        # each subprocess is supposed disconnect connect again: https://jira.mongodb.org/browse/PYTHON-2090
        try:
            from mongoengine import disconnect
            disconnect()
        except Exception:
            pass

        infrastructure.setup()
