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

"""
All APIs are served by one Flask app (:py:mod:`nomad.api.app`) under different paths.
"""

from flask import Flask, jsonify, url_for
from flask_restplus import Api, fields
from flask_cors import CORS
from werkzeug.exceptions import HTTPException
from werkzeug.wsgi import DispatcherMiddleware
import os.path
import inspect
from datetime import datetime
import pytz

from nomad import config, utils

base_path = config.services.api_base_path
""" Provides the root path of the nomad APIs. """


@property  # type: ignore
def specs_url(self):
    """
    Fixes issue where swagger-ui makes a call to swagger.json over HTTP.
    This can ONLY be used on servers that actually use HTTPS.  On servers that use HTTP,
    this code should not be used at all.
    """
    return url_for(self.endpoint('specs'), _external=True, _scheme='https')


if config.services.https:
    Api.specs_url = specs_url


app = Flask(
    __name__,
    static_url_path='/docs',
    static_folder=os.path.abspath(os.path.join(os.path.dirname(__file__), '../../docs/.build/html')))
""" The Flask app that serves all APIs. """

app.config.APPLICATION_ROOT = base_path  # type: ignore
app.config.RESTPLUS_MASK_HEADER = False  # type: ignore
app.config.RESTPLUS_MASK_SWAGGER = False  # type: ignore
app.config.SWAGGER_UI_OPERATION_ID = True  # type: ignore
app.config.SWAGGER_UI_REQUEST_DURATION = True  # type: ignore


def api_base_path_response(env, resp):
    resp('200 OK', [('Content-Type', 'text/plain')])
    return [
        ('Development nomad api server. Api is served under %s/.' %
            config.services.api_base_path).encode('utf-8')]


app.wsgi_app = DispatcherMiddleware(  # type: ignore
    api_base_path_response, {config.services.api_base_path: app.wsgi_app})


CORS(app)

api = Api(
    app, version='1.0', title='nomad@FAIRDI API',
    description='Official API for nomad@FAIRDI services.',
    validate=True)
""" Provides the flask restplust api instance """


@app.errorhandler(Exception)
@api.errorhandler
def handle(error: Exception):
    status_code = getattr(error, 'code', 500)
    if not isinstance(status_code, int):
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
        utils.get_logger(__name__).error('internal server error', exc_info=error)
    return response


@app.route('/alive')
def alive():
    """ Simply endpoint to utilize kubernetes liveness/readiness probing. """
    return "I am, alive!"


def with_logger(func):
    """
    Decorator for endpoint implementations that provides a pre configured logger and
    automatically logs errors on all 500 responses.
    """
    signature = inspect.signature(func)
    has_logger = 'logger' in signature.parameters
    wrapper_signature = signature.replace(parameters=tuple(
        param for param in signature.parameters.values()
        if param.name != 'logger'
    ))

    def wrapper(*args, **kwargs):
        if has_logger:
            args = inspect.getcallargs(wrapper, *args, **kwargs)
            logger_args = {
                k: v for k, v in args.items()
                if k in ['upload_id', 'calc_id']}
            logger = utils.get_logger(__name__, **logger_args)
            args.update(logger=logger)
        try:
            return func(**args)
        except HTTPException as e:
            if getattr(e, 'code', None) == 500:
                logger.error('Internal server error', exc_info=e)
            raise e
        except Exception as e:
            logger.error('Internal server error', exc_info=e)
            raise e

    wrapper.__signature__ = wrapper_signature
    return wrapper


class RFC3339DateTime(fields.DateTime):

    def format(self, value):
        if isinstance(value, datetime):
            return super().format(value.replace(tzinfo=pytz.utc))
        else:
            return str(value)


rfc3339DateTime = RFC3339DateTime()
