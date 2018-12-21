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

from flask import Flask, jsonify
from flask_restplus import Api
from flask_cors import CORS
from werkzeug.exceptions import HTTPException
import os.path

from nomad import config

base_path = config.services.api_base_path
""" Provides the root path of the nomad APIs. """

app = Flask(
    __name__,
    static_url_path='%s/docs' % base_path,
    static_folder=os.path.abspath(os.path.join(os.path.dirname(__file__), '../../docs/.build/html')))
""" The Flask app that serves all APIs. """

CORS(app)

api = Api(app)
""" Provides the flask restful api instance """


@app.errorhandler(HTTPException)
def handle(error):
    status_code = getattr(error, 'code', 500)
    name = getattr(error, 'name', 'Internal Server Error')
    description = getattr(error, 'description', None)
    data = dict(
        code=status_code,
        name=name,
        description=description)
    data.update(getattr(error, 'data', []))
    response = jsonify(data)
    response.status_code = status_code
    return response
