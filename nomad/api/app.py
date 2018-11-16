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

from flask import Flask, g
from flask_restful import Api, abort
from flask_cors import CORS
from flask_httpauth import HTTPBasicAuth
import os.path

from nomad import config, infrastructure
from nomad.user import User
from nomad.processing import Upload

base_path = config.services.api_base_path

app = Flask(
    __name__,
    static_url_path='%s/docs' % base_path,
    static_folder=os.path.abspath(os.path.join(os.path.dirname(__file__), '../docs/.build/html')))
CORS(app)

app.config['SECRET_KEY'] = config.services.api_secret

auth = HTTPBasicAuth()
api = Api(app)


@app.before_first_request
def setup():
    infrastructure.setup()


@auth.verify_password
def verify_password(username_or_token, password):
    # first try to authenticate by token
    g.user = User.verify_auth_token(username_or_token)
    if not g.user:
        # try to authenticate with username/password
        try:
            g.user = User.verify_user_password(username_or_token, password)
        except Exception:
            return False

    if not g.user:
        return True  # anonymous access

    return True


def login_really_required(func):
    @auth.login_required
    def wrapper(*args, **kwargs):
        if g.user is None:
            abort(401, message='Anonymous access is forbidden, authorization required')
        else:
            return func(*args, **kwargs)
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    return wrapper


@app.route('/api/token')
@login_really_required
def get_auth_token():
    assert False, 'All authorization is none via NOMAD-coe repository GUI'
    # TODO all authorization is done via NOMAD-coe repository GUI
    # token = g.user.generate_auth_token(600)
    # return jsonify({'token': token.decode('ascii'), 'duration': 600})


@app.route('%s/admin/<string:operation>' % base_path, methods=['POST'])
def call_admin_operation(operation):
    if operation == 'repair_uploads':
        Upload.repair_all()
    if operation == 'reset':
        infrastructure.reset()
    else:
        abort(400, message='Unknown operation %s' % operation)

    return 'done', 200
