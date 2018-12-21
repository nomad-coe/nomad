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
Endpoints can use *flask_httpauth* based authentication either with basic HTTP
authentication or access tokens. Currently the authentication is validated against
users and sessions in the NOMAD-coe repository postgres db.

.. autodata:: base_path

There are two authentication "schemes" to authenticate users. First we use
HTTP Basic Authentication (username, password), which also works with username=token,
password=''. Second, there is a curstom HTTP header 'X-Token' that can be used to
give a token. The first precedes the second. The used tokens are given and stored
by the NOMAD-coe repository GUI.

Authenticated user information is available via FLASK's build in flask.g.user object.
It is set to None, if no user information is available.

There are two decorators for FLASK API endpoints that can be used if endpoints require
authenticated user information for authorization or otherwise.

.. autofunction:: login_if_available
.. autofunction:: login_really_required
"""

from flask import g, request
from flask_restplus import abort
from flask_httpauth import HTTPBasicAuth

from nomad import config
from nomad.coe_repo import User

from .app import app, base_path

app.config['SECRET_KEY'] = config.services.api_secret
auth = HTTPBasicAuth()


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


def login_if_available(func):
    """
    A decorator for API endpoint implementations that might authenticate users, but
    provide limited functionality even without users.
    """
    @auth.login_required
    def wrapper(*args, **kwargs):
        # TODO the cutom X-Token based authentication should be replaced by a real
        # Authentication header based token authentication
        if not g.user and 'X-Token' in request.headers:
            token = request.headers['X-Token']
            g.user = User.verify_auth_token(token)
            if not g.user:
                abort(401, message='Provided access token is not valid or does not exist.')

        return func(*args, **kwargs)

    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    return wrapper


def login_really_required(func):
    """
    A decorator for API endpoint implementations that forces user authentication on
    endpoints.
    """
    @login_if_available
    def wrapper(*args, **kwargs):
        if g.user is None:
            abort(401, message='Anonymous access is forbidden, authorization required')
        else:
            return func(*args, **kwargs)
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    return wrapper


@app.route('%s/token' % base_path)
@login_really_required
def get_auth_token():
    """
    Get a token for authenticated users. This is currently disabled and all authentication
    matters are solved by the NOMAD-coe repository GUI.

    .. :quickref: Get a token to authenticate the user in follow up requests.

    :resheader Content-Type: application/json
    :status 200: calc successfully retrieved
    :returns: an authentication token that is valid for 10 minutes.
    """
    assert False, 'All authorization is none via NOMAD-coe repository GUI'
    # TODO all authorization is done via NOMAD-coe repository GUI
    # token = g.user.generate_auth_token(600)
    # return jsonify({'token': token.decode('ascii'), 'duration': 600})
