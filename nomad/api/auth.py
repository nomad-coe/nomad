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

from flask import g, request, make_response
from flask_restplus import abort, Resource
from flask_httpauth import HTTPBasicAuth

from nomad import config, processing, files, utils, coe_repo
from nomad.coe_repo import User, LoginException

from .app import app, api

app.config['SECRET_KEY'] = config.services.api_secret
auth = HTTPBasicAuth()


# Authentication scheme definitions, for swagger only.
api.authorizations = {
    'HTTP Basic': {
        'type': 'basic'
    },
    'X-Token': {
        'type': 'apiKey',
        'in': 'header',
        'name': 'X-Token'
    }
}


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
    @api.response(401, 'Not authorized, some data require authentication and authorization')
    @api.doc(security=list(api.authorizations.keys()))
    @auth.login_required
    def wrapper(*args, **kwargs):
        # TODO the cutom X-Token based authentication should be replaced by a real
        # Authentication header based token authentication
        if not g.user and 'X-Token' in request.headers:
            token = request.headers['X-Token']
            g.user = User.verify_auth_token(token)
            if not g.user:
                abort(401, message='Not authorized, some data require authentication and authorization')

        return func(*args, **kwargs)

    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    return wrapper


def login_really_required(func):
    """
    A decorator for API endpoint implementations that forces user authentication on
    endpoints.
    """
    @api.response(401, 'Authentication required or not authorized to access requested data')
    @api.doc(security=list(api.authorizations.keys()))
    @login_if_available
    def wrapper(*args, **kwargs):
        if g.user is None:
            abort(401, message='Authentication required or not authorized to access requested data')
        else:
            return func(*args, **kwargs)
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    return wrapper


ns = api.namespace(
    'auth',
    description='Authentication related endpoints.')


@ns.route('/token')
class TokenResource(Resource):
    @api.doc('get_token')
    @api.response(200, 'Token send', headers={'Content-Type': 'text/plain; charset=utf-8'})
    @login_really_required
    def get(self):
        """
        Get the access token for the authenticated user.

        You can use basic authentication to access this endpoint and receive a
        token for further api access. This token will expire at some point and presents
        a more secure method of authentication.
        """
        try:
            response = make_response(g.user.get_auth_token().decode('utf-8'))
            response.headers['Content-Type'] = 'text/plain; charset=utf-8'
            return response
        except LoginException:
            abort(
                401,
                message='You are not propertly logged in at the NOMAD coe repository, '
                        'there is no token for you.')


def create_authorization_predicate(upload_hash, calc_hash=None):
    """
    Returns a predicate that determines if the logged in user has the authorization
    to access the given upload and calculation.
    """
    def func():
        if g.user is None:
            # guest users don't have authorized access to anything
            return False

        # look in repository
        upload = coe_repo.Upload.from_upload_hash(upload_hash)
        if upload is not None:
            return upload.user_id == g.user.user_id

        # look in staging
        staging_upload = processing.Upload.get(upload_hash)
        if staging_upload is not None:
            return str(g.user.user_id) == str(staging_upload.user_id)

        # There are no db entries for the given resource
        if files.UploadFiles.get(upload_hash) is not None:
            logger = utils.get_logger(__name__, upload_hash=upload_hash, calc_hash=calc_hash)
            logger.error('Upload files without respective db entry')

        raise KeyError
    return func
