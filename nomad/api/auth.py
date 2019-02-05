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
from flask_restplus import abort, Resource, fields
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
    if username_or_token is None or username_or_token == '':
        g.user = None
        return True

    if password is None or password == '':
        g.user = User.verify_auth_token(username_or_token)
        return g.user is not None
    else:
        try:
            g.user = User.verify_user_password(username_or_token, password)
        except Exception as e:
            utils.get_logger(__name__).error('could not verify password', exc_info=e)
            return False

        return g.user is not None


@auth.error_handler
def auth_error_handler():
    abort(401, 'Could not authenticate user, bad credentials')


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


user_model = api.model('User', {
    'first_name': fields.String(description='The user\'s first name'),
    'last_name': fields.String(description='The user\'s last name'),
    'email': fields.String(description='Guess what, the user\'s email'),
    'affiliation': fields.String(description='The user\'s affiliation'),
    'token': fields.String(
        description='The access token that authenticates the user with the API. '
        'User the HTTP header "X-Token" to provide it in API requests.')
})


@ns.route('/user')
class UserResource(Resource):
    @api.doc('get_user')
    @api.marshal_with(user_model, skip_none=True, code=200, description='User data send')
    @login_really_required
    def get(self):
        """
        Get user information including a long term access token for the authenticated user.

        You can use basic authentication to access this endpoint and receive a
        token for further api access. This token will expire at some point and presents
        a more secure method of authentication.
        """
        try:
            return g.user
        except LoginException:
            abort(
                401,
                message='User not logged in, provide credentials via Basic HTTP authentication.')


token_model = api.model('Token', {
    'user': fields.Nested(user_model),
    'token': fields.String(description='The short term token to sign URLs'),
    'experies_at': fields.DateTime(desription='The time when the token expires')
})


signature_token_argument = dict(
    name='token', type=str, help='Token that signs the URL and authenticates the user',
    location='args')


@ns.route('/token')
class TokenResource(Resource):
    @api.doc('get_token')
    @api.marshal_with(token_model, skip_none=True, code=200, description='Token send')
    @login_really_required
    def get(self):
        """
        Generates a short (10s) term JWT token that can be used to authenticate the user in
        URLs towards most API get request, e.g. for file downloads on the
        raw or archive api endpoints. Use the token query parameter to sign URLs.
        """
        token, expires_at = g.user.get_signature_token()
        return {
            'user': g.user,
            'token': token,
            'expires_at': expires_at.isoformat()
        }


def with_signature_token(func):
    """
    A decorator for API endpoint implementations that validates signed URLs.
    """
    @api.response(401, 'Invalid or expired signature token')
    def wrapper(*args, **kwargs):
        token = request.args.get('token', None)
        if token is not None:
            try:
                g.user = coe_repo.User.verify_signature_token(token)
            except LoginException:
                abort(401, 'Invalid or expired signature token')

        return func(*args, **kwargs)
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    return wrapper


def create_authorization_predicate(upload_id, calc_id=None):
    """
    Returns a predicate that determines if the logged in user has the authorization
    to access the given upload and calculation.
    """
    def func():
        if g.user is None:
            # guest users don't have authorized access to anything
            return False

        # look in repository
        upload = coe_repo.Upload.from_upload_id(upload_id)
        if upload is not None:
            return upload.user_id == g.user.user_id

        # look in staging
        staging_upload = processing.Upload.get(upload_id)
        if staging_upload is not None:
            return str(g.user.user_id) == str(staging_upload.user_id)

        # There are no db entries for the given resource
        if files.UploadFiles.get(upload_id) is not None:
            logger = utils.get_logger(__name__, upload_id=upload_id, calc_id=calc_id)
            logger.error('Upload files without respective db entry')

        raise KeyError
    return func
