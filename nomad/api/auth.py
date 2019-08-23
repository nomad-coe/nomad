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
The API is protected with *keycloak* and *OpenIDConnect*. All API endpoints that require
or support authentication accept OIDC bearer tokens via HTTP header (``Authentication``,
recommended), query (``access_token``), or form parameter (``access_token``). These
token can be acquired from the NOMAD keycloak server or through the ``/auth`` endpoint
that also supports HTTP Basic authentication and passes the given credentials to
keycloak.

Authenticated user information is available via FLASK's build in flask.g.user object.
It is set to None, if no user information is available.

There are three decorators for FLASK API endpoints that can be used to protect
endpoints that require or support authentication.

.. autofunction:: login_if_available
.. autofunction:: login_really_required
.. autofunction:: admin_login_required
"""

from flask import g, request
from flask_restplus import abort, Resource, fields
import functools
import jwt
import datetime

from nomad import config, processing, files, utils, infrastructure, datamodel

from .app import api, RFC3339DateTime


def login_if_available(token_only: bool = True):
    """
    A decorator for API endpoint implementations that might authenticate users, but
    provide limited functionality even without users.
    """
    def decorator(func):
        @functools.wraps(func)
        @api.response(401, 'Not authorized, some data require authentication and authorization')
        @api.doc(security=list('OpenIDConnect Bearer Token'))
        def wrapper(*args, **kwargs):
            user_or_error = infrastructure.keycloak.authorize_flask(token_only)
            if user_or_error is None:
                pass
            elif isinstance(user_or_error, datamodel.User):
                g.user = user_or_error
            else:
                abort(401, message=user_or_error)

            return func(*args, **kwargs)

        return wrapper

    return decorator


def login_really_required(token_only: bool = True):
    """
    A decorator for API endpoint implementations that forces user authentication on
    endpoints.
    """
    def decorator(func):
        @functools.wraps(func)
        @api.response(401, 'Not authorized, this endpoint requires authorization')
        @login_if_available(token_only)
        def wrapper(*args, **kwargs):
            if g.user is None:
                abort(401, 'Not authorized, this endpoint requires authorization')

            return func(*args, **kwargs)

        return wrapper

    return decorator


def admin_login_required(func):
    """
    A decorator for API endpoint implementations that should only work for the admin user.
    """
    @functools.wraps(func)
    @api.response(401, 'Authentication required or not authorized as admin user. Only admin can access this endpoint.')
    @login_really_required
    def wrapper(*args, **kwargs):
        if not g.user.is_admin:
            abort(401, message='Only the admin user use this endpoint')

        return func(*args, **kwargs)

    return wrapper


ns = api.namespace(
    'auth',
    description='Authentication related endpoints.')


user_model = api.model('User', {
    'user_id': fields.Integer(description='The id to use in the repo db, make sure it does not already exist.'),
    'first_name': fields.String(description='The user\'s first name'),
    'last_name': fields.String(description='The user\'s last name'),
    'email': fields.String(description='Guess what, the user\'s email'),
    'affiliation': fields.Nested(model=api.model('Affiliation', {
        'name': fields.String(description='The name of the affiliation', default='not given'),
        'address': fields.String(description='The address of the affiliation', default='not given')})),
    'password': fields.String(description='The bcrypt 2y-indented password for initial and changed password'),
    'token': fields.String(
        description='The access token that authenticates the user with the API. '
        'User the HTTP header "X-Token" to provide it in API requests.'),
    'created': RFC3339DateTime(description='The create date for the user.')
})


@ns.route('/')
class AuthResource(Resource):
    @api.doc('get_token')
    @api.marshal_with(user_model, skip_none=True, code=200, description='User info send')
    @login_really_required(token_only=False)
    def get(self):
        return g.user


token_model = api.model('Token', {
    'user': fields.Nested(user_model),
    'token': fields.String(description='The short term token to sign URLs'),
    'expiries_at': RFC3339DateTime(desription='The time when the token expires')
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
        expires_at = datetime.datetime.utcnow() + datetime.timedelta(seconds=10)
        token = jwt.encode(
            dict(user=g.user.user_id, exp=expires_at),
            config.services.api_secret, 'HS256').decode('utf-8')

        return {
            'user': g.user,
            'token': token,
            'expires_at': expires_at.isoformat()
        }


def with_signature_token(func):
    """
    A decorator for API endpoint implementations that validates signed URLs.
    """
    @functools.wraps(func)
    @api.response(401, 'Invalid or expired signature token')
    def wrapper(*args, **kwargs):
        token = request.args.get('token', None)
        if token is not None:
            try:
                decoded = jwt.decode(token, config.services.api_secret, algorithms=['HS256'])
                user = datamodel.User.get(decoded['user'])
                if user is None:
                    abort(401, 'User for token does not exist')
                else:
                    g.user = user
            except KeyError:
                abort(401, 'Token with invalid/unexpected payload')
            except jwt.ExpiredSignatureError:
                abort(401, 'Expired token')
            except jwt.InvalidTokenError:
                abort(401, 'Invalid token')

        return func(*args, **kwargs)

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
        elif g.user.is_admin:
            # the admin user does have authorization to access everything
            return True

        # look in mongodb
        processing.Upload.get(upload_id).user_id == g.user.user_id

        # There are no db entries for the given resource
        if files.UploadFiles.get(upload_id) is not None:
            logger = utils.get_logger(__name__, upload_id=upload_id, calc_id=calc_id)
            logger.error('Upload files without respective db entry')

        raise KeyError
    return func
