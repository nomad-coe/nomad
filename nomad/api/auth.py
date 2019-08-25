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

.. autofunction:: authenticate
"""
from flask import g, request
from flask_restplus import abort, Resource, fields
import functools
import jwt
import datetime
import hmac
import hashlib
import uuid

from nomad import config, processing, utils, infrastructure, datamodel

from .app import api, RFC3339DateTime


# Authentication scheme definitions, for swagger
api.authorizations = {
    'HTTP Basic Authentication': {
        'type': 'basic'
    },
    'OpenIDConnect Bearer Token': {
        'type': 'apiKey',
        'in': 'header',
        'name': 'Authorization'
    },
    'NOMAD upload token': {
        'type': 'apiKey',
        'in': 'query',
        'name': 'token'
    },
    'NOMAD signature': {
        'type': 'apiKey',
        'in': 'query',
        'name': 'signature_token'
    }
}


def generate_upload_token(user):
    """
    Generates a short user authenticating token based on its keycloak UUID.
    It can be used to authenticate users in less security relevant but short curl commands.

    It uses the users UUID as urlsafe base64 encoded payload with a HMACSHA1 signature.
    """
    payload = uuid.UUID(user.user_id).bytes
    signature = hmac.new(
        bytes(config.services.api_secret, 'utf-8'),
        msg=payload,
        digestmod=hashlib.sha1)

    return '%s.%s' % (
        utils.base64_encode(payload),
        utils.base64_encode(signature.digest()))


def verify_upload_token(token) -> str:
    """
    Verifies the upload token generated with :func:`generate_upload_token`.

    Returns: The user UUID or None if the toke could not be verified.
    """
    payload, signature = token.split('.')
    payload = utils.base64_decode(payload)
    signature = utils.base64_decode(signature)

    compare = hmac.new(
        bytes(config.services.api_secret, 'utf-8'),
        msg=payload,
        digestmod=hashlib.sha1)

    if signature != compare.digest():
        return None

    return str(uuid.UUID(bytes=payload))


def authenticate(
        basic: bool = False, upload_token: bool = False, signature_token: bool = False,
        required: bool = False, admin_only: bool = False):
    """
    A decorator to protect API endpoints with authentication. Uses keycloak access
    token to authenticate users. Other methods might apply. Will abort with 401
    if necessary.

    Arguments:
        basic: Also allow Basic HTTP authentication
        upload_token: Also allow upload_token
        signature_token: Also allow signed urls
        required: Authentication is required
        admin_only: Only the admin user is allowed to use the endpoint.
    """
    methods = ['OpenIDConnect Bearer Token']
    if basic:
        methods.append('HTTP Basic Authentication')
    if upload_token:
        methods.append('NOMAD upload token')
    if signature_token:
        methods.append('NOMAD signature')

    def decorator(func):
        @functools.wraps(func)
        @api.response(401, 'Not authorized, some data require authentication and authorization')
        @api.doc(security=methods)
        def wrapper(*args, **kwargs):
            g.user = None

            if upload_token and 'token' in request.args:
                token = request.args['token']
                user_id = verify_upload_token(token)
                if user_id is not None:
                    g.user = infrastructure.keycloak.get_user(user_id)

            elif signature_token and 'signature_token' in request.args:
                token = request.args.get('signature_token', None)
                if token is not None:
                    try:
                        decoded = jwt.decode(token, config.services.api_secret, algorithms=['HS256'])
                        user = datamodel.User.get(decoded['user'])
                        if user is None:
                            abort(401, 'User for the given signature does not exist')
                        else:
                            g.user = user
                    except KeyError:
                        abort(401, 'Token with invalid/unexpected payload')
                    except jwt.ExpiredSignatureError:
                        abort(401, 'Expired token')
                    except jwt.InvalidTokenError:
                        abort(401, 'Invalid token')

            elif 'token' in request.args:
                abort(401, 'Queram param token not supported for this endpoint')

            user_or_error = infrastructure.keycloak.authorize_flask(basic=basic)
            if user_or_error is not None:
                if isinstance(user_or_error, datamodel.User):
                    g.user = user_or_error
                else:
                    abort(401, message=user_or_error)

            if required and g.user is None:
                abort(401, message='Authentication is required for this endpoint')
            if admin_only and (g.user is None or not g.user.is_admin):
                abort(401, message='Only the admin user is allowed to use this endpoint')

            return func(*args, **kwargs)

        return wrapper

    return decorator


ns = api.namespace(
    'auth',
    description='Authentication related endpoints.')


user_model = api.model('User', {
    'user_id': fields.String(description='The users UUID.'),
    'name': fields.String('The publically visible user name.'),
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
    @api.doc('get_user')
    @api.marshal_with(user_model, skip_none=True, code=200, description='User info send')
    @authenticate(required=True, basic=True)
    def get(self):
        return g.user


token_model = api.model('Token', {
    'user': fields.Nested(user_model, skip_none=True),
    'token': fields.String(description='The short term token to sign URLs'),
    'expiries_at': RFC3339DateTime(desription='The time when the token expires')
})


@ns.route('/token')
class TokenResource(Resource):
    @api.doc('get_token')
    @api.marshal_with(token_model, skip_none=True, code=200, description='Token send')
    @authenticate(required=True)
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
            'expires_at': expires_at.isoformat(),
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

        # look in mongo
        try:
            upload = processing.Upload.get(upload_id)
            return g.user.user_id == upload.user_id

        except KeyError as e:
            logger = utils.get_logger(__name__, upload_id=upload_id, calc_id=calc_id)
            logger.error('Upload files without respective db entry')
            raise e

    return func
