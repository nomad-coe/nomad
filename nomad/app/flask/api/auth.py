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

'''
The API is protected with *keycloak* and *OpenIDConnect*. All API endpoints that require
or support authentication accept OIDC bearer tokens via HTTP header (``Authentication``).
These token can be acquired from the NOMAD keycloak server or through the ``/auth`` endpoint
that also supports HTTP Basic authentication and passes the given credentials to
keycloak. For GUI's it is recommended to accquire an access token through the regular OIDC
login flow.

Authenticated user information is available via FLASK's build in flask.g.user object.
It is set to None, if no user information is available. To protect endpoints use the following
decorator.

.. autofunction:: authenticate

To allow authentification with signed urls, use this decorator:

.. autofunction:: with_signature_token
'''
from flask import g, request
from flask_restplus import abort, Resource, fields
import functools
import jwt
import datetime
import hmac
import hashlib
import uuid

from nomad import config, processing, utils, infrastructure, datamodel
from nomad.metainfo.flask_extension import generate_flask_restplus_model

from .api import api


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


def authenticate(
        basic: bool = False, upload_token: bool = False, signature_token: bool = False,
        required: bool = False, admin_only: bool = False):
    '''
    A decorator to protect API endpoints with authentication. Uses keycloak access
    token to authenticate users. Other methods might apply. Will abort with 401
    if necessary.

    Arguments:
        basic: Also allow Basic HTTP authentication
        upload_token: Also allow upload_token
        signature_token: Also allow signed urls
        required: Authentication is required
        admin_only: Only the admin user is allowed to use the endpoint.
    '''
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
                try:
                    token = request.args['token']
                    payload, signature = token.split('.')
                    payload = utils.base64_decode(payload)
                    signature = utils.base64_decode(signature)

                    compare = hmac.new(
                        bytes(config.services.api_secret, 'utf-8'),
                        msg=payload,
                        digestmod=hashlib.sha1)

                    if signature != compare.digest():
                        return None

                    user_id = str(uuid.UUID(bytes=payload))
                    g.user = infrastructure.keycloak.get_user(user_id)
                except KeyError:
                    abort(401, 'Invalid token')

            elif signature_token and 'signature_token' in request.args:
                token = request.args.get('signature_token', None)
                try:
                    decoded = jwt.decode(token, config.services.api_secret, algorithms=['HS256'])
                    g.user = datamodel.User.get(user_id=decoded['user'])
                except KeyError:
                    abort(401, 'Token with invalid/unexpected payload')
                except jwt.ExpiredSignatureError:
                    abort(401, 'Expired token')
                except jwt.InvalidTokenError:
                    abort(401, 'Invalid token')

            elif 'token' in request.args:
                abort(401, 'Query param token not supported for this endpoint')

            elif 'signature_token' in request.args:
                abort(401, 'Query param signature_token not supported for this endpoint')

            else:
                try:
                    g.user, g.oidc_access_token = infrastructure.keycloak.auth(request.headers, allow_basic=basic)
                except infrastructure.KeycloakError as e:
                    abort(401, message=str(e))

            if config.oasis.allowed_users is not None:
                if g.user is None:
                    abort(401, message='Authentication is required for this Oasis')
                if g.user.email not in config.oasis.allowed_users:
                    abort(401, message='You are not authorized to access this Oasis')

            if required and g.user is None:
                abort(401, message='Authentication is required for this endpoint')
            if admin_only and (g.user is None or not g.user.is_admin):
                abort(401, message='Only the admin user is allowed to use this endpoint')

            return func(*args, **kwargs)

        return wrapper

    return decorator


def generate_upload_token(user):
    payload = uuid.UUID(user.user_id).bytes
    signature = hmac.new(
        bytes(config.services.api_secret, 'utf-8'),
        msg=payload,
        digestmod=hashlib.sha1)

    return '%s.%s' % (
        utils.base64_encode(payload),
        utils.base64_encode(signature.digest()))


ns = api.namespace(
    'auth',
    description='Authentication related endpoints.')


auth_model = api.model('Auth', {
    'access_token': fields.String(description='The OIDC access token'),
    'upload_token': fields.String(description='A short token for human readable upload URLs'),
    'signature_token': fields.String(description='A short term token to sign URLs')
})


@ns.route('/')
class AuthResource(Resource):
    @api.doc('get_auth')
    @api.marshal_with(auth_model, skip_none=True, code=200, description='Auth info send')
    @authenticate(required=True, basic=True)
    def get(self):
        '''
        Provides authentication information. This endpoint requires authentification.
        Like all endpoints the OIDC access token based authentification. In additional,
        basic HTTP authentification can be used. This allows to login and acquire an
        access token.

        The response contains a short (10s) term JWT token that can be used to sign
        URLs with a ``signature_token`` query parameter, e.g. for file downloads on the
        raw or archive api endpoints; a short ``upload_token`` that is used in
        ``curl`` command line based uploads; and the OIDC JWT access token.
        '''

        def signature_token():
            expires_at = datetime.datetime.utcnow() + datetime.timedelta(seconds=10)
            return jwt.encode(
                dict(user=g.user.user_id, exp=expires_at),
                config.services.api_secret, 'HS256').decode('utf-8')

        try:
            return {
                'upload_token': generate_upload_token(g.user),
                'signature_token': signature_token(),
                'access_token': g.oidc_access_token
            }

        except KeyError:
            abort(401, 'The authenticated user does not exist')


user_model = generate_flask_restplus_model(api, datamodel.User.m_def)
users_model = api.model('UsersModel', {
    'users': fields.Nested(user_model, skip_none=True)
})


users_parser = api.parser()
users_parser.add_argument(
    'query', default='',
    help='Only return users that contain this string in their names, usernames, or emails.')


@ns.route('/users')
class UsersResource(Resource):
    @api.doc('get_users')
    @api.marshal_with(users_model, code=200, description='User suggestions send')
    @api.expect(users_parser, validate=True)
    def get(self):
        ''' Get existing users. '''
        args = users_parser.parse_args()

        return dict(users=infrastructure.keycloak.search_user(args.get('query')))

    @api.doc('invite_user')
    @api.marshal_with(user_model, code=200, skip_none=True, description='User invited')
    @api.expect(user_model, validate=True)
    def put(self):
        ''' Invite a new user. '''
        if config.keycloak.oasis:
            abort(400, 'User invide does not work this NOMAD OASIS')

        json_data = request.get_json()
        try:
            user = datamodel.User.m_from_dict(json_data)
        except Exception as e:
            abort(400, 'Invalid user data: %s' % str(e))

        if user.email is None:
            abort(400, 'Invalid user data: email is required')

        try:
            error = infrastructure.keycloak.add_user(user, invite=True)
        except KeyError as e:
            abort(400, 'Invalid user data: %s' % str(e))

        if error is not None:
            abort(400, 'Could not invite user: %s' % error)

        return datamodel.User.get(username=user.username), 200


def with_signature_token(func):
    '''
    A decorator for API endpoint implementations that validates signed URLs. Token to
    sign URLs can be retrieved via the ``/auth`` endpoint.
    '''
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
    '''
    Returns a predicate that determines if the logged in user has the authorization
    to access the given upload and calculation.
    '''
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
            if g.user.user_id == upload.user_id:
                return True

            if calc_id is not None:
                try:
                    calc = processing.Calc.get(calc_id)
                except KeyError:
                    return False
                return g.user.user_id in calc.metadata.get('shared_with', [])

            return False

        except KeyError as e:
            logger = utils.get_logger(__name__, upload_id=upload_id, calc_id=calc_id)
            logger.error('Upload files without respective db entry')
            raise e

    return func
