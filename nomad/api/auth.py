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

from typing import Tuple
from flask import g, request
from flask_restplus import abort, Resource, fields
from datetime import datetime
import functools
import basicauth

from nomad import config, processing, files, utils, coe_repo, infrastructure
from nomad.coe_repo import LoginException

from .app import api, RFC3339DateTime, oidc


class User:
    """
    A data class that holds all information for a single user. This can be the logged in
    and authenticated user, or other users (i.e. co-authors, etc.).
    """
    def __init__(
            self, email, name=None, first_name='', last_name='', affiliation=None,
            created: datetime = None, **kwargs):
        assert email is not None, 'Users must have an email, it is used as unique id'

        self.email = email

        first_name = kwargs.get('firstName', first_name)
        last_name = kwargs.get('lastName', last_name)
        name = kwargs.get('username', name)
        created_timestamp = kwargs.get('createdTimestamp', None)

        if len(last_name) > 0 and len(first_name) > 0:
            name = '%s, %s' % (last_name, first_name)
        elif len(last_name) != 0:
            name = last_name
        elif len(first_name) != 0:
            name = first_name
        elif name is None:
            name = 'unnamed user'

        self.name = name

        if created is not None:
            self.created = None
        elif created_timestamp is not None:
            self.created = datetime.fromtimestamp(created_timestamp)
        else:
            self.created = None

        # TODO affliation


def _validate_token(require_token: bool = True, **kwargs) -> Tuple[bool, str]:
    """
    Uses OIDC to check if the request carries token based authentication and if
    this authentication is valid.

    Returns: A tuple with bool and potential error message
    """
    token = None
    if 'Authorization' in request.headers and request.headers['Authorization'].startswith('Bearer '):
        token = request.headers['Authorization'].split(None, 1)[1].strip()
    if 'access_token' in request.form:
        token = request.form['access_token']
    elif 'access_token' in request.args:
        token = request.args['access_token']

    validity = oidc.validate_token(token, **kwargs)

    if validity:
        g.oidc_id_token = g.oidc_token_info

    return (validity is True) or (not require_token), validity


def _get_user():
    """
    Retrieves OIDC user info and populate the global flask ``g.user`` variable.
    """
    if g.oidc_id_token:
        try:
            g.user = User(**oidc.user_getinfo([
                'email', 'firstName', 'lastName', 'username', 'createdTimestamp']))
        except Exception as e:
            ## TODO logging
            raise e
    else:
        g.user = None


def login_if_available(func):
    """
    A decorator for API endpoint implementations that might authenticate users, but
    provide limited functionality even without users.
    """
    @functools.wraps(func)
    @api.response(401, 'Not authorized, some data require authentication and authorization')
    @api.doc(security=list('OpenIDConnect Bearer Token'))
    def wrapper(*args, **kwargs):
        valid, msg = _validate_token(require_token=False)
        if valid:
            _get_user()
            return func(*args, **kwargs)
        else:
            abort(401, message=msg)

    return wrapper


def login_really_required(func):
    """
    A decorator for API endpoint implementations that forces user authentication on
    endpoints.
    """
    @functools.wraps(func)
    @api.response(401, 'Not authorized, this endpoint required authorization')
    @api.doc(security=list('OpenIDConnect Bearer Token'))
    def wrapper(*args, **kwargs):
        valid, msg = _validate_token(require_token=True)
        if valid:
            _get_user()
            return func(*args, **kwargs)
        else:
            abort(401, message=msg)

    return wrapper


def admin_login_required(func):
    """
    A decorator for API endpoint implementations that should only work for the admin user.
    """
    @functools.wraps(func)
    @api.response(401, 'Authentication required or not authorized as admin user. Only admin can access this endpoint.')
    @api.doc(security=list('OpenIDConnect Bearer Token'))
    @oidc.accept_token(require_token=True)
    def wrapper(*args, **kwargs):
        if oidc.user_getfield('email') == config.keycloak.adminEmail:
            return func(*args, **kwargs)
        else:
            abort(401, message='Only the admin user can perform reset.')

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
    @login_if_available
    def get(self):
        if g.user is not None:
            return g.user

        if 'Authorization' in request.headers and request.headers['Authorization'].startswith('Basic '):
            try:
                username, password = basicauth.decode(request.headers['Authorization'])
                token = infrastructure.keycloak_oidc_client.token(username=username, password=password)
                validity = oidc.validate_token(token['access_token'])
            except Exception as e:
                # TODO logging
                abort(401, message='Could not authenticate Basic auth: %s' % str(e))

            if validity is not True:
                abort(401, message=validity)
            else:
                g.oidc_id_token = g.oidc_token_info
                _get_user()
        else:
            abort(401, message='Authentication credentials found in your request')

        if g.user is None:
            abort(401, message='User not authenticated')

        return g.user


@ns.route('/user')
class UserResource(Resource):
    @api.doc('create_user')
    @api.expect(user_model)
    @api.response(400, 'Invalid user data')
    @api.marshal_with(user_model, skip_none=True, code=200, description='User created')
    @admin_login_required
    def put(self):
        """
        Creates a new user account. Currently only the admin user is allows. The
        NOMAD-CoE repository GUI should be used to create user accounts for now.
        Passwords have to be encrypted by the client with bcrypt and 2y indent.
        """
        data = request.get_json()
        if data is None:
            data = {}

        for required_key in ['last_name', 'first_name', 'password', 'email']:
            if required_key not in data:
                abort(400, message='The %s is missing' % required_key)

        if 'user_id' in data:
            if coe_repo.User.from_user_id(data['user_id']) is not None:
                abort(400, 'User with given user_id %d already exists.' % data['user_id'])

        user = coe_repo.User.create_user(
            email=data['email'], password=data.get('password', None), crypted=True,
            first_name=data['first_name'], last_name=data['last_name'],
            created=data.get('created', datetime.utcnow()),
            affiliation=data.get('affiliation', None), token=data.get('token', None),
            user_id=data.get('user_id', None))

        return user, 200


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
        elif g.user.user_id == 0:
            # the admin user does have authorization to access everything
            return True

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
