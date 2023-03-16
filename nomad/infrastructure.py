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
This module provides function to establish connections to the database, searchengine, etc.
infrastructure services. Usually everything is setup at once with :func:`setup`. This
is run once for each *api* and *worker* process. Individual functions for partial setups
exist to facilitate testing, aspects of :py:mod:`nomad.cli`, etc.
'''

import os.path
import os
import shutil
from elasticsearch_dsl import connections
from mongoengine import connect, disconnect
from mongoengine.connection import ConnectionFailure
import smtplib
from email.mime.text import MIMEText
from keycloak import KeycloakOpenID, KeycloakAdmin
from keycloak.exceptions import KeycloakAuthenticationError, KeycloakGetError
import json
import jwt
from datetime import datetime
import re
import unidecode

from nomad import config, utils
from nomad.utils.structlogging import get_logger

# The metainfo is defined and used during imports. This is problematic.
# We import all parsers very early in the infrastructure setup. This will populate
# the metainfo with parser specific definitions, before the metainfo might be used.
from nomad.parsing import parsers  # pylint: disable=unused-import

# TODO put somemore thought into warnings
import warnings
warnings.filterwarnings("ignore")

logger = get_logger(__name__)

elastic_client = None
''' The elastic search client. '''

mongo_client = None
''' The pymongo mongodb client. '''


def setup():
    '''
    Uses the current configuration (nomad/config.py and environment) to setup all the
    infrastructure services (repository db, mongo, elastic search) and logging.
    Will create client instances for the databases and has to be called before they
    can be used.
    '''
    setup_files()
    setup_mongo()
    setup_elastic()


def setup_files():
    for directory in [config.fs.public, config.fs.staging, config.fs.tmp]:
        if not os.path.exists(directory):
            os.makedirs(directory)


def setup_mongo(client=False):
    ''' Creates connection to mongodb. '''
    global mongo_client
    try:
        mongo_client = connect(db=config.mongo.db_name, host=config.mongo.host, port=config.mongo.port)
    except ConnectionFailure:
        disconnect()
        mongo_client = connect(db=config.mongo.db_name, host=config.mongo.host, port=config.mongo.port)

    logger.info('setup mongo connection')
    return mongo_client


def setup_elastic():
    ''' Creates connection to elastic search. '''
    global elastic_client
    elastic_client = connections.create_connection(
        hosts=['%s:%d' % (config.elastic.host, config.elastic.port)],
        timeout=config.elastic.timeout, max_retries=10, retry_on_timeout=True)
    logger.info('setup elastic connection')
    from nomad.metainfo.elasticsearch_extension import create_indices as create_v1_indices
    create_v1_indices()
    logger.info('initialized v1 elastic indices')

    return elastic_client


class KeycloakError(Exception): pass


class Keycloak():
    '''
    A class that encapsulates all keycloak related functions for easier mocking and
    configuration
    '''
    def __init__(self):
        self.__oidc_client = None
        self.__public_keys = None

    @property
    def _oidc_client(self):
        if self.__oidc_client is None:
            self.__oidc_client = KeycloakOpenID(
                server_url=config.keycloak.server_url,
                client_id=config.keycloak.client_id,
                realm_name=config.keycloak.realm_name,
                client_secret_key=config.keycloak.client_secret)

        return self.__oidc_client

    @property
    def _public_keys(self):
        if self.__public_keys is None:
            try:
                jwks = self._oidc_client.certs()
                self.__public_keys = {}
                for jwk in jwks['keys']:
                    kid = jwk['kid']
                    self.__public_keys[kid] = jwt.algorithms.RSAAlgorithm.from_jwk(
                        json.dumps(jwk))
            except Exception as e:
                self.__public_keys = None
                raise e

        return self.__public_keys

    def refresh_token(self, access_token: str, refresh_token: str, **kwargs) -> str:
        return self._oidc_client.refresh_token(refresh_token)

    def basicauth(self, username: str, password: str) -> str:
        '''
        Performs basic authentication and returns an access token.

        Raises:
            KeycloakError
        '''
        try:
            token_info = self._oidc_client.token(username=username, password=password)
        except KeycloakAuthenticationError as e:
            raise KeycloakError(e)
        except Exception as e:
            logger.error('cannot perform basicauth', exc_info=e)
            raise e

        return token_info['access_token']

    def decode_access_token(self, access_token: str) -> dict:
        try:
            kid = jwt.get_unverified_header(access_token)['kid']
            key = keycloak._public_keys.get(kid)
            if key is None:
                logger.error('The user provided keycloak public key does not exist. Does the UI use the right realm?')
                raise KeycloakError(utils.strip('''
                    Could not validate credentials.
                    The user provided keycloak public key does not exist.
                    Does the UI use the right realm?'''))

            issuer = f'{config.keycloak.public_server_url.rstrip("/")}/realms/{config.keycloak.realm_name}'
            options = dict(verify_aud=False, verify_exp=True, verify_iss=True)
            return jwt.decode(
                access_token, key=key, algorithms=['RS256'], options=options,
                issuer=issuer)
        except jwt.InvalidTokenError:
            raise KeycloakError('Could not validate credentials. The given token is invalid.')

    def tokenauth(self, access_token: str) -> object:
        '''
        Authenticates the given access_token

        Returns:
            The user
        '''
        try:
            payload = self.decode_access_token(access_token)

            user_id: str = payload.get('sub')
            if user_id is None:
                raise KeycloakError(utils.strip('''
                    Could not validate credentials.
                    The given token does not contain a user_id.'''))

            from nomad import datamodel
            return datamodel.User(
                user_id=user_id,
                username=payload.get('preferred_username', None),
                email=payload.get('email', None),
                first_name=payload.get('given_name', None),
                last_name=payload.get('family_name', None))

        except Exception as e:
            logger.error('cannot perform tokenauth', exc_info=e)
            raise e


keycloak = Keycloak()


class UserManagement():
    def add_user(self, user, bcrypt_password=None, invite=False):
        '''
        Adds the given :class:`nomad.datamodel.User` instance to the configured keycloak
        realm using the keycloak admin API.
        '''
        raise NotImplementedError()

    def search_user(self, query: str):
        raise NotImplementedError()

    def get_user(self, user_id: str = None, username: str = None, email: str = None):
        '''
        Retrives all available information about a user from the local keycloak admin
        interface or the central NOMAD installation. This can be used to retrieve
        complete user information, because the info solely gathered from tokens is generally
        incomplete.
        '''
        raise NotImplementedError()


class OasisUserManagement(UserManagement):
    def __init__(self, users_api_url: str = None):
        if users_api_url:
            self._users_api_url = users_api_url
        else:
            self._users_api_url = f'{config.oasis.central_nomad_deployment_url}/v1/users'

    def add_user(self, user, bcrypt_password=None, invite=False):
        raise NotImplementedError(
            'Adding a user is not possible for an Oasis using the central user management.')

    def __user_from_api_user(self, api_user):
        from nomad import datamodel
        del api_user['is_admin']
        del api_user['is_oasis_admin']
        return datamodel.User.m_from_dict(api_user)

    def search_user(self, query: str):
        import requests
        response = requests.get(self._users_api_url, params=dict(prefix=query))
        if response.status_code != 200:
            raise KeycloakError('Could not request central nomad\'s user management.')

        return list(self.__user_from_api_user(user) for user in response.json()['data'])

    def get_user(self, user_id: str = None, username: str = None, email: str = None):
        import requests

        kwargs = {}
        if user_id:
            kwargs['user_id'] = user_id
        elif username:
            kwargs['username'] = username
        elif email:
            kwargs['email'] = email
        else:
            return None

        response = requests.get(self._users_api_url, params=kwargs)
        if response.status_code != 200:
            raise KeycloakError('Could not request central nomad\'s user management.')

        data = response.json()
        if len(data['data']) == 0:
            return None

        return self.__user_from_api_user(data['data'][0])


class KeycloakUserManagement(UserManagement):
    def __init__(self):
        self.__admin_client = None

    def __create_username(self, user):
        if user.first_name is not None and user.last_name is not None:
            user.username = '%s%s' % (user.first_name[:1], user.last_name)
        elif user.last_name is not None:
            user.username = user.last_name
        elif '@' in user.username:
            user.username = user.username.split('@')[0]

        user.username = unidecode.unidecode(user.username.lower())
        user.username = re.sub(r'[^0-9a-zA-Z_\-\.]+', '', user.username)

        index = 1
        try:
            while self.get_user(username=user.username):
                user.username += '%d' % index
                index += 1
        except KeyError:
            pass

    def add_user(self, user, bcrypt_password=None, invite=False):
        from nomad import datamodel
        if not isinstance(user, datamodel.User):
            if 'user_id' not in user:
                user['user_id'] = 'not set'

            if 'password' in user:
                bcrypt_password = user.pop('password')

            created = user.get('created', None)
            if created is not None and not isinstance(created, datetime):
                user['created'] = datetime.fromtimestamp(created / 1000)

            user = datamodel.User(**user)

        if user.username is None or not re.match(r'^[a-zA-Z0-9_\-\.]+$', user.username):
            self.__create_username(user)

        keycloak_user = dict(
            id=user.user_id if user.user_id != 'not set' else None,
            email=user.email,
            username=user.username,
            firstName=user.first_name,
            lastName=user.last_name,
            attributes=dict(
                repo_user_id=user.repo_user_id,
                affiliation=user.affiliation if user.affiliation is not None else '',
                affiliation_address=user.affiliation_address if user.affiliation_address is not None else ''),
            createdTimestamp=user.created.timestamp() * 1000 if user.created is not None else None,
            enabled=True,
            emailVerified=True)

        if invite:
            keycloak_user['requiredActions'] = [
                'UPDATE_PASSWORD', 'UPDATE_PROFILE', 'VERIFY_EMAIL']

        if bcrypt_password is not None:
            keycloak_user['credentials'] = [dict(
                type='password',
                hashedSaltedValue=bcrypt_password,
                algorithm='bcrypt')]

        keycloak_user = {
            key: value for key, value in keycloak_user.items()
            if value is not None}

        if user.user_id != 'not_set':
            try:
                self._admin_client.get_user(user.user_id)
                return 'User %s with given id already exists' % user.email
            except KeycloakGetError:
                pass

        if self._admin_client.get_user_id(user.email) is not None:
            return 'User with email %s already exists' % user.email

        try:
            self._admin_client.create_user(keycloak_user)
        except KeycloakGetError as e:
            try:
                return json.loads(e.response_body)['errorMessage']
            except Exception:
                return str(e)
        except Exception as e:
            return str(e)

        if invite:
            try:
                user = self.get_user(username=user.username)
                self._admin_client.send_verify_email(user_id=user.user_id)
            except Exception as e:
                logger.error('could not send verify email', exc_info=e)

        return None

    def __user_from_keycloak_user(self, keycloak_user):
        from nomad import datamodel

        kwargs = {key: value[0] for key, value in keycloak_user.get('attributes', {}).items()}
        oasis_admin = kwargs.pop('is_oasis_admin', None) is not None
        return datamodel.User(
            m_ignore_additional_keys=True,
            user_id=keycloak_user['id'],
            email=keycloak_user.get('email'),
            username=keycloak_user.get('username'),
            first_name=keycloak_user.get('firstName'),
            last_name=keycloak_user.get('lastName'),
            is_oasis_admin=oasis_admin,
            created=datetime.fromtimestamp(keycloak_user['createdTimestamp'] / 1000),
            **kwargs)

    def search_user(self, query: str):
        kwargs = {}
        if query is not None:
            kwargs['query'] = dict(search=query, max=1000)
        else:
            kwargs['query'] = dict(max=1000)
        try:
            keycloak_results = self._admin_client.get_users(**kwargs)
        except Exception as e:
            logger.error('Could not retrieve users from keycloak', exc_info=e)
            raise e

        return [
            self.__user_from_keycloak_user(keycloak_user)
            for keycloak_user in keycloak_results]

    def get_user(self, user_id: str = None, username: str = None, email: str = None):
        if username is not None and user_id is None:
            with utils.lnr(logger, 'Could not use keycloak admin client'):
                user_id = self._admin_client.get_user_id(username)

            if user_id is None:
                raise KeyError('User with username %s does not exist' % username)

        if email is not None and user_id is None:
            with utils.lnr(logger, 'Could not use keycloak admin client'):
                users = self._admin_client.get_users(query=dict(email=email))

            if len(users) > 0:
                user_id = users[0]['id']

            if user_id is None:
                raise KeyError('User with email %s does not exist' % email)

        assert user_id is not None, 'Could not determine user from given kwargs'

        try:
            keycloak_user = self._admin_client.get_user(user_id)

        except Exception as e:
            if str(getattr(e, 'response_code', 404)) == '404':
                raise KeyError('User does not exist')

            logger.error('Could not retrieve user from keycloak', exc_info=e)
            raise e

        return self.__user_from_keycloak_user(keycloak_user)

    @property
    def _admin_client(self):
        if True:  # TODO (self.__admin_client is None:), client becomes unusable after 60s
            self.__admin_client = KeycloakAdmin(
                server_url=config.keycloak.server_url,
                username=config.keycloak.username,
                password=config.keycloak.password,
                realm_name=config.keycloak.realm_name,
                verify=True)
            self.__admin_client.realm_name = config.keycloak.realm_name

        return self.__admin_client


user_management: UserManagement
if config.oasis.uses_central_user_management:
    user_management = OasisUserManagement()
else:
    user_management = KeycloakUserManagement()


def reset(remove: bool):
    '''
    Resets the databases mongo, elastic/entries, and all files. Be careful.
    In contrast to :func:`remove`, it will only remove the contents of dbs and indicies.
    This function just attempts to remove everything, there is no exception handling
    or any warranty it will succeed.

    Args:
        remove: Do not try to recreate empty databases, remove entirely.
    '''
    try:
        if not mongo_client:
            setup_mongo()
        mongo_client.drop_database(config.mongo.db_name)
        logger.info('mongodb resetted')
    except Exception as e:
        logger.error('exception reset mongodb', exc_info=e)

    try:
        from nomad.metainfo.elasticsearch_extension import create_indices, delete_indices

        if not elastic_client:
            setup_elastic()

        delete_indices()

        if not remove:
            create_indices()

        logger.info('elastic index resetted')
    except Exception as e:
        logger.error('exception resetting elastic', exc_info=e)

    try:
        shutil.rmtree(config.fs.staging, ignore_errors=True)
        shutil.rmtree(config.fs.public, ignore_errors=True)

        # delete tmp without the folder
        if os.path.isdir(config.fs.tmp):
            for sub_path in os.listdir(config.fs.tmp):
                path = os.path.join(config.fs.tmp, sub_path)
                try:
                    if os.path.isfile(path):
                        os.unlink(path)
                    elif os.path.isdir(path): shutil.rmtree(path, ignore_errors=True)
                except Exception:
                    pass

        logger.info('files resetted')
    except Exception as e:
        logger.error('exception deleting files', exc_info=e)


def send_mail(name: str, email: str, message: str, subject: str):
    """Used to programmatically send mails.

    Args:
        name: The email recipient name.
        email: The email recipient address.
        messsage: The email body.
        subject: The subject line.
    """
    if not config.mail.enabled:
        return

    logger = utils.get_logger(__name__)
    server = smtplib.SMTP(config.mail.host, config.mail.port)

    if config.mail.port == 995:
        try:
            server.starttls()
        except Exception as e:
            logger.warning('Could not use TTS', exc_info=e)

    if config.mail.with_login:
        try:
            server.login(config.mail.user, config.mail.password)
        except Exception as e:
            logger.warning('Could not log into mail server', exc_info=e)

    msg = MIMEText(message)
    msg['Subject'] = subject
    msg['To'] = name
    msg['From'] = config.mail.from_address
    to_addrs = [email]

    if config.mail.cc_address is not None:
        msg['Cc'] = 'The nomad team <%s>' % config.mail.cc_address
        to_addrs.append(config.mail.cc_address)

    try:
        server.send_message(msg, from_addr=config.mail.from_address, to_addrs=to_addrs)
    except Exception as e:
        logger.error('Could not send email', exc_info=e)

    server.quit()
