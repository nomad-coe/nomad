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
This module provides function to establish connections to the database, searchengine, etc.
infrastructure services. Usually everything is setup at once with :func:`setup`. This
is run once for each *api* and *worker* process. Individual functions for partial setups
exist to facilitate testing, :py:mod:`nomad.migration`, aspects of :py:mod:`nomad.cli`, etc.
"""

from typing import Union
import os.path
import shutil
from elasticsearch.exceptions import RequestError
from elasticsearch_dsl import connections
from mongoengine import connect
import smtplib
from email.mime.text import MIMEText
from keycloak import KeycloakOpenID, KeycloakAdmin
from flask_oidc import OpenIDConnect
import json
from flask import g, request
import basicauth

from nomad import config, utils

logger = None

elastic_client = None
""" The elastic search client. """

mongo_client = None
""" The pymongo mongodb client. """


def setup():
    """
    Uses the current configuration (nomad/config.py and environment) to setup all the
    infrastructure services (repository db, mongo, elastic search) and logging.
    Will create client instances for the databases and has to be called before they
    can be used.
    """
    global elastic_client
    setup_logging()
    setup_mongo()
    setup_elastic()


def setup_logging():
    utils.configure_logging()

    global logger
    logger = utils.get_logger(__name__)

    logger.info(
        'setup logging',
        logstash=config.logstash.enabled,
        logstash_host=config.logstash.host,
        logstash_port=config.logstash.tcp_port,
        logstash_level=config.logstash.level)


def setup_mongo():
    """ Creates connection to mongodb. """
    global mongo_client
    mongo_client = connect(db=config.mongo.db_name, host=config.mongo.host, port=config.mongo.port)
    logger.info('setup mongo connection')
    return mongo_client


def setup_elastic():
    """ Creates connection to elastic search. """
    global elastic_client
    elastic_client = connections.create_connection(
        hosts=['%s:%d' % (config.elastic.host, config.elastic.port)],
        timeout=60, max_retries=10, retry_on_timeout=True)
    logger.info('setup elastic connection')

    try:
        from nomad.search import Entry
        Entry.init(index=config.elastic.index_name)
        Entry._index._name = config.elastic.index_name
        logger.info('initialized elastic index', index_name=config.elastic.index_name)
    except RequestError as e:
        if e.status_code == 400 and 'resource_already_exists_exception' in e.error:
            # happens if two services try this at the same time
            pass
        else:
            raise e

    return elastic_client


class Keycloak():
    """
    A class that encapsulates all keycloak related functions for easier mocking and
    configuration
    """
    def __init__(self):
        self._flask_oidc = None
        self.__oidc_client = None
        self.__admin_client = None

    def configure_flask(self, app):
        oidc_issuer_url = '%s/realms/%s' % (config.keycloak.server_url.rstrip('/'), config.keycloak.realm_name)
        oidc_client_secrets = dict(
            client_id=config.keycloak.client_id,
            client_secret=config.keycloak.client_secret_key,
            issuer=oidc_issuer_url,
            auth_uri='%s/protocol/openid-connect/auth' % oidc_issuer_url,
            token_uri='%s/protocol/openid-connect/token' % oidc_issuer_url,
            userinfo_uri='%s/protocol/openid-connect/userinfo' % oidc_issuer_url,
            token_introspection_uri='%s/protocol/openid-connect/token/introspect' % oidc_issuer_url,
            redirect_uris=['http://localhost/fairdi/nomad/latest'])
        oidc_client_secrets_file = os.path.join(config.fs.tmp, 'oidc_client_secrets')
        with open(oidc_client_secrets_file, 'wt') as f:
            json.dump(dict(web=oidc_client_secrets), f)
        app.config.update(dict(
            SECRET_KEY=config.services.api_secret,
            OIDC_CLIENT_SECRETS=oidc_client_secrets_file,
            OIDC_OPENID_REALM=config.keycloak.realm_name))

        self._flask_oidc = OpenIDConnect(app)

    @property
    def _oidc_client(self):
        if self.__oidc_client is None:
            self.__oidc_client = KeycloakOpenID(
                server_url=config.keycloak.server_url,
                client_id=config.keycloak.client_id,
                realm_name=config.keycloak.realm_name,
                client_secret_key=config.keycloak.client_secret_key)

        return self.__oidc_client

    def authorize_flask(self, token_only: bool = True) -> Union[str, object]:
        token = None
        if 'Authorization' in request.headers and request.headers['Authorization'].startswith('Bearer '):
            token = request.headers['Authorization'].split(None, 1)[1].strip()
        elif 'access_token' in request.form:
            token = request.form['access_token']
        elif 'access_token' in request.args:
            token = request.args['access_token']
        elif 'Authorization' in request.headers and request.headers['Authorization'].startswith('Basic '):
            if token_only:
                return 'Basic authentication not allowed, use Bearer token instead'

            try:
                username, password = basicauth.decode(request.headers['Authorization'])
                token_info = self._oidc_client.token(username=username, password=password)
                token = token_info['access_token']
            except Exception as e:
                # TODO logging
                return 'Could not authenticate Basic auth: %s' % str(e)

        if token is not None:
            validity = self._flask_oidc.validate_token(token)

            if validity is not True:
                return validity

            else:
                g.oidc_id_token = g.oidc_token_info
                return self.get_user()

        else:
            return None

    def get_user(self, user_id: str = None, email: str = None) -> object:
        from nomad import datamodel

        if email is not None:
            try:
                user_id = self._admin_client.get_user_id(email)
            except Exception:
                raise KeyError('User does not exist')

        if user_id is None and g.oidc_id_token is not None and self._flask_oidc is not None:
            try:
                return datamodel.User(token=g.oidc_id_token, **self._flask_oidc.user_getinfo([
                    'email', 'firstName', 'lastName', 'username', 'createdTimestamp']))
            except Exception as e:
                # TODO logging
                raise e

        assert user_id is not None, 'Could not determine user from given kwargs'

        try:
            keycloak_user = self._admin_client.get_user(user_id)
        except Exception:
            raise KeyError('User does not exist')

        return datamodel.User(**keycloak_user)

    @property
    def _admin_client(self):
        if self.__admin_client is None:
            self.__admin_client = KeycloakAdmin(
                server_url=config.keycloak.server_url,
                username=config.keycloak.username,
                password=config.keycloak.password,
                realm_name='master',
                verify=True)
            self.__admin_client.realm_name = config.keycloak.realm_name

        return self.__admin_client


keycloak = Keycloak()


def reset():
    """
    Resets the databases mongo, elastic/calcs, and all files. Be careful.
    In contrast to :func:`remove`, it will only remove the contents of dbs and indicies.
    This function just attempts to remove everything, there is no exception handling
    or any warranty it will succeed.
    """
    try:
        if not mongo_client:
            setup_mongo()
        mongo_client.drop_database(config.mongo.db_name)
        logger.info('mongodb resetted')
    except Exception as e:
        logger.error('exception reset mongodb', exc_info=e)

    try:
        if not elastic_client:
            setup_elastic()
        elastic_client.indices.delete(index=config.elastic.index_name)
        from nomad.search import Entry
        Entry.init(index=config.elastic.index_name)
        logger.info('elastic index resetted')
    except Exception as e:
        logger.error('exception resetting elastic', exc_info=e)

    try:
        shutil.rmtree(config.fs.staging, ignore_errors=True)
        shutil.rmtree(config.fs.public, ignore_errors=True)
        # delete tmp without the folder
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


def remove():
    """
    Removes the databases mongo, elastic, and all files. Be careful.
    This function just attempts to remove everything, there is no exception handling
    or any warranty it will succeed.
    """
    try:
        if not mongo_client:
            setup_mongo()
        mongo_client.drop_database(config.mongo.db_name)
        logger.info('mongodb deleted')
    except Exception as e:
        logger.error('exception deleting mongodb', exc_info=e)

    try:
        if not elastic_client:
            setup_elastic()
        elastic_client.indices.delete(index=config.elastic.index_name)
        logger.info('elastic index')
    except Exception as e:
        logger.error('exception deleting elastic', exc_info=e)

    logger.info('reset files')
    try:
        shutil.rmtree(config.fs.staging, ignore_errors=True)
        shutil.rmtree(config.fs.public, ignore_errors=True)
        shutil.rmtree(config.fs.tmp, ignore_errors=True)
    except Exception as e:
        logger.error('exception deleting files', exc_info=e)


def send_mail(name: str, email: str, message: str, subject: str):
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
    msg['From'] = 'The NOMAD team <%s>' % config.mail.from_address
    msg['To'] = name
    to_addrs = [email]

    if config.mail.cc_address is not None:
        msg['Cc'] = 'The NOMAD team <%s>' % config.mail.cc_address
        to_addrs.append(config.mail.cc_address)

    try:
        server.send_message(msg, from_addr=config.mail.from_address, to_addrs=to_addrs)
    except Exception as e:
        logger.error('Could not send email', exc_info=e)

    server.quit()
