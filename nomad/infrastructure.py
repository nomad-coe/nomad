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
exist to facilitate testing, :py:mod:`nomad.migration`, aspects of :py:mod:`nomad.client`, etc.
"""

import os.path
import shutil
from contextlib import contextmanager
import psycopg2
import psycopg2.extensions
from sqlalchemy import create_engine
from sqlalchemy.orm import Session
from elasticsearch.exceptions import RequestError
from elasticsearch_dsl import connections
from mongoengine import connect
from passlib.hash import bcrypt
import smtplib
from email.mime.text import MIMEText

from nomad import config, utils

logger = None

elastic_client = None
""" The elastic search client. """

mongo_client = None
""" The pymongo mongodb client. """

repository_db = None
""" The repository postgres db sqlalchemy session. """
repository_db_conn = None
""" The repository postgres db sqlalchemy connection. """


def setup():
    """
    Uses the current configuration (nomad/config.py and environemnt) to setup all the
    infrastructure services (repository db, mongo, elastic search) and logging.
    Will create client instances for the databases and has to be called before they
    can be used.
    """
    global elastic_client
    setup_logging()
    setup_mongo()
    setup_elastic()
    setup_repository_db(readonly=False)


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
        hosts=['%s:%d' % (config.elastic.host, config.elastic.port)])
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


def setup_repository_db(**kwargs):
    """ Creates a connection and stores it in the module variables. """
    repo_args = dict(readonly=False)
    repo_args.update(kwargs)
    connection, db = sqlalchemy_repository_db(**kwargs)

    global repository_db
    global repository_db_conn

    repository_db_conn, repository_db = connection, db
    logger.info('setup repository db connection')

    return repository_db_conn, repository_db


def sqlalchemy_repository_db(exists: bool = False, readonly: bool = True, **kwargs):
    """
    Returns SQLAlchemy connection and session for the given db parameters.

    Arguments:
        exists: Set to False to check and ensure db and schema existence
        readonly: Set to False for a write enabled connection
        **kwargs: Overwrite `config.repository_db` parameters
    """
    dbname = kwargs.get('dbname', config.repository_db.dbname)
    db_exists = exists
    if not db_exists:
        try:
            with repository_db_connection(dbname=dbname):
                logger.info('repository db postgres database already exists')
                db_exists = True
        except psycopg2.OperationalError as e:
            if not ('database "%s" does not exist' % dbname) in str(e):
                raise e

    if not db_exists:
        logger.info('repository db postgres database does not exist')
        try:
            with repository_db_connection(dbname='postgres', with_trans=False) as con:
                with con.cursor() as cursor:
                    cursor.execute("CREATE DATABASE %s  ;" % dbname)
                logger.info('repository db postgres database created')
        except Exception as e:
            logger.info('could not create repository db postgres database', exc_info=e)
            raise e

    # ensure that the schema exists
    schema_exists = exists
    if not schema_exists:
        with repository_db_connection(dbname=dbname) as conn:
            with conn.cursor() as cur:
                cur.execute(
                    "select exists(select * from information_schema.tables "
                    "where table_name='users')")
                schema_exists = cur.fetchone()[0]
        if not schema_exists:
            logger.info('repository db postgres schema does not exists')
            reset_repository_db_schema(dbname=dbname)
        else:
            logger.info('repository db postgres schema already exists')

    # set the admin user password
    if not exists:
        with repository_db_connection(dbname=dbname) as conn:
            with conn.cursor() as cur:
                cur.execute(
                    "UPDATE public.users SET password='%s' WHERE user_id=0;" %
                    bcrypt.encrypt(config.services.admin_password, ident='2y'))

    def no_flush():
        pass

    params = config.repository_db._asdict()
    params.update(**kwargs)
    url = 'postgresql://%s:%s@%s:%d/%s' % utils.to_tuple(params, 'user', 'password', 'host', 'port', 'dbname')
    engine = create_engine(url, echo=False)

    repository_db_conn = engine.connect()
    repository_db = Session(bind=repository_db_conn, autocommit=True)
    if readonly:
        repository_db.flush = no_flush

    return repository_db_conn, repository_db


def set_pid_prefix(prefix=7000000, target_db=None):
    if target_db is None:
        target_db = repository_db

    target_db.begin()
    target_db.execute('ALTER SEQUENCE calculations_calc_id_seq RESTART WITH %d' % prefix)
    target_db.commit()
    logger.info('set pid prefix', pid_prefix=prefix)


def reset(repo_content_only: bool = False):
    """
    Resets the databases mongo, elastic/calcs, repository db and all files. Be careful.
    In contrast to :func:`remove`, it will only remove the contents of dbs and indicies.
    This function just attempts to remove everything, there is no exception handling
    or any warranty it will succeed.

    Arguments:
        repo_content_only: True will only remove the calc/upload data from the repo db.
            But still reset all other dbs.
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
        if repo_content_only:
            reset_repository_db_content()
        else:
            reset_repository_db()
        logger.info('repository db resetted')
    except Exception as e:
        logger.error('exception resetting repository db', exc_info=e)

    logger.info('reset files')
    try:
        shutil.rmtree(config.fs.objects, ignore_errors=True)
        shutil.rmtree(config.fs.tmp, ignore_errors=True)
    except Exception as e:
        logger.error('exception deleting files', exc_info=e)


def remove():
    """
    Removes the databases mongo, elastic, repository db, and all files. Be careful.
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

    try:
        if repository_db is not None:
            repository_db.expunge_all()
            repository_db.invalidate()
        if repository_db_conn is not None:
            repository_db_conn.close()
        with repository_db_connection(dbname='postgres', with_trans=False) as con:
            with con.cursor() as cur:
                cur.execute('DROP DATABASE IF EXISTS %s' % config.repository_db.dbname)
        logger.info('repository db deleted')
    except Exception as e:
        logger.error('exception deleting repository db', exc_info=e)

    logger.info('reset files')
    try:
        shutil.rmtree(config.fs.objects, ignore_errors=True)
        shutil.rmtree(config.fs.tmp, ignore_errors=True)
    except Exception as e:
        logger.error('exception deleting files', exc_info=e)


@contextmanager
def repository_db_connection(dbname=None, with_trans=True):
    """ Contextmanager for a psycopg2 session for the NOMAD-coe repository postgresdb """
    repository_db_dict = config.repository_db._asdict()
    if dbname is not None:
        repository_db_dict.update(dbname=dbname)
    conn_str = "host='%s' port=%d dbname='%s' user='%s' password='%s'" % (
        repository_db_dict['host'],
        repository_db_dict['port'],
        repository_db_dict['dbname'],
        repository_db_dict['user'],
        repository_db_dict['password'])

    conn = psycopg2.connect(conn_str)
    if not with_trans:
        conn.set_isolation_level(psycopg2.extensions.ISOLATION_LEVEL_AUTOCOMMIT)
    try:
        yield conn
    except Exception as e:
        logger.error('Unhandled exception within repository db connection.', exc_info=e)
        conn.rollback()
        conn.close()
        raise e

    conn.commit()
    conn.close()


def reset_repository_db():
    """ Drops the existing NOMAD-coe repository postgres schema and creates a new minimal one. """
    global repository_db
    global repository_db_conn

    # invalidate and close all connections and sessions
    if repository_db is not None:
        repository_db.expunge_all()
        repository_db.invalidate()
        repository_db.close_all()
    if repository_db_conn is not None:
        repository_db_conn.close()
        repository_db_conn.engine.dispose()

    # perform the reset
    reset_repository_db_schema()

    # try tp repair existing db connections
    if repository_db is not None:
        new_connection, repository_db = setup_repository_db(exists=False)
        repository_db.bind = new_connection
        repository_db_conn = new_connection


def reset_repository_db_schema(**kwargs):
    with repository_db_connection(with_trans=False, **kwargs) as conn:
        with conn.cursor() as cur:
            cur.execute("DROP SCHEMA IF EXISTS public CASCADE;")

            cur.execute(
                "CREATE SCHEMA public;"
                "GRANT ALL ON SCHEMA public TO postgres;"
                "GRANT ALL ON SCHEMA public TO public;")
            sql_file = os.path.join(os.path.dirname(__file__), 'empty_repository_db.sql')
            cur.execute(open(sql_file, 'r').read())
    logger.info('(re-)created repository db postgres schema')


def reset_repository_db_content():
    tables = [
        'calcsets',
        'calculations',
        'citations',
        'coauthorships',
        'codefamilies',
        'codeversions',
        'doi_mapping',
        'metadata',
        'metadata_citations',
        'ownerships',
        'shareships',
        'spacegroups',
        'struct_ratios',
        'tags',
        'topics',
        'uploads'
    ]
    with repository_db_connection(with_trans=False) as conn:
        with conn.cursor() as cur:
            cur.execute('TRUNCATE %s CASCADE;' % ', '.join(tables))

    logger.info('removed repository db content')


def send_mail(name: str, email: str, message: str, subject: str):
    if config.mail.host is None or config.mail.host.strip() == '':
        return

    logger = utils.get_logger(__name__)
    server = smtplib.SMTP(config.mail.host, config.mail.port)

    if config.mail.port == 995:
        try:
            server.starttls()
        except Exception as e:
            logger.warning('Could use TTS', exc_info=e)

    if config.mail.user is not None:
        try:
            server.login("youremailusername", "password")
        except Exception as e:
            logger.warning('Could not log into mail server', exc_info=e)

    msg = MIMEText(message)
    msg['Subject'] = subject
    msg['From'] = 'nomad@fairdi webmaster'
    msg['To'] = name

    try:
        server.send_message(msg, config.mail.from_address, email)
    except Exception as e:
        logger.error('Could send email', exc_info=e)

    server.quit()
