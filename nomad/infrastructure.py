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
infrastructure services.
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

from nomad import config, utils

logger = utils.get_logger(__name__)

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
    setup_repository_db()


def setup_logging():
    utils.configure_logging()
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


def setup_elastic():
    """ Creates connection to elastic search. """
    global elastic_client
    elastic_client = connections.create_connection(
        hosts=['%s:%d' % (config.elastic.host, config.elastic.port)])
    logger.info('setup elastic connection')

    try:
        from nomad.repo import RepoCalc
        RepoCalc.init()
        from nomad.search import Entry
        Entry.init()
    except RequestError as e:
        if e.status_code == 400 and 'resource_already_exists_exception' in e.error:
            pass  # happens if two services try this at the same time
        else:
            raise e
    else:
        logger.info('init elastic index')


def setup_repository_db():
    """
    Makes sure that a minimal NOMAD-coe repository postgres db exists.
    Returns:
        An sqlalchemy session for the NOMAD-coe repository postgres db.
    """
    # ensure that the database exists
    exists = False
    try:
        with repository_db_connection():
            logger.info('repository db postgres database already exists')
            exists = True
    except psycopg2.OperationalError as e:
        if not ('database "%s" does not exist' % config.repository_db.dbname) in str(e):
            raise e
    if not exists:
        logger.info('repository db postgres database does not exist')
        try:
            with repository_db_connection(dbname='postgres', with_trans=False) as con:
                with con.cursor() as cursor:
                    cursor.execute("CREATE DATABASE %s  ;" % config.repository_db.dbname)
                logger.info('repository db postgres database created')
        except Exception as e:
            logger.info('could not create repository db postgres database', exc_info=e)
            raise e

    # ensure that the schema exists
    with repository_db_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(
                "select exists(select * from information_schema.tables "
                "where table_name='users')")
            exists = cur.fetchone()[0]
    if not exists:
        logger.info('repository db postgres schema does not exists')
        reset_repository_db()
    else:
        logger.info('repository db postgres schema already exists')

    # set the admin user password
    with repository_db_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(
                "UPDATE public.users SET password='%s' WHERE user_id=1;" %
                bcrypt.encrypt(config.services.admin_password, ident='2y'))

    global repository_db
    global repository_db_conn

    url = 'postgresql://%s:%s@%s:%d/%s' % (
        config.repository_db.user,
        config.repository_db.password,
        config.repository_db.host,
        config.repository_db.port,
        config.repository_db.dbname)
    engine = create_engine(url, echo=False)

    repository_db_conn = engine.connect()
    repository_db = Session(bind=repository_db_conn, autocommit=True)
    logger.info('setup repository db')


def reset():
    """
    Resets the databases mongo, elastic/calcs, repository db and all files. Be careful.
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
        from nomad.repo import RepoCalc
        RepoCalc.init()
        logger.info('elastic index resetted')
    except Exception as e:
        logger.error('exception resetting elastic', exc_info=e)

    try:
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
        return

    conn.commit()
    conn.close()


def reset_repository_db():
    """ Drops the existing NOMAD-coe repository postgres schema and creates a new minimal one. """
    old_repository_db = repository_db

    # invalidate and close all connections and sessions
    if repository_db is not None:
        repository_db.expunge_all()
        repository_db.invalidate()
    if repository_db_conn is not None:
        repository_db_conn.close()

    # perform the reset
    with repository_db_connection(with_trans=False) as conn:
        with conn.cursor() as cur:
            cur.execute("DROP SCHEMA IF EXISTS public CASCADE;")

            cur.execute(
                "CREATE SCHEMA public;"
                "GRANT ALL ON SCHEMA public TO postgres;"
                "GRANT ALL ON SCHEMA public TO public;")
            sql_file = os.path.join(os.path.dirname(__file__), 'empty_repository_db.sql')
            cur.execute(open(sql_file, 'r').read())
            logger.info('(re-)created repository db postgres schema')

    # try tp repair existing db connections
    if old_repository_db is not None:
        setup_repository_db()
        old_repository_db.bind = repository_db_conn
