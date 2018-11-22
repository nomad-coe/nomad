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
from sqlalchemy import create_engine
from sqlalchemy.orm import Session
from elasticsearch.exceptions import RequestError
from elasticsearch_dsl import connections
from mongoengine import connect

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
    mongo_client = connect(db=config.mongo.users_db, host=config.mongo.host, port=config.mongo.port)
    logger.info('setup mongo connection')


def setup_elastic():
    """ Creates connection to elastic search. """
    global elastic_client
    elastic_client = connections.create_connection(hosts=[config.elastic.host])
    logger.info('setup elastic connection')

    try:
        from nomad.repo import RepoCalc
        RepoCalc.init()
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
    # ensure that the schema exists
    with repository_db_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(
                "select exists(select * from information_schema.tables "
                "where table_name='users')")
            exists = cur.fetchone()[0]

    if not exists:
        reset_repository_db()

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
    """ Resets the databases mongo, elastic/calcs, and repository db. Be careful. """
    logger.info('reset mongodb')
    mongo_client.drop_database(config.mongo.users_db)

    logger.info('reset elastic search')
    elastic_client.indices.delete(index=config.elastic.calc_index)
    from nomad.repo import RepoCalc
    RepoCalc.init()

    logger.info('reset repository db')
    reset_repository_db()

    logger.info('reset files')
    shutil.rmtree(config.fs.objects, ignore_errors=True)
    shutil.rmtree(config.fs.tmp, ignore_errors=True)


@contextmanager
def repository_db_connection():
    """ Contextmanager for a psycopg2 session for the NOMAD-coe repository postgresdb """
    conn_str = "host='%s' port=%d dbname='%s' user='%s' password='%s'" % (
        config.repository_db.host,
        config.repository_db.port,
        config.repository_db.dbname,
        config.repository_db.user,
        config.repository_db.password)

    conn = psycopg2.connect(conn_str)
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
    """ Drops the existing NOMAD-coe repository postgres db and creates a new minimal one. """
    with repository_db_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(
                "DROP SCHEMA public CASCADE;"
                "CREATE SCHEMA public;"
                "GRANT ALL ON SCHEMA public TO postgres;"
                "GRANT ALL ON SCHEMA public TO public;")
            sql_file = os.path.join(os.path.dirname(__file__), 'empty_repository_db.sql')
            cur.execute(open(sql_file, 'r').read())


if __name__ == '__main__':
    reset_repository_db()
