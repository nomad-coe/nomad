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

from mongoengine import connect
from elasticsearch_dsl import connections
from elasticsearch.exceptions import RequestError

from nomad import config, utils

logger = utils.get_logger(__name__)

elastic_client = None
""" The elastic search client. """


def setup():
    """ Creates connections to mongodb and elastic search. """
    global elastic_client
    setup_mongo()
    setup_elastic()

    from nomad import user
    user.ensure_test_users()


def setup_mongo():
    """ Creates connection to mongodb. """
    connect(db=config.mongo.users_db, host=config.mongo.host, port=config.mongo.port)
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
