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

from flask_restplus import abort
from flask import request
from elasticsearch_dsl import Q
from cachetools import cached, TTLCache

from nomad import search

from .api import api


ns = api.namespace('v1', description='The version v1 API namespace with all OPTiMaDe endpoints.')


# TODO replace with decorator that filters response_fields
def base_request_args():
    if request.args.get('response_format', 'json') != 'json':
        abort(400, 'Response format is not supported.')

    properties_str = request.args.get('response_fields', None)
    if properties_str is not None:
        return properties_str.split(',')
    return None


def base_search_request():
    ''' Creates a search request for all public and optimade enabled data. '''
    return search.SearchRequest().owner('public', None).search_parameter('processed', True).query(
        Q('exists', field='dft.optimade.elements'))  # TODO use the elastic annotations when done


@cached(TTLCache(maxsize=1, ttl=60 * 60))
def nentries():
    ''' Gives the overall number of public calculations. '''
    return search.SearchRequest().owner(owner_type='public').execute()['total']
