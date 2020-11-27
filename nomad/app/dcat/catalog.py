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

from flask_restplus import Resource, reqparse, fields
from flask import Response
from elasticsearch_dsl import Q

from nomad import search

from .api import api
from .mapping import Mapping

ns = api.namespace('catalog', description='The API for DCAT catalog.')

iso8601 = fields.DateTime(dt_format='iso8601')

arg_parser = reqparse.RequestParser()
arg_parser.add_argument('format', type=str, choices=[
    'xml',
    'n3',
    'turtle',
    'nt',
    'pretty-xml',
    'trig'])
arg_parser.add_argument('after', type=str)
arg_parser.add_argument(
    'modified_since', type=lambda x: iso8601.parse(x),
    help='A yyyy-MM-ddTHH:mm:ss (RFC3339) maximum entry time (e.g. upload time)')


@ns.route('/')
class Catalog(Resource):
    @api.doc('get_dcat_datasets')
    @api.expect(arg_parser)
    @api.produces(['application/xml'])
    @api.response(404, 'There is no entry with the given id.')
    @api.response(401, 'This entry is not publically accessible.')
    @api.response(200, 'Data send', headers={'Content-Type': 'application/xml'})
    def get(self):
        ''' Returns a page of DCAT datasets. '''

        args = arg_parser.parse_args()
        format_ = args.get('format')
        if format_ is None:
            format_ = 'xml'

        modified_since = args.get('modified_since', None)
        modified_since = iso8601.parse(modified_since) if modified_since is not None else None
        after = args.get('after', None)

        search_request = search.SearchRequest().owner('public')
        if modified_since is not None:
            modified_clause = Q('range', upload_time=dict(gte=modified_since))
            modified_clause |= Q('range', last_edit=dict(gte=modified_since))
            modified_clause |= Q('range', last_processing=dict(gte=modified_since))
            search_request.q &= modified_clause

        es_search = search_request._search.query(search_request.q)
        if after is not None:
            es_search = es_search.extra(search_after=[after], sort='calc_id')
        es_response = es_search.execute()

        mapping = Mapping()
        mapping.map_catalog(es_response.hits)

        return Response(
            mapping.g.serialize(format=format_).decode('utf-8'), 200,
            {'Content-Type': 'application/xml'})
