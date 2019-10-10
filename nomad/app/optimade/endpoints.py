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

from flask_restplus import Resource, abort
from flask import request

from nomad import search
from nomad.metainfo.optimade import OptimadeStructureEntry

from .api import api, url
from .models import json_api_single_response_model, entry_listing_endpoint_parser, Meta, \
    Links, CalculationDataObject, single_entry_endpoint_parser, base_endpoint_parser
from .filterparser import parse_filter, FilterException


# TODO replace with decorator that filters response_fields
def base_request_args():
    if request.args.get('response_format', 'json') != 'json':
        abort(400, 'Response format is not supported.')

    properties_str = request.args.get('request_fields', None)
    if properties_str is not None:
        return properties_str.split(',')
    return None


@api.route('/calculations')
class CalculationList(Resource):
    @api.doc('list_calculations')
    @api.response(400, 'Invalid requests, e.g. bad parameter.')
    @api.expect(entry_listing_endpoint_parser, validate=True)
    @api.marshal_with(json_api_single_response_model, skip_none=True, code=200)
    def get(self):
        """ Retrieve a list of calculations that match the given Optimade filter expression. """
        request_fields = base_request_args()

        try:
            filter = request.args.get('filter', None)
            page_limit = int(request.args.get('page_limit', 10))
            page_number = int(request.args.get('page_number', 1))
            sort = request.args.get('sort', 'chemical_formula_reduced'),

        except Exception:
            abort(400, message='bad parameter types')  # TODO Specific json API error handling

        search_request = search.SearchRequest().owner('all', None)
        if filter is not None:
            try:
                search_request.query(parse_filter(filter))
            except FilterException as e:
                abort(400, message='Could not parse filter expression: %s' % str(e))

        result = search_request.execute_paginated(
            page=page_number,
            per_page=page_limit)
        # order_by='optimade.%s' % sort)  # TODO map the Optimade property

        available = result['pagination']['total']
       
        print(result['results'][0]['optimade'].keys())
        raise 

        return dict(
            meta=Meta(
                query=request.url,
                returned=len(result['results']),
                available=available,
                last_id=result['results'][-1]['calc_id'] if available > 0 else None),
            links=Links(
                'calculations',
                available=available,
                page_number=page_number,
                page_limit=page_limit,
                sort=sort, filter=filter),
            data=[CalculationDataObject(d, request_fields=request_fields) for d in result['results']]
        ), 200


@api.route('/calculations/<string:id>')
class Calculation(Resource):
    @api.doc('retrieve_calculation')
    @api.response(400, 'Invalid requests, e.g. bad parameter.')
    @api.response(404, 'Id does not exist.')
    @api.expect(single_entry_endpoint_parser, validate=True)
    @api.marshal_with(json_api_single_response_model, skip_none=True, code=200)
    def get(self, id: str):
        """ Retrieve a single calculation for the given id. """
        request_fields = base_request_args()
        search_request = search.SearchRequest().owner('all', None).search_parameters(calc_id=id)

        result = search_request.execute_paginated(
            page=1,
            per_page=1)

        available = result['pagination']['total']
        if available == 0:
            abort(404, 'The calculation with id %s does not exist' % id)

        print('================', result['results'][0])
        raise
        return dict(
            meta=Meta(query=request.url, returned=1),
            data=CalculationDataObject(result['results'][0], request_fields=request_fields)
        ), 200


@api.route('/info/calculation')
class CalculationInfo(Resource):
    @api.doc('calculations_info')
    @api.response(400, 'Invalid requests, e.g. bad parameter.')
    @api.expect(base_endpoint_parser, validate=True)
    @api.marshal_with(json_api_single_response_model, skip_none=True, code=200)
    def get(self):
        """ Returns information relating to the API implementation- """
        base_request_args()

        result = {
            'type': 'info',
            'id': 'calculation',
            'attributes': {
                'description': 'A calculations entry.',
                # TODO non optimade, nomad specific properties
                'properties': {
                    attr.name: dict(description=attr.description)
                    for attr in OptimadeStructureEntry.m_def.attributes.values()
                },
                'formats': ['json'],
                'output_fields_by_format': {
                    'json': OptimadeStructureEntry.m_def.attributes.keys()
                }
            }
        }

        return dict(
            meta=Meta(query=request.url, returned=1),
            data=result
        ), 200


@api.route('/info')
class Info(Resource):
    @api.doc('info')
    @api.response(400, 'Invalid requests, e.g. bad parameter.')
    @api.expect(base_endpoint_parser, validate=True)
    @api.marshal_with(json_api_single_response_model, skip_none=True, code=200)
    def get(self):
        """ Returns information relating to the API implementation- """
        base_request_args()

        result = {
            'type': 'info',
            'id': '/',
            'attributes': {
                'api_version': '0.10.0',
                'available_api_versions': [{
                    'url': url(),
                    'version': '0.10.0'
                }],
                'formats': ['json'],
                'entry_types_by_format': {
                    'json': ['calculations', 'info']
                },
                'available_endpoints': ['calculations', 'info'],
                'is_index': False
            }
        }

        return dict(
            meta=Meta(query=request.url, returned=1),
            data=result
        ), 200
