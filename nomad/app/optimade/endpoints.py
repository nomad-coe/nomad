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
from elasticsearch_dsl import Q

from nomad import search
from nomad.metainfo.optimade import OptimadeEntry

from .api import api, url
from .models import json_api_single_response_model, entry_listing_endpoint_parser, Meta, \
    Links, CalculationDataObject, single_entry_endpoint_parser, base_endpoint_parser, \
    json_api_info_response_model, json_api_list_response_model
from .filterparser import parse_filter, FilterException


ns = api.namespace('', description='The (only) API namespace with all OPTiMaDe endpoints.')


# TODO replace with decorator that filters response_fields
def base_request_args():
    if request.args.get('response_format', 'json') != 'json':
        abort(400, 'Response format is not supported.')

    properties_str = request.args.get('request_fields', None)
    if properties_str is not None:
        return properties_str.split(',')
    return None


def base_search_request():
    """ Creates a search request for all public and optimade enabled data. """
    return search.SearchRequest().owner('all', None).query(
        Q('exists', field='optimade.nelements'))  # TODO use the elastic annotations when done


@ns.route('/calculations')
class CalculationList(Resource):
    @api.doc('list_calculations')
    @api.response(400, 'Invalid requests, e.g. bad parameter.')
    @api.expect(entry_listing_endpoint_parser, validate=True)
    @api.marshal_with(json_api_list_response_model, skip_none=True, code=200)
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

        search_request = base_search_request().include('calc_id')

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
        results = search.to_calc_with_metadata(result['results'])
        assert len(results) == len(result['results']), 'Mongodb and elasticsearch are not consistent'

        return dict(
            meta=Meta(
                query=request.url,
                returned=len(results),
                available=available,
                last_id=results[-1].calc_id if available > 0 else None),
            links=Links(
                'calculations',
                available=available,
                page_number=page_number,
                page_limit=page_limit,
                sort=sort, filter=filter),
            data=[CalculationDataObject(d, request_fields=request_fields) for d in results]
        ), 200


@ns.route('/calculations/<string:id>')
class Calculation(Resource):
    @api.doc('retrieve_calculation')
    @api.response(400, 'Invalid requests, e.g. bad parameter.')
    @api.response(404, 'Id does not exist.')
    @api.expect(single_entry_endpoint_parser, validate=True)
    @api.marshal_with(json_api_single_response_model, skip_none=True, code=200)
    def get(self, id: str):
        """ Retrieve a single calculation for the given id. """
        request_fields = base_request_args()
        search_request = base_search_request().search_parameters(calc_id=id)

        result = search_request.execute_paginated(
            page=1,
            per_page=1)

        available = result['pagination']['total']
        results = search.to_calc_with_metadata(result['results'])
        assert len(results) == len(result['results']), 'Mongodb and elasticsearch are not consistent'

        if available == 0:
            abort(404, 'The calculation with id %s does not exist' % id)

        return dict(
            meta=Meta(query=request.url, returned=1),
            data=CalculationDataObject(results[0], request_fields=request_fields)
        ), 200


@ns.route('/info/calculations')
class CalculationInfo(Resource):
    @api.doc('calculations_info')
    @api.response(400, 'Invalid requests, e.g. bad parameter.')
    @api.expect(base_endpoint_parser, validate=True)
    @api.marshal_with(json_api_info_response_model, skip_none=True, code=200)
    def get(self):
        """ Returns information relating to the API implementation- """
        base_request_args()

        result = {
            'description': 'a calculation entry',
            'properties': {
                attr.name: dict(description=attr.description)
                for attr in OptimadeEntry.m_def.all_properties.values()},
            'formats': ['json'],
            'output_fields_by_format': {
                'json': OptimadeEntry.m_def.all_properties.keys()}
        }

        return dict(
            meta=Meta(query=request.url, returned=1),
            data=result
        ), 200


@ns.route('/info')
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
