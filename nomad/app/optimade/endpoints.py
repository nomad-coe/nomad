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

from typing import List, Dict, Any
from flask_restplus import Resource, abort
from flask import request
from elasticsearch_dsl import Q

from nomad import search, files, datamodel, config
from nomad.datamodel import OptimadeEntry

from .api import api, url, base_request_args
from .models import json_api_single_response_model, entry_listing_endpoint_parser, Meta, \
    Links as LinksModel, single_entry_endpoint_parser, base_endpoint_parser, \
    json_api_info_response_model, json_api_list_response_model, EntryDataObject, \
    ToplevelLinks, get_entry_properties, json_api_structure_response_model, \
    json_api_structures_response_model
from .filterparser import parse_filter, FilterException

ns = api.namespace('v0', description='The version v0 API namespace with all OPTiMaDe endpoints.')


def base_search_request():
    ''' Creates a search request for all public and optimade enabled data. '''
    return search.SearchRequest().owner('all', None).search_parameter('processed', True).query(
        Q('exists', field='dft.optimade.elements'))  # TODO use the elastic annotations when done


def to_calc_with_metadata(results: List[Dict[str, Any]]):
    ''' Translates search results into :class:`EntryMetadata` objects read from archive. '''

    upload_files_cache: Dict[str, files.UploadFiles] = {}

    def transform(result):
        calc_id, upload_id = result['calc_id'], result['upload_id']
        upload_files = upload_files_cache.get(upload_id)

        if upload_files is None:
            upload_files = files.UploadFiles.get(upload_id)
            upload_files_cache[upload_id] = upload_files

        archive = upload_files.read_archive(calc_id)  # , access='public')
        metadata = archive[calc_id]['section_metadata'].to_dict()
        return datamodel.EntryMetadata.m_from_dict(metadata)

    result = [transform(result) for result in results]

    for upload_files in upload_files_cache.values():
        upload_files.close()

    return result


# TODO the Entry/ListEntry endpoints for References, Calculations, Structures should
# reuse more code.
# Calculations are identical to structures. Not sure if this is what the optimade
# specification intends.
@ns.route('/calculations')
class CalculationList(Resource):
    @api.doc('list_calculations')
    @api.response(400, 'Invalid requests, e.g. bad parameter.')
    @api.expect(entry_listing_endpoint_parser, validate=True)
    @api.marshal_with(json_api_list_response_model, skip_none=True, code=200)
    def get(self):
        ''' Returns a list of calculations that match the given optimade filter expression. '''
        request_fields = base_request_args()

        try:
            filter = request.args.get('filter', None)
            page_limit = int(request.args.get('page_limit', 10))
            page_number = int(request.args.get('page_number', 1))
            sort = request.args.get('sort', 'chemical_formula_reduced')

        except Exception:
            abort(400, message='bad parameter types')  # TODO Specific json API error handling

        search_request = base_search_request().include('calc_id', 'upload_id')

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
        results = to_calc_with_metadata(result['results'])
        assert len(results) == len(result['results']), 'archive and elasticsearch are not consistent'

        return dict(
            meta=Meta(
                query=request.url,
                returned=len(results),
                available=available,
                last_id=results[-1].calc_id if available > 0 else None),
            links=LinksModel(
                'calculations',
                available=available,
                page_number=page_number,
                page_limit=page_limit,
                sort=sort, filter=filter),
            data=[EntryDataObject(d, optimade_type='calculations', request_fields=request_fields) for d in results]
        ), 200


@ns.route('/calculations/<string:id>')
class Calculation(Resource):
    @api.doc('retrieve_calculation')
    @api.response(400, 'Invalid requests, e.g. bad parameter.')
    @api.response(404, 'Id does not exist.')
    @api.expect(single_entry_endpoint_parser, validate=True)
    @api.marshal_with(json_api_single_response_model, skip_none=True, code=200)
    def get(self, id: str):
        ''' Retrieve a single calculation for the given id '''
        request_fields = base_request_args()
        search_request = base_search_request().search_parameters(calc_id=id)

        result = search_request.execute_paginated(
            page=1,
            per_page=1)

        available = result['pagination']['total']
        results = to_calc_with_metadata(result['results'])
        assert len(results) == len(result['results']), 'Mongodb and elasticsearch are not consistent'

        if available == 0:
            abort(404, 'The calculation with id %s does not exist' % id)

        return dict(
            meta=Meta(query=request.url, returned=1),
            data=EntryDataObject(results[0], optimade_type='calculations', request_fields=request_fields)
        ), 200


@ns.route('/info/calculations')
class CalculationInfo(Resource):
    @api.doc('calculations_info')
    @api.response(400, 'Invalid requests, e.g. bad parameter.')
    @api.expect(base_endpoint_parser, validate=True)
    @api.marshal_with(json_api_info_response_model, skip_none=True, code=200)
    def get(self):
        ''' Returns information about the calculation endpoint implementation '''
        base_request_args()

        result = {
            'description': 'a calculation entry',
            'properties': get_entry_properties(),
            'formats': ['json'],
            'output_fields_by_format': {
                'json': list(OptimadeEntry.m_def.all_properties.keys())}
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
        ''' Returns information about this optimade implementation '''
        base_request_args()

        result = {
            'type': 'info',
            'id': '/',
            'attributes': {
                'api_version': '0.10.1',
                'available_api_versions': [{
                    'url': url(),
                    'version': '0.10.1'
                }],
                'formats': ['json'],
                'entry_types_by_format': {
                    'json': ['structures', 'calculations', 'info']
                },
                'available_endpoints': ['structures', 'calculations', 'info'],
                'is_index': False
            }
        }

        return dict(
            meta=Meta(query=request.url, returned=1),
            data=result
        ), 200


def execute_search(**kwargs):
    filter = kwargs.get('filter')
    page_number = kwargs.get('page_number')
    page_limit = kwargs.get('page_limit')
    search_request = base_search_request().include('calc_id', 'upload_id')

    if filter is not None:
        try:
            search_request.query(parse_filter(filter))
        except FilterException as e:
            abort(400, message='Could not parse filter expression: %s' % str(e))

    result = search_request.execute_paginated(
        page=page_number,
        per_page=page_limit)
    # order_by='optimade.%s' % sort)  # TODO map the Optimade property

    return result


# TODO This does not return reference
# TODO This also needs a single entry endpoint?
# TODO This also needs an info endpoint
# @ns.route('/references')
# class References(Resource):
#     @api.doc('references')
#     @api.response(400, 'Invalid requests, e.g. bad parameter.')
#     @api.response(422, 'Validation error')
#     @api.expect(entry_listing_endpoint_parser, validate=True)
#     @api.marshal_with(json_api_references_response_model, skip_none=True, code=200)
#     def get(self):
#         ''' Returns references for the structures that match the given optimade filter expression'''
#         try:
#             filter = request.args.get('filter', None)
#             page_limit = int(request.args.get('page_limit', 10))
#             page_number = int(request.args.get('page_number', 1))
#             sort = request.args.get('sort', 'chemical_formula_reduced'),

#         except Exception:
#             abort(400, message='bad parameter types')  # TODO Specific json API error handling

#         result = execute_search(
#             filter=filter, page_limit=page_limit, page_number=page_number, sort=sort)
#         available = result['pagination']['total']
#         results = to_calc_with_metadata(result['results'])
#         assert len(results) == len(result['results']), 'Mongodb and elasticsearch are not consistent'

#         # TODO References are about returning user provided references to paper or web resources.
#         # The ReferenceObject does not have this kind of information.
#         # TODO Why is TopLevelLinks different from LinksModel. Any what is "TopLevel" about it.
#         return dict(
#             meta=Meta(
#                 query=request.url,
#                 returned=len(results),
#                 available=available,
#                 last_id=results[-1].calc_id if available > 0 else None),
#             links=ToplevelLinks(
#                 'structures',
#                 available=available,
#                 page_number=page_number,
#                 page_limit=page_limit,
#                 sort=sort, filter=filter),
#             data=[ReferenceObject(d) for d in results]
#         ), 200


@ns.route('/links')
class Links(Resource):
    @api.doc('links')
    @api.response(400, 'Invalid requests, e.g. bad parameter.')
    @api.expect(base_endpoint_parser, validate=True)
    @api.marshal_with(json_api_list_response_model, skip_none=True, code=200)
    def get(self):
        ''' Returns information about related optimade databases '''
        base_request_args()

        result = [
            {
                "type": "parent",
                "id": "index",
                "attributes": {
                    "name": config.meta.name,
                    "description": config.meta.description,
                    "base_url": {
                        "href": url(version=None, prefix='index'),
                    },
                    "homepage": config.meta.homepage
                }
            }
        ]

        return dict(
            meta=Meta(query=request.url, returned=1),
            data=result
        ), 200


@ns.route('/structures')
class StructureList(Resource):
    @api.doc('structures')
    @api.response(400, 'Invalid requests, e.g. bad parameter.')
    @api.response(422, 'Validation error')
    @api.expect(entry_listing_endpoint_parser, validate=True)
    @api.marshal_with(json_api_structures_response_model, skip_none=True, code=200)
    def get(self):
        ''' Retrieve the structures that match the given optimade filter expression '''
        request_fields = base_request_args()

        try:
            filter = request.args.get('filter', None)
            page_limit = int(request.args.get('page_limit', 10))
            page_number = int(request.args.get('page_number', 1))
            sort = request.args.get('sort', 'chemical_formula_reduced'),

        except Exception:
            abort(400, message='bad parameter types')  # TODO Specific json API error handling

        result = execute_search(
            filter=filter, page_limit=page_limit, page_number=page_number, sort=sort)
        available = result['pagination']['total']
        results = to_calc_with_metadata(result['results'])
        assert len(results) == len(result['results']), 'Mongodb and elasticsearch are not consistent'

        return dict(
            meta=Meta(
                query=request.url,
                returned=len(results),
                available=available,
                last_id=results[-1].calc_id if available > 0 else None),
            links=ToplevelLinks(
                'structures',
                available=available,
                page_number=page_number,
                page_limit=page_limit,
                sort=sort, filter=filter
            ),
            data=[EntryDataObject(d, optimade_type='structures', request_fields=request_fields) for d in results]
        ), 200


@ns.route('/structures/<string:id>')
class Structure(Resource):
    @api.doc('structure')
    @api.response(400, 'Invalid requests, e.g. bad parameter.')
    @api.response(404, 'Id does not exist.')
    @api.expect(single_entry_endpoint_parser, validate=True)
    @api.marshal_with(json_api_structure_response_model, skip_none=True, code=200)
    def get(self, id: str):
        ''' Retrieve a single structure for the given id '''
        request_fields = base_request_args()
        search_request = base_search_request().search_parameters(calc_id=id)

        result = search_request.execute_paginated(
            page=1,
            per_page=1)

        available = result['pagination']['total']
        results = to_calc_with_metadata(result['results'])
        assert len(results) == len(result['results']), 'Mongodb and elasticsearch are not consistent'

        if available == 0:
            abort(404, 'The calculation with id %s does not exist' % id)

        return dict(
            meta=Meta(query=request.url, returned=1),
            data=EntryDataObject(results[0], optimade_type='structures', request_fields=request_fields)
        ), 200


@ns.route('/info/structures')
class StructuresInfo(Resource):
    @api.doc('structures_info')
    @api.response(400, 'Invalid requests, e.g. bad parameter.')
    @api.expect(base_endpoint_parser, validate=True)
    @api.marshal_with(json_api_info_response_model, skip_none=True, code=200)
    def get(self):
        ''' Returns information about the structures endpoint implementation '''
        base_request_args()

        result = {
            'description': 'a structure entry',
            'properties': get_entry_properties(),
            'formats': ['json'],
            'output_fields_by_format': {
                'json': list(OptimadeEntry.m_def.all_properties.keys())}
        }

        return dict(
            meta=Meta(query=request.url, returned=1),
            data=result
        ), 200
