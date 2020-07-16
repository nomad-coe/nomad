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

from flask_restplus import Resource
from flask import request

from nomad import config

from .api import api, url
from .common import base_request_args
from .models import json_api_single_response_model, base_endpoint_parser, json_api_single_response_model, Meta, json_api_list_response_model

ns = api.namespace('index/v1', description='This is the OPTiMaDe index for NOMAD\' implementations.')


@ns.route('/info')
class Info(Resource):
    @api.doc('index_info')
    @api.response(400, 'Invalid requests, e.g. bad parameter.')
    @api.expect(base_endpoint_parser, validate=True)
    @api.marshal_with(json_api_single_response_model, skip_none=True, code=200)
    def get(self):
        ''' Returns information relating to the API implementation- '''
        base_request_args()

        result = {
            'type': 'info',
            'id': '/',
            'attributes': {
                'api_version': '1.0.0',
                'available_api_versions': [{
                    'url': url(prefix='index'),
                    'version': '1.0.0'
                }],
                'formats': ['json'],
                'entry_types_by_format': {
                    'json': ['links', 'info']
                },
                'available_endpoints': ['links', 'info'],
                'is_index': True
            },
            'relationships': {
                'default': {
                    'data': {
                        'id': 'v1',
                        'type': 'links'
                    }
                }
            }
        }

        return dict(
            meta=Meta(query=request.url, returned=1),
            data=result
        ), 200


@ns.route('/links')
class Links(Resource):
    @api.doc('index_links')
    @api.response(400, 'Invalid requests, e.g. bad parameter.')
    @api.expect(base_endpoint_parser, validate=True)
    @api.marshal_with(json_api_list_response_model, skip_none=True, code=200)
    def get(self):
        ''' Returns information relating to the API implementation- '''
        base_request_args()

        result = [
            {
                "type": "child",
                "id": "v1",
                "attributes": {
                    "name": config.meta.name,
                    "description": config.meta.description,
                    "base_url": {
                        "href": url(version=None),
                    },
                    "homepage": config.meta.homepage
                }
            }
        ]

        return dict(
            meta=Meta(query=request.url, returned=1),
            data=result
        ), 200
