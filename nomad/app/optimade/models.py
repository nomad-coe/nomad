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
All the API flask restplus models.
"""

from typing import Dict, Any, Set
from flask_restplus import fields
import datetime
import math

from nomad import config
from nomad.app.utils import RFC3339DateTime

from .api import api, base_url, url


# TODO error/warning objects

json_api_meta_object_model = api.model('JsonApiMetaObject', {
    'query': fields.Nested(model=api.model('JsonApiQuery', {
        'representation': fields.String(
            required=True,
            description='A string with the part of the URL following the base URL')
    }, description='Information on the query that was requested.')),

    'api_version': fields.String(
        required=True,
        description='A string containing the version of the API implementation.'),

    'time_stamp': RFC3339DateTime(
        required=True,
        description='A timestamp containing the date and time at which the query was executed.'),

    'data_returned': fields.Integer(
        required=True,
        description='An integer containing the number of data objects returned for the query.'),

    'more_data_available': fields.Boolean(
        required=True,
        description='False if all data for this query has been returned, and true if not.'),

    'provider': fields.Nested(
        required=True, skip_none=True,
        description='Information on the database provider of the implementation.',
        model=api.model('JsonApiProvider', {
            'name': fields.String(
                required=True,
                description='A short name for the database provider.'),
            'description': fields.String(
                required=True,
                description='A longer description of the database provider.'),
            'prefix': fields.String(
                required=True,
                description='Database-provider-specific prefix.'),
            'homepage': fields.String(
                required=False,
                description='Homepage of the database provider'),
            'index_base_url': fields.String(
                required=False,
                description='Base URL for the index meta-database.')
        })),

    'data_available': fields.Integer(
        required=False,
        description=('An integer containing the total number of data objects available in '
                     'the database.')),

    'last_id': fields.String(
        required=False,
        description='A string containing the last ID returned'),

    'response_message': fields.String(
        required=False,
        description='Response string from the server.'),

    'implementation': fields.Nested(
        required=False, skip_none=True,
        description='Server implementation details.',
        model=api.model('JsonApiImplementation', {
            'name': fields.String(
                description='Name of the implementation'),
            'version': fields.String(
                description='Version string of the current implementation.'),
            'source_url': fields.String(
                description=' URL of the implementation source, either downloadable archive or version control system.'),
            'maintainer': fields.Nested(
                skip_none=True,
                description='Details about the maintainer of the implementation',
                model=api.model('JsonApiMaintainer', {
                    'email': fields.String()
                })
            )
        }))
})


class Meta():

    def __init__(self, query: str, returned: int, available: int = None, last_id: str = None):
        self.query = dict(representation=query)
        self.api_version = '0.10.0'
        self.time_stamp = datetime.datetime.now()
        self.data_returned = returned
        self.more_data_available = available > returned if available is not None else False
        self.provider = dict(
            name='NOMAD',
            description='The NOvel MAterials Discovery project and database.',
            prefix='nomad',
            homepage='https//nomad-coe.eu',
            index_base_url=base_url
        )

        self.data_available = available
        self.last_id = last_id
        self.implementation = dict(
            name='nomad@fairdi',
            version=config.version,
            source_url='https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR',
            maintainer=dict(email='markus.scheidgen@physik.hu-berlin.de'))


json_api_links_model = api.model('JsonApiLinks', {
    'base_url': fields.String(
        description='The base URL of the implementation'),

    'next': fields.String(
        description=('A link to fetch the next set of results. When the current response '
                     'is the last page of data, this field MUST be either omitted or null.')),

    'prev': fields.String(
        description=('The previous page of data. null or omitted when the current response '
                     'is the first page of data.')),

    'last': fields.String(
        description='The last page of data.'),

    'first': fields.String(
        description='The first page of data.')

})


def Links(endpoint: str, available: int, page_number: int, page_limit: int, **kwargs):
    last_page = math.ceil(available / page_limit)

    rest = dict(page_limit=page_limit)
    rest.update(**{key: value for key, value in kwargs.items() if value is not None})

    result = dict(
        base_url=url(),
        first=url(endpoint, page_number=1, **rest),
        last=url(endpoint, page_number=last_page, **rest))

    if page_number > 1:
        result['prev'] = url(endpoint, page_number=page_number - 1, **rest)
    if page_number * page_limit < available:
        result['next'] = url(endpoint, page_number=page_number + 1, **rest)

    return result


json_api_response_model = api.model('JsonApiResponse', {
    'links': fields.Nested(
        required=False,
        description='Links object with pagination links.',
        skip_none=True,
        model=json_api_links_model),

    'meta': fields.Nested(
        required=True,
        description='JSON API meta object.',
        model=json_api_meta_object_model),

    'included': fields.List(
        fields.Arbitrary(),
        required=False,
        description=('A list of JSON API resource objects related to the primary data '
                     'contained in data. Responses that contain related resources under '
                     'included are known as compound documents in the JSON API.'))
})

json_api_data_object_model = api.model('JsonApiDataObject', {
    'type': fields.String(
        description='The type of the object [structure or calculations].'),

    'id': fields.String(
        description='The id of the object.'),

    'immutable_id': fields.String(
        description='The entries immutable id.'),

    'last_modified': RFC3339DateTime(
        description='Date and time representing when the entry was last modified.'),

    'attributes': fields.Raw(
        description='A dictionary, containing key-value pairs representing the entries properties')

    # further optional fields: links, meta, relationships
})


class CalculationDataObject:
    def __init__(self, search_entry_dict: Dict[str, Any], request_fields: Set[str] = None):
        attrs = search_entry_dict

        attrs = {
            key: value for key, value in search_entry_dict['optimade'].items()
            if key is not 'elements_ratios' and (request_fields is None or key in request_fields)
        }
        if request_fields is None or 'elements_ratios' in request_fields:
            attrs['elements_ratios'] = [
                d['elements_ratios'] for d in search_entry_dict['optimade']['elements_ratios']
            ]
        attrs.update(**{
            '_nomad_%s' % key: value for key, value in search_entry_dict.items()
            if key != 'optimade' and (request_fields is None or '_nomad_%s' % key in request_fields)
        })

        self.type = 'calculation'
        self.id = search_entry_dict['calc_id']
        self.immutable_id = search_entry_dict['calc_id']
        self.last_modified = search_entry_dict.get(
            'last_processing', search_entry_dict.get('upload_time', None))
        self.attributes = attrs


class Property:
    @staticmethod
    def from_nomad_to_optimade(name: str):
        if name.startswith('optimade.'):
            return name[9:]
        else:
            return '_nomad_%s' % name

    @staticmethod
    def from_optimade_to_nomad(name: str):
        if name.startswith('_nomad_'):
            return name[7:]
        else:
            return 'optimade.%s' % name


json_api_single_response_model = api.inherit(
    'JsonApiSingleResponse', json_api_response_model, {
        'data': fields.Nested(
            model=json_api_data_object_model,
            required=True,
            description=('The returned response object.'))
    })

json_api_list_response_model = api.inherit(
    'JsonApiSingleResponse', json_api_response_model, {
        'data': fields.List(
            fields.Nested(json_api_data_object_model),
            required=True,
            description=('The list of returned response objects.'))
    })


base_endpoint_parser = api.parser()
base_endpoint_parser.add_argument(
    'response_format', type=str,
    help=('The output format requested. Defaults to the format string "json", which '
          'specifies the standard output format described in this specification.'))
base_endpoint_parser.add_argument(
    'email_address', type=str,
    help=('An email address of the user making the request. The email SHOULD be that of a '
          'person and not an automatic system.'))
base_endpoint_parser.add_argument(
    'response_fields', action='split', type=str,
    help=('A comma-delimited set of fields to be provided in the output. If provided, only '
          'these fields MUST be returned and no others.'))

entry_listing_endpoint_parser = base_endpoint_parser.copy()
entry_listing_endpoint_parser.add_argument(
    'filter', type=str, help='An optimade filter string.')
entry_listing_endpoint_parser.add_argument(
    'page_limit', type=int,
    help='Sets a numerical limit on the number of entries returned.')
entry_listing_endpoint_parser.add_argument(
    'page_number', type=int,
    help='Sets the page number to return.')
entry_listing_endpoint_parser.add_argument(
    'sort', type=str,
    help='Name of the property to sort the results by.')

single_entry_endpoint_parser = base_endpoint_parser.copy()
