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

'''
All the API flask restplus models.
'''

from typing import Set, List, Dict, Any
from flask_restplus import fields
import datetime
import math
from cachetools import cached

from nomad import config, datamodel, files
from nomad.app.common import RFC3339DateTime
from nomad.datamodel import EntryMetadata
from nomad.datamodel.dft import DFTMetadata
from nomad.datamodel.optimade import OptimadeEntry

from .api import api, url


# TODO error/warning objects

json_api_meta_object_model = api.model('MetaObject', {
    'query': fields.Nested(model=api.model('Query', {
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
        model=api.model('Provider', {
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
        model=api.model('Implementation', {
            'name': fields.String(
                description='Name of the implementation'),
            'version': fields.String(
                description='Version string of the current implementation.'),
            'source_url': fields.String(
                description=' URL of the implementation source, either downloadable archive or version control system.'),
            'maintainer': fields.Nested(
                skip_none=True,
                description='Details about the maintainer of the implementation',
                model=api.model('Maintainer', {
                    'email': fields.String()
                })
            )
        }))
})


class Meta():

    def __init__(
            self, query: str, returned: int, available: int = None, last_id: str = None):

        self.query = dict(representation=query)
        self.api_version = '1.0.0'
        self.time_stamp = datetime.datetime.now()
        self.data_returned = returned
        self.more_data_available = available > returned if available is not None else False
        self.provider = dict(
            name=config.meta.name,
            description=config.meta.description,
            prefix='nomad',
            homepage=config.meta.homepage,
            index_base_url=url(version=None, prefix='index')
        )

        self.data_available = available
        self.last_id = last_id
        self.implementation = dict(
            name='nomad@fairdi',
            version=config.meta.version,
            source_url=config.meta.source_url,
            maintainer=dict(email=config.meta.maintainer_email))


class Links:
    def __init__(self, endpoint: str, returned: int, page_number: int, page_limit: int, **kwargs):
        last_page = math.ceil(returned / page_limit)

        rest = dict(page_limit=page_limit)
        rest.update(**{key: value for key, value in kwargs.items() if value is not None})
        self.self = url()
        self.related = None
        self.first = url(endpoint, page_number=1, **rest)
        self.last = url(endpoint, page_number=last_page, **rest)
        self.prev = url(endpoint, page_number=max((page_number - 1, 1)), **rest)
        self.next = url(endpoint, page_number=min((page_number + 1, last_page)), **rest)


json_api_links_model = api.model('ApiLinks', {
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


json_api_response_model = api.model('Response', {
    'links': fields.Nested(
        required=False,
        description='Links object with pagination links.',
        skip_none=True,
        model=json_api_links_model),

    'meta': fields.Nested(
        required=True, skip_none=True,
        description='JSON API meta object.',
        model=json_api_meta_object_model),

    'included': fields.List(
        fields.Arbitrary(),
        required=False, skip_none=True,
        description=('A list of JSON API resource objects related to the primary data '
                     'contained in data. Responses that contain related resources under '
                     'included are known as compound documents in the JSON API.'))
})

json_api_data_object_model = api.model('DataObject', {
    'type': fields.String(
        description='The type of the object [structure or calculations].'),

    'id': fields.String(
        description='The id of the object.'),

    'attributes': fields.Raw(
        description='A dictionary, containing key-value pairs representing the entries properties'),

    'relationships': fields.Raw(
        description='In accordance with section Relationships, all entry types MAY use relationships to describe relations to other entries.')
})


json_api_info_model = api.model('CalculationInfo', {
    'description': fields.String(
        description='Description of the entry'),

    'properties': fields.Raw(
        description=('A dictionary describing queryable properties for this '
                     'entry type, where each key is a property name')),

    'formats': fields.List(
        fields.String(),
        required=True,
        description='List of output formats available for this type of entry'),

    'output_fields_by_format': fields.Raw(
        description=('Dictionary of available output fields for this entry'
                     'type, where the keys are the values of the formats list'
                     'and the values are the keys of the properties dictionary'))

})


@cached({})
def get_entry_properties(include_optimade: bool = True):
    properties = {
        attr.name: dict(description=attr.description)
        for attr in OptimadeEntry.m_def.all_properties.values()
        if include_optimade}

    def add_nmd_properties(prefix, section_cls):
        for quantity in section_cls.m_def.all_quantities.values():
            name = prefix + quantity.name
            properties[name] = dict(description=quantity.description)

    add_nmd_properties('_nmd_', EntryMetadata)
    add_nmd_properties('_nmd_dft_', DFTMetadata)

    return properties


class EntryDataObject:
    def __init__(self, calc: EntryMetadata, optimade_type: str, response_fields: Set[str] = None):

        def include(key):
            if optimade_type == 'calculations' and key not in ['immutable_id', 'last_modified']:
                return False

            if response_fields is None or (key in response_fields):
                return True

            return False

        attrs = {key: value for key, value in calc.dft.optimade.m_to_dict().items() if include(key)}
        if include('immutable_id'):
            attrs['immutable_id'] = calc.calc_id
        if include('last_modified'):
            attrs['last_modified'] = calc.last_processing if calc.last_processing is not None else calc.upload_time

        if response_fields is not None:
            for request_field in response_fields:
                if not request_field.startswith('_nmd_'):
                    continue

                try:
                    if request_field.startswith('_nmd_dft_'):
                        attrs[request_field] = getattr(calc.dft, request_field[9:])
                    else:
                        attrs[request_field] = getattr(calc, request_field[5:])
                except AttributeError:
                    # if unknown properties where provided, we will ignore them
                    pass

        self.type = optimade_type
        self.id = calc.calc_id
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
    'SingleResponse', json_api_response_model, {
        'data': fields.Nested(
            model=json_api_data_object_model,
            required=True, skip_none=True,
            description=('The returned response object.'))
    })

json_api_list_response_model = api.inherit(
    'ListResponse', json_api_response_model, {
        'data': fields.List(
            fields.Nested(json_api_data_object_model, skip_none=True),
            required=True,
            description=('The list of returned response objects.'))
    })

json_api_info_response_model = api.inherit(
    'InfoResponse', json_api_response_model, {
        'data': fields.Nested(
            model=json_api_info_model,
            required=True,
            description=('The returned response object.'))
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
