#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

from flask import Blueprint, Response
from flask_restplus import Api, reqparse
import urllib.parse
from rdflib import Graph

from nomad import config

blueprint = Blueprint('dcat', __name__)

base_url = config.api_url(api='dcat')


def url(*args, **kwargs):
    ''' Returns the full dcat api url for the given path (args) and query (kwargs) parameters. '''
    url = f'{base_url.rstrip("/")}/{"/".join(args).lstrip("/")}'

    if len(kwargs) > 0:
        return f'{url}?{urllib.parse.urlencode(kwargs)}'
    else:
        return url


api = Api(
    blueprint,
    version='1.0', title='NOMAD\'s API for servicing dcat resources',
    description='NOMAD\'s API for serving dcat resources',
    validate=True)


# For some unknown reason it is necessary for each fr api to have a handler.
# Otherwise the global app error handler won't be called.
@api.errorhandler(Exception)
def errorhandler(error):
    '''When an internal server error is caused by an unexpected exception.'''
    return str(error)


arg_parser = reqparse.RequestParser()
arg_parser.add_argument('format', type=str, choices=[
    'xml',
    'n3',
    'turtle',
    'nt',
    'pretty-xml',
    'trig'])


def rdf_respose(g: Graph) -> Response:
    args = arg_parser.parse_args()
    format_ = args.get('format')
    if format_ is None:
        format_ = 'pretty-xml'
    content_type = 'application/xml' if format in ['xml', 'pretty-xml'] else 'text/%s' % format_
    return Response(
        g.serialize(format=format_).decode('utf-8'), 200,
        {'Content-Type': content_type})
