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
from flask_restplus import Resource, abort
from elasticsearch.exceptions import NotFoundError

from nomad import search
from nomad.app.flask.api.auth import authenticate

from .api import api, arg_parser, rdf_respose, response_types
from .mapping import Mapping

ns = api.namespace('datasets', description='The API for DCAT datasets.')


@ns.route('/<string:entry_id>')
class Dataset(Resource):
    @api.doc('get_dcat_dataset')
    @api.expect(arg_parser)
    @api.representation('application/xml')
    @api.produces(response_types)
    @api.response(404, 'There is no entry with the given id.')
    @api.response(401, 'This entry is not publically accessible.')
    @api.response(200, 'Data send', headers={'Content-Type': 'application/xml'})
    @authenticate()
    def get(self, entry_id):
        ''' Returns a DCAT dataset for a given NOMAD entry id. '''
        try:
            entry = search.entry_document.get(entry_id)
        except NotFoundError:
            abort(404, message='There is no calculation with id %s' % entry_id)

        if entry.with_embargo or not entry.published:
            abort(401, message='Not authorized to access %s' % entry_id)

        mapping = Mapping()
        mapping.map_entry(entry)
        return rdf_respose(mapping.g)
