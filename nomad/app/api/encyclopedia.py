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
The encyclopedia API of the nomad@FAIRDI APIs.
'''

from flask_restplus import Resource, abort, fields
from elasticsearch_dsl import Search

from .api import api
from .auth import authenticate

from nomad import config

ns = api.namespace('encyclopedia', description='Access encyclopedia metadata.')
search = Search(index=config.elastic.index_name)


@ns.route('/materials/<string:material_id>')
class EncMaterialResource(Resource):
    @api.response(404, 'The material does not exist')
    @api.response(401, 'Not authorized to access the material')
    @api.response(200, 'Metadata send', fields.Raw)
    @api.doc('get_enc_material')
    @authenticate()
    def get(self, material_id):
        """Used to retrive basic information related to the specified material.
        """
        def add_result(result, key, function, default=""):
            try:
                value = function()
            except Exception:
                value = default
            result[key] = value

        # Find the first entry with this material id and take information from
        # there. In principle all other entries should have the same
        # information.
        s = search.query('term', encyclopedia__material__material_id=material_id)
        response = s.execute()

        if len(response) == 0:
            abort(404, message='There is no material {}'.format(material_id))

        entry = response[0]

        # Create result JSON
        result = {}
        result["material_id"] = material_id
        add_result(result, "bravais_lattice", lambda: entry.encyclopedia.material.bulk.bravais_lattice, ""),
        add_result(result, "crystal_system", lambda: entry.encyclopedia.material.bulk.crystal_system, "")
        add_result(result, "formula", lambda: entry.encyclopedia.material.formula, "")
        add_result(result, "formula_reduced", lambda: entry.encyclopedia.material.formula_reduced, "")
        add_result(result, "material_name", lambda: entry.encyclopedia.material.material_name, "")
        add_result(result, "point_group", lambda: entry.encyclopedia.material.bulk.point_group, "")
        add_result(result, "space_group", lambda: entry.encyclopedia.material.bulk.space_group_number, "")
        add_result(result, "structure_type", lambda: entry.encyclopedia.material.bulk.structure_type, "")
        add_result(result, "system_type", lambda: entry.encyclopedia.material.material_type, "")

        return result, 200
