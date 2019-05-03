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
API endpoint that deliver backend configuration details.
"""

from flask_restplus import Resource, fields

from nomad import config, parsing, normalizing, datamodel, gitinfo

from .app import api


ns = api.namespace('info', description='Access to nomad configuration details.')

domain_quantity_model = api.model('DomainQuantity', {
    'name': fields.String,
    'description': fields.String,
    'multi': fields.Boolean,
    'order_default': fields.Boolean
})

metainfo_model = api.model('Metainfo', {
    'all_package': fields.String(description='Name of the metainfo package that references all available packages, i.e. the complete metainfo.'),
    'root_sections': fields.List(fields.String, description='Name of the topmost section, e.g. section_run for computational material science data.')
})

domain_model = api.model('Domain', {
    'name': fields.String,
    'quantities': fields.List(fields.Nested(model=domain_quantity_model)),
    'aggregations_names': fields.List(fields.String),
    'metrics_names': fields.List(fields.String),
    'metainfo': fields.Nested(model=metainfo_model)
})

git_info_model = api.model('GitInfo', {
    'ref': fields.String,
    'version': fields.String,
    'log': fields.String
})

info_model = api.model('Info', {
    'parsers': fields.List(fields.String),
    'normalizers': fields.List(fields.String),
    'domain': fields.Nested(model=domain_model),
    'version': fields.String,
    'release': fields.String,
    'git': fields.Nested(model=git_info_model)
})


@ns.route('/')
class InfoResource(Resource):
    @api.doc('get_info')
    @api.marshal_with(info_model, skip_none=True, code=200, description='Info send')
    def get(self):
        """ Return information about the nomad backend and its configuration. """
        return {
            'parsers': [key[8:] for key in parsing.parser_dict.keys()],
            'normalizers': [normalizer.__name__ for normalizer in normalizing.normalizers],
            'domain': {
                'name': datamodel.Domain.instance.name,
                'quantities': [quantity for quantity in datamodel.Domain.instance.quantities.values()],
                'metrics_names': datamodel.Domain.instance.metrics_names,
                'aggregations_names': datamodel.Domain.instance.aggregations_names,
                'metainfo': {
                    'all_package': datamodel.Domain.instance.metainfo_all_package,
                    'root_sections': datamodel.Domain.instance.root_sections
                }
            },
            'version': config.version,
            'release': config.release,
            'git': {
                'ref': gitinfo.ref,
                'version': gitinfo.version,
                'log': gitinfo.log
            }
        }, 200
