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
API endpoint that deliver backend configuration details.
'''

from typing import Dict, Any
from flask_restplus import Resource, fields
from datetime import datetime

from nomad import config, normalizing, datamodel, gitinfo, search
from nomad.parsing import parsers, MatchingParser

from .api import api


ns = api.namespace('info', description='Access to nomad configuration details.')

metainfo_model = api.model('Metainfo', {
    'all_package': fields.String(description='Name of the metainfo package that references all available packages, i.e. the complete metainfo.'),
    'root_section': fields.String(description='Name of the topmost section, e.g. section_run for computational material science data.')
})

domain_model = api.model('Domain', {
    'name': fields.String,
    'metainfo': fields.Nested(model=metainfo_model)
})

git_info_model = api.model('GitInfo', {
    'ref': fields.String,
    'version': fields.String,
    'commit': fields.String,
    'log': fields.String
})

statistics_info_model = api.model('StatisticsInfo', {
    'n_entries': fields.Integer(description='Number of entries in NOMAD'),
    'n_uploads': fields.Integer(description='Number of uploads in NOMAD'),
    'n_quantities': fields.Integer(description='Accumulated number of quantities over all entries in the Archive'),
    'n_calculations': fields.Integer(description='Accumulated number of calculations, e.g. total energy calculations in the Archive'),
    # TODO
    # 'raw_file_size': fields.Integer(description='Total amount of raw files in TB'),
    # 'archive_file_size': fields.Integer(description='Total amount of binary archive data in TB')
})

code_info_model = api.model('CodeInfo', {
    'code_name': fields.String(description='Name of the code or input format', allow_null=True),
    'code_homepage': fields.String(description='Homepage of the code or input format', allow_null=True)
}, allow_null=True, skip_none=True)

info_model = api.model('Info', {
    'parsers': fields.List(fields.String),
    'metainfo_packages': fields.List(fields.String),
    'codes': fields.List(fields.Nested(code_info_model)),
    'normalizers': fields.List(fields.String),
    'domains': fields.List(fields.Nested(model=domain_model)),
    'statistics': fields.Nested(model=statistics_info_model, description='General NOMAD statistics'),
    'search_quantities': fields.Raw(),
    'version': fields.String,
    'release': fields.String,
    'git': fields.Nested(model=git_info_model),
    'oasis': fields.Boolean
})


_statistics: Dict[str, Any] = None


def statistics():
    global _statistics
    if _statistics is None or datetime.now().timestamp() - _statistics.get('timestamp', 0) > 3600 * 24:
        _statistics = dict(timestamp=datetime.now().timestamp())
        _statistics.update(
            **search.SearchRequest().global_statistics().execute()['global_statistics'])

    return _statistics


@ns.route('/')
class InfoResource(Resource):
    @api.doc('get_info')
    @api.marshal_with(info_model, skip_none=True, code=200, description='Info send')
    def get(self):
        ''' Return information about the nomad backend and its configuration. '''
        codes_dict = {}
        for parser in parsers.parser_dict.values():
            if isinstance(parser, MatchingParser) and parser.domain == 'dft':
                code_name = parser.code_name
                if code_name in codes_dict:
                    continue
                codes_dict[code_name] = dict(code_name=code_name, code_homepage=parser.code_homepage)
        codes = sorted(list(codes_dict.values()), key=lambda code_info: code_info['code_name'].lower())

        return {
            'parsers': [
                key[key.index('/') + 1:]
                for key in parsers.parser_dict.keys()],
            'metainfo_packages': ['general', 'general.experimental', 'common', 'public'] + sorted([
                key[key.index('/') + 1:]
                for key in parsers.parser_dict.keys()]),
            'codes': codes,
            'normalizers': [normalizer.__name__ for normalizer in normalizing.normalizers],
            'statistics': statistics(),
            'domains': [
                {
                    'name': domain_name,
                    'metainfo': {
                        'all_package': domain['metainfo_all_package'],
                        'root_section': domain['root_section']
                    }
                }
                for domain_name, domain in datamodel.domains.items()
            ],
            'search_quantities': {
                s.qualified_name: {
                    'name': s.qualified_name,
                    'description': s.description,
                    'many': s.many
                }
                for s in search.search_quantities.values()
                if 'optimade' not in s.qualified_name
            },
            'version': config.meta.version,
            'release': config.meta.release,
            'git': {
                'ref': gitinfo.ref,
                'version': gitinfo.version,
                'commit': gitinfo.commit,
                'log': gitinfo.log
            },
            'oasis': config.keycloak.oasis
        }, 200
