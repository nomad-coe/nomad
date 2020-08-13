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
The archive API of the nomad@FAIRDI APIs. This API is about serving processed
(parsed and normalized) calculation data in nomad's *meta-info* format.
'''

from flask_restplus import abort, Resource
import importlib

from nomad.metainfo.legacy import python_package_mapping, LegacyMetainfoEnvironment
from nomad.metainfo import Package
from nomad.parsing.parsers import parsers

from .api import api


ns = api.namespace(
    'metainfo',
    description='Access the NOMAD Metainfo (i.e. the archive\'s schema/definitions).')


@ns.route('/')
class AllMetainfoResource(Resource):
    @api.doc('get_all_metainfo')
    @api.response(200, 'Metainfo send')
    def get(self):
        '''
        Returns all metainfo packages.
        '''
        # Ensure all metainfo is loaded
        for parser in parsers:
            _ = parser.metainfo_env

        return {
            key: value.m_to_dict()
            for key, value in Package.registry.items()}


@ns.route('/<string:metainfo_package_name>')
class MetainfoResource(Resource):
    @api.doc('get_metainfo')
    @api.response(404, 'Package (e.g. code, parser, converter) does not exist')
    @api.response(200, 'Metainfo send')
    def get(self, metainfo_package_name):
        '''
        Get a JSON representation of the NOMAD Metainfo.

        You can get the metainfo for 'common', and parser/code metainfo packages.
        Parser/code packages constain the necessary definitions that the respective
        parser/code might use. 'Common' contains all non specific general definitions.
        Other required packages might also be returned, e.g. a parser might organize its
        definitions in multiple packages.
        '''
        package = metainfo_package_name
        if package.endswith('.json'):
            package = package[:-5]

        try:
            try:
                python_module = importlib.import_module(package)
            except ImportError:
                python_package_name, _ = python_package_mapping(package)
                python_module = importlib.import_module(python_package_name)

            metainfo_package = getattr(python_module, 'm_package')
        except (ImportError, KeyError, FileNotFoundError, AttributeError):
            abort(404, message='Metainfo package %s does not exist.' % package)

        result = {}
        result[metainfo_package.name] = metainfo_package.m_to_dict()
        for dependency in metainfo_package.dependencies:
            result[dependency.name] = dependency.m_to_dict()

        return result


@ns.route('/legacy/<string:metainfo_package_name>')
class LegacyMetainfoResource(Resource):
    @api.doc('get_legacy_metainfo')
    @api.response(404, 'Package (e.g. code, parser, converter) does not exist')
    @api.response(200, 'Metainfo send')
    def get(self, metainfo_package_name):
        '''
        Get a JSON representation of the NOMAD Metainfo in its old legacy JSON format.

        You can get the metainfo for 'common', and parser/code metainfo packages.
        Parser/code packages constain the necessary definitions that the respective
        parser/code might use. 'Common' contains all non specific general definitions.
        Other required packages might also be returned, e.g. a parser might organize its
        definitions in multiple packages.
        '''
        try:
            metainfo = LegacyMetainfoEnvironment.from_legacy_package_path(metainfo_package_name)
        except (ImportError, KeyError, FileNotFoundError, AttributeError):
            abort(404, message='Metainfo package %s does not exist.' % metainfo_package_name)

        if isinstance(metainfo, LegacyMetainfoEnvironment):
            return metainfo.to_legacy_dict(metainfo.packages)
        else:
            abort(404, message='Metainfo package %s is not a legacy package.' % metainfo_package_name)
