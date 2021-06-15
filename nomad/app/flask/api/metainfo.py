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

'''
The archive API of the nomad@FAIRDI APIs. This API is about serving processed
(parsed and normalized) calculation data in nomad's *meta-info* format.
'''

from flask_restplus import abort, Resource
import importlib

from nomad.metainfo import Package

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
        from nomad.parsing.parsers import parsers  # pylint: disable=unused-import
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
        Get a JSON representation of a NOMAD Metainfo package. The package name is
        the qualified Python name of the respective module that contains the definitions.
        Examples are `nomad.datamodel.metainfo.common_dft` or `vaspparser.metainfo`.
        If the desired package depends on other packages, these will also be contain in
        the results.
        '''
        package = metainfo_package_name

        try:
            python_module = importlib.import_module(package)
            metainfo_package = getattr(python_module, 'm_package')
        except (ImportError, KeyError, FileNotFoundError, AttributeError):
            abort(404, message='Metainfo package %s does not exist.' % package)

        result = {}
        result[metainfo_package.name] = metainfo_package.m_to_dict()
        for dependency in metainfo_package.dependencies:
            result[dependency.name] = dependency.m_to_dict()

        return result
