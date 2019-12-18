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
The archive API of the nomad@FAIRDI APIs. This API is about serving processed
(parsed and normalized) calculation data in nomad's *meta-info* format.
"""

from typing import Dict, Any
from io import BytesIO
import os.path
from flask import send_file
from flask_restplus import abort, Resource
import json
import importlib
import contextlib

import nomad_meta_info

from nomad.files import UploadFiles, Restricted
from nomad import utils, search

from .auth import authenticate, create_authorization_predicate
from .api import api
from .repo import search_request_parser, add_query
from .common import calc_route, streamed_zipfile

ns = api.namespace(
    'archive',
    description='Access archive data and archive processing logs.')


@calc_route(ns, '/logs')
class ArchiveCalcLogResource(Resource):
    @api.doc('get_archive_logs')
    @api.response(404, 'The upload or calculation does not exist')
    @api.response(401, 'Not authorized to access the data.')
    @api.response(200, 'Archive data send', headers={'Content-Type': 'application/plain'})
    @authenticate(signature_token=True)
    def get(self, upload_id, calc_id):
        """
        Get calculation processing log.

        Calcs are references via *upload_id*, *calc_id* pairs.
        """
        archive_id = '%s/%s' % (upload_id, calc_id)

        upload_files = UploadFiles.get(
            upload_id, is_authorized=create_authorization_predicate(upload_id, calc_id))

        if upload_files is None:
            abort(404, message='Upload %s does not exist.' % upload_id)

        try:
            return send_file(
                upload_files.archive_log_file(calc_id, 'rb'),
                mimetype='text/plain',
                as_attachment=True,
                cache_timeout=0,
                attachment_filename='%s.log' % archive_id)
        except Restricted:
            abort(401, message='Not authorized to access %s/%s.' % (upload_id, calc_id))
        except KeyError:
            abort(404, message='Calculation %s does not exist.' % archive_id)


@calc_route(ns)
class ArchiveCalcResource(Resource):
    @api.doc('get_archive_calc')
    @api.response(404, 'The upload or calculation does not exist')
    @api.response(401, 'Not authorized to access the data.')
    @api.response(200, 'Archive data send')
    @authenticate(signature_token=True)
    def get(self, upload_id, calc_id):
        """
        Get calculation data in archive form.

        Calcs are references via *upload_id*, *calc_id* pairs.
        """
        archive_id = '%s/%s' % (upload_id, calc_id)

        upload_file = UploadFiles.get(
            upload_id, is_authorized=create_authorization_predicate(upload_id, calc_id))

        if upload_file is None:
            abort(404, message='Archive %s does not exist.' % upload_id)

        try:
            return send_file(
                upload_file.archive_file(calc_id, 'rb'),
                mimetype='application/json',
                as_attachment=True,
                cache_timeout=0,
                attachment_filename='%s.json' % archive_id)
        except Restricted:
            abort(401, message='Not authorized to access %s/%s.' % (upload_id, calc_id))
        except KeyError:
            abort(404, message='Calculation %s does not exist.' % archive_id)


archives_from_query_parser = search_request_parser.copy()
archives_from_query_parser.add_argument(
    name='compress', type=bool, help='Use compression on .zip files, default is not.',
    location='args')


@ns.route('/query')
class ArchiveQueryResource(Resource):
    manifest_quantities = ['upload_id', 'calc_id', 'external_id', 'raw_id', 'pid', 'calc_hash']

    @api.doc('archives_from_query')
    @api.response(400, 'Invalid requests, e.g. wrong owner type or bad search parameters')
    @api.expect(archives_from_query_parser, validate=True)
    @api.response(200, 'File(s) send', headers={'Content-Type': 'application/zip'})
    @authenticate(signature_token=True)
    def get(self):
        """
        Get calculation data in archive form from all query results.

        See ``/repo`` endpoint for documentation on the search
        parameters.

        Zip files are streamed; instead of 401 errors, the zip file will just not contain
        any files that the user is not authorized to access.

        The zip file will contain a ``manifest.json`` with the repository meta data.
        """
        try:
            args = archives_from_query_parser.parse_args()
            compress = args.get('compress', False)
        except Exception:
            abort(400, message='bad parameter types')

        search_request = search.SearchRequest()
        add_query(search_request, search_request_parser.parse_args())

        calcs = search_request.execute_scan(order_by='upload_id')

        def generator():
            manifest = {}
            for entry in calcs:
                upload_id = entry['upload_id']
                calc_id = entry['calc_id']
                upload_files = UploadFiles.get(
                    upload_id, create_authorization_predicate(upload_id))
                if upload_files is None:
                    utils.get_logger(__name__).error('upload files do not exist', upload_id=upload_id)
                    continue

                if hasattr(upload_files, 'zipfile_cache'):
                    zipfile_cache = upload_files.zipfile_cache()
                else:
                    zipfile_cache = contextlib.suppress()

                with zipfile_cache:
                    yield (
                        '%s.%s' % (calc_id, upload_files._archive_ext), calc_id,
                        lambda calc_id: upload_files.archive_file(calc_id, 'rb'),
                        lambda calc_id: upload_files.archive_file_size(calc_id))

                manifest[calc_id] = {
                    key: entry[key]
                    for key in ArchiveQueryResource.manifest_quantities
                    if entry.get(key) is not None
                }

            try:
                manifest_contents = json.dumps(manifest).encode('utf-8')
            except Exception as e:
                manifest_contents = json.dumps(
                    dict(error='Could not create the manifest: %s' % (e))).encode('utf-8')
                utils.get_logger(__name__).error(
                    'could not create raw query manifest', exc_info=e)

            yield (
                'manifest.json', 'manifest',
                lambda *args: BytesIO(manifest_contents),
                lambda *args: len(manifest_contents))

        return streamed_zipfile(
            generator(), zipfile_name='nomad_archive.zip', compress=compress)


@ns.route('/metainfo/<string:metainfo_package_name>')
@api.doc(params=dict(metainfo_package_name='The name of the metainfo package.'))
class MetainfoResource(Resource):
    @api.doc('get_metainfo')
    @api.response(404, 'The metainfo does not exist')
    @api.response(200, 'Metainfo data send')
    def get(self, metainfo_package_name):
        """
        Get a metainfo definition file.
        """
        try:
            return load_metainfo(metainfo_package_name), 200
        except FileNotFoundError:
            parser_prefix = metainfo_package_name[:-len('.nomadmetainfo.json')]

            try:
                return load_metainfo(dict(
                    parser='%sparser' % parser_prefix,
                    path='%s.nomadmetainfo.json' % parser_prefix)), 200
            except FileNotFoundError:
                abort(404, message='The metainfo %s does not exist.' % metainfo_package_name)


metainfo_main_path = os.path.dirname(os.path.abspath(nomad_meta_info.__file__))


def load_metainfo(
        package_name_or_dependency: str, dependency_source: str = None,
        loaded_packages: Dict[str, Any] = None) -> Dict[str, Any]:
    """
    Loads the given metainfo package and all its dependencies. Returns a dict with
    all loaded package_names and respective packages.

    Arguments:
        package_name_or_dependency: The name of the package, or a nomadmetainfo dependency object.
        dependency_source: The path of the metainfo that uses this function to load a relative dependency.
        loaded_packages: Give a dict and the function will added freshly loaded packages
            to it and return it.
    """
    if loaded_packages is None:
        loaded_packages = {}

    if isinstance(package_name_or_dependency, str):
        package_name = package_name_or_dependency
        metainfo_path = os.path.join(metainfo_main_path, package_name)
    else:
        dependency = package_name_or_dependency
        if 'relativePath' in dependency:
            if dependency_source is None:
                raise Exception(
                    'Can only load relative dependency from within another metainfo package')

            metainfo_path = os.path.join(
                os.path.dirname(dependency_source), dependency['relativePath'])

        elif 'metainfoPath' in dependency:
            metainfo_path = os.path.join(metainfo_main_path, dependency['metainfoPath'])

        elif 'parser' in dependency:
            parser = dependency['parser']
            path = dependency['path']
            try:
                parser_module = importlib.import_module(parser).__file__
            except Exception:
                raise Exception('Parser not installed %s for metainfo path %s' % (parser, metainfo_path))

            parser_directory = os.path.dirname(parser_module)
            metainfo_path = os.path.join(parser_directory, path)

        else:
            raise Exception('Invalid dependency type in metainfo package %s' % metainfo_path)

        package_name = os.path.basename(metainfo_path)

    package_name = os.path.basename(package_name)

    if package_name in loaded_packages:
        return loaded_packages

    with open(metainfo_path, 'rt') as f:
        metainfo_json = json.load(f)

    loaded_packages[package_name] = metainfo_json

    for dependency in metainfo_json.get('dependencies', []):
        load_metainfo(dependency, dependency_source=metainfo_path, loaded_packages=loaded_packages)

    return loaded_packages
