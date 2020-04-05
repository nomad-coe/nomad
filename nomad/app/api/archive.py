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

from typing import Dict, Any
from io import BytesIO
import os.path
from flask import request, g
from flask_restplus import abort, Resource, fields
import json
import orjson
import importlib
import urllib.parse

import metainfo

from nomad.files import UploadFiles, Restricted
from nomad.archive import query_archive
from nomad import search, config
from nomad.app import common

from .auth import authenticate, create_authorization_predicate
from .api import api
from .common import calc_route, streamed_zipfile, search_model, add_search_parameters, apply_search_parameters, query_model


ns = api.namespace(
    'archive',
    description='Access archive data and archive processing logs.')


@calc_route(ns, '/logs')
class ArchiveCalcLogResource(Resource):
    @api.doc('get_archive_logs')
    @api.response(404, 'The upload or calculation does not exist')
    @api.response(401, 'Not authorized to access the data.')
    @api.response(200, 'Archive data send')
    @authenticate(signature_token=True)
    def get(self, upload_id, calc_id):
        '''
        Get calculation processing log.

        Calcs are references via *upload_id*, *calc_id* pairs.
        '''
        archive_id = '%s/%s' % (upload_id, calc_id)

        upload_files = UploadFiles.get(
            upload_id, is_authorized=create_authorization_predicate(upload_id, calc_id))

        if upload_files is None:
            abort(404, message='Upload %s does not exist.' % upload_id)

        try:
            with upload_files.read_archive(calc_id) as archive:
                return [entry.to_dict() for entry in archive[calc_id]['processing_logs']]

        except Restricted:
            abort(401, message='Not authorized to access %s/%s.' % (upload_id, calc_id))
        except KeyError:
            abort(404, message='Calculation %s does not exist.' % archive_id)


@calc_route(ns)
class ArchiveCalcResource(Resource):
    @api.doc('get_archive_calc')
    @api.response(404, 'The upload or calculation does not exist')
    @api.response(401, 'Not authorized to access the data.')
    @api.response(200, 'Archive data send', headers={'Content-Type': 'application/json'})
    @authenticate(signature_token=True)
    def get(self, upload_id, calc_id):
        '''
        Get calculation data in archive form.

        Calcs are references via *upload_id*, *calc_id* pairs.
        '''
        archive_id = '%s/%s' % (upload_id, calc_id)

        upload_files = UploadFiles.get(
            upload_id, is_authorized=create_authorization_predicate(upload_id, calc_id))

        if upload_files is None:
            abort(404, message='Archive %s does not exist.' % upload_id)

        try:
            with upload_files.read_archive(calc_id) as archive:
                return {
                    key: value
                    for key, value in archive[calc_id].to_dict().items()
                    if key != 'processing_logs'}

        except Restricted:
            abort(401, message='Not authorized to access %s/%s.' % (upload_id, calc_id))
        except KeyError:
            abort(404, message='Calculation %s does not exist.' % archive_id)


_archive_download_parser = api.parser()
add_search_parameters(_archive_download_parser)
_archive_download_parser.add_argument(
    name='compress', type=bool, help='Use compression on .zip files, default is not.',
    location='args')


@ns.route('/download')
class ArchiveDownloadResource(Resource):
    manifest_quantities = ['upload_id', 'calc_id', 'external_id', 'raw_id', 'pid', 'calc_hash']

    @api.doc('archive_download')
    @api.response(400, 'Invalid requests, e.g. wrong owner type or bad search parameters')
    @api.expect(_archive_download_parser, validate=True)
    @api.response(200, 'File(s) send', headers={'Content-Type': 'application/zip'})
    @authenticate(signature_token=True)
    def get(self):
        '''
        Get calculation data in archive form from all query results.

        See ``/repo`` endpoint for documentation on the search
        parameters.

        Zip files are streamed; instead of 401 errors, the zip file will just not contain
        any files that the user is not authorized to access.

        The zip file will contain a ``manifest.json`` with the repository meta data.
        '''
        try:
            args = _archive_download_parser.parse_args()
            compress = args.get('compress', False)
        except Exception:
            abort(400, message='bad parameter types')

        search_request = search.SearchRequest()
        apply_search_parameters(search_request, args)
        search_request.include('calc_id', 'upload_id', 'mainfile')

        calcs = search_request.execute_scan(
            order_by='upload_id',
            size=config.services.download_scan_size,
            scroll=config.services.download_scan_timeout)

        def generator():
            try:
                manifest = {}
                upload_files = None

                for entry in calcs:
                    upload_id = entry['upload_id']
                    calc_id = entry['calc_id']
                    if upload_files is None or upload_files.upload_id != upload_id:
                        if upload_files is not None:
                            upload_files.close()

                        upload_files = UploadFiles.get(
                            upload_id, create_authorization_predicate(upload_id))

                        if upload_files is None:
                            common.logger.error('upload files do not exist', upload_id=upload_id)
                            continue

                    with upload_files.read_archive(calc_id) as archive:
                        f = BytesIO(orjson.dumps(
                            archive[calc_id].to_dict(),
                            option=orjson.OPT_INDENT_2 | orjson.OPT_NON_STR_KEYS))

                        yield (
                            '%s.%s' % (calc_id, 'json'), calc_id,
                            lambda calc_id: f,
                            lambda calc_id: f.getbuffer().nbytes)

                    manifest[calc_id] = {
                        key: entry[key]
                        for key in ArchiveDownloadResource.manifest_quantities
                        if entry.get(key) is not None
                    }

                if upload_files is not None:
                    upload_files.close()

                try:
                    manifest_contents = json.dumps(manifest).encode('utf-8')
                except Exception as e:
                    manifest_contents = json.dumps(
                        dict(error='Could not create the manifest: %s' % (e))).encode('utf-8')
                    common.logger.error(
                        'could not create raw query manifest', exc_info=e)

                yield (
                    'manifest.json', 'manifest',
                    lambda *args: BytesIO(manifest_contents),
                    lambda *args: len(manifest_contents))

            except Exception as e:
                common.logger.warning(
                    'unexpected error while streaming raw data from query',
                    exc_info=e,
                    query=urllib.parse.urlencode(request.args, doseq=True))

        return streamed_zipfile(
            generator(), zipfile_name='nomad_archive.zip', compress=compress)


_archive_query_model = api.inherit('ArchiveSearch', search_model, {
    'query': fields.Nested(query_model, description='The query used to find the requested entries.', skip_none=True),
    'query_schema': fields.Raw(description='The query schema that defines what archive data to retrive.')
})


@ns.route('/query')
class ArchiveQueryResource(Resource):
    @api.doc('post_archive_query')
    @api.response(400, 'Invalid requests, e.g. wrong owner type or bad search parameters')
    @api.response(401, 'Not authorized to access the data.')
    @api.response(404, 'The upload or calculation does not exist')
    @api.expect(_archive_query_model)
    @api.marshal_with(_archive_query_model, skip_none=True, code=200, description='Archive search results sent')
    @authenticate()
    def post(self):
        '''
        Post a query schema and return it filled with archive data.

        See ``/repo`` endpoint for documentation on the search
        parameters.

        The actual data are in results and a supplementary python code (curl) to
        execute search is in python (curl).
        '''
        try:
            data_in = request.get_json()
            scroll = data_in.get('scroll', None)
            if scroll:
                scroll_id = scroll.get('scroll_id')
                scroll = True

            pagination = data_in.get('pagination', {})
            page = pagination.get('page', 1)
            per_page = pagination.get('per_page', 10 if not scroll else 1000)

            query = data_in.get('query', {})
            query_schema = data_in.get('query_schema', '*')

        except Exception:
            abort(400, message='bad parameter types')

        if not (page >= 1 and per_page > 0):
            abort(400, message='invalid pagination')

        search_request = search.SearchRequest()
        if g.user is not None:
            search_request.owner('all', user_id=g.user.user_id)
        else:
            search_request.owner('all')

        apply_search_parameters(search_request, query)
        search_request.include('calc_id', 'upload_id', 'with_embargo', 'parser_name')

        try:
            if scroll:
                results = search_request.execute_scrolled(
                    scroll_id=scroll_id, size=per_page, order_by='upload_id')
                results['scroll']['scroll'] = True

            else:
                results = search_request.execute_paginated(
                    per_page=per_page, page=page, order_by='upload_id')

        except search.ScrollIdNotFound:
            abort(400, 'The given scroll_id does not exist.')
        except KeyError as e:
            import traceback
            traceback.print_exc()
            abort(400, str(e))

        data = []
        calcs = results['results']
        upload_files = None
        current_upload_id = None
        for entry in calcs:
            upload_id = entry['upload_id']
            calc_id = entry['calc_id']
            if upload_files is None or current_upload_id != upload_id:
                if upload_files is not None:
                    upload_files.close()

                upload_files = UploadFiles.get(upload_id, create_authorization_predicate(upload_id))

                if upload_files is None:
                    return []

                current_upload_id = upload_id

            if entry['with_embargo']:
                access = 'restricted'
            else:
                access = 'public'

            try:
                with upload_files.read_archive(calc_id, access) as archive:
                    data.append({
                        'calc_id': calc_id,
                        'parser_name': entry['parser_name'],
                        'archive': query_archive(
                            archive, {calc_id: query_schema})[calc_id]
                    })

            except Restricted:
                # optimize and not access restricted for same upload again
                pass

        # assign archive data to results
        results['results'] = data

        return results, 200


@ns.route('/metainfo/<string:metainfo_package_name>')
@api.doc(params=dict(metainfo_package_name='The name of the metainfo package.'))
class MetainfoResource(Resource):
    @api.doc('get_metainfo')
    @api.response(404, 'The metainfo does not exist')
    @api.response(200, 'Metainfo data send')
    def get(self, metainfo_package_name):
        '''
        Get a metainfo definition file.
        '''
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


metainfo_main_path = os.path.dirname(os.path.abspath(metainfo.__file__))


def load_metainfo(
        package_name_or_dependency: str, dependency_source: str = None,
        loaded_packages: Dict[str, Any] = None) -> Dict[str, Any]:
    '''
    Loads the given metainfo package and all its dependencies. Returns a dict with
    all loaded package_names and respective packages.

    Arguments:
        package_name_or_dependency: The name of the package, or a nomadmetainfo dependency object.
        dependency_source: The path of the metainfo that uses this function to load a relative dependency.
        loaded_packages: Give a dict and the function will added freshly loaded packages
            to it and return it.
    '''
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
