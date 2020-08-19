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
The mirror API of the nomad@FAIRDI APIs. Allows to export upload metadata.
'''

from flask import request, send_file
from flask_restplus import Resource, abort, fields

from nomad import processing as proc
from nomad.datamodel import Dataset
from nomad.doi import DOI
from nomad.files import PublicUploadFiles
from nomad import infrastructure

from .api import api
from .auth import authenticate
from .common import upload_route

ns = api.namespace('mirror', description='Export upload (and all calc) metadata.')


mirror_upload_model = api.model('MirrorUpload', {
    'upload_id': fields.String(description='The id of the exported upload'),
    'upload': fields.String(description='The upload metadata as mongoengine json string'),
    'calcs': fields.List(fields.Raw, description='All upload calculation metadata as mongo SON'),
    'datasets': fields.Raw(description='All upload datasets as dict id->mongo SON'),
    'dois': fields.Raw(description='All upload dois as dict id->mongo SON'),
    'upload_files_path': fields.String(description='The path to the local uploads file folder')
})

mirror_query_model = api.model('MirrorQuery', {
    'query': fields.Raw(
        description='Mongoengine query that is used to search for uploads to mirror.')
})

_Dataset = Dataset.m_def.a_mongo.mongo_cls


@ns.route('/')
class MirrorUploadsResource(Resource):
    @api.doc('get_uploads_mirror')
    @api.response(400, 'Bad query')
    @api.marshal_with(
        mirror_upload_model, skip_none=True, code=200, as_list=True,
        description='Uploads exported')
    @api.expect(mirror_query_model)
    @authenticate(admin_only=True)
    def post(self):
        json_data = request.get_json()
        if json_data is None:
            query = {}
        else:
            query = json_data.get('query', {})

        try:
            return [
                {
                    'upload_id': upload.upload_id,
                    'upload_files_path': upload.upload_files.os_path
                }
                for upload in proc.Upload.objects(**query)], 200
        except Exception as e:
            abort(400, message='Could not query mongodb: %s' % str(e))


def _upload_data(upload_id, upload_json, calcs_col, datasets_col, dois_col):
    calcs = []
    datasets = {}
    dois = {}
    for calc in calcs_col.find(dict(upload_id=upload_id)):
        calcs.append(calc)
        for dataset in calc['metadata'].get('datasets', []):
            if dataset not in datasets:
                datasets[dataset] = datasets_col.find_one(dict(_id=dataset))
                doi = datasets[dataset].get('doi', None)
                if doi is not None:
                    doi_obj = dois_col.find_one(dict(_id=doi))
                    if doi_obj is not None:
                        dois[doi] = doi_obj

    return {
        'upload_id': upload_id,
        'upload': upload_json,
        'calcs': calcs,
        'datasets': datasets,
        'dois': dois
    }


@upload_route(ns)
class MirrorUploadResource(Resource):
    @api.response(400, 'Not available for the given upload, e.g. upload not published.')
    @api.response(404, 'The upload does not exist')
    @api.marshal_with(mirror_upload_model, skip_none=True, code=200, description='Upload exported')
    @api.doc('get_upload_mirror')
    @authenticate(admin_only=True)
    def get(self, upload_id):
        '''
        Export upload (and all calc) metadata for mirrors.
        '''
        try:
            upload = proc.Upload.get(upload_id)
        except KeyError:
            abort(404, message='Upload with id %s does not exist.' % upload_id)

        if upload.tasks_running or upload.process_running:
            abort(400, message='Only non processing uploads can be exported')

        upload_data = _upload_data(
            upload.upload_id,
            upload.to_json(),
            calcs_col=proc.Calc._get_collection(),
            datasets_col=_Dataset._get_collection(),
            dois_col=DOI._get_collection())
        upload_data.update(upload_files_path=upload.upload_files.os_path)
        return upload_data, 200


_mirror_files_parser = api.parser()
_mirror_files_parser.add_argument(
    'prefix', type='str', help='File to download archive or raw', location='args')


@upload_route(ns, '/files')
class MirrorFilesResource(Resource):
    @api.doc('download_files_mirror')
    @api.expect(_mirror_files_parser, validate=True)
    @api.response(400, 'Invalid requests, e.g. wrong owner type or bad search parameters')
    @api.response(404, 'The upload or calculation does not exist')
    @authenticate(admin_only=True)
    def get(self, upload_id):
        '''
        Download archive and raw files for mirrors
        '''
        try:
            args = request.args
            prefix = args.get('prefix')
            assert prefix in ('archive', 'raw')

        except Exception:
            abort(400, message='bad parameter types')

        try:
            upload_files = PublicUploadFiles(upload_id)

            if prefix == 'raw':
                ext = 'plain'
                ending = 'zip'

            elif prefix == 'archive':
                ext = 'msg'
                ending = 'msg'

            elif prefix == 'legacy-archive':
                ext = 'json'
                ending = 'zip'

            else:
                abort(400, message='Unsupported prefix.')

            fileobj = upload_files._file_object(prefix, 'public', ext, ending)
            if not fileobj.exists():
                raise KeyError

            return send_file(
                open(fileobj.os_path, 'rb'),
                mimetype='application/zip',
                as_attachment=True,
                cache_timeout=0,
                attachment_filename=fileobj.os_path)

        except KeyError:
            abort(404, message='Upload %d does not exist' % upload_id)


@ns.route('/users')
class MirrorUsersResource(Resource):
    @api.doc('downlod_users_mirror')
    @api.response(400, 'Unsuccessful userlist query')
    @authenticate(admin_only=True)
    def get(self):
        '''
        Download user list for mirrors
        '''
        try:
            users = infrastructure.keycloak.search_user(query='')
            result = dict()
            for user in users:
                credentials = user.m_to_dict()
                credentials.pop('email', None)
                result[user.username] = credentials

            return result, 200

        except Exception:
            abort(400, message='Failed to fetch users')
