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
The mirror API of the nomad@FAIRDI APIs. Allows to export upload metadata.
"""

from flask import request
from flask_restplus import Resource, abort, fields

from nomad import processing as proc

from .app import api
from .auth import authenticate
from .common import upload_route

ns = api.namespace('mirror', description='Export upload (and all calc) metadata.')


mirror_upload_model = api.model('MirrorUpload', {
    'upload_id': fields.String(description='The id of the exported upload'),
    'upload': fields.String(description='The upload metadata as mongoengine json string'),
    'calcs': fields.List(fields.String, description='All upload calculation metadata as mongoengine json strings'),
    'upload_files_path': fields.String(description='The path to the local uploads file folder')
})

mirror_query_model = api.model('MirrorQuery', {
    'query': fields.Raw(
        description='Mongoengine query that is used to search for uploads to mirror.')
})


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


@upload_route(ns)
class MirrorUploadResource(Resource):
    @api.response(400, 'Not available for the given upload, e.g. upload not published.')
    @api.response(404, 'The upload does not exist')
    @api.marshal_with(mirror_upload_model, skip_none=True, code=200, description='Upload exported')
    @api.doc('get_upload_mirror')
    @authenticate(admin_only=True)
    def get(self, upload_id):
        """
        Export upload (and all calc) metadata for mirrors.
        """
        try:
            upload = proc.Upload.get(upload_id)
        except KeyError:
            abort(404, message='Upload with id %s does not exist.' % upload_id)

        if not upload.published:
            abort(400, message='Only published uploads can be exported')

        return {
            'upload_id': upload.upload_id,
            'upload': upload.to_json(),
            'calcs': [calc.to_json() for calc in proc.Calc.objects(upload_id=upload.upload_id)],
            'upload_files_path': upload.upload_files.os_path
        }, 200
