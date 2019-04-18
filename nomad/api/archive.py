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

import os.path

from flask import send_file
from flask_restplus import abort, Resource

import nomad_meta_info
from nomadcore.local_meta_info import load_metainfo

from nomad.files import UploadFiles, Restricted

from .app import api
from .auth import login_if_available, create_authorization_predicate, \
    signature_token_argument, with_signature_token
from .common import calc_route

ns = api.namespace(
    'archive',
    description='Access archive data and archive processing logs.')


archive_file_request_parser = api.parser()
archive_file_request_parser.add_argument(**signature_token_argument)


@calc_route(ns, '/logs')
class ArchiveCalcLogResource(Resource):
    @api.doc('get_archive_logs')
    @api.response(404, 'The upload or calculation does not exist')
    @api.response(401, 'Not authorized to access the data.')
    @api.response(200, 'Archive data send', headers={'Content-Type': 'application/plain'})
    @api.expect(archive_file_request_parser, validate=True)
    @login_if_available
    @with_signature_token
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
    @api.expect(archive_file_request_parser, validate=True)
    @login_if_available
    @with_signature_token
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
                attachment_filename='%s.json' % archive_id)
        except Restricted:
            abort(401, message='Not authorized to access %s/%s.' % (upload_id, calc_id))
        except KeyError:
            abort(404, message='Calculation %s does not exist.' % archive_id)


@ns.route('/metainfo/<string:metainfo_path>')
@api.doc(params=dict(nomad_metainfo_path='A path or metainfo definition file name.'))
class MetainfoResource(Resource):
    @api.doc('get_metainfo')
    @api.response(404, 'The metainfo does not exist')
    @api.response(200, 'Metainfo data send')
    def get(self, metainfo_path):
        """
        Get a metainfo definition file.
        """
        try:
            file_dir = os.path.dirname(os.path.abspath(nomad_meta_info.__file__))
            metainfo_file = os.path.normpath(os.path.join(file_dir, metainfo_path.strip()))

            meta_info, _ = load_metainfo(metainfo_file)
            return meta_info.toJsonList(False), 200
        except FileNotFoundError:
            abort(404, message='The metainfo %s does not exist.' % metainfo_path)
