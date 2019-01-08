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

from nomad import config
from nomad.files import ArchiveFile, ArchiveLogFile
from nomad.utils import get_logger

from .app import api
from .auth import login_if_available
from .common import calc_route

ns = api.namespace(
    'archive',
    description='Access archive data and archive processing logs.')


@calc_route(ns, '/logs')
class ArchiveCalcLogResource(Resource):
    @api.doc('get_archive_logs')
    @api.response(404, 'The upload or calculation does not exist')
    @api.response(200, 'Archive data send', headers={'Content-Type': 'application/plain'})
    @login_if_available
    def get(self, upload_hash, calc_hash):
        """
        Get calculation processing log.

        Calcs are references via *upload_hash*, *calc_hash* pairs.
        """
        archive_id = '%s/%s' % (upload_hash, calc_hash)

        try:
            archive = ArchiveLogFile(archive_id)
            if not archive.exists():
                raise FileNotFoundError()

            archive_path = archive.os_path

            rv = send_file(
                archive_path,
                mimetype='text/plain',
                as_attachment=True,
                attachment_filename=os.path.basename(archive_path))

            return rv
        except FileNotFoundError:
            abort(404, message='Archive/calculation %s does not exist.' % archive_id)
        except Exception as e:
            logger = get_logger(
                __name__, endpoint='logs', action='get',
                upload_hash=upload_hash, calc_hash=calc_hash)
            logger.error('Exception on accessing calc proc log', exc_info=e)
            abort(500, message='Could not accessing the logs.')


@calc_route(ns)
class ArchiveCalcResource(Resource):
    @api.doc('get_archive_calc')
    @api.response(404, 'The upload or calculation does not exist')
    @api.response(200, 'Archive data send')
    @login_if_available
    def get(self, upload_hash, calc_hash):
        """
        Get calculation data in archive form.

        Calcs are references via *upload_hash*, *calc_hash* pairs.
        """
        archive_id = '%s/%s' % (upload_hash, calc_hash)

        try:
            archive = ArchiveFile(archive_id)
            if not archive.exists():
                raise FileNotFoundError()

            archive_path = archive.os_path

            rv = send_file(
                archive_path,
                mimetype='application/json',
                as_attachment=True,
                attachment_filename=os.path.basename(archive_path))

            if config.files.compress_archive:
                rv.headers['Content-Encoding'] = 'gzip'

            return rv
        except FileNotFoundError:
            abort(404, message='Archive %s does not exist.' % archive_id)
        except Exception as e:
            logger = get_logger(
                __name__, endpoint='archive', action='get',
                upload_hash=upload_hash, calc_hash=calc_hash)
            logger.error('Exception on accessing archive', exc_info=e)
            abort(500, message='Could not accessing the archive.')


@ns.route('/metainfo/<string:metainfo_path>')
@api.doc(params=dict(metainfo_path='A path or metainfo definition file name.'))
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
            meta_info_path = os.path.normpath(os.path.join(file_dir, metainfo_path.strip()))

            rv = send_file(
                meta_info_path,
                mimetype='application/json',
                as_attachment=True,
                attachment_filename=os.path.basename(metainfo_path))

            return rv
        except FileNotFoundError:
            abort(404, message='The metainfo %s does not exist.' % metainfo_path)
        except Exception as e:
            logger = get_logger(
                __name__, endpoint='metainfo', action='get', metainfo_path=metainfo_path)
            logger.error('Exception on accessing metainfo', exc_info=e)
            abort(500, message='Could not accessing the metainfo.')
