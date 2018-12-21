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
from flask_restplus import abort

from nomad import config
from nomad.files import ArchiveFile, ArchiveLogFile
from nomad.utils import get_logger

from .app import app, base_path


@app.route('%s/logs/<string:upload_hash>/<string:calc_hash>' % base_path, methods=['GET'])
def get_calc_proc_log(upload_hash, calc_hash):
    """
    Get calculation processing log. Calcs are references via *upload_hash*, *calc_hash*
    pairs.

    .. :quickref: archive; Get calculation processing logs.

    **Example request**:

    .. sourcecode:: http

        GET /nomad/api/logs/W36aqCzAKxOCfIiMFsBJh3nHPb4a/7ddvtfRfZAvc3Crr7jOJ8UH0T34I HTTP/1.1
        Accept: application/json

    :param string upload_hash: the hash of the upload (from uploaded file contents)
    :param string calc_hash: the hash of the calculation (from mainfile)
    :resheader Content-Type: application/json
    :status 200: calc successfully retrieved
    :status 404: calc with given hashes does not exist
    :returns: the log data, a line by line sequence of structured logs
    """
    archive_id = '%s/%s' % (upload_hash, calc_hash)

    try:
        archive = ArchiveLogFile(archive_id)
        if not archive.exists():
            raise FileNotFoundError()

        archive_path = archive.os_path

        rv = send_file(
            archive_path,
            mimetype='application/text',
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


@app.route('%s/archive/<string:upload_hash>/<string:calc_hash>' % base_path, methods=['GET'])
def get_calc(upload_hash, calc_hash):
    """
    Get calculation data in archive form. Calcs are references via *upload_hash*, *calc_hash*
    pairs.

    .. :quickref: archive; Get calculation data in archive form.

    **Example request**:

    .. sourcecode:: http

        GET /nomad/api/archive/W36aqCzAKxOCfIiMFsBJh3nHPb4a/7ddvtfRfZAvc3Crr7jOJ8UH0T34I HTTP/1.1
        Accept: application/json

    :param string upload_hash: the hash of the upload (from uploaded file contents)
    :param string calc_hash: the hash of the calculation (from mainfile)
    :resheader Content-Type: application/json
    :status 200: calc successfully retrieved
    :status 404: calc with given hashes does not exist
    :returns: the metainfo formated JSON data of the requested calculation
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
