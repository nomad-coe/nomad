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
The raw API of the nomad@FAIRDI APIs. Can be used to retrieve raw calculation files.
"""

# TODO implement restrictions based on user, permissions, and upload/calc metadata

import os.path
from zipfile import ZIP_DEFLATED, ZIP_STORED

import zipstream
from flask import Response, request, send_file
from flask_restplus import abort, Resource, fields
from werkzeug.exceptions import HTTPException

from nomad.utils import get_logger
from nomad.uploads import UploadFiles

from .app import api
from .auth import login_if_available

ns = api.namespace('raw', description='Downloading raw data files.')


def fix_file_paths(path):
    """ Removed the leading data from file paths that where given in mainfile uris. """
    # TODO, mainfile URI's should change or this implementation should change
    return path[5:]


raw_file_compress_argument = dict(
    name='compress', type=bool, help='Use compression on .zip files, default is not.',
    location='args')
raw_file_from_path_parser = api.parser()
raw_file_from_path_parser.add_argument(**raw_file_compress_argument)


@ns.route('/<string:upload_hash>/<path:path>')
@api.doc(params={
    'upload_hash': 'The unique hash for the requested upload.',
    'path': 'The path to a file or directory.'
})
@api.header('Content-Type', 'application/gz')
class RawFileFromPathResource(Resource):
    @api.doc('get')
    @api.response(404, 'The upload or path does not exist')
    @api.response(200, 'File(s) send', headers={'Content-Type': 'application/gz'})
    @api.expect(raw_file_from_path_parser, validate=True)
    @login_if_available
    def get(self, upload_hash: str, path: str):
        """
        Get a single raw calculation file or whole directory from a given upload.

        If the given path points to a file, the file is provided. If the given path
        points to an directory, the directory and all contents is provided as .zip file.
        """
        upload_filepath = fix_file_paths(path)

        try:
            upload_files = UploadFiles.get(upload_hash)
        except KeyError:
            abort(404, message='The upload with hash %s does not exist.' % upload_hash)

        if upload_filepath[-1:] == '*':
            upload_filepath = upload_filepath[0:-1]
            files = list(upload_files.raw_file_manifest(path_prefix=upload_filepath))
            if len(files) == 0:
                abort(404, message='There are no files for %s.' % upload_filepath)
            else:
                compress = request.args.get('compress', None) is not None
                return respond_to_get_raw_files(upload_hash, files, compress)

        try:
            with upload_files.raw_file(upload_filepath) as f:
                rv = send_file(
                    f,
                    mimetype='application/octet-stream',
                    as_attachment=True,
                    attachment_filename=os.path.basename(upload_filepath))
                return rv
        except KeyError:
            files = list(file for file in upload_files.raw_file_manifest(upload_filepath))
            if len(files) == 0:
                abort(404, message='The file %s does not exist.' % upload_filepath)
            else:
                abort(404, message='The file %s does not exist, but there are files with matching paths' % upload_filepath, files=files)
        except HTTPException as e:
            raise e
        except Exception as e:
            logger = get_logger(
                __name__, endpoint='raw', action='get',
                upload_hash=upload_hash, upload_filepath=upload_filepath)
            logger.error('Exception on accessing raw data', exc_info=e)
            abort(500, message='Could not accessing the raw data.')


raw_files_request_model = api.model('RawFilesRequest', {
    'files': fields.List(
        fields.String, default=[], description='List of files to download.'),
    'compress': fields.Boolean(
        default=False,
        description='Enable compression, default is not compression.')
})

raw_files_request_parser = api.parser()
raw_files_request_parser.add_argument(**raw_file_compress_argument)
raw_files_request_parser.add_argument(
    'files', required=True, type=str, help='Comma separated list of files to download.', location='args')


@ns.route('/<string:upload_hash>')
@api.doc(params={
    'upload_hash': 'The unique hash for the requested upload.'
})
class RawFilesResource(Resource):
    @api.doc('get_files')
    @api.response(404, 'The upload or path does not exist')
    @api.response(200, 'File(s) send', headers={'Content-Type': 'application/gz'})
    @api.expect(raw_files_request_model, validate=True)
    @login_if_available
    def post(self, upload_hash):
        """ Download multiple raw calculation files. """
        json_data = request.get_json()
        compress = json_data.get('compress', False)
        files = [fix_file_paths(file.strip()) for file in json_data['files']]

        return respond_to_get_raw_files(upload_hash, files, compress)

    @api.doc('get_files_alternate')
    @api.response(404, 'The upload or path does not exist')
    @api.response(200, 'File(s) send', headers={'Content-Type': 'application/gz'})
    @api.expect(raw_files_request_parser, validate=True)
    @login_if_available
    def get(self, upload_hash):
        """ Download multiple raw calculation files. """
        files_str = request.args.get('files', None)
        compress = request.args.get('compress', 'false') == 'true'

        if files_str is None:
            abort(400, message="No files argument given.")
        files = [fix_file_paths(file.strip()) for file in files_str.split(',')]

        return respond_to_get_raw_files(upload_hash, files, compress)


def respond_to_get_raw_files(upload_hash, files, compress=False):
    logger = get_logger(__name__, endpoint='raw', action='get files', upload_hash=upload_hash)

    try:
        upload_file = UploadFiles.get(upload_hash)
    except KeyError:
        abort(404, message='The upload with hash %s does not exist.' % upload_hash)

    def generator():
        """ Stream a zip file with all files using zipstream. """
        def iterator():
            """ Replace the directory based iter of zipstream with an iter over all given files. """
            try:
                for filename in files:
                    # Write a file to the zipstream.
                    try:
                        with upload_file.raw_file(filename) as f:
                            def iter_content():
                                while True:
                                    data = f.read(100000)
                                    if not data:
                                        break
                                    yield data

                            yield dict(arcname=filename, iterable=iter_content())
                    except KeyError as e:
                        # files that are not found, will not be returned
                        pass

            except Exception as e:
                logger.error('Exception while accessing files.', exc_info=e)

        compression = ZIP_DEFLATED if compress else ZIP_STORED
        zip_stream = zipstream.ZipFile(mode='w', compression=compression, allowZip64=True)
        zip_stream.paths_to_write = iterator()

        for chunk in zip_stream:
            yield chunk

    response = Response(generator(), mimetype='application/zip')
    response.headers['Content-Disposition'] = 'attachment; filename={}'.format('%s.zip' % upload_hash)
    return response
