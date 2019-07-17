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

import os.path
from zipfile import ZIP_DEFLATED, ZIP_STORED

import zipstream
from flask import Response, request, send_file, stream_with_context
from flask_restplus import abort, Resource, fields

from nomad.files import UploadFiles, Restricted

from .app import api
from .auth import login_if_available, create_authorization_predicate, \
    signature_token_argument, with_signature_token

ns = api.namespace('raw', description='Downloading raw data files.')

raw_file_list_model = api.model('RawFileList', {
    'upload_id': fields.String(description='The id of the upload.'),
    'directory': fields.String(description='The path to the directory in the upload.'),
    'contents': fields.List(
        fields.Nested(model=api.model('RawFileListContents', {
            'file': fields.String(description='The file name'),
            'size': fields.Integer(description='The file size in bytes')
        })))})

raw_file_compress_argument = dict(
    name='compress', type=bool, help='Use compression on .zip files, default is not.',
    location='args')
raw_file_from_path_parser = api.parser()
raw_file_from_path_parser.add_argument(**raw_file_compress_argument)
raw_file_from_path_parser.add_argument(**signature_token_argument)


@ns.route('/list/<string:upload_id>/<path:directory>')
@api.doc(params={
    'upload_id': 'The unique id for the requested upload.',
    'directory': 'The directory in the upload with the desired contents.'
})
@api.header('Content-Type', 'application/json')
class RawFileList(Resource):
    @api.doc('get')
    @api.response(404, 'The upload or path does not exist')
    @api.response(401, 'Not authorized to access the data.')
    @api.response(200, 'File(s) send', headers={'Content-Type': 'application/json'})
    @api.marshal_with(raw_file_list_model, skip_none=True, code=200, description='File list send')
    @login_if_available
    @with_signature_token
    def get(self, upload_id: str, directory: str):
        """
        Get the contents of the given directory for the given upload.

        If the path points to a file a single entry is returned. If the path
        points to a directory, information on all files in the directory are returned.
        """

        upload_files = UploadFiles.get(upload_id, create_authorization_predicate(upload_id))
        if upload_files is None:
            abort(404, message='The upload with id %s does not exist.' % upload_id)

        files = upload_files.raw_file_list(directory=directory)
        if len(files) == 0:
            abort(404, message='There are no files for %s.' % directory)
        else:
            return {
                'upload_id': upload_id,
                'directory': directory,
                'contents': [dict(file=file, size=size) for file, size in files]
            }


@ns.route('/<string:upload_id>/<path:path>')
@api.doc(params={
    'upload_id': 'The unique id for the requested upload.',
    'path': 'The path to a file or directory.'
})
@api.header('Content-Type', 'application/gz')
class RawFileFromPathResource(Resource):
    @api.doc('get')
    @api.response(404, 'The upload or path does not exist')
    @api.response(401, 'Not authorized to access the data.')
    @api.response(200, 'File(s) send', headers={'Content-Type': 'application/gz'})
    @api.expect(raw_file_from_path_parser, validate=True)
    @login_if_available
    @with_signature_token
    def get(self, upload_id: str, path: str):
        """
        Get a single raw calculation file or whole directory from a given upload.

        If the given path points to a file, the file is provided. If the given path
        points to an directory, the directory and all contents is provided as .zip file.
        Zip files are streamed; instead of 401 errors, the zip file will just not contain
        any files that the user is not authorized to access.
        """
        upload_filepath = path

        # TODO find a better way to all access to certain files
        if os.path.basename(path).endswith('.png'):
            def authorization_predicate(*args, **kwargs):
                return True
        else:
            authorization_predicate = create_authorization_predicate(upload_id)

        upload_files = UploadFiles.get(upload_id, authorization_predicate)
        if upload_files is None:
            abort(404, message='The upload with id %s does not exist.' % upload_id)

        if upload_filepath[-1:] == '*':
            upload_filepath = upload_filepath[0:-1]
            files = list(upload_files.raw_file_manifest(path_prefix=upload_filepath))
            if len(files) == 0:
                abort(404, message='There are no files for %s.' % upload_filepath)
            else:
                compress = request.args.get('compress', None) is not None
                return respond_to_get_raw_files(upload_id, files, compress)

        try:
            return send_file(
                upload_files.raw_file(upload_filepath, 'br'),
                mimetype='application/octet-stream',
                as_attachment=True,
                attachment_filename=os.path.basename(upload_filepath))
        except Restricted:
            abort(401, message='Not authorized to access upload %s.' % upload_id)
        except KeyError:
            files = list(file for file in upload_files.raw_file_manifest(upload_filepath))
            if len(files) == 0:
                abort(404, message='The file %s does not exist.' % upload_filepath)
            else:
                abort(404, message='The file %s does not exist, but there are files with matching paths' % upload_filepath, files=files)


raw_files_request_model = api.model('RawFilesRequest', {
    'files': fields.List(
        fields.String, default=[], description='List of files to download.'),
    'compress': fields.Boolean(
        default=False,
        description='Enable compression, default is not compression.')
})

raw_files_request_parser = api.parser()
raw_files_request_parser.add_argument(
    'files', required=True, type=str, help='Comma separated list of files to download.', location='args')
raw_files_request_parser.add_argument(**raw_file_compress_argument)
raw_file_from_path_parser.add_argument(**signature_token_argument)


@ns.route('/<string:upload_id>')
@api.doc(params={
    'upload_id': 'The unique id for the requested upload.'
})
class RawFilesResource(Resource):
    @api.doc('get_files')
    @api.response(404, 'The upload or path does not exist')
    @api.response(200, 'File(s) send', headers={'Content-Type': 'application/gz'})
    @api.expect(raw_files_request_model, validate=True)
    @login_if_available
    def post(self, upload_id):
        """
        Download multiple raw calculation files in a .zip file.
        Zip files are streamed; instead of 401 errors, the zip file will just not contain
        any files that the user is not authorized to access.
        """
        json_data = request.get_json()
        compress = json_data.get('compress', False)
        files = [file.strip() for file in json_data['files']]

        return respond_to_get_raw_files(upload_id, files, compress)

    @api.doc('get_files_alternate')
    @api.response(404, 'The upload or path does not exist')
    @api.response(200, 'File(s) send', headers={'Content-Type': 'application/gz'})
    @api.expect(raw_files_request_parser, validate=True)
    @login_if_available
    @with_signature_token
    def get(self, upload_id):
        """
        Download multiple raw calculation files.
        Download multiple raw calculation files in a .zip file.
        Zip files are streamed; instead of 401 errors, the zip file will just not contain
        any files that the user is not authorized to access.
        """
        files_str = request.args.get('files', None)
        compress = request.args.get('compress', 'false') == 'true'

        if files_str is None:
            abort(400, message="No files argument given.")
        files = [file.strip() for file in files_str.split(',')]

        return respond_to_get_raw_files(upload_id, files, compress)


def respond_to_get_raw_files(upload_id, files, compress=False):
    upload_files = UploadFiles.get(
        upload_id, create_authorization_predicate(upload_id))
    if upload_files is None:
        abort(404, message='The upload with id %s does not exist.' % upload_id)

    def generator():
        """ Stream a zip file with all files using zipstream. """
        def iterator():
            """ Replace the directory based iter of zipstream with an iter over all given files. """
            for filename in files:
                # Write a file to the zipstream.
                try:
                    with upload_files.raw_file(filename, 'rb') as f:
                        def iter_content():
                            while True:
                                data = f.read(100000)
                                if not data:
                                    break
                                yield data

                        yield dict(arcname=filename, iterable=iter_content())
                except KeyError:
                    # files that are not found, will not be returned
                    pass
                except Restricted:
                    # due to the streaming nature, we cannot raise 401 here
                    # we just leave it out in the download
                    pass

        compression = ZIP_DEFLATED if compress else ZIP_STORED
        zip_stream = zipstream.ZipFile(mode='w', compression=compression, allowZip64=True)
        zip_stream.paths_to_write = iterator()

        for chunk in zip_stream:
            yield chunk

    response = Response(stream_with_context(generator()), mimetype='application/zip')
    response.headers['Content-Disposition'] = 'attachment; filename={}'.format('%s.zip' % upload_id)
    return response
