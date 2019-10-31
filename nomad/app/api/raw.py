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

from typing import IO, Any, Union, Iterable, Tuple, Set, List
import os.path
import zipstream
from flask import Response, request, send_file, stream_with_context
from flask_restplus import abort, Resource, fields
import magic
import sys
import contextlib
import fnmatch
import json

from nomad import search, utils
from nomad.files import UploadFiles, Restricted
from nomad.processing import Calc

from .api import api
from .auth import login_if_available, create_authorization_predicate, \
    signature_token_argument, with_signature_token
from .repo import search_request_parser, add_query

if sys.version_info >= (3, 7):
    import zipfile
else:
    import zipfile37 as zipfile


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
raw_file_from_path_parser.add_argument(
    name='length', type=int, help='Download only x bytes from the given file.',
    location='args')
raw_file_from_path_parser.add_argument(
    name='offset', type=int, help='Start downloading a file\' content from the given offset.',
    location='args')


class FileView:
    """
    File-like wrapper that restricts the contents to a portion of the file.
    Arguments:
        f: the file-like
        offset: the offset
        length: the amount of bytes
    """
    def __init__(self, f, offset, length):
        self.f = f
        self.f_offset = offset
        self.offset = 0
        self.length = length

    def seek(self, offset, whence=0):
        if whence == os.SEEK_SET:
            self.offset = offset
        elif whence == os.SEEK_CUR:
            self.offset += offset
        elif whence == os.SEEK_END:
            self.offset = self.length + offset
        else:
            # Other values of whence should raise an IOError
            return self.f.seek(offset, whence)
        return self.f.seek(self.offset + self.f_offset, os.SEEK_SET)

    def tell(self):
        return self.offset

    def read(self, size=-1):
        self.seek(self.offset)
        if size < 0:
            size = self.length - self.offset
        size = max(0, min(size, self.length - self.offset))
        self.offset += size
        return self.f.read(size)


def get_raw_file_from_upload_path(
        upload_files, upload_filepath, authorization_predicate, mainfile: str = None):
    """
    Helper method used by func:`RawFileFromUploadPathResource.get` and
    func:`RawFileFromCalcPathResource.get`.
    """
    upload_filepath = upload_filepath.rstrip('/')

    if upload_filepath[-1:] == '*':
        upload_filepath = upload_filepath[0:-1]
        wildcarded_files = list(upload_files.raw_file_manifest(path_prefix=upload_filepath))
        if len(wildcarded_files) == 0:
            abort(404, message='There are no files for %s.' % upload_filepath)
        else:
            compress = request.args.get('compress', None) is not None
            return respond_to_get_raw_files(upload_files.upload_id, wildcarded_files, compress)

    try:
        with upload_files.raw_file(upload_filepath, 'br') as raw_file:
            buffer = raw_file.read(2048)
        mime_type = magic.from_buffer(buffer, mime=True)

        try:
            offset = int(request.args.get('offset', 0))
            length = int(request.args.get('length', -1))
        except Exception:
            abort(400, message='bad parameter types')

        if offset < 0:
            abort(400, message='bad offset, length values')
        if offset > 0 and length <= 0:
            abort(400, message='bad offset, length values')

        raw_file = upload_files.raw_file(upload_filepath, 'br')
        raw_file_view: Union[FileView, IO[Any]] = None
        if length > 0:
            raw_file_view = FileView(raw_file, offset, length)
        else:
            raw_file_view = raw_file

        return send_file(
            raw_file_view,
            mimetype=mime_type,
            as_attachment=True,
            cache_timeout=0,
            attachment_filename=os.path.basename(upload_filepath))
    except Restricted:
        abort(401, message='Not authorized to access all files in %s.' % upload_files.upload_id)
    except KeyError:
        directory_files = upload_files.raw_file_list(upload_filepath)
        if len(directory_files) == 0:
            abort(404, message='There is nothing to be found at %s.' % upload_filepath)

        contents = sorted([dict(name=name, size=size) for name, size in directory_files], key=lambda x: '' if x['name'] == mainfile else x['name'])
        # if mainfile is not None:
        #     contents = [mainfile] + [content for content in contents if content['name'] != mainfile]
        return {
            'upload_id': upload_files.upload_id,
            'directory': upload_filepath,
            'contents': contents
        }, 200


@ns.route('/<string:upload_id>/<path:path>')
@api.doc(params={
    'upload_id': 'The unique id for the requested upload.',
    'path': 'The path to a file or directory with optional wildcard.'
})
class RawFileFromUploadPathResource(Resource):
    @api.doc('get')
    @api.response(404, 'The upload or path does not exist')
    @api.response(401, 'Not authorized to access the requested files.')
    @api.response(200, 'File(s) send')
    @api.expect(raw_file_from_path_parser, validate=True)
    @login_if_available
    @with_signature_token
    def get(self, upload_id: str, path: str):
        """ Get a single raw calculation file, directory contents, or whole directory sub-tree
        from a given upload.

        The 'upload_id' parameter needs to identify an existing upload.

        If the upload
        is not yet published or contains requested data with embargo, proper authentication
        is required. This can be done via HTTP headers as usual. But, if you need to
        access files via plain URLs (e.g. for curl, download link, etc.), URLs for
        this endpoint can be token signed (see also /auth/token). For unpublished
        uploads, authentication is required regardless. For (partially) embargoed data,
        multi file downloads work, but will not contain any embargoed data.

        If the given path points to a file, the file is provided with the appropriate
        Content-Type header. A 401 is returned for staging, embargo files with unsigned
        or wrongly signed URLs. When accessing a file, the additional query parameters 'length'
        and 'offset' can be used to partially download a file's content.

        If the given path points to a directory, the content (names, sizes, type) is returned
        as a json body. Only visible items (depending on authenticated user, token) are
        returned.

        If the given path ends with the '*' wildcard character, all upload contents that
        match the given path at the start, will be returned as a .zip file body.
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

        return get_raw_file_from_upload_path(upload_files, upload_filepath, authorization_predicate)


@ns.route('/calc/<string:upload_id>/<string:calc_id>/<path:path>')
@api.doc(params={
    'upload_id': 'The unique id for the requested calc\'s upload.',
    'calc_id': 'The unique calc id for the requested calc',
    'path': 'The path to a file or directory with optional wildcard.'
})
class RawFileFromCalcPathResource(Resource):
    @api.doc('get_file_from_calc')
    @api.response(404, 'The upload or path does not exist')
    @api.response(401, 'Not authorized to access the requested files.')
    @api.response(200, 'File(s) send')
    @api.expect(raw_file_from_path_parser, validate=True)
    @login_if_available
    @with_signature_token
    def get(self, upload_id: str, calc_id: str, path: str):
        """ Get a single raw calculation file, calculation contents, or all files for a
        given calculation.

        The 'upload_id' parameter needs to identify an existing upload.
        The 'calc_id' parameter needs to identify a calculation within in the upload.

        This endpoint behaves exactly like /raw/<upload_id>/<path>, but the path is
        now relative to the calculation and not the upload.
        """
        calc_filepath = path if path is not None else ''
        authorization_predicate = create_authorization_predicate(upload_id)
        upload_files = UploadFiles.get(upload_id, authorization_predicate)
        if upload_files is None:
            abort(404, message='The upload with id %s does not exist.' % upload_id)

        calc = Calc.get(calc_id)
        if calc is None:
            abort(404, message='The calc with id %s does not exist.' % calc_id)
        if calc.upload_id != upload_id:
            abort(404, message='The calc with id %s is not part of the upload with id %s.' % (calc_id, upload_id))

        upload_filepath = os.path.join(os.path.dirname(calc.mainfile), calc_filepath)
        return get_raw_file_from_upload_path(
            upload_files, upload_filepath, authorization_predicate,
            mainfile=os.path.basename(calc.mainfile))


@ns.route('/calc/<string:upload_id>/<string:calc_id>/')
class RawFileFromCalcEmptyPathResource(RawFileFromCalcPathResource):
    @api.doc('get_file_list_from_calc')
    @api.response(404, 'The upload or path does not exist')
    @api.response(401, 'Not authorized to access the requested files.')
    @api.response(200, 'File(s) send')
    @api.expect(raw_file_from_path_parser, validate=True)
    @login_if_available
    @with_signature_token
    def get(self, upload_id: str, calc_id: str):
        """ Get calculation contents.

        This is basically /raw/calc/<upload_id>/<calc_id>/<path> with an empty path, since
        having an empty path parameter is not possible.
        """
        return super().get(upload_id, calc_id, None)


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
        """ Download multiple raw calculation files in a .zip file.

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


raw_file_from_query_parser = search_request_parser.copy()
raw_file_from_query_parser.add_argument(
    name='compress', type=bool, help='Use compression on .zip files, default is not.',
    location='args')
raw_file_from_query_parser.add_argument(
    name='strip', type=bool, help='Removes a potential common path prefix from all file paths.',
    location='args')
raw_file_from_query_parser.add_argument(
    name='file_pattern', type=str,
    help=(
        'A wildcard pattern. Only filenames that match this pattern will be in the '
        'download. Multiple patterns will be combined with logical or'),
    location='args', action='append')


@ns.route('/query')
class RawFileQueryResource(Resource):
    manifest_quantities = ['upload_id', 'calc_id', 'external_id', 'raw_id', 'pid', 'calc_hash']

    @api.doc('raw_files_from_query')
    @api.response(400, 'Invalid requests, e.g. wrong owner type or bad search parameters')
    @api.expect(raw_file_from_query_parser, validate=True)
    @api.response(200, 'File(s) send', headers={'Content-Type': 'application/gz'})
    @login_if_available
    def get(self):
        """ Download a .zip file with all raw-files for all entries that match the given
        search parameters.

        See ``/repo`` endpoint for documentation on the search
        parameters.

        Zip files are streamed; instead of 401 errors, the zip file will just not contain
        any files that the user is not authorized to access.

        The zip file will contain a ``manifest.json`` with the repository meta data.
        """
        patterns: List[str] = None
        try:
            args = raw_file_from_query_parser.parse_args()
            compress = args.get('compress', False)
            strip = args.get('strip', False)
            pattern = args.get('file_pattern', None)
            if isinstance(pattern, str):
                patterns = [pattern]
            elif pattern is None:
                patterns = []
            else:
                patterns = pattern
        except Exception:
            abort(400, message='bad parameter types')

        search_request = search.SearchRequest()
        add_query(search_request, search_request_parser)

        def path(entry):
            return '%s/%s' % (entry['upload_id'], entry['mainfile'])

        calcs = sorted(
            [entry for entry in search_request.execute_scan()],
            key=lambda x: x['upload_id'])

        paths = [path(entry) for entry in calcs]
        if strip:
            common_prefix_len = len(utils.common_prefix(paths))
        else:
            common_prefix_len = 0

        def generator():
            for entry in calcs:
                upload_id = entry['upload_id']
                mainfile = entry['mainfile']
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
                    filenames = upload_files.raw_file_manifest(
                        path_prefix=os.path.dirname(mainfile))
                    for filename in filenames:
                        filename_w_upload = os.path.join(upload_files.upload_id, filename)
                        filename_wo_prefix = filename_w_upload[common_prefix_len:]
                        if len(patterns) == 0 or any(
                                fnmatch.fnmatchcase(os.path.basename(filename_wo_prefix), pattern)
                                for pattern in patterns):

                            yield filename_wo_prefix, filename, upload_files

        try:
            manifest = {
                path(entry): {
                    key: entry[key]
                    for key in RawFileQueryResource.manifest_quantities
                    if entry.get(key) is not None
                }
                for entry in calcs
            }
            manifest_contents = json.dumps(manifest)
        except Exception as e:
            manifest_contents = dict(error='Could not create the manifest: %s' % (e))
            utils.get_logger(__name__).error(
                'could not create raw query manifest', exc_info=e)

        return _streamed_zipfile(
            generator(), zipfile_name='nomad_raw_files.zip', compress=compress,
            manifest=manifest_contents)


def respond_to_get_raw_files(upload_id, files, compress=False):
    upload_files = UploadFiles.get(
        upload_id, create_authorization_predicate(upload_id))
    if upload_files is None:
        abort(404, message='The upload with id %s does not exist.' % upload_id)

    # the zipfile cache allows to access many raw-files from public upload files without
    # having to re-open the underlying zip files all the time
    if hasattr(upload_files, 'zipfile_cache'):
        zipfile_cache = upload_files.zipfile_cache()
    else:
        zipfile_cache = contextlib.suppress()

    with zipfile_cache:
        return _streamed_zipfile(
            [(filename, filename, upload_files) for filename in files],
            zipfile_name='%s.zip' % upload_id, compress=compress)


def _streamed_zipfile(
        files: Iterable[Tuple[str, str, UploadFiles]], zipfile_name: str,
        compress: bool = False, manifest: str = None):
    """
    Creates a response that streams the given files as a streamed zip file. Ensures that
    each given file is only streamed once, based on its filename in the resulting zipfile.

    Arguments:
        files: An iterable of tuples with the filename to be used in the resulting zipfile,
            the filename within the upload, the :class:`UploadFiles` that contains
            the file.
        zipfile_name: A name that will be used in the content disposition attachment
            used as an HTTP respone.
        compress: Uses compression. Default is stored only.
        manifest: Add a ``manifest.json`` with the given content.
    """

    streamed_files: Set[str] = set()

    def generator():
        """ Stream a zip file with all files using zipstream. """
        def iterator():
            """
            Replace the directory based iter of zipstream with an iter over all given
            files.
            """
            # first the manifest
            if manifest is not None:
                yield dict(arcname='manifest.json', iterable=(manifest.encode('utf-8'),))

            # now the actual contents
            for zipped_filename, upload_filename, upload_files in files:
                if zipped_filename in streamed_files:
                    continue
                streamed_files.add(zipped_filename)

                # Write a file to the zipstream.
                try:
                    with upload_files.raw_file(upload_filename, 'rb') as f:
                        def iter_content():
                            while True:
                                data = f.read(1024 * 64)
                                if not data:
                                    break
                                yield data

                        yield dict(arcname=zipped_filename, iterable=iter_content())
                except KeyError:
                    # files that are not found, will not be returned
                    pass
                except Restricted:
                    # due to the streaming nature, we cannot raise 401 here
                    # we just leave it out in the download
                    pass

        compression = zipfile.ZIP_DEFLATED if compress else zipfile.ZIP_STORED
        zip_stream = zipstream.ZipFile(mode='w', compression=compression, allowZip64=True)
        zip_stream.paths_to_write = iterator()

        for chunk in zip_stream:
            yield chunk

    response = Response(stream_with_context(generator()), mimetype='application/zip')
    response.headers['Content-Disposition'] = 'attachment; filename={}'.format(zipfile_name)
    return response
