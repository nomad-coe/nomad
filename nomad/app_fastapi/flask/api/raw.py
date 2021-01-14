#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

'''
The raw API of the nomad@FAIRDI APIs. Can be used to retrieve raw calculation files.
'''

from typing import IO, Any, Union, List
import os.path
from flask import request, send_file
from flask_restplus import abort, Resource, fields
import magic
import fnmatch
import gzip
import lzma
import urllib.parse

from nomad import search, utils, config
from nomad.files import UploadFiles, Restricted
from nomad.processing import Calc

from .. import common
from .api import api
from .auth import authenticate, create_authorization_predicate
from .common import streamed_zipfile, add_search_parameters, apply_search_parameters,\
    search_model


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
raw_file_strip_argument = dict(
    name='strip', type=bool, help='Removes a potential common path prefix from all file paths.',
    location='args')

_raw_file_from_path_parser = api.parser()
_raw_file_from_path_parser.add_argument(**raw_file_compress_argument)
_raw_file_from_path_parser.add_argument(**raw_file_strip_argument)
_raw_file_from_path_parser.add_argument(
    name='length', type=int, help='Download only x bytes from the given file.',
    location='args')
_raw_file_from_path_parser.add_argument(
    name='offset', type=int, help='Start downloading a file\' content from the given offset.',
    location='args')
_raw_file_from_path_parser.add_argument(
    name='decompress', type=int, help='Automatically decompress the file if compressed. Only supports .gz',
    location='args')


class FileView:
    '''
    File-like wrapper that restricts the contents to a portion of the file.
    Arguments:
        f: the file-like
        offset: the offset
        length: the amount of bytes
    '''
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
    '''
    Helper method used by func:`RawFileFromUploadPathResource.get` and
    func:`RawFileFromCalcPathResource.get`.
    '''
    upload_filepath = upload_filepath.rstrip('/')

    if upload_filepath[-1:] == '*':
        upload_filepath = upload_filepath[0:-1]
        wildcarded_files = list(upload_files.raw_file_manifest(path_prefix=upload_filepath))
        if len(wildcarded_files) == 0:
            abort(404, message='There are no files for %s.' % upload_filepath)
        else:
            compress = request.args.get('compress', None) is not None
            strip = request.args.get('strip', None) is not None
            return respond_to_get_raw_files(
                upload_files.upload_id, wildcarded_files, compress=compress, strip=strip)

    try:
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

        if request.args.get('decompress') is not None:
            if upload_filepath.endswith('.gz'):
                filename = upload_filepath[:3]
                upload_filepath = filename
                raw_file = gzip.GzipFile(filename=filename, mode='rb', fileobj=raw_file)
            if upload_filepath.endswith('.xz'):
                filename = upload_filepath[:3]
                upload_filepath = filename
                raw_file = lzma.open(filename=raw_file, mode='rb')

        buffer = raw_file.read(2048)
        raw_file.seek(0)
        mime_type = magic.from_buffer(buffer, mime=True)

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
    @api.expect(_raw_file_from_path_parser, validate=True)
    @authenticate(signature_token=True)
    def get(self, upload_id: str, path: str):
        ''' Get a single raw calculation file, directory contents, or whole directory sub-tree
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
        '''
        # TODO this is a quick fix, since swagger cannot deal with not encoded path parameters
        if path is not None:
            path = urllib.parse.unquote(path)

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
    @api.expect(_raw_file_from_path_parser, validate=True)
    @authenticate(signature_token=True)
    def get(self, upload_id: str, calc_id: str, path: str):
        ''' Get a single raw calculation file, calculation contents, or all files for a
        given calculation.

        The 'upload_id' parameter needs to identify an existing upload.
        The 'calc_id' parameter needs to identify a calculation within in the upload.

        This endpoint behaves exactly like /raw/<upload_id>/<path>, but the path is
        now relative to the calculation and not the upload.
        '''
        # TODO this is a quick fix, since swagger cannot deal with not encoded path parameters
        if path is not None:
            path = urllib.parse.unquote(path)

        calc_filepath = path if path is not None else ''
        authorization_predicate = create_authorization_predicate(upload_id, calc_id=calc_id)
        upload_files = UploadFiles.get(upload_id, authorization_predicate)
        if upload_files is None:
            abort(404, message='The upload with id %s does not exist.' % upload_id)

        try:
            calc = Calc.get(calc_id)
        except KeyError:
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
    @api.expect(_raw_file_from_path_parser, validate=True)
    @authenticate(signature_token=True)
    def get(self, upload_id: str, calc_id: str):
        ''' Get calculation contents.

        This is basically /raw/calc/<upload_id>/<calc_id>/<path> with an empty path, since
        having an empty path parameter is not possible.
        '''
        return super().get(upload_id, calc_id, None)


_raw_files_request_model = api.model('RawFilesRequest', {
    'files': fields.List(
        fields.String, default=[], description='List of files to download.'),
    'compress': fields.Boolean(
        default=False,
        description='Enable compression, default is not compression.')
})

_raw_files_request_parser = api.parser()
_raw_files_request_parser.add_argument(
    'files', required=True, type=str, help='Comma separated list of files to download.', location='args')
_raw_files_request_parser.add_argument(
    'prefix', type=str, help='A common prefix that is prepend to all files', location='args')
_raw_files_request_parser.add_argument(**raw_file_strip_argument)
_raw_files_request_parser.add_argument(**raw_file_compress_argument)


@ns.route('/<string:upload_id>')
@api.doc(params={
    'upload_id': 'The unique id for the requested upload.'
})
class RawFilesResource(Resource):
    @api.doc('get_files')
    @api.response(404, 'The upload or path does not exist')
    @api.response(200, 'File(s) send', headers={'Content-Type': 'application/zip'})
    @api.expect(_raw_files_request_model, validate=True)
    @authenticate()
    def post(self, upload_id):
        ''' Download multiple raw calculation files in a .zip file.

        Zip files are streamed; instead of 401 errors, the zip file will just not contain
        any files that the user is not authorized to access.
        '''
        json_data = request.get_json()
        compress = json_data.get('compress', False)
        files = [file.strip() for file in json_data['files']]

        return respond_to_get_raw_files(upload_id, files, compress)

    @api.doc('get_files_alternate')
    @api.response(404, 'The upload or path does not exist')
    @api.response(200, 'File(s) send', headers={'Content-Type': 'application/zip'})
    @api.expect(_raw_files_request_parser, validate=True)
    @authenticate(signature_token=True)
    def get(self, upload_id):
        '''
        Download multiple raw calculation files.
        Download multiple raw calculation files in a .zip file.
        Zip files are streamed; instead of 401 errors, the zip file will just not contain
        any files that the user is not authorized to access.
        '''
        args = _raw_files_request_parser.parse_args()

        files_str = args.get('files')
        if files_str is None:
            abort(400, message="No files argument given.")

        prefix = args.get('prefix')
        if prefix is not None:
            files = [os.path.join(prefix, file.strip()) for file in files_str.split(',')]
        else:
            files = [file.strip() for file in files_str.split(',')]

        compress = args.get('compress', False)
        strip = args.get('strip', False)

        return respond_to_get_raw_files(upload_id, files, compress=compress, strip=strip)


_raw_file_from_query_parser = api.parser()
add_search_parameters(_raw_file_from_query_parser)
_raw_file_from_query_parser.add_argument(
    name='compress', type=bool, help='Use compression on .zip files, default is not.',
    location='args')
_raw_file_from_query_parser.add_argument(**raw_file_strip_argument)
_raw_file_from_query_parser.add_argument(
    name='file_pattern', type=str,
    help=(
        'A wildcard pattern. Only filenames that match this pattern will be in the '
        'download. Multiple patterns will be combined with logical or'),
    location='args', action='append')


_raw_file_from_query_model = api.inherit('ArchiveSearch', search_model, {
    'compress': fields.Boolean(description='Use compression on .zip files, default is not.', default=False),
    'strip': fields.Boolean(description='Removes a potential common path prefix from all file paths.', default=False),
    'file_pattern': fields.List(fields.String, description=(
        'A wildcard pattern. Only filenames that match this pattern will be in the '
        'download. Multiple patterns will be combined with logical or'), allow_null=True, skip_none=True)
})


@ns.route('/query')
class RawFileQueryResource(Resource):
    manifest_quantities = ['upload_id', 'calc_id', 'external_id', 'raw_id', 'pid', 'calc_hash']

    @api.doc('raw_files_from_query')
    @api.response(400, 'Invalid requests, e.g. wrong owner type or bad search parameters')
    @api.expect(_raw_file_from_query_parser, validate=True)
    @api.response(200, 'File(s) send', headers={'Content-Type': 'application/zip'})
    @authenticate(signature_token=True)
    def get(self):
        ''' Download a .zip file with all raw-files for all entries that match the given
        search parameters.

        See ``/repo`` endpoint for documentation on the search
        parameters.

        Zip files are streamed; instead of 401 errors, the zip file will just not contain
        any files that the user is not authorized to access.

        The zip file will contain a ``manifest.json`` with the repository meta data.
        '''
        logger = common.logger.bind(query=urllib.parse.urlencode(request.args, doseq=True))
        try:
            args = _raw_file_from_query_parser.parse_args()
        except Exception:
            abort(400, message='could not parse request arguments')

        search_request = search.SearchRequest()
        apply_search_parameters(search_request, _raw_file_from_query_parser.parse_args())

        return respond_to_raw_files_query(search_request, args, logger)

    @api.doc('post_raw_files_from_query')
    @api.expect(_raw_file_from_query_model)
    @api.response(400, 'Invalid requests, e.g. wrong owner type or bad search parameters')
    @api.response(200, 'File(s) send', headers={'Content-Type': 'application/zip'})
    @authenticate(signature_token=True)
    def post(self):
        ''' Download a .zip file with all raw-files for all entries that match the given
        search parameters.

        See ``/repo`` endpoint for documentation on the search
        parameters.

        Zip files are streamed; instead of 401 errors, the zip file will just not contain
        any files that the user is not authorized to access.

        The zip file will contain a ``manifest.json`` with the repository meta data.
        '''
        try:
            post_data = request.get_json()
            query = post_data.get('query', {})
            query_expression = {key: val for key, val in query.items() if '$' in key}
        except Exception:
            abort(400, message='could not parse request body')

        logger = common.logger.bind(query=urllib.parse.urlencode(query, doseq=True))

        search_request = search.SearchRequest()
        apply_search_parameters(search_request, query)

        if query_expression:
            try:
                search_request.query_expression(query_expression)
            except AssertionError as e:
                abort(400, str(e))

        return respond_to_raw_files_query(search_request, post_data, logger)


def respond_to_raw_files_query(search_request, args, logger):
    patterns: List[str] = None
    try:
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

    search_request.include('calc_id', 'upload_id', 'mainfile')

    def path(entry):
        return '%s/%s' % (entry['upload_id'], entry['mainfile'])

    calcs = search_request.execute_scan(
        order_by='upload_id',
        size=config.services.download_scan_size,
        scroll=config.services.download_scan_timeout)

    if strip:
        if search_request.execute()['total'] > config.raw_file_strip_cutoff:
            abort(400, 'The requested download has to many files for using "strip".')
        calcs = list(calcs)
        paths = [path(entry) for entry in calcs]
        common_prefix_len = len(utils.common_prefix(paths))
    else:
        common_prefix_len = 0

    def generator():
        try:
            directories = set()
            upload_files = None
            streamed, skipped = 0, 0

            for entry in calcs:
                manifest = {
                    path(entry): {
                        key: entry[key]
                        for key in RawFileQueryResource.manifest_quantities
                        if entry.get(key) is not None
                    }
                }

                upload_id = entry['upload_id']
                mainfile = entry['mainfile']
                if upload_files is None or upload_files.upload_id != upload_id:
                    logger.info('opening next upload for raw file streaming', upload_id=upload_id)
                    if upload_files is not None:
                        upload_files.close()

                    upload_files = UploadFiles.get(upload_id)

                    if upload_files is None:
                        logger.error('upload files do not exist', upload_id=upload_id)
                        continue

                    def open_file(upload_filename):
                        return upload_files.raw_file(upload_filename, 'rb')

                upload_files._is_authorized = create_authorization_predicate(
                    upload_id=upload_id, calc_id=entry['calc_id'])
                directory = os.path.dirname(mainfile)
                directory_w_upload = os.path.join(upload_files.upload_id, directory)
                if directory_w_upload not in directories:
                    streamed += 1
                    directories.add(directory_w_upload)
                    for filename, file_size in upload_files.raw_file_list(directory=directory):
                        filename = os.path.join(directory, filename)
                        filename_w_upload = os.path.join(upload_files.upload_id, filename)
                        filename_wo_prefix = filename_w_upload[common_prefix_len:]
                        if len(patterns) == 0 or any(
                                fnmatch.fnmatchcase(os.path.basename(filename_wo_prefix), pattern)
                                for pattern in patterns):
                            yield (
                                filename_wo_prefix, filename, manifest, open_file,
                                lambda *args, **kwargs: file_size)
                else:
                    skipped += 1

                if (streamed + skipped) % 10000 == 0:
                    logger.info('streaming raw files', streamed=streamed, skipped=skipped)

            if upload_files is not None:
                upload_files.close()

        except Exception as e:
            logger.error('unexpected error while streaming raw data from query', exc_info=e)

    logger.info('start streaming raw files')
    return streamed_zipfile(
        generator(), zipfile_name='nomad_raw_files.zip', compress=compress, manifest=dict())


def respond_to_get_raw_files(upload_id, files, compress=False, strip=False):
    upload_files = UploadFiles.get(
        upload_id, create_authorization_predicate(upload_id))
    if upload_files is None:
        abort(404, message='The upload with id %s does not exist.' % upload_id)

    if strip:
        common_prefix_len = len(utils.common_prefix(files))
    else:
        common_prefix_len = 0

    try:
        return streamed_zipfile(
            [(
                filename[common_prefix_len:], filename, None,
                lambda upload_filename: upload_files.raw_file(upload_filename, 'rb'),
                lambda upload_filename: upload_files.raw_file_size(upload_filename)
            ) for filename in files],
            zipfile_name='%s.zip' % upload_id, compress=compress)
    finally:
        upload_files.close()
