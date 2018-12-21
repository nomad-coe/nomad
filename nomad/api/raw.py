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
from flask_restplus import abort
from werkzeug.exceptions import HTTPException

from nomad.files import RepositoryFile
from nomad.utils import get_logger

from .app import app, base_path


def fix_file_paths(path):
    """ Removed the leading data from file paths that where given in mainfile uris. """
    # TODO, mainfile URI's should change or this implementation should change
    return path[5:]


@app.route('%s/raw/<string:upload_hash>/<path:upload_filepath>' % base_path, methods=['GET'])
def get_raw_file(upload_hash, upload_filepath):
    """
    Get a single raw calculation file from a given upload (or many files via wildcard).

    .. :quickref: raw; Get single raw calculation file.

    **Example request**:

    .. sourcecode:: http

        GET /nomad/api/raw/W36aqCzAKxOCfIiMFsBJh3nHPb4a/Si/si.out HTTP/1.1
        Accept: application/gz

    :param string upload_hash: the hash based identifier of the upload
    :param path upload_filepath: the path to the desired file within the upload;
        can also contain a wildcard * at the end to denote all files with path as prefix
    :qparam compress: any value to use compression for wildcard downloads, default is no compression
    :resheader Content-Type: application/gz
    :status 200: calc raw data successfully retrieved
    :status 404: upload with given hash does not exist or the given file does not exist
    :returns: the gzipped raw data in the body or a zip file when wildcard was used
    """
    upload_filepath = fix_file_paths(upload_filepath)

    repository_file = RepositoryFile(upload_hash)
    if not repository_file.exists():
        abort(404, message='The upload with hash %s does not exist.' % upload_hash)

    if upload_filepath[-1:] == '*':
        upload_filepath = upload_filepath[0:-1]
        files = list(
            file for file in repository_file.manifest
            if file.startswith(upload_filepath))
        if len(files) == 0:
            abort(404, message='There are no files for %s.' % upload_filepath)
        else:
            compress = request.args.get('compress', None) is not None
            return respond_to_get_raw_files(upload_hash, files, compress)

    try:
        the_file = repository_file.get_file(upload_filepath)
        with the_file.open() as f:
            rv = send_file(
                f,
                mimetype='application/octet-stream',
                as_attachment=True,
                attachment_filename=os.path.basename(upload_filepath))
            return rv
    except KeyError:
        files = list(file for file in repository_file.manifest if file.startswith(upload_filepath))
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


@app.route('%s/raw/<string:upload_hash>' % base_path, methods=['GET'])
def get_raw_files(upload_hash):
    """
    Get multiple raw calculation files.

    .. :quickref: raw; Get multiple raw calculation files.

    **Example request**:

    .. sourcecode:: http

        GET /nomad/api/raw/W36aqCzAKxOCfIiMFsBJh3nHPb4a?files=Si/si.out,Si/aux.txt HTTP/1.1
        Accept: application/gz

    :param string upload_hash: the hash based identifier of the upload
    :qparam string files: a comma separated list of file path
    :qparam compress: any value to use compression, default is no compression
    :resheader Content-Type: application/json
    :status 200: calc raw data successfully retrieved
    :status 404: calc with given hash does not exist or one of the given files does not exist
    :returns: a streamed .zip archive with the raw data
    """
    files_str = request.args.get('files', None)
    compress = request.args.get('compress', None) is not None

    if files_str is None:
        abort(400, message="No files argument given.")
    files = [fix_file_paths(file.strip()) for file in files_str.split(',')]

    return respond_to_get_raw_files(upload_hash, files, compress)


@app.route('%s/raw/<string:upload_hash>' % base_path, methods=['POST'])
def get_raw_files_post(upload_hash):
    """
    Get multiple raw calculation files.

    .. :quickref: raw; Get multiple raw calculation files.

    **Example request**:

    .. sourcecode:: http

        POST /nomad/api/raw/W36aqCzAKxOCfIiMFsBJh3nHPb4a HTTP/1.1
        Accept: application/gz
        Content-Type: application/json

        {
            "files": ["Si/si.out", "Si/aux.txt"]
        }

    :param string upload_hash: the hash based identifier of the upload
    :jsonparam files: a comma separated list of file paths
    :jsonparam compress: boolean to enable compression (true), default is not compression (false)
    :resheader Content-Type: application/json
    :status 200: calc raw data successfully retrieved
    :status 404: calc with given hash does not exist or one of the given files does not exist
    :returns: a streamed .zip archive with the raw data
    """
    json_data = request.get_json()
    if json_data is None:
        json_data = {}

    if 'files' not in json_data:
        abort(400, message='No files given, use key "files" in json body to provide file paths.')
    compress = json_data.get('compress', False)
    if not isinstance(compress, bool):
        abort(400, message='Compress value %s is not a bool.' % str(compress))
    files = [fix_file_paths(file.strip()) for file in json_data['files']]

    return respond_to_get_raw_files(upload_hash, files, compress)


def respond_to_get_raw_files(upload_hash, files, compress=False):
    logger = get_logger(__name__, endpoint='raw', action='get files', upload_hash=upload_hash)

    repository_file = RepositoryFile(upload_hash)
    if not repository_file.exists():
        abort(404, message='The upload with hash %s does not exist.' % upload_hash)

    def generator():
        """ Stream a zip file with all files using zipstream. """
        def iterator():
            """ Replace the directory based iter of zipstream with an iter over all given files. """
            try:
                with repository_file.zipped_container.zip_file() as zf:
                    for filename in files:
                        # Write a file to the zipstream.
                        try:
                            with zf.open(repository_file.zipped_container.get_zip_path(filename)) as f:
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
