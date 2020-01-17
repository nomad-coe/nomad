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
Common data, variables, decorators, models used throughout the API.
"""
from typing import Callable, IO, Set, Tuple, Iterable
from flask_restplus import fields
import zipstream
from flask import stream_with_context, Response
import sys

from nomad.app.utils import RFC3339DateTime
from nomad.files import Restricted

from .api import api


if sys.version_info >= (3, 7):
    import zipfile
else:
    import zipfile37 as zipfile


metadata_model = api.model('MetaData', {
    'with_embargo': fields.Boolean(default=False, description='Data with embargo is only visible to the upload until the embargo period ended.'),
    'comment': fields.String(description='The comment are shown in the repository for each calculation.'),
    'references': fields.List(fields.String, descriptions='References allow to link calculations to external source, e.g. URLs.'),
    'coauthors': fields.List(fields.String, description='A list of co-authors given by user_id.'),
    'shared_with': fields.List(fields.String, description='A list of users to share calculations with given by user_id.'),
    '_upload_time': RFC3339DateTime(description='Overrride the upload time.'),
    '_uploader': fields.String(description='Override the uploader with the given user id.'),
    'datasets': fields.List(fields.String, description='A list of dataset ids.')
})

pagination_model = api.model('Pagination', {
    'total': fields.Integer(description='Number of total elements.'),
    'page': fields.Integer(description='Number of the current page, starting with 0.'),
    'per_page': fields.Integer(description='Number of elements per page.')
})
""" Model used in responses with pagination. """


pagination_request_parser = api.parser()
""" Parser used for requests with pagination. """

pagination_request_parser.add_argument(
    'page', type=int, help='The page, starting with 1.', location='args')
pagination_request_parser.add_argument(
    'per_page', type=int, help='Desired calcs per page.', location='args')
pagination_request_parser.add_argument(
    'order_by', type=str, help='The field to sort by.', location='args')
pagination_request_parser.add_argument(
    'order', type=int, help='Use -1 for decending and 1 for acending order.', location='args')


def calc_route(ns, prefix: str = ''):
    """ A resource decorator for /<upload>/<calc> based routes. """
    def decorator(func):
        ns.route('%s/<string:upload_id>/<string:calc_id>' % prefix)(
            api.doc(params={
                'upload_id': 'The unique id for the requested upload.',
                'calc_id': 'The unique id for the requested calculation.'
            })(func)
        )
    return decorator


def upload_route(ns, prefix: str = ''):
    """ A resource decorator for /<upload> based routes. """
    def decorator(func):
        ns.route('%s/<string:upload_id>' % prefix)(
            api.doc(params={
                'upload_id': 'The unique id for the requested upload.'
            })(func)
        )
    return decorator


def streamed_zipfile(
        files: Iterable[Tuple[str, str, Callable[[str], IO], Callable[[str], int]]],
        zipfile_name: str, compress: bool = False):
    """
    Creates a response that streams the given files as a streamed zip file. Ensures that
    each given file is only streamed once, based on its filename in the resulting zipfile.

    Arguments:
        files: An iterable of tuples with the filename to be used in the resulting zipfile,
            an file id within the upload, a callable that gives an binary IO object for the
            file id, and a callable that gives the file size for the file id.
        zipfile_name: A name that will be used in the content disposition attachment
            used as an HTTP respone.
        compress: Uses compression. Default is stored only.
    """

    streamed_files: Set[str] = set()

    def generator():
        """ Stream a zip file with all files using zipstream. """
        def iterator():
            """
            Replace the directory based iter of zipstream with an iter over all given
            files.
            """
            # the actual contents
            for zipped_filename, file_id, open_io, file_size in files:
                if zipped_filename in streamed_files:
                    continue
                streamed_files.add(zipped_filename)

                # Write a file to the zipstream.
                try:
                    f = open_io(file_id)
                    try:
                        def iter_content():
                            while True:
                                data = f.read(1024 * 64)
                                if not data:
                                    break
                                yield data

                        yield dict(
                            arcname=zipped_filename, iterable=iter_content(),
                            buffer_size=file_size(file_id))
                    finally:
                        f.close()
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


def build_snippet(args, base_url):
    str_code = 'import requests\n'
    str_code += 'from urllib.parse import urlencode\n'
    str_code += '\n\n'
    str_code += 'def query_repository(args, base_url):\n'
    str_code += '    url = "%s?%s" % (base_url, urlencode(args))\n'
    str_code += '    response = requests.get(url)\n'
    str_code += '    if response.status_code != 200:\n'
    str_code += '        raise Exception("nomad return status %d" % response.status_code)\n'
    str_code += '    return response.json()\n'
    str_code += '\n\n'
    str_code += 'args = {'
    for key, val in args.items():
        if val is None:
            continue
        if isinstance(val, str):
            str_code += '"%s": "%s", ' % (key, val)
        else:
            str_code += '"%s": %s, ' % (key, val)
    str_code += '}\n'
    str_code += 'base_url = "%s"\n' % base_url
    str_code += 'JSON_DATA = query_repository(args, base_url)\n'

    return str_code
