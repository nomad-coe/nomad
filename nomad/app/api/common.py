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
from typing import Callable, IO, Set, Tuple, Iterable, Dict, Any
from flask_restplus import fields
import zipstream
from flask import stream_with_context, Response, g, abort
from urllib.parse import urlencode

import sys
import os.path

from nomad import search, config
from nomad.app.optimade import filterparser
from nomad.app.utils import RFC3339DateTime, rfc3339DateTime
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

search_model = api.model('Search', {
    'pagination': fields.Nested(pagination_model, skip_none=True),
    'scroll': fields.Nested(allow_null=True, skip_none=True, model=api.model('Scroll', {
        'total': fields.Integer(description='The total amount of hits for the search.'),
        'scroll_id': fields.String(allow_null=True, description='The scroll_id that can be used to retrieve the next page.'),
        'size': fields.Integer(help='The size of the returned scroll page.')})),
    'results': fields.List(fields.Raw, description=(
        'A list of search results. Each result is a dict with quantitie names as key and '
        'values as values')),
})


def add_pagination_parameters(request_parser):
    """ Add pagination parameters to Flask querystring parser. """
    request_parser.add_argument(
        'page', type=int, help='The page, starting with 1.', location='args')
    request_parser.add_argument(
        'per_page', type=int, help='Desired calcs per page.', location='args')
    request_parser.add_argument(
        'order_by', type=str, help='The field to sort by.', location='args')
    request_parser.add_argument(
        'order', type=int, help='Use -1 for decending and 1 for acending order.', location='args')


request_parser = api.parser()
add_pagination_parameters(request_parser)
pagination_request_parser = request_parser.copy()


def add_scroll_parameters(request_parser):
    """ Add scroll parameters to Flask querystring parser. """
    request_parser.add_argument(
        'scroll', type=bool, help='Enable scrolling')
    request_parser.add_argument(
        'scroll_id', type=str, help='The id of the current scrolling window to use.')


def add_search_parameters(request_parser):
    """ Add search parameters to Flask querystring parser. """
    # more search parameters
    request_parser.add_argument(
        'owner', type=str,
        help='Specify which calcs to return: ``all``, ``public``, ``user``, ``staging``, default is ``all``')
    request_parser.add_argument(
        'from_time', type=lambda x: rfc3339DateTime.parse(x),
        help='A yyyy-MM-ddTHH:mm:ss (RFC3339) minimum entry time (e.g. upload time)')
    request_parser.add_argument(
        'until_time', type=lambda x: rfc3339DateTime.parse(x),
        help='A yyyy-MM-ddTHH:mm:ss (RFC3339) maximum entry time (e.g. upload time)')

    # main search parameters
    for quantity in search.quantities.values():
        request_parser.add_argument(
            quantity.name, help=quantity.description,
            action=quantity.argparse_action if quantity.multi else None)


def apply_search_parameters(search_request: search.SearchRequest, args: Dict[str, Any]):
    """
    Help that adds query relevant request args to the given SearchRequest.
    """
    args = {key: value for key, value in args.items() if value is not None}

    # owner
    owner = args.get('owner', 'all')
    try:
        search_request.owner(
            owner,
            g.user.user_id if g.user is not None else None)
    except ValueError as e:
        abort(401, getattr(e, 'message', 'Invalid owner parameter: %s' % owner))
    except Exception as e:
        abort(400, getattr(e, 'message', 'Invalid owner parameter'))

    # time range
    from_time_str = args.get('from_time', None)
    until_time_str = args.get('until_time', None)

    try:
        from_time = rfc3339DateTime.parse(from_time_str) if from_time_str is not None else None
        until_time = rfc3339DateTime.parse(until_time_str) if until_time_str is not None else None
        search_request.time_range(start=from_time, end=until_time)
    except Exception:
        abort(400, message='bad datetime format')

    # optimade
    try:
        optimade = args.get('optimade', None)
        if optimade is not None:
            q = filterparser.parse_filter(optimade)
            search_request.query(q)
    except filterparser.FilterException:
        abort(400, message='could not parse optimade query')

    # search parameter
    search_request.search_parameters(**{
        key: value for key, value in args.items()
        if key not in ['optimade'] and key in search.quantities})


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


def query_api_url(*args, query_string: Dict[str, Any] = None):
    """
    Creates a API URL.
    Arguments:
        *args: URL path segments after the API base URL
        query_string: A dict with query string parameters
    """
    url = os.path.join(config.api_url(False), *args)
    if query_string is not None:
        url = '%s?%s' % (url, urlencode(query_string))

    return url


def query_api_python(*args, **kwargs):
    """
    Creates a string of python code to execute a search query to the repository using
    the requests library.
    """
    url = query_api_url(*args, **kwargs)
    return '''import requests
response = requests.get("{}")
data = response.json()'''.format(url)


def query_api_curl(*args, **kwargs):
    """
    Creates a string of curl command to execute a search query to the repository.
    """
    url = query_api_url(*args, **kwargs)
    return 'curl -X GET %s -H  "accept: application/json" --output "nomad.json"' % url
