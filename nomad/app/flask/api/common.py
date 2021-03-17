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
Common data, variables, decorators, models used throughout the API.
'''
from typing import Callable, IO, Set, Tuple, Iterable, Dict, Any
from flask_restplus import fields
from flask import request, make_response
import zipstream
from flask import stream_with_context, Response, g, abort
from urllib.parse import urlencode
import pprint
import io
import json

import sys
import os.path
import gzip
from functools import wraps

from nomad import search, config, datamodel, utils
from nomad.app.optimade import filterparser
from nomad.files import Restricted

from ..common import RFC3339DateTime, rfc3339DateTime
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
    'per_page': fields.Integer(description='Number of elements per page.'),
    'order_by': fields.String(description='Sorting criterion.'),
    'order': fields.Integer(description='Sorting order -1 for descending, 1 for asceding.')
})
''' Model used in responses with pagination. '''

scroll_model = api.model('Scroll', {
    'scroll': fields.Boolean(default=False, description='Flag if scrolling is enables.'),
    'total': fields.Integer(default=0, description='The total amount of hits for the search.'),
    'scroll_id': fields.String(default=None, allow_null=True, description='The scroll_id that can be used to retrieve the next page.'),
    'size': fields.Integer(default=0, help='The size of the returned scroll page.')})
''' Model used in responses with scroll. '''

aggregation_model = api.model('Aggregation', {
    'after': fields.String(description='The after key for the current request.', allow_null=True),
    'total': fields.Integer(default=0, description='The total amount of hits for the search.'),
    'per_page': fields.Integer(default=0, help='The size of the requested page.', allow_null=True)})
''' Model used in responses with id aggregation. '''

query_model_fields = {
    qualified_name: quantity.flask_field
    for qualified_name, quantity in search.search_quantities.items()}

query_model_fields.update(**{
    'owner': fields.String(description='The group the calculations belong to.', allow_null=True, skip_none=True),
    'domain': fields.String(description='Specify the domain to search in: %s, default is ``%s``' % (
        ', '.join(['``%s``' % domain for domain in datamodel.domains]), config.meta.default_domain)),
    'from_time': fields.Raw(description='The minimum entry time.', allow_null=True, skip_none=True),
    'until_time': fields.Raw(description='The maximum entry time.', allow_null=True, skip_none=True)
})

query_model_fields.update(**{
    '$and': fields.List(fields.Raw, description=(
        'List of queries which must be present in search results.')),
    '$or': fields.List(fields.Raw, description=(
        'List of queries which should be present in search results.')),
    '$not': fields.List(fields.Raw, description=(
        'List of queries which must not be present in search results.')),
    '$lt': fields.Raw(description=(
        'Dict of quantiy name: value such that search results should have values '
        'less than value.')),
    '$lte': fields.Raw(description=(
        'Dict of quantiy name: value such that search results should have values '
        'less than or equal to value')),
    '$gt': fields.Raw(description=(
        'Dict of quantiy name: value such that search results should have values '
        'greater than value')),
    '$gte': fields.Raw(description=(
        'Dict of quantiy name: value such that search results should have values '
        'greater than or equal to value')),
})

query_model = api.model('Query', query_model_fields)

search_model_fields = {
    'query': fields.Nested(query_model, allow_null=True, skip_none=True),
    'pagination': fields.Nested(pagination_model, allow_null=True, skip_none=True),
    'scroll': fields.Nested(scroll_model, allow_null=True, skip_none=True),
    'aggregation': fields.Nested(aggregation_model, allow_null=True),
    'results': fields.List(fields.Raw(allow_null=True, skip_none=True), description=(
        'A list of search results. Each result is a dict with quantitie names as key and '
        'values as values'), allow_null=True, skip_none=True),
    'code': fields.Nested(api.model('Code', {
        'repo_url': fields.String(description=(
            'An encoded URL for the search query on the repo api.')),
        'python': fields.String(description=(
            'A piece of python code snippet which can be executed to reproduce the api result.')),
        'curl': fields.String(description=(
            'A curl command which can be executed to reproduce the api result.')),
        'clientlib': fields.String(description=(
            'A piece of python code which uses NOMAD\'s client library to access the archive.'))
    }), allow_null=True, skip_none=True)}

search_model = api.model('Search', search_model_fields)


def add_pagination_parameters(request_parser):
    ''' Add pagination parameters to Flask querystring parser. '''
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
    ''' Add scroll parameters to Flask querystring parser. '''
    request_parser.add_argument(
        'scroll', type=bool, help='Enable scrolling')
    request_parser.add_argument(
        'scroll_id', type=str, help='The id of the current scrolling window to use.')


def add_search_parameters(request_parser):
    ''' Add search parameters to Flask querystring parser. '''
    # more search parameters
    request_parser.add_argument(
        'domain', type=str,
        help='Specify the domain to search in: %s, default is ``%s``' % (
            ', '.join(['``%s``' % domain for domain in datamodel.domains]),
            config.meta.default_domain))
    request_parser.add_argument(
        'owner', type=str,
        help='Specify which calcs to return: ``visible``, ``public``, ``all``, ``user``, ``staging``, default is ``visible``')
    request_parser.add_argument(
        'from_time', type=lambda x: rfc3339DateTime.parse(x),
        help='A yyyy-MM-ddTHH:mm:ss (RFC3339) minimum entry time (e.g. upload time)')
    request_parser.add_argument(
        'until_time', type=lambda x: rfc3339DateTime.parse(x),
        help='A yyyy-MM-ddTHH:mm:ss (RFC3339) maximum entry time (e.g. upload time)')
    request_parser.add_argument(
        'dft.optimade', type=str,
        help='A search query in the optimade filter language.')
    request_parser.add_argument(
        'query', type=str,
        help='A json serialized structured search query (as used in POST reuquests).')

    # main search parameters
    for qualified_name, quantity in search.search_quantities.items():
        request_parser.add_argument(
            qualified_name, help=quantity.description, action=quantity.argparse_action)


_search_quantities = set(search.search_quantities.keys())


def apply_search_parameters(search_request: search.SearchRequest, args: Dict[str, Any]):
    '''
    Help that adds query relevant request args to the given SearchRequest.
    '''
    args = {key: value for key, value in args.items() if value is not None}

    # domain
    domain = args.get('domain')
    if domain is not None:
        search_request.domain(domain=domain)

    # owner
    owner = args.get('owner', 'visible')
    try:
        search_request.owner(
            owner,
            g.user.user_id if g.user is not None else None)
    except search.AuthenticationRequiredError as e:
        abort(401, str(e))
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
        abort(400, 'bad datetime format')

    # optimade
    try:
        optimade = args.get('dft.optimade', None)
        if optimade is not None:
            q = filterparser.parse_filter(
                optimade, nomad_properties=domain, without_prefix=True)
            search_request.query(q)
    except filterparser.FilterException as e:
        abort(400, 'Could not parse optimade query: %s' % str(e))

    # search expression
    query_str = args.get('query', None)
    if query_str is not None:
        try:
            query = json.loads(query_str)
        except Exception as e:
            abort(400, 'Could not JSON parse query expression: %s' % str(e))

        try:
            search_request.query_expression(query)
        except Exception as e:
            abort(400, 'Invalid query expression: %s' % str(e))

    # search parameter
    search_request.search_parameters(**{
        key: value for key, value in args.items()
        if key in _search_quantities})


def calc_route(ns, prefix: str = ''):
    ''' A resource decorator for /<upload>/<calc> based routes. '''
    def decorator(func):
        ns.route('%s/<string:upload_id>/<string:calc_id>' % prefix)(
            api.doc(params={
                'upload_id': 'The unique id for the requested upload.',
                'calc_id': 'The unique id for the requested calculation.'
            })(func)
        )
    return decorator


def upload_route(ns, prefix: str = ''):
    ''' A resource decorator for /<upload> based routes. '''
    def decorator(func):
        ns.route('%s/<string:upload_id>' % prefix)(
            api.doc(params={
                'upload_id': 'The unique id for the requested upload.'
            })(func)
        )
    return decorator


def streamed_zipfile(
        files: Iterable[Tuple[str, str, dict, Callable[[str], IO], Callable[[str], int]]],
        zipfile_name: str, compress: bool = False, manifest: dict = None):
    '''
    Creates a response that streams the given files as a streamed zip file. Ensures that
    each given file is only streamed once, based on its filename in the resulting zipfile.

    Arguments:
        files: An iterable of tuples with the filename to be used in the resulting zipfile,
            an file id within the upload, an optional manifest, a callable that gives an
            binary IO object for the file id, and a callable that gives the file size for
            the file id.
        zipfile_name: A name that will be used in the content disposition attachment
            used as an HTTP respone.
        compress: Uses compression. Default is stored only.
        manifest: The dict contents of the manifest file. Will be extended if the files
            provide manifest information.
    '''

    streamed_files: Set[str] = set()

    def generator():
        ''' Stream a zip file with all files using zipstream. '''
        def iterator():
            '''
            Replace the directory based iter of zipstream with an iter over all given
            files.
            '''
            collected_manifest = manifest
            # the actual contents
            for zipped_filename, file_id, manifest_part, open_io, file_size in files:
                if manifest_part is not None:
                    if collected_manifest is None:
                        collected_manifest = {}
                    collected_manifest.update(manifest_part)

                if zipped_filename in streamed_files:
                    continue

                # Write a file to the zipstream.
                try:
                    f = open_io(file_id)
                    try:
                        def iter_content():
                            while True:
                                data = f.read(1024 * 64)
                                if len(data) == 0:
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
                except Exception as e:
                    utils.get_logger('__name__').error(
                        'unexpected exception while streaming zipfile', exc_info=e)
                    manifest.setdefault('errors', []).append(
                        'unexpected exception while streaming zipfile: %s' % str(e))

            try:
                if collected_manifest is not None:
                    manifest_content = json.dumps(collected_manifest).encode('utf-8')
                    yield dict(
                        arcname='manifest.json', iterable=[manifest_content],
                        buffer_size=len(manifest_content))
            except Exception as e:
                utils.get_logger('__name__').error(
                    'could not stream zipfile manifest', exc_info=e)

        compression = zipfile.ZIP_DEFLATED if compress else zipfile.ZIP_STORED
        zip_stream = zipstream.ZipFile(mode='w', compression=compression, allowZip64=True)
        zip_stream.paths_to_write = iterator()

        for chunk in zip_stream:
            yield chunk

    response = Response(stream_with_context(generator()), mimetype='application/zip')
    response.headers['Content-Disposition'] = 'attachment; filename={}'.format(zipfile_name)
    return response


def _query_api_url(*args, query: Dict[str, Any] = None):
    '''
    Creates a API URL.
    Arguments:
        *args: URL path segments after the API base URL
        query_string: A dict with query string parameters
    '''
    url = os.path.join(config.api_url(False), *args)
    if query is not None and len(query) > 0:
        url = '%s?%s' % (url, urlencode(query, doseq=True))

    return url


def query_api_python(query):
    '''
    Creates a string of python code to execute a search query to the repository using
    the requests library.
    '''
    query = _filter_api_query(query)
    url = _query_api_url('archive', 'query')
    return '''import requests
response = requests.post('{}', json={{
    'query': {{
{}
    }}
}})
data = response.json()'''.format(
        url,
        ',\n'.join([
            '        \'%s\': %s' % (key, pprint.pformat(value, compact=True))
            for key, value in query.items()])
    )


def _filter_api_query(query):
    def normalize_value(key, value):
        quantity = search.search_quantities.get(key)
        if quantity.many:
            if not isinstance(value, list):
                return [value]
            else:
                return value
        elif isinstance(value, list) and len(value) == 1:
            return value[0]

        return value

    result = {
        key: normalize_value(key, value) for key, value in query.items()
        if key in search.search_quantities and (key != 'domain' or value != config.meta.default_domain)
    }

    for key in ['dft.optimade']:
        if key in query:
            result[key] = query[key]

    return result


def query_api_repo_url(query):
    '''
    Creates an encoded URL string access a search query on the repo api.
    '''
    query = dict(query)
    for to_delete in ['per_page', 'page', 'exclude']:
        if to_delete in query:
            del(query[to_delete])
    for key, value in dict(order_by=['upload_time'], order=['-1'], domain=['dft'], owner=['public']).items():
        if key in query and query[key] == value:
            del(query[key])
    return _query_api_url('repo', query=query)


def query_api_clientlib(**kwargs):
    '''
    Creates a string of python code to execute a search query on the archive using
    the client library.
    '''
    query = _filter_api_query(kwargs)

    out = io.StringIO()
    out.write('from nomad import client, config\n')
    out.write('config.client.url = \'%s\'\n' % config.api_url(ssl=False))
    out.write('results = client.query_archive(query={%s' % ('' if len(kwargs) == 0 else '\n'))
    out.write(',\n'.join([
        '    \'%s\': %s' % (key, pprint.pformat(value, compact=True))
        for key, value in query.items()]))
    out.write('})\n')
    out.write('print(results)\n')

    return out.getvalue()


def query_api_curl(query):
    '''
    Creates a string of curl command to execute a search query and download the respective
    archives in a .zip file.
    '''
    url = _query_api_url('archive', 'download', query=query)
    return 'curl "%s" --output nomad.zip' % url


def enable_gzip(level: int = 1, min_size: int = 1024):
    """
    Args:
        level: The gzip compression level from 1-9
        min_size: The minimum response size in bytes for which the compression
            will be enabled.
    """
    def inner(function):
        @wraps(function)
        def wrapper(*args, **kwargs):
            response = make_response(function(*args, **kwargs))
            if response.status_code == 200:
                accept_encoding = request.headers["Accept-Encoding"]
                content_length = int(response.headers["Content-Length"])
                if "gzip" in accept_encoding and content_length >= min_size:
                    data = response.data
                    data = gzip.compress(data, level)
                    response.data = data
                    response.headers['Content-Length'] = len(data)
                    response.headers["Content-Encoding"] = "gzip"
                    return response
            return response
        return wrapper
    return inner
