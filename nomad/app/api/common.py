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

from flask_restplus import fields

from nomad.app.utils import RFC3339DateTime

from .api import api


dataset_model = api.model('DataSet', {
    'id': fields.Integer(required=True, description='The repository db dataset id'),
    '_doi': fields.String(description='The DOI of the dataset'),
    '_name': fields.String(description='The unique dataset name')
})

metadata_model = api.model('MetaData', {
    'with_embargo': fields.Boolean(default=False, description='Data with embargo is only visible to the upload until the embargo period ended.'),
    'comment': fields.String(description='The comment are shown in the repository for each calculation.'),
    'references': fields.List(fields.String, descriptions='References allow to link calculations to external source, e.g. URLs.'),
    'coauthors': fields.List(fields.String, description='A list of co-authors given by user_id.'),
    'shared_with': fields.List(fields.String, description='A list of users to share calculations with given by user_id.'),
    '_upload_time': RFC3339DateTime(description='Overrride the upload time.'),
    '_uploader': fields.String(description='Override the uploader with the given user id.'),
    'datasets': fields.List(fields.Nested(model=dataset_model, skip_none=True), description='A list of datasets.')
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
