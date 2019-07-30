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
The repository API of the nomad@FAIRDI APIs. Currently allows to resolve repository
meta-data.
"""

from flask_restplus import Resource, abort, fields
from flask import request, g
from elasticsearch_dsl import Q
from elasticsearch.exceptions import NotFoundError
import datetime

from nomad import search

from .app import api, rfc3339DateTime
from .auth import login_if_available
from .common import pagination_model, pagination_request_parser, calc_route

ns = api.namespace('repo', description='Access repository metadata.')


@calc_route(ns)
class RepoCalcResource(Resource):
    @api.response(404, 'The upload or calculation does not exist')
    @api.response(401, 'Not authorized to access the calculation')
    @api.response(200, 'Metadata send', fields.Raw)
    @api.doc('get_repo_calc')
    @login_if_available
    def get(self, upload_id, calc_id):
        """
        Get calculation metadata in repository form.

        Repository metadata only entails the quantities shown in the repository.
        Calcs are references via *upload_id*, *calc_id* pairs.
        """
        try:
            calc = search.Entry.get(calc_id)
        except NotFoundError:
            abort(404, message='There is no calculation %s/%s' % (upload_id, calc_id))

        if calc.with_embargo or not calc.published:
            if g.user is None:
                abort(401, message='Not logged in to access %s/%s.' % (upload_id, calc_id))

            is_owner = g.user.user_id == 0
            if not is_owner:
                for owner in calc.owners:
                    # At somepoint ids will be emails (strings) anyways.
                    # Right now it is hard to make sure that both are either str or int.
                    if str(owner.user_id) == str(g.user.user_id):
                        is_owner = True
                        break
            if not is_owner:
                abort(401, message='Not authorized to access %s/%s.' % (upload_id, calc_id))

        return calc.to_dict(), 200


repo_calcs_model = api.model('RepoCalculations', {
    'pagination': fields.Nested(pagination_model, allow_null=True),
    'scroll': fields.Nested(allow_null=True, skip_none=True, model=api.model('Scroll', {
        'total': fields.Integer(description='The total amount of hits for the search.'),
        'scroll_id': fields.String(allow_null=True, description='The scroll_id that can be used to retrieve the next page.'),
        'size': fields.Integer(help='The size of the returned scroll page.')})),
    'results': fields.List(fields.Raw, description=(
        'A list of search results. Each result is a dict with quantitie names as key and '
        'values as values')),
    'quantities': fields.Raw(description=(
        'A dict with all aggregations. Each aggregation is dictionary with a metrics dict as '
        'value and quantity value as key. The metrics are code runs(calcs), %s. '
        'There is a pseudo quantity "total" with a single value "all" that contains the metrics over all results. ' %
        ', '.join(search.metrics_names)))
})


def add_common_parameters(request_parser):
    request_parser.add_argument(
        'owner', type=str,
        help='Specify which calcs to return: ``all``, ``public``, ``user``, ``staging``, default is ``all``')
    request_parser.add_argument(
        'from_time', type=lambda x: rfc3339DateTime.parse(x),
        help='A yyyy-MM-ddTHH:mm:ss (RFC3339) minimum entry time (e.g. upload time)')
    request_parser.add_argument(
        'until_time', type=lambda x: rfc3339DateTime.parse(x),
        help='A yyyy-MM-ddTHH:mm:ss (RFC3339) maximum entry time (e.g. upload time)')

    for quantity in search.search_quantities.values():
        request_parser.add_argument(
            quantity.name, type=str, help=quantity.description,
            action='append' if quantity.multi else None)


repo_request_parser = pagination_request_parser.copy()
add_common_parameters(repo_request_parser)
repo_request_parser.add_argument(
    'scroll', type=bool, help='Enable scrolling')
repo_request_parser.add_argument(
    'scroll_id', type=str, help='The id of the current scrolling window to use.')
repo_request_parser.add_argument(
    'metrics', type=str, help=(
        'Metrics to aggregate over all quantities and their values as comma separated list. '
        'Possible values are %s.' % ', '.join(search.metrics_names)))


def create_owner_query():
    owner = request.args.get('owner', 'all')

    # TODO this should be removed after migration
    # if owner == 'migrated':
    #     q = Q('term', published=True) & Q('term', with_embargo=False)
    #     if g.user is not None:
    #         q = q | Q('term', owners__user_id=g.user.user_id)
    #     q = q & ~Q('term', **{'uploader.user_id': 1})  # pylint: disable=invalid-unary-operand-type
    if owner == 'all':
        q = Q('term', published=True) & Q('term', with_embargo=False)
        if g.user is not None:
            q = q | Q('term', owners__user_id=g.user.user_id)
    elif owner == 'public':
        q = Q('term', published=True) & Q('term', with_embargo=False)
    elif owner == 'user':
        if g.user is None:
            abort(401, message='Authentication required for owner value user.')

        q = Q('term', owners__user_id=g.user.user_id)
    elif owner == 'staging':
        if g.user is None:
            abort(401, message='Authentication required for owner value user.')
        q = Q('term', published=False) & Q('term', owners__user_id=g.user.user_id)
    elif owner == 'admin':
        if g.user is None or not g.user.is_admin:
            abort(401, message='This can only be used by the admin user.')
        q = None
    else:
        abort(400, message='Invalid owner value. Valid values are all|user|staging, default is all')

    return q


def create_search_parameters():
    """ Helper that creates a request.args dict with isolated search parameters """
    search_parameters = dict(**request.args)
    search_parameters.pop('owner', None)
    search_parameters.pop('scroll', None)
    search_parameters.pop('scroll_id', None)
    search_parameters.pop('per_page', None)
    search_parameters.pop('page', None)
    search_parameters.pop('order', None)
    search_parameters.pop('order_by', None)
    search_parameters.pop('metrics', None)
    search_parameters.pop('from_time', None)
    search_parameters.pop('until_time', None)
    search_parameters.pop('size', None)
    search_parameters.pop('after', None)

    return search_parameters


@ns.route('/')
class RepoCalcsResource(Resource):
    @api.doc('search')
    @api.response(400, 'Invalid requests, e.g. wrong owner type or bad search parameters')
    @api.expect(repo_request_parser, validate=True)
    @api.marshal_with(repo_calcs_model, skip_none=True, code=200, description='Search results send')
    @login_if_available
    def get(self):
        """
        Search for calculations in the repository form, paginated.

        The ``owner`` parameter determines the overall entries to search through.
        Possible values are: ``all`` (show all entries visible to the current user), ``public``
        (show all publically visible entries), ``user`` (show all user entries, requires login),
        ``staging`` (show all user entries in staging area, requires login).

        You can use the various quantities to search/filter for. For some of the
        indexed quantities this endpoint returns aggregation information. This means
        you will be given a list of all possible values and the number of entries
        that have the certain value. You can also use these aggregations on an empty
        search to determine the possible values.

        The pagination parameters allows determine which page to return via the
        ``page`` and ``per_page`` parameters. Pagination however, is limited to the first
        100k (depending on ES configuration) hits.

        An alternative to pagination is to use ``scroll`` and ``scroll_id``. With ``scroll``
        you will get a ``scroll_id`` on the first request. Each call with ``scroll`` and
        the respective ``scroll_id`` will return the next ``per_page`` (here the default is 1000)
        results. Scroll however, ignores ordering and does not return aggregations.
        The scroll view used in the background will stay alive for 1 minute between requests.
        If the given ``scroll_id`` is not available anymore, a HTTP 400 is raised.

        The search will return aggregations on a predefined set of quantities. Aggregations
        will tell you what quantity values exist and how many entries match those values.

        Ordering is determined by ``order_by`` and ``order`` parameters.
        """

        try:
            scroll = bool(request.args.get('scroll', False))
            scroll_id = request.args.get('scroll_id', None)
            page = int(request.args.get('page', 1))
            per_page = int(request.args.get('per_page', 10 if not scroll else 1000))
            order = int(request.args.get('order', -1))
            metrics_str = request.args.get('metrics', '')

            from_time = rfc3339DateTime.parse(request.args.get('from_time', '2000-01-01'))
            until_time_str = request.args.get('until_time', None)
            until_time = rfc3339DateTime.parse(until_time_str) if until_time_str is not None else datetime.datetime.utcnow()
            time_range = (from_time, until_time)

            metrics = [
                metric for metric in metrics_str.split(',')
                if metric in search.metrics_names]
        except Exception:
            abort(400, message='bad parameter types')

        order_by = request.args.get('order_by', 'formula')

        try:
            assert page >= 1
            assert per_page > 0
        except AssertionError:
            abort(400, message='invalid pagination')

        if order not in [-1, 1]:
            abort(400, message='invalid pagination')

        q = create_owner_query()

        # TODO this should be removed after migration
        without_currupted_mainfile = ~Q('term', code_name='currupted mainfile')  # pylint: disable=invalid-unary-operand-type
        q = q & without_currupted_mainfile if q is not None else without_currupted_mainfile

        search_parameters = create_search_parameters()

        try:
            if scroll:
                results = search.scroll_search(
                    q=q, scroll_id=scroll_id, size=per_page, search_parameters=search_parameters)

            else:
                results = search.metrics_search(
                    q=q, per_page=per_page, page=page, order=order, order_by=order_by,
                    time_range=time_range, metrics_to_use=metrics, search_parameters=search_parameters)

                # TODO just a work around to make things prettier
                quantities = results['quantities']
                if 'code_name' in quantities and 'currupted mainfile' in quantities['code_name']:
                    del(quantities['code_name']['currupted mainfile'])

            return results, 200
        except search.ScrollIdNotFound:
            abort(400, 'The given scroll_id does not exist.')
        except KeyError as e:
            abort(400, str(e))


repo_quantity_values_model = api.model('RepoQuantityValues', {
    'quantities': fields.Raw(description='''
        A dict with the requested quantity as single key.
        The value is a dictionary with 'after' and 'values' keys.
        The 'values' key holds a dict with actual values as keys and their entry count
        as values (i.e. number of entries with that value). ''')
})

repo_quantity_search_request_parser = api.parser()
add_common_parameters(repo_quantity_search_request_parser)
repo_quantity_search_request_parser.add_argument(
    'after', type=str, help='The after value to use for "scrolling".')
repo_request_parser.add_argument(
    'size', type=int, help='The max size of the returned values.')


@ns.route('/<string:quantity>')
class RepoQuantityResource(Resource):
    @api.doc('quantity_search')
    @api.response(400, 'Invalid requests, e.g. wrong owner type, bad quantity, bad search parameters')
    @api.expect(repo_quantity_search_request_parser, validate=True)
    @api.marshal_with(repo_quantity_values_model, skip_none=True, code=200, description='Search results send')
    @login_if_available
    def get(self, quantity: str):
        """
        Retrieve quantity values from entries matching the search.

        You can use the various quantities to search/filter for. For some of the
        indexed quantities this endpoint returns aggregation information. This means
        you will be given a list of all possible values and the number of entries
        that have the certain value. You can also use these aggregations on an empty
        search to determine the possible values.

        There is no ordering and no pagination. Instead there is an 'after' key based
        scrolling. The result will contain an 'after' value, that can be specified
        for the next request. You can use the 'size' and 'after' parameters accordingly.

        The result will contain a 'quantities' key with the given quantity and the
        respective values (upto 'size' many). For the rest of the values use the
        'after' parameter accordingly.
        """

        try:
            after = request.args.get('after', None)
            size = int(request.args.get('size', 100))

            from_time = rfc3339DateTime.parse(request.args.get('from_time', '2000-01-01'))
            until_time_str = request.args.get('until_time', None)
            until_time = rfc3339DateTime.parse(until_time_str) if until_time_str is not None else datetime.datetime.utcnow()
            time_range = (from_time, until_time)
        except Exception:
            abort(400, message='bad parameter types')

        try:
            assert size >= 0
        except AssertionError:
            abort(400, message='invalid size')

        q = create_owner_query()
        search_parameters = create_search_parameters()

        try:
            results = search.quantity_search(
                q=q, time_range=time_range, search_parameters=search_parameters,
                quantities={quantity: after}, size=size, with_entries=False)

            return results, 200
        except KeyError as e:
            import traceback
            traceback.print_exc()
            abort(400, 'Given quantity does not exist: %s' % str(e))
