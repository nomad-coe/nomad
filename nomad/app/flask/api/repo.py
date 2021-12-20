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
The repository API of the nomad@FAIRDI APIs. Currently allows to resolve repository
meta-data.
'''

from typing import List, Dict, Any
from flask_restplus import Resource, abort, fields
from flask import request, g
from elasticsearch_dsl import Q
from elasticsearch.exceptions import NotFoundError
import elasticsearch.helpers
from datetime import datetime

from nomad import search, utils, datamodel, processing as proc, infrastructure, files, metainfo
from nomad.datamodel import Dataset, User, EditableUserMetadata

from .. import common
from ..common import RFC3339DateTime, DotKeyNested
from .api import api
from .auth import authenticate
from .common import search_model, calc_route, add_pagination_parameters,\
    add_scroll_parameters, add_search_parameters, apply_search_parameters,\
    query_api_python, query_api_curl, query_api_clientlib, query_api_repo_url, \
    _search_quantities

ns = api.namespace('repo', description='Access repository metadata.')


@calc_route(ns)
class RepoCalcResource(Resource):
    @api.response(404, 'The upload or calculation does not exist')
    @api.response(401, 'Not authorized to access the calculation')
    @api.response(200, 'Metadata send', fields.Raw)
    @api.doc('get_repo_calc')
    @authenticate()
    def get(self, upload_id, calc_id):
        '''
        Get calculation metadata in repository form.

        Repository metadata only entails the quantities shown in the repository.
        Calcs are references via *upload_id*, *calc_id* pairs.
        '''
        try:
            calc = search.entry_document.get(calc_id)
        except NotFoundError:
            abort(404, message='There is no calculation %s/%s' % (upload_id, calc_id))

        if calc.with_embargo or not calc.published:
            if g.user is None:
                abort(401, message='Not logged in to access %s/%s.' % (upload_id, calc_id))

            if not (any(g.user.user_id == user.user_id for user in calc.owners) or g.user.is_admin):
                abort(401, message='Not authorized to access %s/%s.' % (upload_id, calc_id))

        result = search._es_to_entry_dict(calc, required=None)
        result['code'] = {
            'python': query_api_python(dict(upload_id=upload_id, calc_id=calc_id)),
            'curl': query_api_curl(dict(upload_id=upload_id, calc_id=calc_id)),
            'clientlib': query_api_clientlib(upload_id=[upload_id], calc_id=[calc_id])
        }

        return result, 200


_search_request_parser = api.parser()
add_pagination_parameters(_search_request_parser)
add_scroll_parameters(_search_request_parser)
add_search_parameters(_search_request_parser)
_search_request_parser.add_argument(
    'date_histogram', type=bool, help='Add an additional aggregation over the upload time')
_search_request_parser.add_argument(
    'interval', type=str, help='Interval to use for upload time aggregation.')
_search_request_parser.add_argument(
    'metrics', type=str, action='append', help=(
        'Metrics to aggregate over all quantities and their values as comma separated list. '
        'Possible values are %s.' % ', '.join(search.metrics.keys())))
_search_request_parser.add_argument(
    'statistics', type=str, action='append', help=(
        'Quantities for which to aggregate values and their metrics.'))
_search_request_parser.add_argument(
    'exclude', type=str, action='split', help='Excludes the given keys in the returned data.')
for group_name in search.groups:
    _search_request_parser.add_argument(
        group_name, type=bool, help=('Return %s group data.' % group_name))
    _search_request_parser.add_argument(
        '%s_after' % group_name, type=str,
        help='The last %s id of the last scroll window for the %s group' % (group_name, group_name))

_repo_calcs_model_fields = {
    'statistics': fields.Raw(description=(
        'A dict with all statistics. Each statistic is dictionary with a metrics dict as '
        'value and quantity value as key. The possible metrics are code runs(calcs), %s. '
        'There is a pseudo quantity "total" with a single value "all" that contains the '
        ' metrics over all results. ' % ', '.join(search.metrics.keys())))}

for group_name in search.groups:
    _repo_calcs_model_fields[group_name] = (DotKeyNested if '.' in group_name else fields.Nested)(api.model('RepoGroup', {
        'after': fields.String(description='The after value that can be used to retrieve the next %s.' % group_name),
        'values': fields.Raw(description='A dict with %s as key. The values are dicts with "total" and "examples" keys.' % group_name)
    }), skip_none=True)

for qualified_name, quantity in search.search_quantities.items():
    _repo_calcs_model_fields[qualified_name] = fields.Raw(
        description=quantity.description, allow_null=True, skip_none=True)

_repo_calcs_model_fields.update(**{
    'date_histogram': fields.Boolean(default=False, description='Add an additional aggregation over the upload time', allow_null=True, skip_none=True),
    'interval': fields.String(description='Interval to use for upload time aggregation.', allow_null=True, skip_none=True),
    'metrics': fields.List(fields.String, description=(
        'Metrics to aggregate over all quantities and their values as comma separated list. '
        'Possible values are %s.' % ', '.join(search.metrics.keys())), allow_null=True, skip_none=True),
    'statistics_required': fields.List(fields.String, description='Quantities for which to aggregate values and their metrics.', allow_null=True, skip_none=True),
    'exclude': fields.List(fields.String, description='Excludes the given keys in the returned data.', allow_null=True, skip_none=True)
})

_repo_calcs_model = api.inherit('RepoCalculations', search_model, _repo_calcs_model_fields)


@ns.route('/')
class RepoCalcsResource(Resource):
    @api.doc('search')
    @api.response(400, 'Invalid requests, e.g. wrong owner type or bad search parameters')
    @api.expect(_search_request_parser, validate=True)
    @api.marshal_with(_repo_calcs_model, skip_none=True, code=200, description='Search results send')
    @authenticate()
    def get(self):
        '''
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

        Ordering is determined by ``order_by`` and ``order`` parameters. Default is
        ``upload_time`` in decending order.
        '''

        try:
            parsed_args = _search_request_parser.parse_args()
            args = {
                key: value for key, value in parsed_args.items()
                if value is not None}

            scroll = args.get('scroll', False)
            scroll_id = args.get('scroll_id', None)
            page = args.get('page', 1)
            per_page = args.get('per_page', 10 if not scroll else 1000)
            order = args.get('order', -1)
            order_by = args.get('order_by', 'upload_time')

            date_histogram = args.get('date_histogram', False)
            interval = args.get('interval', '1M')
            metrics: List[str] = request.args.getlist('metrics')
            statistics = args.get('statistics', [])
        except Exception as e:
            abort(400, message='bad parameters: %s' % str(e))

        for metric in metrics:
            if metric not in search.metrics:
                abort(400, message='there is no metric %s' % metric)

        search_request = search.SearchRequest()
        apply_search_parameters(search_request, args)
        if date_histogram:
            search_request.date_histogram(interval=interval, metrics_to_use=metrics)

        try:
            assert page >= 1
            assert per_page >= 0
        except AssertionError:
            abort(400, message='invalid pagination')

        if order not in [-1, 1]:
            abort(400, message='invalid pagination')

        if len(statistics) > 0:
            search_request.statistics(statistics, metrics_to_use=metrics)

        group_metrics = [
            group_quantity.metric_name
            for group_name, group_quantity in search.groups.items()
            if args.get(group_name, False)]
        total_metrics = metrics + group_metrics
        if len(total_metrics) > 0:
            search_request.totals(metrics_to_use=total_metrics)

        if 'exclude' in parsed_args:
            excludes = parsed_args['exclude']
            if excludes is not None:
                search_request.exclude(*excludes)

        try:
            if scroll:
                results = search_request.execute_scrolled(scroll_id=scroll_id, size=per_page)

            else:
                for group_name, group_quantity in search.groups.items():
                    if args.get(group_name, False):
                        kwargs: Dict[str, Any] = {}
                        if group_name == 'uploads_grouped':
                            kwargs.update(order_by='upload_time', order='desc')
                        search_request.quantity(
                            group_quantity.qualified_name, size=per_page, examples=1,
                            after=request.args.get('%s_after' % group_name, None),
                            **kwargs)

                results = search_request.execute_paginated(
                    per_page=per_page, page=page, order=order, order_by=order_by)

                # TODO just a work around to make things prettier
                if 'statistics' in results:
                    statistics = results['statistics']
                    if 'code_name' in statistics and 'currupted mainfile' in statistics['code_name']:
                        del(statistics['code_name']['currupted mainfile'])

                if 'quantities' in results:
                    quantities = results.pop('quantities')

                for group_name, group_quantity in search.groups.items():
                    if args.get(group_name, False):
                        results[group_name] = quantities[group_quantity.qualified_name]

            # build python code/curl snippet
            code_args = request.args.to_dict(flat=False)
            if 'statistics' in code_args:
                del(code_args['statistics'])
            results['code'] = {
                'repo_url': query_api_repo_url(code_args),
                'curl': query_api_curl(code_args),
                'python': query_api_python(code_args),
                'clientlib': query_api_clientlib(**code_args)
            }

            return results, 200
        except search.ScrollIdNotFound:
            abort(400, 'The given scroll_id does not exist.')
        except KeyError as e:
            import traceback
            traceback.print_exc()
            abort(400, str(e))

    @api.doc('post_search')
    @api.response(400, 'Invalid requests, e.g. wrong owner type or bad search parameters')
    @api.expect(_repo_calcs_model)
    @api.marshal_with(_repo_calcs_model, skip_none=True, code=200, description='Search results send')
    @authenticate()
    def post(self):
        '''
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

        Ordering is determined by ``order_by`` and ``order`` parameters. Default is
        ``upload_time`` in decending order.
        '''
        try:
            data_in = request.get_json()

            Scroll = data_in.get('scroll', {})
            scroll = Scroll.get('scroll', False)
            scroll_id = Scroll.get('scroll_id', None)
            pagination = data_in.get('pagination', {})
            page = pagination.get('page', 1)
            per_page = pagination.get('per_page', 10 if not scroll else 1000)
            order = pagination.get('order', -1)
            order_by = pagination.get('order_by', 'upload_time')

            date_histogram = data_in.get('date_histogram', False)
            interval = data_in.get('interval', '1M')
            metrics: List[str] = data_in.get('metrics', [])
            statistics = data_in.get('statistics_required', [])

            query = data_in.get('query', {})
            query_expression = {key: val for key, val in query.items() if '$' in key}
        except Exception as e:
            abort(400, message='bad parameters: %s' % str(e))

        for metric in metrics:
            if metric not in search.metrics:
                abort(400, message='there is no metric %s' % metric)

        search_request = search.SearchRequest()
        apply_search_parameters(search_request, query)
        if date_histogram:
            search_request.date_histogram(interval=interval, metrics_to_use=metrics)

        if query_expression:
            try:
                search_request.query_expression(query_expression)
            except AssertionError as e:
                abort(400, str(e))

        try:
            assert page >= 1
            assert per_page >= 0
        except AssertionError:
            abort(400, message='invalid pagination')

        if order not in [-1, 1]:
            abort(400, message='invalid pagination')

        if len(statistics) > 0:
            search_request.statistics(statistics, metrics_to_use=metrics)

        group_metrics = [
            group_quantity.metric_name
            for group_name, group_quantity in search.groups.items()
            if group_name in data_in]
        total_metrics = metrics + group_metrics
        if len(total_metrics) > 0:
            search_request.totals(metrics_to_use=total_metrics)

        if 'exclude' in data_in:
            excludes = data_in['exclude']
            if excludes is not None:
                search_request.exclude(*excludes)

        try:
            if scroll:
                results = search_request.execute_scrolled(scroll_id=scroll_id, size=per_page)

            else:
                for group_name, group_quantity in search.groups.items():
                    if group_name in data_in:
                        kwargs: Dict[str, Any] = {}
                        if group_name == 'uploads_grouped':
                            kwargs.update(order_by='upload_time', order='desc')
                        search_request.quantity(
                            group_quantity.qualified_name, size=per_page, examples=1,
                            after=data_in[group_name].get('after', None),
                            **kwargs)

                results = search_request.execute_paginated(
                    per_page=per_page, page=page, order=order, order_by=order_by)

                # TODO just a work around to make things prettier
                if 'statistics' in results:
                    statistics = results['statistics']
                    if 'code_name' in statistics and 'currupted mainfile' in statistics['code_name']:
                        del(statistics['code_name']['currupted mainfile'])

                if 'quantities' in results:
                    quantities = results.pop('quantities')

                for group_name, group_quantity in search.groups.items():
                    if group_name in data_in:
                        results[group_name] = quantities[group_quantity.qualified_name]

            # build python code/curl snippet
            code_args = dict(data_in)
            if 'statistics' in code_args:
                del(code_args['statistics'])
            results['code'] = {
                'curl': query_api_curl(code_args),
                'python': query_api_python(code_args),
                'clientlib': query_api_clientlib(**code_args)
            }

            return results, 200
        except search.ScrollIdNotFound:
            abort(400, 'The given scroll_id does not exist.')
        except KeyError as e:
            import traceback
            traceback.print_exc()
            abort(400, str(e))


_query_model_parameters = {
    'owner': fields.String(description='Specify which calcs to return: ``all``, ``public``, ``user``, ``staging``, default is ``all``'),
    'from_time': RFC3339DateTime(description='A yyyy-MM-ddTHH:mm:ss (RFC3339) minimum entry time (e.g. upload time)'),
    'until_time': RFC3339DateTime(description='A yyyy-MM-ddTHH:mm:ss (RFC3339) maximum entry time (e.g. upload time)')
}

for qualified_name, quantity in search.search_quantities.items():
    if quantity.many_and == 'append' or quantity.many_or == 'append':
        def field(**kwargs):
            return fields.List(fields.String(**kwargs))
    else:
        field = fields.String
    _query_model_parameters[qualified_name] = field(description=quantity.description)

_repo_query_model = api.model('RepoQuery', _query_model_parameters, skip_none=True)


def repo_edit_action_field(quantity):
    if quantity.is_scalar:
        return fields.Nested(_repo_edit_action_model, description=quantity.description, skip_none=True)
    else:
        return fields.List(
            fields.Nested(_repo_edit_action_model, skip_none=True), description=quantity.description)


_repo_edit_action_model = api.model('RepoEditAction', {
    'value': fields.String(description='The value/values that is set as a string.'),
    'success': fields.Boolean(description='If this can/could be done. Only in API response.'),
    'message': fields.String(descriptin='A message that details the action result. Only in API response.')
})

_repo_edit_model = api.model('RepoEdit', {
    'verify': fields.Boolean(description='If true, no action is performed.'),
    'query': fields.Nested(_repo_query_model, skip_none=True, description='New metadata will be applied to query results.'),
    'actions': fields.Nested(
        api.model('RepoEditActions', {
            quantity.name: repo_edit_action_field(quantity)
            for quantity in EditableUserMetadata.m_def.definitions
            if isinstance(quantity, metainfo.Quantity)
        }), skip_none=True,
        description='Each action specifies a single value (even for multi valued quantities).'),
    'success': fields.Boolean(description='If the overall edit can/could be done. Only in API response.'),
    'message': fields.String(description='A message that details the overall edit result. Only in API response.')
})

_editable_quantities = {
    quantity.name: quantity for quantity in EditableUserMetadata.m_def.definitions}


def edit(parsed_query: Dict[str, Any], mongo_update: Dict[str, Any] = None, re_index=True) -> List[str]:
    # get all calculations that have to change
    with utils.timer(common.logger, 'edit query executed'):
        search_request = search.SearchRequest().include('calc_id', 'upload_id')
        apply_search_parameters(search_request, parsed_query)
        upload_ids = set()
        calc_ids = []

        for hit in search_request.execute_scan():
            calc_ids.append(hit['calc_id'])
            upload_ids.add(hit['upload_id'])

    # perform the update on the mongo db
    with utils.timer(common.logger, 'edit mongo update executed', size=len(calc_ids)):
        if mongo_update is not None:
            n_updated = proc.Calc.objects(calc_id__in=calc_ids).update(multi=True, **mongo_update)
            if n_updated != len(calc_ids):
                common.logger.error('edit repo did not update all entries', payload=mongo_update)

    # re-index the affected entries in elastic search
    with utils.timer(common.logger, 'edit elastic update executed', size=len(calc_ids)):
        if re_index:
            def elastic_updates():
                upload_files_cache: Dict[str, files.UploadFiles] = dict()

                for calc in proc.Calc.objects(calc_id__in=calc_ids):
                    upload_id = calc.upload_id
                    upload_files = upload_files_cache.get(upload_id)
                    if upload_files is None:
                        upload_files = files.UploadFiles.get(upload_id, is_authorized=lambda: True)
                        upload_files_cache[upload_id] = upload_files

                    try:
                        entry_metadata = calc.entry_metadata(upload_files)
                        entry = entry_metadata.a_elastic.create_index_entry().to_dict(include_meta=True)
                        entry['_op_type'] = 'index'

                        yield entry

                    except Exception as e:
                        common.logger.error('edit repo could not create index doc', exc_info=e)

                for upload_files in upload_files_cache.values():
                    upload_files.close()

            _, failed = elasticsearch.helpers.bulk(
                infrastructure.elastic_client, elastic_updates(), stats_only=True)
            search.refresh()
            if failed > 0:
                common.logger.error(
                    'edit repo with failed elastic updates',
                    payload=mongo_update, nfailed=len(failed))

    return list(upload_ids)


def get_uploader_ids(query):
    ''' Get all the uploader from the query, to check coauthers and shared_with for uploaders. '''
    search_request = search.SearchRequest()
    apply_search_parameters(search_request, query)
    search_request.quantity(name='uploader_id')
    return search_request.execute()['quantities']['uploader_id']['values']


@ns.route('/edit')
class EditRepoCalcsResource(Resource):
    @api.doc('edit_repo')
    @api.response(400, 'Invalid requests, e.g. wrong owner type or bad search parameters')
    @api.expect(_repo_edit_model)
    @api.marshal_with(_repo_edit_model, skip_none=True, code=200, description='Edit verified/performed')
    @authenticate()
    def post(self):
        ''' Edit repository metadata. '''

        # basic body parsing and some semantic checks
        json_data = request.get_json()
        if json_data is None:
            json_data = {}
        query = json_data.get('query', {})

        owner = query.get('owner', 'user')
        if owner not in ['user', 'staging']:
            abort(400, 'Not a valid owner for edit %s. Edit can only be performed in user or staging' % owner)
        query['owner'] = owner

        if 'actions' not in json_data:
            abort(400, 'Missing key actions in edit data')
        actions = json_data['actions']
        verify = json_data.get('verify', False)

        # preparing the query of entries that are edited
        parsed_query = {}
        for quantity_name, value in query.items():
            if quantity_name in _search_quantities:
                quantity = search.search_quantities[quantity_name]
                if quantity.many:
                    if not isinstance(value, list):
                        value = value.split(',')
                parsed_query[quantity_name] = value
        parsed_query['owner'] = owner
        parsed_query['domain'] = query.get('domain')

        # checking the edit actions and preparing a mongo update on the fly
        json_data['success'] = True
        mongo_update = {}
        uploader_ids = None
        lift_embargo = False
        has_error = False
        removed_datasets = None

        with utils.timer(common.logger, 'edit verified'):
            for action_quantity_name, quantity_actions in actions.items():
                quantity = _editable_quantities.get(action_quantity_name)
                if quantity is None:
                    abort(400, 'Unknown quantity %s' % action_quantity_name)

                quantity_flask = quantity.m_get_annotations('flask', {})
                if quantity_flask.get('admin_only', False):
                    if not g.user.is_admin():
                        abort(404, 'Only the admin user can set %s' % quantity.name)

                if isinstance(quantity_actions, list) == quantity.is_scalar:
                    abort(400, 'Wrong shape for quantity %s' % action_quantity_name)

                if not isinstance(quantity_actions, list):
                    quantity_actions = [quantity_actions]

                flask_verify = quantity_flask.get('verify', None)
                mongo_key = 'metadata__%s' % quantity.name
                has_error = False
                for action in quantity_actions:
                    action['success'] = True
                    action['message'] = None
                    action_value = action.get('value')
                    action_value = action_value if action_value is None else action_value.strip()

                    if action_quantity_name == 'with_embargo':
                        # ignore the actual value ... just lift the embargo
                        mongo_value = False
                        lift_embargo = True

                        # check if necessary
                        search_request = search.SearchRequest()
                        apply_search_parameters(search_request, parsed_query)
                        search_request.q = search_request.q & Q('term', with_embargo=True)
                        if search_request.execute()['total'] == 0:
                            action['success'] = False
                            has_error = True
                            action['message'] = 'There is no embargo to lift'
                            continue

                    elif action_value is None:
                        mongo_value = None

                    elif action_value == '':
                        mongo_value = None

                    elif flask_verify == datamodel.User:
                        try:
                            mongo_value = User.get(user_id=action_value).user_id
                        except KeyError:
                            action['success'] = False
                            has_error = True
                            action['message'] = 'User does not exist'
                            continue

                        if uploader_ids is None:
                            uploader_ids = get_uploader_ids(parsed_query)
                        if action_value in uploader_ids:
                            action['success'] = False
                            has_error = True
                            action['message'] = 'This user is already an uploader of one entry in the query'
                            continue

                    elif flask_verify == datamodel.Dataset:
                        try:
                            mongo_value = Dataset.m_def.a_mongo.get(
                                user_id=g.user.user_id, name=action_value).dataset_id
                        except KeyError:
                            action['message'] = 'Dataset does not exist and will be created'
                            mongo_value = None
                            if not verify:
                                dataset = Dataset(
                                    dataset_id=utils.create_uuid(), user_id=g.user.user_id,
                                    name=action_value, created=datetime.utcnow())
                                dataset.a_mongo.create()
                                mongo_value = dataset.dataset_id

                    else:
                        mongo_value = action_value

                    # verify if the edit action creates consistent shared_with and with_embargo for
                    # whole uploads
                    if action_quantity_name == 'shared_with' or action_quantity_name == 'with_embargo':
                        search_request = search.SearchRequest()
                        apply_search_parameters(search_request, parsed_query)
                        search_request.quantity('upload_id')
                        uploads = search_request.execute()['quantities']['upload_id']['values']
                        for upload_id, upload_data in uploads.items():
                            search_request = search.SearchRequest().search_parameters(upload_id=upload_id)
                            entries = search_request.execute()['total']
                            if entries != upload_data['total']:
                                action['success'] = False
                                action['message'] = json_data.get('message', '') + (
                                    'Edit would create an upload with inconsistent shared_with or with_embargo. '
                                    'You can only set those for all entries of an upload.')
                                has_error = True
                                continue

                    if len(quantity.shape) == 0:
                        mongo_update[mongo_key] = mongo_value
                    else:
                        mongo_values = mongo_update.setdefault(mongo_key, [])
                        if mongo_value is not None:
                            if mongo_value in mongo_values:
                                action['success'] = False
                                has_error = True
                                action['message'] = 'Duplicate values are not allowed'
                                continue
                            mongo_values.append(mongo_value)

                if len(quantity_actions) == 0 and len(quantity.shape) > 0:
                    mongo_update[mongo_key] = []

                if action_quantity_name == 'datasets':
                    # check if datasets edit is allowed and if datasets have to be removed
                    search_request = search.SearchRequest()
                    apply_search_parameters(search_request, parsed_query)
                    search_request.quantity(name='dataset_id')
                    old_datasets = list(
                        search_request.execute()['quantities']['dataset_id']['values'].keys())

                    removed_datasets = []
                    for dataset_id in old_datasets:
                        if dataset_id not in mongo_update.get(mongo_key, []):
                            removed_datasets.append(dataset_id)

                    doi_ds = Dataset.m_def.a_mongo.objects(
                        dataset_id__in=removed_datasets, doi__ne=None).first()
                    if doi_ds is not None:
                        json_data['success'] = False
                        json_data['message'] = json_data.get('message', '') + \
                            'Edit would remove entries from a dataset with DOI (%s) ' % doi_ds.name
                        has_error = True

        # stop here, if client just wants to verify its actions
        if verify:
            return json_data, 200

        # stop if the action were not ok
        if has_error:
            return json_data, 400

        # perform the change
        mongo_update['metadata__last_edit'] = datetime.utcnow()
        upload_ids = edit(parsed_query, mongo_update, True)

        # lift embargo
        if lift_embargo:
            for upload_id in upload_ids:
                upload = proc.Upload.get(upload_id)
                if upload.published:
                    upload.re_pack()

        # remove potentially empty old datasets
        if removed_datasets is not None:
            for dataset in removed_datasets:
                if proc.Calc.objects(metadata__datasets=dataset).first() is None:
                    Dataset.m_def.a_mongo.objects(dataset_id=dataset).delete()

        return json_data, 200


_repo_quantity_search_request_parser = api.parser()
add_search_parameters(_repo_quantity_search_request_parser)
_repo_quantity_search_request_parser.add_argument(
    'after', type=str, help='The after value to use for "scrolling".')
_repo_quantity_search_request_parser.add_argument(
    'size', type=int, help='The max size of the returned values.')
_repo_quantity_search_request_parser.add_argument(
    'value', type=str, help='A partial value. Only values that include this will be returned')

_repo_quantity_model = api.model('RepoQuantity', {
    'after': fields.String(description='The after value that can be used to retrieve the next set of values.'),
    'values': fields.Raw(description='A dict with values as key. Values are dicts with "total" and "examples" keys.')
})

_repo_quantity_values_model = api.model('RepoQuantityValues', {
    'quantity': fields.Nested(_repo_quantity_model, allow_null=True)
})


@ns.route('/quantity/<string:quantity>')
class RepoQuantityResource(Resource):
    @api.doc('quantity_search')
    @api.response(400, 'Invalid requests, e.g. wrong owner type, bad quantity, bad search parameters')
    @api.expect(_repo_quantity_search_request_parser, validate=True)
    @api.marshal_with(_repo_quantity_values_model, skip_none=True, code=200, description='Search results send')
    @authenticate()
    def get(self, quantity: str):
        '''
        Retrieve quantity values from entries matching the search.

        You can use the various quantities to search/filter for. For some of the
        indexed quantities this endpoint returns aggregation information. This means
        you will be given a list of all possible values and the number of entries
        that have the certain value. You can also use these aggregations on an empty
        search to determine the possible values.

        There is no ordering and no pagination. Instead there is an 'after' key based
        scrolling. The result will contain an 'after' value, that can be specified
        for the next request. You can use the 'size' and 'after' parameters accordingly.

        The result will contain a 'quantity' key with quantity values and the "after"
        value. There will be upto 'size' many values. For the rest of the values use the
        "after" parameter in another request.
        '''

        search_request = search.SearchRequest()
        args = {
            key: value
            for key, value in _repo_quantity_search_request_parser.parse_args().items()
            if value is not None}

        apply_search_parameters(search_request, args)
        after = args.get('after', None)
        size = args.get('size', 100)

        try:
            assert size >= 0
        except AssertionError:
            abort(400, message='invalid size')

        try:
            search_request.quantity(quantity, size=size, after=after)
            results = search_request.execute()
            quantities = results.pop('quantities')
            results['quantity'] = quantities[quantity]

            return results, 200
        except KeyError as e:
            import traceback
            traceback.print_exc()
            abort(400, 'Given quantity does not exist: %s' % str(e))


_repo_suggestions_search_request_parser = api.parser()
add_search_parameters(_repo_suggestions_search_request_parser)
_repo_suggestions_search_request_parser.add_argument(
    'size', type=int, help='The max size of the returned values.')
_repo_suggestions_search_request_parser.add_argument(
    'include', type=str, help='A substring that all values need to include.')

_repo_suggestions_model = api.model('RepoSuggestionsValues', {
    'suggestions': fields.List(fields.String, description='A list with the suggested values.')
})


@ns.route('/suggestions/<string:quantity>')
class RepoSuggestionsResource(Resource):
    @api.doc('suggestions_search')
    @api.response(400, 'Invalid requests, e.g. wrong owner type, bad quantity, bad search parameters')
    @api.expect(_repo_suggestions_search_request_parser, validate=True)
    @api.marshal_with(_repo_suggestions_model, skip_none=True, code=200, description='Suggestions send')
    @authenticate()
    def get(self, quantity: str):
        '''
        Retrieve the top values for the given quantity from entries matching the search.
        Values can be filtered by to include a given value.

        There is no ordering, no pagination, and no scroll interface.

        The result will contain a 'suggestions' key with values. There will be upto 'size' many values.
        '''

        search_request = search.SearchRequest()
        args = {
            key: value
            for key, value in _repo_suggestions_search_request_parser.parse_args().items()
            if value is not None}

        apply_search_parameters(search_request, args)
        size = args.get('size', 20)
        include = args.get('include', None)

        try:
            assert size >= 0
        except AssertionError:
            abort(400, message='invalid size')

        try:
            search_request.statistic(quantity, size=size, include=include, order=dict(_key='desc'))
            results = search_request.execute()
            values = {
                value: metric['code_runs']
                for value, metric in results['statistics'][quantity].items()
                if metric['code_runs'] > 0}
            results['suggestions'] = sorted(
                values.keys(), key=lambda value: values[value], reverse=True)

            return results, 200
        except KeyError as e:
            import traceback
            traceback.print_exc()
            abort(400, 'Given quantity does not exist: %s' % str(e))


_repo_quantities_search_request_parser = api.parser()
add_search_parameters(_repo_quantities_search_request_parser)
_repo_quantities_search_request_parser.add_argument(
    'quantities', type=str, action='append',
    help='The quantities to retrieve values from')
_repo_quantities_search_request_parser.add_argument(
    'size', type=int, help='The max size of the returned values.')

_repo_quantities_model = api.model('RepoQuantitiesResponse', {
    'quantities': fields.Nested(api.model('RepoQuantities', {
        quantity: fields.List(fields.Nested(_repo_quantity_model))
        for quantity in search.search_quantities
    }))
})


@ns.route('/quantities')
class RepoQuantitiesResource(Resource):
    @api.doc('quantities_search')
    @api.response(400, 'Invalid requests, e.g. wrong owner type, bad quantity, bad search parameters')
    @api.expect(_repo_quantities_search_request_parser, validate=True)
    @api.marshal_with(_repo_quantities_model, skip_none=True, code=200, description='Search results send')
    @authenticate()
    def get(self):
        '''
        Retrieve quantity values for multiple quantities at once.

        You can use the various quantities to search/filter for. For some of the
        indexed quantities this endpoint returns aggregation information. This means
        you will be given a list of all possible values and the number of entries
        that have the certain value. You can also use these aggregations on an empty
        search to determine the possible values.

        There is no ordering and no pagination and not after key based scrolling. Instead
        there is an 'after' key based scrolling.

        The result will contain a 'quantities' key with a dict of quantity names and the
        retrieved values as values.
        '''

        search_request = search.SearchRequest()
        args = {
            key: value
            for key, value in _repo_quantities_search_request_parser.parse_args().items()
            if value is not None}

        apply_search_parameters(search_request, args)
        quantities = args.get('quantities', [])
        size = args.get('size', 5)

        try:
            assert size >= 0
        except AssertionError:
            abort(400, message='invalid size')

        for quantity in quantities:
            try:
                search_request.quantity(quantity, size=size)
            except KeyError as e:
                import traceback
                traceback.print_exc()
                abort(400, 'Given quantity does not exist: %s' % str(e))

        return search_request.execute(), 200


_repo_calc_id_model = api.model('RepoCalculationId', {
    'upload_id': fields.String(), 'calc_id': fields.String()
})


@ns.route('/pid/<path:pid>')
class RepoPidResource(Resource):
    @api.doc('resolve_pid')
    @api.response(404, 'Entry with PID does not exist')
    @api.marshal_with(_repo_calc_id_model, skip_none=True, code=200, description='Entry resolved')
    @authenticate()
    def get(self, pid: str):
        if '/' in pid:
            prefix, pid = pid.split('/')
            if prefix != '21.11132':
                abort(400, 'Wrong PID format')
            try:
                pid_int = utils.decode_handle_id(pid)
            except ValueError:
                abort(400, 'Wrong PID format')
        else:
            try:
                pid_int = int(pid)
            except ValueError:
                abort(400, 'Wrong PID format')

        search_request = search.SearchRequest().include('upload_id', 'calc_id')

        if g.user is not None:
            search_request.owner('all', user_id=g.user.user_id)
        else:
            search_request.owner('all')

        search_request.search_parameter('pid', pid_int)

        results = list(search_request.execute_scan())
        total = len(results)

        if total == 0:
            abort(404, 'Entry with PID %s does not exist' % pid)

        if total > 1:
            common.logger.error('Two entries for the same pid', pid=pid_int)

        result = results[0]
        return dict(
            upload_id=result['upload_id'],
            calc_id=result['calc_id'])
