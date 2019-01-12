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

from elasticsearch.exceptions import NotFoundError
from flask import g, request
from flask_restplus import Resource, abort, fields

from nomad.repo import RepoCalc

from .app import api
from .auth import login_if_available
from .common import pagination_model, pagination_request_parser, calc_route

ns = api.namespace('repo', description='Access repository metadata, edit user metadata.')


@calc_route(ns)
class RepoCalcResource(Resource):
    @api.response(404, 'The upload or calculation does not exist')
    @api.response(200, 'Metadata send')
    @api.doc('get_repo_calc')
    def get(self, upload_id, calc_hash):
        """
        Get calculation metadata in repository form.

        Repository metadata only entails the quanties shown in the repository.
        This is basically the elastic search index entry for the
        requested calculations. Calcs are references via *upload_id*, *calc_hash*
        pairs.
        """
        try:
            return RepoCalc.get(id='%s/%s' % (upload_id, calc_hash)).json_dict, 200
        except NotFoundError:
            abort(404, message='There is no calculation for %s/%s' % (upload_id, calc_hash))
        except Exception as e:
            abort(500, message=str(e))


repo_calcs_model = api.model('RepoCalculations', {
    'pagination': fields.Nested(pagination_model),
    'results': fields.List(fields.Raw)
})

repo_request_parser = pagination_request_parser.copy()
repo_request_parser.add_argument(
    'owner', type=str,
    help='Specify which calcs to return: ``all``, ``user``, ``staging``, default is ``all``')


@ns.route('/')
class RepoCalcsResource(Resource):
    @api.doc('get_calcs')
    @api.response(400, 'Invalid requests, e.g. wrong owner type')
    @api.expect(repo_request_parser, validate=True)
    @api.marshal_with(repo_calcs_model, skip_none=True, code=200, description='Metadata send')
    @login_if_available
    def get(self):
        """
        Get *'all'* calculations in repository from, paginated.
        """
        page = int(request.args.get('page', 1))
        per_page = int(request.args.get('per_page', 10))
        owner = request.args.get('owner', 'all')

        try:
            assert page >= 1
            assert per_page > 0
        except AssertionError:
            abort(400, message='invalid pagination')

        if owner == 'all':
            search = RepoCalc.search().query('match_all')
        elif owner == 'user':
            if g.user is None:
                abort(401, message='Authentication required for owner value user.')
            search = RepoCalc.search().query('match_all')
            search = search.filter('term', user_id=str(g.user.user_id))
        elif owner == 'staging':
            if g.user is None:
                abort(401, message='Authentication required for owner value user.')
            search = RepoCalc.search().query('match_all')
            search = search.filter('term', user_id=str(g.user.user_id)).filter('term', staging=True)
        else:
            abort(400, message='Invalid owner value. Valid values are all|user|staging, default is all')

        search = search[(page - 1) * per_page: page * per_page]
        return {
            'pagination': {
                'total': search.count(),
                'page': page,
                'per_page': per_page
            },
            'results': [result.json_dict for result in search]
        }, 200
