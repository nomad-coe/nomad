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

from nomad.files import UploadFiles, Restricted

from .app import api
from .auth import login_if_available, create_authorization_predicate
from .common import pagination_model, pagination_request_parser, calc_route

ns = api.namespace('repo', description='Access repository metadata.')


@calc_route(ns)
class RepoCalcResource(Resource):
    @api.response(404, 'The upload or calculation does not exist')
    @api.response(401, 'Not authorized to access the calculation')
    @api.response(200, 'Metadata send')
    @api.doc('get_repo_calc')
    @login_if_available
    def get(self, upload_id, calc_id):
        """
        Get calculation metadata in repository form.

        Repository metadata only entails the quanties shown in the repository.
        Calcs are references via *upload_id*, *calc_id* pairs.
        """
        # TODO use elastic search instead of the files
        upload_files = UploadFiles.get(upload_id, create_authorization_predicate(upload_id, calc_id))
        if upload_files is None:
            abort(404, message='There is no upload %s' % upload_id)

        try:
            return upload_files.metadata.get(calc_id), 200
        except Restricted:
            abort(401, message='Not authorized to access %s/%s.' % (upload_id, calc_id))
        except KeyError:
            abort(404, message='There is no calculation for %s/%s' % (upload_id, calc_id))


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

        This is currently not implemented!
        """
        return dict(pagination=dict(total=0, page=1, per_page=10), results=[]), 200
        # page = int(request.args.get('page', 1))
        # per_page = int(request.args.get('per_page', 10))
        # owner = request.args.get('owner', 'all')

        # try:
        #     assert page >= 1
        #     assert per_page > 0
        # except AssertionError:
        #     abort(400, message='invalid pagination')

        # if owner == 'all':
        #     search = RepoCalc.search().query('match_all')
        # elif owner == 'user':
        #     if g.user is None:
        #         abort(401, message='Authentication required for owner value user.')
        #     search = RepoCalc.search().query('match_all')
        #     search = search.filter('term', user_id=str(g.user.user_id))
        # elif owner == 'staging':
        #     if g.user is None:
        #         abort(401, message='Authentication required for owner value user.')
        #     search = RepoCalc.search().query('match_all')
        #     search = search.filter('term', user_id=str(g.user.user_id)).filter('term', staging=True)
        # else:
        #     abort(400, message='Invalid owner value. Valid values are all|user|staging, default is all')

        # search = search[(page - 1) * per_page: page * per_page]
        # return {
        #     'pagination': {
        #         'total': search.count(),
        #         'page': page,
        #         'per_page': per_page
        #     },
        #     'results': [result.json_dict for result in search]
        # }, 200
