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

from flask import g
from flask_restplus import abort, Resource

from nomad import infrastructure, config

from .app import api
from .auth import login_really_required


ns = api.namespace('admin', description='Administrative operations')


@ns.route('/<string:operation>')
@api.doc(params={'operation': 'The operation to perform.'})
class AdminOperationsResource(Resource):
    # TODO in production this requires authorization
    @api.response(200, 'Operation performed')
    @api.response(404, 'Operation does not exist')
    @api.response(400, 'Operation not available/disabled')
    @login_really_required
    def post(self, operation):
        """
        Allows to perform administrative operations on the nomad services.

        The possible operations are ``reset`` and ``remove``.

        The ``reset`` operation will attempt to clear the contents of all databased and
        indices.

        The ``remove``operation will attempt to remove all databases. Expect the
        api to stop functioning after this request.

        Reset and remove can be disabled.
        """
        if g.user.email != 'admin':
            abort(401, message='Only the admin user can perform this operation.')

        if operation == 'reset':
            if config.services.disable_reset:
                abort(400, message='Operation is disabled')
            infrastructure.reset()
        elif operation == 'remove':
            if config.services.disable_reset:
                abort(400, message='Operation is disabled')
            infrastructure.remove()
        else:
            abort(404, message='Unknown operation %s' % operation)

        return dict(messager='Operation %s performed.' % operation), 200
