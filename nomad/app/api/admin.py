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

from flask import request
from flask_restplus import abort, Resource, fields

from nomad import infrastructure, config

from .api import api
from .auth import admin_login_required


ns = api.namespace('admin', description='Administrative operations')


@ns.route('/reset')
class AdminResetResource(Resource):
    @api.doc('exec_reset_command')
    @api.response(200, 'Reset performed')
    @api.response(400, 'Reset not available/disabled')
    @admin_login_required
    def post(self):
        """
        The ``reset`` command will attempt to clear the contents of all databased and
        indices.

        Nomad can be configured to disable reset and the command might not be available.
        """
        if config.services.disable_reset:
            abort(400, message='Operation is disabled')

        infrastructure.reset(repo_content_only=True)

        return dict(messager='Reset performed.'), 200


@ns.route('/remove')
class AdminRemoveResource(Resource):
    @api.doc('exec_remove_command')
    @api.response(200, 'Remove performed')
    @api.response(400, 'Remove not available/disabled')
    @admin_login_required
    def post(self):
        """
        The ``remove``command will attempt to remove all databases. Expect the
        api to stop functioning after this request.

        Nomad can be configured to disable remove and the command might not be available.
        """

        if config.services.disable_reset:
            abort(400, message='Operation is disabled')

        infrastructure.remove()

        return dict(messager='Remove performed.'), 200


pidprefix_model = api.model('PidPrefix', {
    'prefix': fields.Integer(description='The prefix. All new calculations will get an id that is greater.', required=True)
})


# TODO remove after migration
@ns.route('/pidprefix')
class AdminPidPrefixResource(Resource):
    @api.doc('exec_pidprefix_command')
    @api.response(200, 'Pid prefix set')
    @api.response(400, 'Bad pid prefix data')
    @api.expect(pidprefix_model)
    @admin_login_required
    def post(self):
        """
        The ``pidprefix``command will set the pid counter to the given value.

        This might be useful while migrating data with old pids.
        """

        infrastructure.set_pid_prefix(**request.get_json())

        return dict(messager='PID prefix set.'), 200
