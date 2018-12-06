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

from flask_restful import abort

from nomad import infrastructure
from nomad.processing import Upload

from .app import app, base_path


# TODO in production this requires authorization
@app.route('%s/admin/<string:operation>' % base_path, methods=['POST'])
def call_admin_operation(operation):
    """
    Allows to perform administrative operations on the nomad services. The possible
    operations are *repair_uploads*
    (cleans incomplete or otherwise unexpectedly failed uploads), *reset* (clears all
    databases and resets nomad).

    .. :quickref: Allows to perform administrative operations on the nomad services.

    :param string operation: the operation to perform
    :status 400: unknown operation
    :status 200: operation successfully started
    :returns: an authentication token that is valid for 10 minutes.
    """
    if operation == 'repair_uploads':
        Upload.repair_all()
    if operation == 'reset':
        infrastructure.reset()
    else:
        abort(400, message='Unknown operation %s' % operation)

    return 'done', 200
