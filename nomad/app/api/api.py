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

from flask import Blueprint
from flask_restplus import Api


blueprint = Blueprint('api', __name__)

api = Api(
    blueprint,
    version='1.0', title='NOMAD API',
    description='Official NOMAD API',
    validate=True)
''' Provides the flask restplus api instance for the regular NOMAD api'''


# For some unknown reason it is necessary for each fr api to have a handler.
# Otherwise the global app error handler won't be called.
@api.errorhandler(Exception)
def errorhandler(error):
    '''When an internal server error is caused by an unexpected exception.'''
    return str(error)
