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
import urllib.parse

from nomad import config

blueprint = Blueprint('optimade', __name__)

base_url = 'https://%s/%s/optimade' % (
    config.services.api_host.strip('/'),
    config.services.api_base_path.strip('/'))


def url(endpoint: str = None, version='v1', prefix=None, **kwargs):
    ''' Returns the full optimade api url (for a given endpoint) including query parameters. '''
    if endpoint is not None:
        url = '/' + endpoint
    else:
        url = ''

    if version is not None:
        url = '/' + version + url

    if prefix is not None:
        url = '/' + prefix + url

    url = base_url + url

    if len(kwargs) > 0:
        return '%s?%s' % (url, urllib.parse.urlencode(kwargs))
    else:
        return url


api = Api(
    blueprint,
    version='1.0', title='NOMAD\'s OPTiMaDe API implementation',
    description='NOMAD\'s OPTiMaDe API implementation, version 1.0.0.',
    validate=True)
''' Provides the flask restplust api instance for the optimade api'''


# For some unknown reason it is necessary for each fr api to have a handler.
# Otherwise the global app error handler won't be called.
@api.errorhandler(Exception)
def errorhandler(error):
    '''When an internal server error is caused by an unexpected exception.'''
    return str(error)
