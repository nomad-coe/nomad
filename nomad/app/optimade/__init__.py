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

from .filterparser import parse_filter

"""
The optimade implementation of NOMAD.
"""

blueprint = Blueprint('optimade', __name__)

api = Api(
    blueprint,
    version='1.0', title='NOMAD optimade PI',
    description='The NOMAD optimade API',
    validate=True)
""" Provides the flask restplust api instance for the optimade api"""
