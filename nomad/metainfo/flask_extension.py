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

from flask_restplus import fields

from nomad.app.flask.common import RFC3339DateTime

from .metainfo import Section, Quantity, Datetime, Capitalized, MEnum, JSON


def field(quantity: Quantity):
    ''' Returns a flask restplus field with quantity type and shape. '''
    field = None
    if quantity.type == int:
        field = fields.Integer
    elif quantity.type == float:
        field = fields.Float
    elif quantity.type == str or quantity.type == Capitalized:
        field = fields.String
    elif quantity.type == bool:
        field = fields.Boolean
    elif quantity.type == Datetime:
        field = RFC3339DateTime
    elif isinstance(quantity.type, MEnum):
        field = fields.String
    elif quantity.type == JSON:
        field = fields.Arbitrary
    else:
        raise NotImplementedError

    result = field(description=quantity.description)

    if len(quantity.shape) == 0:
        return result
    elif len(quantity.shape) == 1:
        return fields.List(result)
    else:
        raise NotImplementedError


def generate_flask_restplus_model(api, section_def: Section):
    return api.model(section_def.name, {
        name: field(quantity)
        for name, quantity in section_def.all_quantities.items()
    })
