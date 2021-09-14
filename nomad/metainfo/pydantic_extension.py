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

'''
Allows to create pydantic models from section definitions.
'''

from typing import cast
from pydantic import create_model, Field, BaseConfig
from datetime import datetime

from .metainfo import DefinitionAnnotation, Definition, Section, Quantity, Datetime, MEnum, Capitalized, JSON


class _OrmConfig(BaseConfig):
    orm_mode = True


class PydanticModel(DefinitionAnnotation):
    '''
    This annotation class can be used to extend metainfo sections. It will create a
    pydantic model from the section definition. Its a SectionAnnotation and allows
    to create pydantic model instances from section instances.

    Attributes:
        model: The pydantic model that represents the section defintion.
    '''
    def __init__(self):
        self.model = None

    def to_pydantic(self, section):
        ''' Returns the pydantic model instance for the given section. '''
        return self.model.from_orm(section)

    def init_annotation(self, definition: Definition):
        section_definition = cast(Section, definition)
        name = section_definition.name

        def create_field(quantity: Quantity):
            pydantic_type: type = None
            if quantity.type == Datetime:
                pydantic_type = datetime
            elif isinstance(quantity.type, MEnum):
                pydantic_type = str
            elif quantity.type == Capitalized:
                pydantic_type = str
            elif quantity.type == JSON:
                pydantic_type = dict
            else:
                pydantic_type = quantity.type

            return pydantic_type, Field(quantity.default, description=quantity.description)

        fields = {
            name: create_field(quantity)
            for name, quantity in section_definition.all_quantities.items()
        }

        self.model = create_model(name, __config__=_OrmConfig, **fields)
