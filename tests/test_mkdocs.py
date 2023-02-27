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

from enum import Enum
import pytest
from typing import Union, List, Dict, Optional
from typing_extensions import Literal, Annotated  # type: ignore
from pydantic import BaseModel, Field
from nomad.config.models import WidgetHistogram, WidgetTerms
from nomad.mkdocs import (
    get_field_type_info,
    get_field_description,
    get_field_default,
    get_field_options,
    get_field_deprecated
)


class MyStrEnum(str, Enum):
    TEST = 'test'


class MyIntEnum(int, Enum):
    TEST = 1


@pytest.mark.parametrize(
    'type_, name, classes',
    [
        pytest.param(str, 'str', {str}, id='str'),
        pytest.param(int, 'int', {int}, id='int'),
        pytest.param(float, 'float', {float}, id='float'),
        pytest.param(list, 'list', {list}, id='list'),
        pytest.param(dict, 'dict', {dict}, id='dict'),
        pytest.param(set, 'set', {set}, id='set'),
        pytest.param(tuple, 'tuple', {tuple}, id='tuple'),
        pytest.param(WidgetHistogram, 'WidgetHistogram', {WidgetHistogram}, id='basemodel'),
        pytest.param(Enum, 'Enum', {Enum}, id='class'),
        pytest.param(Optional[WidgetHistogram], 'WidgetHistogram', {WidgetHistogram}, id='optional-ignored'),
        pytest.param(Union[str, WidgetHistogram], 'Union[str, WidgetHistogram]', {Union, str, WidgetHistogram}, id='union'),
        pytest.param(List[Union[str, WidgetHistogram]], 'List[Union[str, WidgetHistogram]]', {list, Union, str, WidgetHistogram}, id='list-with-union'),
        pytest.param(Dict[str, WidgetHistogram], 'Dict[str, WidgetHistogram]', {dict, str, WidgetHistogram}, id='dict'),
        pytest.param(Literal["test"], 'str', {Literal}, id='literal-not-shown-str'),
        pytest.param(Literal[1], 'int', {Literal}, id='literal-not-shown-int'),
        pytest.param(MyStrEnum, 'str', {MyStrEnum}, id='enum-string'),
        pytest.param(MyIntEnum, 'int', {MyIntEnum}, id='enum-int'),
        pytest.param(
            List[Annotated[Union[WidgetTerms, WidgetHistogram], Field(discriminator="type")]],  # type: ignore
            'List[Union[WidgetTerms, WidgetHistogram]]',
            {list, Union, WidgetTerms, WidgetHistogram},
            id='annotated-ignored'
        ),
    ]
)
def test_field_type_info(type_, name, classes):

    class Test(BaseModel):
        a: type_ = Field()

    name_found, classes_found = get_field_type_info(Test.__fields__['a'])
    assert name_found == name
    assert classes_found == classes


@pytest.mark.parametrize(
    'description',
    [
        pytest.param(None, id='no-description'),
        pytest.param('This is a test description.', id='string-description'),
    ]
)
def test_field_description(description):

    class Test(BaseModel):
        a: str = Field(description=description)

    description_found = get_field_description(Test.__fields__['a'])
    assert description_found == description


@pytest.mark.parametrize(
    'default, default_str',
    [
        pytest.param(None, None, id='no-default'),
        pytest.param('test', '`test`', id='str-default'),
        pytest.param(1, '`1`', id='int-default'),
        pytest.param({'test': 'test'}, 'Complex object, default value not displayed.', id='complex-default'),
    ]
)
def test_field_default(default, default_str):

    class Test(BaseModel):
        a: str = Field(default)

    default_found = get_field_default(Test.__fields__['a'])
    assert default_found == default_str


@pytest.mark.parametrize(
    'type_, options',
    [
        pytest.param(str, {}, id='no-options'),
        pytest.param(MyStrEnum, {'test': None}, id='str-options'),
        pytest.param(MyIntEnum, {'1': None}, id='int-options'),
    ]
)
def test_field_options(type_, options):

    class Test(BaseModel):
        a: type_ = Field()

    options_found = get_field_options(Test.__fields__['a'])
    assert len(options_found) == len(options)
    for key in options_found:
        assert options_found[key] == options[key]


@pytest.mark.parametrize('deprecated', [True, False])
def test_field_deprecated(deprecated):

    class Test(BaseModel):
        a: str = Field(deprecated=deprecated)

    deprecated_found = get_field_deprecated(Test.__fields__['a'])
    assert deprecated_found == deprecated
