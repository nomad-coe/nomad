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

import logging
from typing import List, Dict, Tuple, Any, Optional, Union, cast, TypeVar
from pydantic import BaseModel, Field, root_validator, Extra  # pylint: disable=unused-import

ConfigBaseModelBound = TypeVar('ConfigBaseModelBound', bound='ConfigBaseModel')


class ConfigBaseModel(BaseModel, extra=Extra.ignore):
    """Customized base class that logs a warning when extra fields are specified."""

    def customize(
        self: ConfigBaseModelBound,
        custom_settings: Union[ConfigBaseModelBound, Dict[str, Any]],
    ) -> ConfigBaseModelBound:
        """
        Returns a new config object, created by taking a copy of the current config and
        updating it with the settings defined in `custom_settings`. The `custom_settings` can
        be a NomadSettings or a dictionary (in the latter case it must not contain any new keys
        (keys not defined in this NomadSettings). If it does, an exception will be raised.
        """

        rv = self.copy(deep=True)

        if custom_settings:
            if isinstance(custom_settings, BaseModel):
                for field_name in custom_settings.__fields__.keys():
                    try:
                        setattr(rv, field_name, getattr(custom_settings, field_name))
                    except Exception:
                        raise AssertionError(f'Invalid setting: {field_name}')
            elif isinstance(custom_settings, dict):
                for key, value in custom_settings.items():
                    if value is None:
                        continue
                    try:
                        setattr(rv, key, value)
                    except Exception:
                        raise AssertionError(f'Invalid setting: ({key}: {value})')

        return cast(ConfigBaseModelBound, rv)

    @root_validator(pre=True)
    def __print_extra_field__(cls, values):  # pylint: disable=no-self-argument
        extra_fields = values.keys() - cls.__fields__.keys()

        def list_items(items):
            return ', '.join([f'"{x}"' for x in items])

        if extra_fields:
            logger = logging.getLogger(__name__)
            logger.warning(
                f'The following unsupported keys were found in your configuration, '
                f'e.g. nomad.yaml: {list_items(extra_fields)}.'
            )

        return values


class OptionsBase(ConfigBaseModel):
    """The most basic model for defining the availability of different options."""

    include: Optional[List[str]] = Field(
        description="""
        List of included options. If not explicitly defined, all of the options will
        be included by default.
    """
    )
    exclude: Optional[List[str]] = Field(
        description="""
        List of excluded options. Has higher precedence than include.
    """
    )

    def filter(self, value: str) -> bool:
        """Determines is a value fitting this specification."""
        included = not self.include or value in self.include or '*' in self.include
        excluded = self.exclude and (value in self.exclude or '*' in self.exclude)

        return included and not excluded


class OptionsGlob(ConfigBaseModel):
    """Controls the availability of different options with the possibility of
    using glob/wildcard syntax.
    """

    include: Optional[List[str]] = Field(
        description="""
        List of included options. Supports glob/wildcard syntax.
    """
    )
    exclude: Optional[List[str]] = Field(
        description="""
        List of excluded options. Supports glob/wildcard syntax. Has higher precedence than include.
    """
    )


class Options(OptionsBase):
    """Common configuration class used for enabling/disabling certain
    elements and defining the configuration of each element.
    """

    options: Optional[Dict[str, Any]] = Field(
        {}, description='Contains the available options.'
    )

    def filtered_keys(self) -> List[str]:
        """Returns a list of keys that fullfill the include/exclude
        requirements.
        """
        if self.include is None or '*' in self.include:
            include = list(self.options.keys())
        else:
            include = self.include
        if self.exclude is not None and '*' in self.exclude:
            return []
        else:
            exclude = self.exclude or []
        return [key for key in include if key not in exclude]

    def filtered_values(self) -> List[Any]:
        """Returns a list of values that fullfill the include/exclude
        requirements.
        """
        return [
            self.options[key] for key in self.filtered_keys() if key in self.options
        ]

    def filtered_items(self) -> List[Tuple[str, Any]]:
        """Returns a list of key/value pairs that fullfill the include/exclude
        requirements.
        """
        return [
            (key, self.options[key])
            for key in self.filtered_keys()
            if key in self.options
        ]


class OptionsSingle(Options):
    """Represents options where one value can be selected."""

    selected: str = Field(description='Selected option.')


class OptionsMulti(Options):
    """Represents options where multiple values can be selected."""

    selected: List[str] = Field(description='Selected options.')
