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

import importlib

from nomad.config import config
from nomad.config.models.plugins import Schema

from . import (
    annotations,
)  # Should be imported first to register the annotations before they are used
from .simulation import m_env
from .downloads import m_package
from .eln.labfolder import m_package
from .eln.openbis import m_package
from .plot import m_package


class SchemaInterface:
    def __init__(self, path) -> None:
        self._module = None
        self._path = path

    @property
    def module(self):
        if self._module is None:
            self._module = importlib.import_module(self._path)
        return self._module

    def __getattr__(self, name: str):
        return getattr(self.module, name, None)


simulationworkflowschema, runschema = None, None
for plugin in config.plugins.filtered_values():
    if isinstance(plugin, Schema):
        if plugin.name == 'simulationworkflowschema':
            simulationworkflowschema = SchemaInterface(plugin.python_package)
        elif plugin.name == 'runschema':
            runschema = SchemaInterface(plugin.python_package)
        else:
            importlib.import_module(plugin.python_package)

SCHEMA_IMPORT_ERROR = 'Schema not defined.'
