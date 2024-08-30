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
#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD.
# See https://nomad-lab.eu for further info.
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
from nomad.utils import get_logger
from nomad.metainfo import Environment
from .run import Run
from .calculation import Calculation
from .method import Method
from .system import System
from . import run
from . import method
from . import calculation
from . import system
from . import workflow
from . import legacy_workflows

# Disable warning as it is always imported.
# get_logger(__name__).warning('Schema is deprecated, use plugins.')

m_env = Environment()
m_env.m_add_sub_section(Environment.packages, run.m_package)
m_env.m_add_sub_section(Environment.packages, method.m_package)
m_env.m_add_sub_section(Environment.packages, calculation.m_package)
m_env.m_add_sub_section(Environment.packages, system.m_package)
m_env.m_add_sub_section(Environment.packages, workflow.m_package)  # noqa
m_env.m_add_sub_section(Environment.packages, legacy_workflows.m_package)  # noqa
