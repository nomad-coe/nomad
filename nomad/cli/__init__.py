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

'''
Command line interface (CLI) for nomad. Provides a group/sub-command structure, think git,
that offers various functionality to the command line user.

Use it from the command line with ``nomad --help`` or ``python -m nomad.cli --help`` to learn
more.

The CLI uses lazy_import for lazy loading modules. This has some limitations. You will
break lazy loading if an ``from x import y`` is used in the cli code. You will also
have to add imports via :func:`nomad.cli.lazy_import.lazy_module` before importing
them.
'''
from nomad.cli import lazy_import

lazy_import.lazy_module('click')
lazy_import.lazy_module('logging')
lazy_import.lazy_module('os')
lazy_import.lazy_module('typing')
lazy_import.lazy_module('json')
lazy_import.lazy_module('sys')
lazy_import.lazy_module('nomad.config')
lazy_import.lazy_module('nomad.infrastructure')
lazy_import.lazy_module('nomad.utils')
lazy_import.lazy_module('nomad.parsing')
lazy_import.lazy_module('nomad.parsing.parsers')
lazy_import.lazy_module('nomad.normalizing')
lazy_import.lazy_module('nomad.datamodel')
lazy_import.lazy_module('nomad.search')
lazy_import.lazy_module('nomad.metainfo')
lazy_import.lazy_module('nomad.atomutils')
lazy_import.lazy_module('nomad.processing')
lazy_import.lazy_module('nomad.client')

from . import dev, admin, parse, client  # noqa
from .cli import run_cli, cli  # noqa
