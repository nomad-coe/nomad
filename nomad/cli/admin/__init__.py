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

from nomad.cli import lazy_import

lazy_import.lazy_module('typing')
lazy_import.lazy_module('click')
lazy_import.lazy_module('asyncio')
lazy_import.lazy_module('concurrent.futures')
lazy_import.lazy_module('datetime')
lazy_import.lazy_module('time')
lazy_import.lazy_module('elasticsearch_dsl')
lazy_import.lazy_module('elasticsearch')
lazy_import.lazy_module('os')
lazy_import.lazy_module('shutil')
lazy_import.lazy_module('tabulate')
lazy_import.lazy_module('sys')
lazy_import.lazy_module('io')
lazy_import.lazy_module('re')
lazy_import.lazy_module('uuid')
lazy_import.lazy_module('json')
lazy_import.lazy_module('threading')
lazy_import.lazy_module('numpy')
lazy_import.lazy_module('requests')
lazy_import.lazy_module('pymongo')
lazy_import.lazy_module('mongoengine')
lazy_import.lazy_module('ase')
lazy_import.lazy_module('bs4')
lazy_import.lazy_module('matid')
lazy_import.lazy_module('matid.symmetry.symmetryanalyzer')
lazy_import.lazy_module('matid.utils.segfault_protect')
lazy_import.lazy_module('nomad.atomutils')
lazy_import.lazy_module('nomad.normalizing')
lazy_import.lazy_module('nomad.processing')
lazy_import.lazy_module('nomad.search')
lazy_import.lazy_module('nomad.datamodel')
lazy_import.lazy_module('nomad.infrastructure')
lazy_import.lazy_module('nomad.utils')
lazy_import.lazy_module('nomad.config')
lazy_import.lazy_module('nomad.files')
lazy_import.lazy_module('nomad.archive')

from . import admin, uploads, entries, run, clean, users, migrate  # noqa
