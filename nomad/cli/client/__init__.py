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

from nomad.cli import lazy_import

lazy_import.lazy_module('sys')
lazy_import.lazy_module('os')
lazy_import.lazy_module('io')
lazy_import.lazy_module('re')
lazy_import.lazy_module('shutil')
lazy_import.lazy_module('requests')
lazy_import.lazy_module('click')
lazy_import.lazy_module('typing')
lazy_import.lazy_module('datetime')
lazy_import.lazy_module('subprocess')
lazy_import.lazy_module('tarfile')
lazy_import.lazy_module('threading')
lazy_import.lazy_module('json')
lazy_import.lazy_module('bravado.exception')
lazy_import.lazy_module('bravado.requests_client')
lazy_import.lazy_module('bravado.client')
lazy_import.lazy_module('urllib.parse')
lazy_import.lazy_module('keycloak')
lazy_import.lazy_module('time')
lazy_import.lazy_module('matplotlib.scale')
lazy_import.lazy_module('matplotlib.pyplot')
lazy_import.lazy_module('matplotlib.ticker')
lazy_import.lazy_module('matplotlib.transforms')
lazy_import.lazy_module('numpy')
lazy_import.lazy_module('nomad.config')
lazy_import.lazy_module('nomad.utils')
lazy_import.lazy_module('nomad.processing')
lazy_import.lazy_module('nomad.files')
lazy_import.lazy_module('nomad.search')
lazy_import.lazy_module('nomad.datamodel')
lazy_import.lazy_module('nomad.parsing')
lazy_import.lazy_module('nomad.parsing.parsers')
lazy_import.lazy_module('nomad.infrastructure')
lazy_import.lazy_module('nomad.doi')
lazy_import.lazy_module('nomad.client')

from . import local, upload, integrationtests, mirror, statistics, update_database  # noqa
from .client import create_client  # noqa
from .upload import stream_upload_with_client  # noqa
