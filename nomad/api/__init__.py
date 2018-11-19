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

"""
This module comprises the nomad@FAIRDI APIs.

The different APIs are upload, repository (raw data and search), and archive.

There is a separate documentation for the API endpoints from a client perspective.

.. autodata:: app

.. automodule:: nomad.api.app
.. automodule:: nomad.api.upload
.. automodule:: nomad.api.repository
.. automodule:: nomad.api.archive
"""

from .app import app
from . import upload, repository, archive
