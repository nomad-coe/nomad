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
Swagger/bravado based python client library for the API and various usefull shell commands.
"""

from . import local, migration, upload, integrationtests, parse
from .main import cli, create_client
from .upload import stream_upload_with_client


def run_cli():
    cli()  # pylint: disable=E1120
