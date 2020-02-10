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

from structlog import BoundLogger
from flask_restplus import fields
from datetime import datetime
import pytz

from nomad import config


logger: BoundLogger = None
""" A logger pre configured with information about the current request. """

base_path = config.services.api_base_path
""" Provides the root path of the nomad APIs. """


class RFC3339DateTime(fields.DateTime):

    def format(self, value):
        if isinstance(value, datetime):
            return super().format(value.replace(tzinfo=pytz.utc))
        else:
            return str(value)


rfc3339DateTime = RFC3339DateTime()
